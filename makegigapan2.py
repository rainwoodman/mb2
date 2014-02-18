import sharedmem
from gaepsi.render import *
from matplotlib import cm
from sys import argv
from gaepsi.snapshot import Snapshot
from gaepsi.gaplot import create_cosmology
import os.path
from gaepsi.compiledbase.geometry import Cubenoid
from sys import argv

execfile(argv[1])
# center = [ 62995.8671875,  40180.12109375, 52267.77734375]
# pos of most massive halo from halo catalog
# unfold from here.

#snapname = '/physics/yfeng1/e5/run/000/output/snapdir_1366/snapshot_1366.%d'
snap = Snapshot(snapname % 0, 'cmugadget')
cosmology = create_cosmology(snap.C)
Nfile = snap.C['Nfiles']
boxsize = snap.C['boxsize']

cub = Cubenoid(M, [0, 0, 0], [boxsize, boxsize, boxsize], 
        center=0, neworigin=0)
ratio = cub.newboxsize[0] / cub.newboxsize[1]
width = int(width * dpi)
length = int(ratio * width)
print 'real length is', 1.0 * length / dpi, 'inch'
print 'kpc per inch is', cub.newboxsize[0] / length * dpi
print 'kpc per dot is', cub.newboxsize[0] / length
print 'pixel size is wxl', width, length
print 'total pixels = ', width * length / 1024 / 1024., 'M'

ccdbuffer = ccdbuffer % (width, length)
imagename = imagename % (width, length)

def main():
    camera = Camera(int(length), int(width))
    camera.fade = False

    camera.lookat(pos=[
        cub.newboxsize[0] * 0.5, 
        cub.newboxsize[1] * 0.5, 
        cub.newboxsize[2]], 
        target=[
        cub.newboxsize[0] * 0.5, 
        cub.newboxsize[1] * 0.5, 
        cub.newboxsize[2] * 0.5],
        up=[0, 1, 0])
    camera.ortho(near=0, far=cub.newboxsize[2], 
            extent=[-cub.newboxsize[0] * 0.5, 
                     cub.newboxsize[0] * 0.5, 
                    -cub.newboxsize[1] * 0.5,
                     cub.newboxsize[1] * 0.5])
    if not os.path.exists(ccdbuffer):
        CCD = numpy.memmap(ccdbuffer, mode='w+', dtype=('f4', 2), shape=camera.shape)
        CCD[...] = 0
        def prefetcher(fid):
            print 'pref', fid
            try:
                snap = Snapshot(snapname % fid, 'cmugadget')
                snap[0, 'pos']
                snap[0, 'ie']
                snap[0, 'ye']
                snap[0, 'sml']
                snap[0, 'mass']
            except:
                pass
            return 
        chunksize = 64
        def prefetch(fid):
            loadnext = [ 
                    sharedmem.background(
                        lambda fidnext: prefetcher(fidnext),
                        fidnext)
                    for fidnext in range(fid, fid + chunksize) if fidnext < Nfile ]
            return loadnext
        loadnext = prefetch(0)
        for fid in range(0, Nfile, chunksize):
            print fid, 'start'
            for thread in loadnext:
                thread.wait()
            del loadnext
            nextsnap = [Snapshot(snapname % f, 'cmugadget')
                    for f in range(fid, fid + chunksize) if f < Nfile]
            loadnext = prefetch(fid + chunksize)
            makegigapan(nextsnap, camera, CCD)
            del nextsnap
            print fid, 'done'
        CCD[..., 0] /= CCD[..., 1]
    else:
        CCD = numpy.memmap(ccdbuffer, mode='r', dtype=('f4', 2))
    CCD = CCD.reshape(-1, 2)


    print 'CCD done'
    C, L = CCD[..., 0], CCD[..., 1]
    C = C.reshape(-1)
    L = L.reshape(-1)

    with sharedmem.MapReduce() as pool:
        chunksize = 1024 * 1024 * 20
        def work(i):
            sl = slice(i, i + chunksize)
            l = nl_(L[sl])
            return l.vmax, l.vmin
        vmax, vmin = numpy.array(pool.map(work, range(0, len(L), chunksize))).T
        
    vmax = numpy.max(vmax) - 0.2
    vmin = numpy.min(vmin)

    print vmax, vmin

    im = numpy.memmap(imagename, 
            mode='w+', dtype=('u1', 3), shape=L.shape)
    with sharedmem.MapReduce() as pool:
        chunksize = 1024 * 1024 * 128
        def work(i):
            sl = slice(i, i + chunksize)
            im[sl] = image(nl_(C[sl], vmin=3.0, vmax=8.0), nl_(L[sl],
                vmax=vmax, vmin=vmin) ** 0.8, cmap=cm.RdBu_r)[..., :3]
        pool.map(work, range(0, len(L), chunksize))
    im.flush()

def makegigapan(snaps, camera, CCD):
    x = []
    y = []
    z = []
    T = []
    sml = []
    mass = []
    print 'reading'
    Len = numpy.array([snap.C['N'][0] for snap in snaps], dtype='intp')
    End = Len.cumsum()
    Start = End.copy()
    Start[1:] = End[:-1]
    Start[0] = 0
    N = Len.sum()

    x = sharedmem.empty(N, dtype='f4')
    y = sharedmem.empty(N, dtype='f4')
    z = sharedmem.empty(N, dtype='f4')
    T = sharedmem.empty(N, dtype='f4')
    sml = sharedmem.empty(N, dtype='f4')
    mass = sharedmem.empty(N, dtype='f4')

    with sharedmem.MapReduce() as pool:
        def work(i):
            sl = slice(Start[i], End[i])
            snap = snaps[i]
            x[sl] = snap[0, 'pos'][:, 0]
            y[sl] = snap[0, 'pos'][:, 1]
            z[sl] = snap[0, 'pos'][:, 2]
            cub.apply(x[sl], y[sl], z[sl])
            ie = snap[0, 'ie']
            ye = snap[0, 'ye']
            T[sl] = cosmology.ie2T(ie=ie, ye=ye, Xh=0.76)
            sml[sl] = snap[0, 'sml']
            mass[sl] = snap[0, 'mass']
        pool.map(work, range(len(snaps)))

    print 'painting'
    paint((x, y, z), T, mass, sml, camera, CCD, normalize=False, direct_write=True)

def makeblackholes(bhfile):
    bhs = numpy.load(bhfile)
    mass = bhs['mass'] # avoid crazy accretion mass 
    x, y, z = numpy.copy(bhs['pos']).T
    cub.apply(x, y, z)
    
    uvt = camera.transform(x, y, z)
    u = uvt[..., 0] 
    v = uvt[..., 1] 
    numpy.savez(bhcat, u=u, v=v, mass=bhs['mass'])


main()
