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
        CCD = numpy.empty(camera.shape, dtype=('f4', 2))
        CCD[...] = 0
        def prefetch(fid):
            snap = Snapshot(snapname % fid, 'cmugadget')
            snap[0, 'pos']
            snap[0, 'ie']
            snap[0, 'ye']
            snap[0, 'sml']
            snap[0, 'mass']
            return 
        nextsnap = prefetch(0)
        for fid in range(Nfile):
            if fid < Nfile - 1:
                fidnext = fid + 1
                loadnext = sharedmem.background(
                        lambda : prefetch(fidnext))
            print fid, 'start'
            nextsnap = Snapshot(snapname % fid, 'cmugadget')
            makegigapan(nextsnap, camera, CCD)
            print fid, 'done'
            if fid < Nfile - 1 :
                loadnext.wait() 
        CCD[..., 0] /= CCD[..., 1]
        CCD.tofile(ccdbuffer)
    else:
        CCD = numpy.fromfile(ccdbuffer, dtype=('f4', 2))\
                .reshape(camera.shape[0], camera.shape[1], 2)


    C, L = CCD[..., 0], CCD[..., 1]
    l = nl_(L)
    vmax = l.vmax - 0.2
    vmin = l.vmin #+ 0.2
    
    im = numpy.memmap(imagename, 
            mode='w+', dtype=('u1', 3), shape=L.shape)
    with sharedmem.MapReduce() as pool:
        def work(i):
            im[i] = image(nl_(C[i], vmin=3.0, vmax=8.0), nl_(L[i],
                vmax=vmax, vmin=vmin) ** 0.8, cmap=cm.rainbow)[..., :3]
        pool.map(work, range(len(L)))
    im.flush()

def makegigapan(snap, camera, CCD):
    x, y, z = sharedmem.copy(snap[0, 'pos']).T
    ie = snap[0, 'ie']
    ye = snap[0, 'ye']
    T = sharedmem.empty_like(ye)
    with sharedmem.MapReduce() as pool:
        chunksize = 1024 * 1024
        def work(i):
            sl = slice(i, i+chunksize)
            cub.apply(x[sl], y[sl], z[sl])
            cosmology.ie2T(ie=ie[sl], ye=ye[sl], Xh=0.76, out=T[sl])
        pool.map(work, range(0, len(x), chunksize))

    sml = snap[0, 'sml']
    mass = snap[0, 'mass']
    paint((x, y, z), T, mass, sml, camera, CCD, normalize=False)

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
