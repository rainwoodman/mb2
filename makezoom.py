from gaepsi.gaplot import *
from lib import *
from matplotlib.figure import Figure
from matplotlib import cm
from sys import argv
from sys import stdin
import os.path
from StringIO import StringIO
# center = [ 62995.8671875,  40180.12109375, 52267.77734375]
# pos of most massive halo from halo catalog
# unfold from here.
M = numpy.matrix("4,0,1; 1,2,-1; 2,1,0")
zooms = [50000, 30000, 20000, 10000, 5000, 3000, 1000]

print 'started'

def makezoom(snap, context, centers, zooms, unfold):
    print context.fids()
    fids = context.fids()

    C = sharedmem.empty((len(centers), len(zooms), 
        context.shape[0], context.shape[1]), dtype='f8')
    L = sharedmem.empty((len(centers), len(zooms), 
        context.shape[0], context.shape[1]), dtype='f8')
    with sharedmem.Pool() as pool:
        def work(fid):
            ct = GaplotContext(context.snapname,
                    context.format, 
                    periodic=True, np=4)
            ct.reshape(context.shape)
            ct.schema('gas', 0, ['mass', 'ie', 'ye', 'sml'])
            with pool.critical:
                gas = ct.read('gas',fids=[fid])
            print 'read', len(gas), 'particles'
            if unfold:
                ct.unfold(M);
            ct.makeT('gas')
            print gas['T'].max(), gas['T'].min()
            for j, c in enumerate(context['centers']['pos']):
                for i, size in enumerate(zooms):
                    ct.view(center=c, size=size)
                    C_, L_ = ct.paint('gas', 'T', 'mass', 'sml')
                    C_[L_ == 0] = 0
                    with pool.critical:
                        C[j, i] += C_ * L_
                        L[j, i] += L_
        pool.map(work, fids)
        C[...] /= L
        return C, L

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('snap', type=SnapDir)
#    parser.add_argument('mode', choices=['make', 'continue'])
    parser.add_argument('prefix')
    parser.add_argument('--unfold', action='store_true', default=False)
    args = parser.parse_args()
    snap = args.snap
    snapid = snap.snapid
    prefix = args.prefix

    centers = numpy.loadtxt(stdin)
#    centers = [[ 37384.87890625,  40500.9375    ,  60814.69921875]]

    context = GaplotContext()
    context.reshape((2000, 2000))
    context.use(snap.snapname, 'cmugadget', periodic=True)
    bh = snap.readbh()
    context['bh'] = bh
    bad = context['bh']['bhmass'] > 2.0
    context['bh']['bhmass'][bad] = 2.0
    centers = numpy.atleast_2d(centers)
    print centers.shape
    print centers
    ctr = Field(numpoints=len(centers))
    ctr['pos'] = numpy.array(centers)
    print ctr['pos']
    context['centers'] = ctr
#    center = bh['pos'][((bh['pos'] - numpy.array(center)[None, :])**2).sum(axis=-1).argmin()].copy()
    #center = bh[0]['pos'].copy()
    if args.unfold:
        context.unfold(M)

    try:
        C = numpy.fromfile('%szooms-C-%03d.f8' % (prefix, snapid)).reshape(
                len(centers), 
                len(zooms),
                context.shape[0], context.shape[1])
        L = numpy.fromfile('%szooms-L-%03d.f8' % (prefix, snapid)).reshape(C.shape)
        print 'using existing files'
    except Exception as e:
        print e
        print 'making '
        C, L = makezoom(snap, context, centers, zooms, args.unfold)
        C.tofile('%szooms-C-%03d.f8' % (prefix, snapid))
        L.tofile('%szooms-L-%03d.f8' % (prefix, snapid))

    def work(j):
        for i, size in enumerate(zooms):
            figure = Figure(dpi=400)
            ax = context.newaxes(figure)
            vmax = nl_(L[:, i]).vmax - 0.2
            vmin = nl_(L[:, i]).vmin + 0.2
            context.view(center=context['centers']['pos'][j], size=size)
            context.imshow(nl_(C[j, i], vmin=3.5, vmax=7.5), nl_(L[j, i], vmax=vmax, vmin=vmin))
            X, Y, B = context.transform('bh', 'bhmass')
            a = (-B).argsort()[:5]
            context.scatter(X[a], Y[a], nl_(B[a], vmax=0, vmin=-4) * size * 0.2, fancy=True, marker='x', color='w')
            #context.scatter(X, Y, size * 0.2, fancy=True, marker='x', color='w')
            context.frame(axis=False, scale={'color':'w'})
            ax.text(0.1, .8, 'z=%0.2g' % context.C['redshift'], color='w',
                    transform=ax.transAxes) 
            context.savefig('zoom-%s-region-%02d-%03d-%04d.png' % (prefix, j, snapid,
                size), dpi=400)
    with sharedmem.Pool() as pool:
        pool.map(work, range(len(centers)))

if __name__ == '__main__':
    main()
