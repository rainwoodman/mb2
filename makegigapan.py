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
#M = numpy.matrix("4,0,1; 1,2,-1; 2,1,0")
M = numpy.matrix("2,1,1  ;1,4,-2 ; 1,0,1") # 16x9 on gigapan, 9x16 othersize
print 'started'

def makegigapan(snap, context):
    fids = context.fids()

    def work(fid):
        ct = GaplotContext(context.snapname,
                context.format, 
                periodic=True)
        ct.schema('gas', 0, ['mass', 'ie', 'ye', 'sml'])
        ct.attach(context.CCD)
        gas = ct.read('gas',fids=[fid])
        print 'read', len(gas), 'particles'
        ct.unfold(M)
        ct.view(center=ct.origin + ct.boxsize * 0.5,
            size=ct.boxsize)
        ct.makeT('gas')
        ct.paint('gas', 'T', 'mass', 'sml', preserve=True)
    map(work, fids)
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('snap', type=SnapDir)
#    parser.add_argument('mode', choices=['make', 'continue'])
    parser.add_argument('imsize', type=float, help='size of image in Gpixels')
    args = parser.parse_args()
    snap = args.snap
    snapid = snap.snapid


    context = GaplotContext()
    context.use(snap.snapname, 'cmugadget', periodic=True)
    bh = snap.readbh()
    context['bh'] = bh
    bad = context['bh']['bhmass'] > 2.0
    context['bh']['bhmass'][bad] = 2.0
    context.unfold(M)
    context.view(center=context.origin + context.boxsize * 0.5,
            size=context.boxsize)

    fac = (args.imsize * 1e9/ numpy.prod(context.boxsize[:2])) ** 0.5
    shape =(
            int(context.boxsize[0] * fac),
            int(context.boxsize[1] * fac))
    print shape
    try:
        CCD = numpy.memmap('gigapan-CCD-%03d-%dx%d.f8' % (snapid, shape[0],
            shape[1]), dtype='f4', mode='r').reshape(
                shape[0], shape[1], 2)
        context.attach(CCD)
        print 'using existing files'
    except Exception as e:
        print e
        print 'making '
        CCD = numpy.memmap('gigapan-CCD-%03d-%dx%d.f8' % (snapid, shape[0],
            shape[1]), mode='w+', dtype='f4',
                shape=(shape[0], shape[1], 2))
        context.attach(CCD)
        makegigapan(snap, context)
        CCD.flush()
    C, L = CCD[..., 0].copy(), CCD[..., 1]
    C /= L
    figure = Figure(dpi=400)
    ax = context.newaxes(figure)
    vmax = nl_(L).vmax - 0.2
    vmin = nl_(L).vmin + 0.2
    print vmax
    image(nl_(C, vmin=3.5, vmax=7.5), nl_(L, vmax=vmax, vmin=vmin))
    im = numpy.memmap('gigapan-%03d-%dx%d.raw' % (snapid, shape[1], shape[0]), 
            mode='w+', dtype=('u1', 3), shape=L.shape)
    for i in range(len(L)):
        im[i] = image(nl_(C[i], vmin=3.0, vmax=8.0), nl_(L[i],
            vmax=vmax, vmin=vmin))[..., :3]

    X, Y, I = context.transform('bh')
    bh = context['bh']
    B = bh['bhmass'][I]
    a = B.argsort()
    n = nl_(bh['bhmass'])
    print X.min(), X.max(), Y.min(), Y.max()
    addspikes(im, X[a], Y[a], nl_(B[a], vmax=n.vmax, vmin=n.vmin) * 10, color='w')
    im.flush()
if __name__ == '__main__':
    main()
