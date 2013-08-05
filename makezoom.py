from gaepsi.gaplot import *
from matplotlib.figure import Figure
from matplotlib import cm
from sys import argv
center = [ 62995.8671875,  40180.12109375, 52267.77734375]

zoom = [100000, 50000, 10000, 5000, 3000, 2500, 1250, 1000, 600, 500, 300, 200,
        100]

print 'started'

prefix = argv[1] #'/physics/yfeng1/e5/run/000/output/'
snapid = int(argv[2])
use(prefix + '/snapdir_%03d/snapshot_%03d.%%d' %(snapid,
    snapid), 'cmugadget', periodic=True)
#bh = numpy.load('../bhcorr/bh_%03d.npz' % snapid)
bh = read('bh')
newbh = Field(numpoints=len(bh))

for component in ['pos', 'bhmass', 'bhmdot']:
    newbh[component] = bh[component].copy()

context['bh'] = newbh
center = bh[0]['pos'].copy()
print context.fids()
fids = context.fids()
M = numpy.matrix("4,0,1; 1,2,-1; 2,1,0")
#gas = read('gas')
#print 'read', len(gas), 'particles'
#makeT('gas')
#print gas['T'].max(), gas['T'].min()
unfold(M, center=center)
print context.boxsize
context.reshape((2000, 2000))
C = sharedmem.empty((len(zoom), context.shape[0], context.shape[1]), dtype='f8')
L = sharedmem.empty((len(zoom), context.shape[0], context.shape[1]), dtype='f8')
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
        ct.unfold(M, center=center);
        ct.makeT('gas')
        print gas['T'].max(), gas['T'].min()
        for i, size in enumerate(zoom):
            ct.view(center=ct.origin + ct.boxsize * 0.5, size=size)
            C_, L_ = ct.paint('gas', 'T', 'mass', 'sml')
            C_[L_ == 0] = 0
            with pool.critical:
                C[i] += C_ * L_
                L[i] += L_
    pool.map(work, fids)

C[...] /= L

C.tofile('zooms-C-%03d.f8' % snapid)
L.tofile('zooms-L-%03d.f8' % snapid)

for i, size in enumerate(zoom):
    view(center=context.origin + context.boxsize * 0.5, size=zoom[i])
    figure = Figure(dpi=200)
    newaxes(figure)
    vmax = nl_(L[i]).vmax - 0.4
    print vmax
    imshow(nl_(C[i], vmin=3.5, vmax=7.5), nl_(L[i], vmax=vmax, vmin='40db'))
    X, Y, B = transform('bh', 'bhmass')
    scatter(X, Y, nl_(B) * size * 0.2, fancy=True, marker='x', color='w')
    frame(axis=False, scale={'color':'w'})
    print_png('zoom-sim-%03d-%04d.png' % (snapid, size))
