from gaepsi.gaplot import *
from matplotlib.figure import Figure
from matplotlib import cm
from sys import argv

print 'started'

prefix = argv[1] #'/physics/yfeng1/e5/run/000/output/'
snapid = int(argv[2])
use(prefix + '/snapdir_%03d/snapshot_%03d.%%d' %(snapid,
    snapid), 'cmugadget', periodic=True)
print context.fids()
fids = context.fids()
M = numpy.matrix("4,0,1; 1,2,-1; 2,1,0")
#gas = read('gas')
#print 'read', len(gas), 'particles'
#makeT('gas')
#print gas['T'].max(), gas['T'].min()
unfold(M)
view(center=context.boxsize / 2.0, size=context.boxsize)
print context.boxsize
context.reshape((2290, 1230))
C = numpy.zeros(context.shape, dtype='f8')
L = numpy.zeros(context.shape, dtype='f8')
with sharedmem.Pool() as pool:
    def work(fid):
        ct = GaplotContext(context.snapname,
                context.format, 
                periodic=True, np=0)
        ct.reshape(context.shape)
        ct.schema('gas', 0, ['mass', 'ie', 'ye', 'sml'])
        gas = ct.read('gas',fids=[fid])
        print 'read', len(gas), 'particles'
        ct.unfold(M);
        ct.makeT('gas')
        print gas['T'].max(), gas['T'].min()
        ct.view(center=ct.boxsize / 2.0, size=ct.boxsize)
        print 'paint'
        C, L = ct.paint('gas', 'T', 'mass', 'sml')
        print 'done'
        return C, L
    def reduce(C_, L_):
        C_[L_ == 0] = 0
        C[...] += C_ * L_
        L[...] += L_
    pool.map(work, fids, reduce=reduce)

C[...] /= L

C.tofile('C-%03d.f8' % snapid)
L.tofile('L-%03d.f8' % snapid)

figure = Figure(dpi=200)
newaxes(figure)
vmax = nl_(L).vmax - 0.4
print vmax
imshow(nl_(C, vmin=3.5, vmax=7.5), nl_(L, vmax=vmax, vmin='40db'))
frame(axis=False, scale={'color':'w'})
print_png('full-sim-%03d.png' % (snapid))



