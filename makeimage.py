from gaepsi.gaplot import *

print 'started'
use('../../snapdir_085/snapshot_085.%d', 'cmugadget', periodic=True)
context.reshape((220000, 123600))
#context.reshape((220000/10, 123600/10))
schema('gas', 0, ['mass', 'ie', 'ye', 'sml'])
gas = read('gas', np=4)
print 'read', len(gas), 'particles'
makeT('gas')
del gas['ie']
del gas['ye']
print 'makeT'

M = numpy.matrix("7 6 6 ; 6 -3 2;  -1 4 1")
print 'unfold'
unfold(M);
view(center=context.boxsize / 2.0, size=context.boxsize)

print 'paint'
C, L = paint('gas', 'T', 'mass', 'sml', dtype='f4')

print 'saving'
C.tofile('C0-085.f4')
L.tofile('L0-085.f4')
print 'done'
