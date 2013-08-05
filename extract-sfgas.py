from gaepsi.snapshot import Snapshot
import sharedmem
import numpy
from sys import argv
import os.path

first = Snapshot(argv[1] % 0, 'cmugadget')

def work(i):
    snap = Snapshot(argv[1] % i, 'cmugadget', template=first)

    print 'in', i
    sfr = snap[0, 'sfr']
    id = snap[0, 'id']
    mask = sfr > 0
    print 'done', i, mask.sum(), sfr.max(), sfr.min()
    return mask.sum()

with sharedmem.Pool() as pool:
    size = pool.map(work, range(first.C['Nfiles']), ordered=True)

print 'phase 1 done'
size = numpy.array(size, 'i8')
total = size.sum()
offset = numpy.concatenate(([0], size.cumsum()))

print 'allocation total pars', total
idall = sharedmem.empty(total, dtype=first.schema['id'].dtype)
sfrall = sharedmem.empty(total, dtype='f4')

def work(i):
    snap = Snapshot(argv[1] % i, 'cmugadget', template=first)
    print 'i', i
    sfr = snap[0, 'sfr']
    id = snap[0, 'id']
    print 'i', i, 'read'
    mask = sfr > 0
    idall[offset[i]:offset[i]+size[i]] = id[mask]
    sfrall[offset[i]:offset[i]+size[i]] = sfr[mask]
    print 'i', i, 'assign'
    print 'i', i, mask.sum()
with sharedmem.Pool() as pool:
    pool.map(work, range(first.C['Nfiles']))
print 'phase 2 done'

arg = sharedmem.argsort(idall)
idall = idall[arg]
sfrall = sfrall[arg]

idall.tofile('sfgasid')
sfrall.tofile('sfgassfr')

