from gaepsi.snapshot import Snapshot
from gaepsi.tools.analyze import HaloCatalog
import numpy
import time
import sharedmem
import os.path
from sys import argv
firstcat = Snapshot(argv[1] % 0, 'cmugadget.GroupTab', idtype='u8', floattype='f4')
Ncat = (firstcat[0, 'length'] > 8500000).nonzero()[0].max()
print Ncat
cat = HaloCatalog(argv[1], 'cmugadget', idtype='u8', floattype='f4', count=Ncat)
print 'halo catalog loaded', Ncat
first = Snapshot(argv[2] % 0, 'cmugadget')

BoxSize = first.C['boxsize']
class Timer:
    def __enter__(self):
        self.start = time.clock()
        return self
    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start

def work(inputname):
    input = Snapshot(inputname, 'cmugadget', template=first)
    id = input[None, 'id']
    t1 = 0
    t2 = 0
    t3 = 0
    #fix this!###
    for i in [Ncat]: # range(Ncat):
        masks = []
        outputname = 'halo-%03d-' % i + os.path.basename(inputname)
        output = Snapshot(outputname , 'cmugadget', create=True, template=first)
        for ptype in range(6):
          with Timer() as timer:
            id = input[ptype, 'id']
            mask = cat.mask(id, groupid=i)
            masks.append(mask) 
            output.C['N'][ptype] = mask.sum()
          t1 += timer.interval
        for blk in input.schema:
          with Timer() as timer:
            junk = input[None, blk]
            for ptype in range(6):
                if not (ptype, blk) in input: continue
                if output.C['N'][ptype] == 0: continue
                output[ptype, blk][:] = input[ptype, blk][masks[ptype]]
          t2 += timer.interval
        with Timer() as timer:
            output.save()
        t3 += timer.interval
    return inputname, t1, t2, t3
def reduce(inputname, t1, t2, t3):
    print inputname, 'times', t1, t2, t3
#with sharedmem.Pool() as pool:
#    pool.map(work, [argv[2] % i for i in
#        numpy.random.permutation(range(first.C['Nfiles']))],
#            reduce=reduce)

for i in [45]: #range(Ncat):
    outputname = 'halo-%03d-' % i + os.path.basename(argv[2])
    Ntot = numpy.zeros(6, 'i8')
    for j in range(first.C['Nfiles']):
        output = Snapshot(outputname % j, 'cmugadget', template=first)
        Ntot += output.C['N']
    print i, Ntot
    for j in range(first.C['Nfiles']):
        output = Snapshot(outputname % j, 'cmugadget', template=first)
        output.C['Ntot'][:] = Ntot
        output.save('header')

