from gaepsi.snapshot import Snapshot
from gaepsi.tools.analyze import HaloCatalog
import numpy
import time
import sharedmem
import os.path
from sys import argv
firstcat = Snapshot(argv[1] % 0, 'cmugadget.GroupTab', idtype='u8', floattype='f4')
Ncat = (firstcat[0, 'length'] > 8500000).nonzero()[0].max()
class Timer:
    def __enter__(self):
        self.start = time.clock()
        return self
    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start

TheBigGuys = range(Ncat)

def outputname(haloid, j):
    outputname = '%03d/' % haloid + os.path.basename(argv[2])
    return outputname % j

def inputname(j):
    return argv[2] % j

for i in TheBigGuys:
    try:
        os.makedirs(os.path.dirname(outputname(i, 0)))
        print 'made dir', os.path.dirname(outputname(i, 0))
    except Exception as e:
        print e
        pass

print Ncat
cat = HaloCatalog(argv[1], 'cmugadget', idtype='u8', floattype='f4', count=Ncat)
print 'halo catalog loaded', Ncat
first = Snapshot(argv[2] % 0, 'cmugadget')

Nfiles = first.C['Nfiles']


def work(j):
    input = Snapshot(inputname(j), 'cmugadget', template=first)
    id = input[None, 'id']
    t1 = 0
    t2 = 0
    t3 = 0
    for i in TheBigGuys:
        masks = []
        output = Snapshot(outputname(i, j), 'cmugadget', create=True, template=first)
        for ptype in range(6):
          with Timer() as timer:
            id = input[ptype, 'id']
            mask = cat.mask(id, groupid=i)
            masks.append(mask) 
            output.C['N'][ptype] = mask.sum()
          t1 += timer.interval
        for blk in input.schema:
          if not (None, blk) in input: continue
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
    return j, t1, t2, t3

def reduce(j, t1, t2, t3):
    print j, Nfiles, 'times', t1, t2, t3

with sharedmem.Pool() as pool:
    pool.map(work, range(Nfiles),
            reduce=reduce)

for i in TheBigGuys:
    Ntot = numpy.zeros(6, 'i8')
    for j in range(Nfiles):
        output = Snapshot(outputname(i, j), 'cmugadget', template=first)
        Ntot += output.C['N']
    print i, Ntot
    for j in range(Nfiles):
        output = Snapshot(outputname(i, j), 'cmugadget', template=first)
        output.C['Ntot'][:] = Ntot
        output.save('header')

