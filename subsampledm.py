from gaepsi.snapshot import Snapshot
import sharedmem
import numpy
from sys import argv
from sys import stdout
import os.path

first = Snapshot(argv[1] % 0, 'cmugadget')
seedtable = numpy.random.randint(1<<21, size=first.C['Nfiles'])

with sharedmem.Pool() as pool:
    def work(i):
        snap = Snapshot(argv[1] % i, 'cmugadget', template=first)
        numpy.random.seed(seedtable[i])
#        print 'in', i
        with pool.critical:
            pos = snap[1, 'pos']
        r = numpy.random.uniform(size=len(pos))
        mask = r < 100.**3 / first.C['Ntot'][1]
        pos = pos[mask]
#        print 'done', i, mask.sum(), sfr.max(), sfr.min()
        with pool.critical:
            numpy.savetxt(stdout, pos, fmt='%0.7g')
            stdout.flush()

    pool.map(work, range(first.C['Nfiles']))

