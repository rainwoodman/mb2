from gaepsi.snapshot import Snapshot
from gaepsi.tools.analyze import HaloCatalog
import numpy
import time
import sharedmem
import os.path
from sys import argv
first = Snapshot(argv[1] % 0, 'cmugadget', idtype='u8')
class Timer:
    def __enter__(self):
        self.start = time.clock()
        return self
    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start

Ntot = first.C['Ntot']
Nout = first.C['Nfiles']
for i in range(Nout):
    output = Snapshot(argv[2] % i, 'cmugadget', template=first,
            create=True)
    output.C['N'] = Ntot * (i + 1) // Nout - Ntot * i // Nout
    output.save()

for ptype in range(6):
    buffer = {}
    j = 0
    for blk in first.schema:
        buffer[blk] = numpy.empty(0, first.schema[blk].dtype)

    for i in range(Nout):
        output = Snapshot(argv[2] % i, 'cmugadget', template=first,
                readonly=False)
        while output.C['N'][ptype] > len(buffer['pos']):
            print 'need to read ', j, output.C['N'][ptype], len(buffer['pos'])
            newinput = Snapshot(argv[1] % j, 'cmugadget', template=first)
            j = j + 1
            for blk in newinput.schema:
                if not (ptype, blk) in newinput: continue
                buffer[blk] = numpy.append(buffer[blk], newinput[ptype, blk],
                        axis=0)
        for blk in first.schema:
            if not (ptype, blk) in output: continue
            output[ptype, blk][:] = buffer[blk][:output.C['N'][ptype]]
            buffer[blk] = buffer[blk][output.C['N'][ptype]:].copy()
        output.save()





    
