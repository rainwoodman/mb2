from sys import argv
from gaepsi.snapshot import Snapshot
import numpy
import os.path
import sharedmem
idsfile = argv[1]
sfgasdir = argv[2]

first = Snapshot(idsfile % 0, 'cmugadget.GroupIDs', idtype='u8',
        floattype='f4')

sfgasid = numpy.fromfile(sfgasdir + '/sfgasid', dtype='u8')
sfgassfr = numpy.fromfile(sfgasdir + '/sfgassfr', dtype='f4')
matchedfile = file(sfgasdir + '/matchedsfr', 'w')
matchedidfile = file(sfgasdir + '/matchedid', 'w')
with sharedmem.Pool(use_threads=False) as pool:
    def work(i):
       snap = Snapshot(idsfile % i, 'cmugadget.GroupIDs', 
               idtype='u8', floattype='f4')
   
       id = snap[0, 'id']
       result = numpy.zeros(len(id), dtype='f4')
       lookedid = numpy.zeros(len(id), dtype='u8')
       ind = sfgasid.searchsorted(id)
       ind.clip(0, len(sfgasid) - 1, out=ind)
   
       found = sfgasid[ind] == id
       result[found] = sfgassfr[ind[found]]
       with pool.ordered:
           print 'start', i, matchedidfile.tell() / 8., len(id), id[0]
           matchedfile.seek(0, 2)
           result.tofile(matchedfile)
           matchedfile.flush()
           matchedfile.seek(0, 2)
           matchedidfile.seek(0, 2)
           id.tofile(matchedidfile)
           matchedidfile.flush()
           matchedidfile.seek(0, 2)
           print i, 'found', found.sum(), matchedidfile.tell() / 8.
       
       return found.sum()

    total_found = numpy.sum(pool.map(work, range(first.C['Nfiles'])))

print total_found, len(sfgasid)
assert total_found <= len(sfgasid)
