import readsubhalo
from sys import argv
import sharedmem
import numpy
def process(snapid):
    snapdir = readsubhalo.SnapDir(snapid, '../')

    sfr = snapdir.load(0, 'sfr')
    chunksize = 64 * 1024
    def work(i):
        return sfr[i:i+chunksize].sum(dtype='f8')
    with sharedmem.MapReduce() as pool:
        sfrsum = numpy.sum(pool.map(work, range(0, len(sfr), chunksize)))

    bhmdot = snapdir.load(5, 'bhmdot').copy()
    # fix the ugly things
    bhmdot[bhmdot > 1e3] = 0
    print snapid, snapdir.redshift, sfrsum, bhmdot.sum(dtype='f8'), len(bhmdot)

for snapid in argv[1:]:
    process(snapid)
