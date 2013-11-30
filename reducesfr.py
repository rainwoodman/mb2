from gaepsi.snapshot import Snapshot
from sys import argv
print len(argv[1:])
sfr = 0
for fn in argv[1:]:
    snap = Snapshot(fn, 'cmugadget')
    assert snap.C['Nfiles'] == len(argv[1:])
    sfr += snap[0, 'sfr'].sum(dtype='f8')
    print sfr
