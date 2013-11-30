import numpy
import os.path
from readsubhalo import *
from readsubhalo import SnapDir as SnapDirReadOnly

def startend(tab, type):
    """ returns the star end of type from a tab's lenbytype"""
    end = tab['lenbytype'][:, type].cumsum()
    start = end.copy()
    start[0] = 0
    start[1:] = end[:-1]
    return start, end


from gaepsi.snapshot import Snapshot
from gaepsi.field import Field
from gaepsi.cosmology import Cosmology

class SnapDir(SnapDirReadOnly):
    def __init__(self, snapid, ROOT='../'):
        SnapDirReadOnly.__init__(self, snapid, ROOT)

        snapid = self.snapid
        self.snapname = ROOT + '/snapdir/snapdir_%03d/snapshot_%03d.%%d' %(snapid, snapid)
        first = Snapshot(self.snapname % 0, 'cmugadget')
        self.first = first
        self.C = first.C
        C = first.C
        self.cosmology = Cosmology(h=C['h'], M=C['OmegaM'], L=C['OmegaL'])
        self.U = self.cosmology.U
        self.schema = first.schema
        self.collectivedir = ROOT + \
        '/isolatedir/isolate_%03d/%%03d' % snapid

        self.tabfile = ROOT + \
        '/groups_and_subgroups/groups_%03d/subhalo_tab_%03d.%%d' %(snapid, snapid)
        if not os.path.isfile(self.tabfile % 0):
            self.tabfile = ROOT + \
        '/groups_and_subgroups/no_collective_halos/groups_%03d/subhalo_tab_%03d.%%d' \
             %(snapid, snapid)

        self.grouptabfile = self.tabfile.replace('subhalo', 'group')
    def readbh(self):
        rawbh = numpy.load(ROOT + '/bhcorr/bh_%03d.npz' % self.snapid)
        bh = Field(numpoints=len(rawbh['bhmass']))
        assert self.C['Ntot'][5] == len(bh)
        for component in rawbh.files:
            bh[component] = rawbh[component].copy()
        return bh
    def opensnap(self, fid):
        return Snapshot(self.snapname %fid, 'cmugadget', template=self.first)

    def collectivetabfile(self, groupid):
        return self.collectivedir % groupid \
            + '/groups_%03d/subhalo_tab_%03d.%%d' % (self.snapid, self.snapid)
    def hascollective(self, groupid, group):
        """ split is at 8.5e6 """
        if self.snapid <= 73: return False
        if group['len'] > 8.8e6:
            return os.path.isfile(self.collectivetabfile(groupid) % 0)
        return False

class Depot:
    def __init__(self, dtype, nfile):
        self.ifile = 0
        self.nfile = int(nfile)
        self.dtype = numpy.dtype(dtype)
        self.pool = numpy.empty(0, dtype)
        self.total_read = 0
    def __len__(self): return self.nfile
    @property
    def size(self): return self.nfile

    def fetch(self, i):
        """ override to read in i-th shipment """
        raise UnimplementedError("method is unimplemented")

    def consume(self, n=None):
        """return n items, if n is None return all items"""
        while n is None or len(self.pool) < n:
            try: 
                self._supply()
            except StopIteration:
                break
        if n is not None:
            if len(self.pool) < n:
                raise ValueError("insufficient supply limit reached")
            rt = self.pool[:n].copy()
            self.pool = self.pool[n:].copy()
        else:
            rt = self.pool
            self.pool = numpy.empty(0, self.dtype)
        self.total_read = self.total_read + len(rt)
        return rt

    def _supply(self):
        if self.ifile == self.nfile:
            raise StopIteration
        self.pool = numpy.append(self.pool, 
                self.fetch(self.ifile), axis=0)
        self.ifile += 1

def wrong_file_or_die(filename, expectedsize):
    if os.path.isfile(filename):
        size = os.path.getsize(filename)
        if size == expectedsize:
            raise Exception("File seems to be right. quit!")
        else:
            print 'will overwrite', filename

