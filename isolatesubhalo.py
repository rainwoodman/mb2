from gaepsi.snapshot import Snapshot
from gaepsi.cosmology import Cosmology
from lib import * 

import numpy
import argparse
from sys import argv
import os.path
ROOT = '/physics/yfeng1/mb2'
class SubHaloParDepot(Depot):
    def __init__(self, tabfile):
        self.tabfile = tabfile
        self.posvelfile = self.tabfile.replace('_tab_', '_posvel_')
        self.idsfile = self.tabfile.replace('_tab_', '_ids_')
        first = Snapshot(tabfile %0, 'cmugadget.SubHaloTab', idtype='u8',
            floattype='f4')
        Depot.__init__(self, pdtype, first.C['Nfiles'])

    def fetch(self, id):
        f = file(self.posvelfile % id , mode='r')
        header = numpy.fromfile(f, dtype='i4', count=9)
        size = header[2]
        rt = numpy.empty(shape=size, dtype=pdtype)
    
        rt['pos'] = numpy.fromfile(f, dtype=('f4', 3), count=size)
        rt['vel'] = numpy.fromfile(f, dtype=('f4', 3), count=size)
        rt['type'] = numpy.fromfile(f, dtype='u1', count=size)
        rt['mass'] = numpy.fromfile(f, dtype='f4', count=size)
        # do not use it
        age = numpy.fromfile(f, dtype='f4', count=size)
        rt['id'] = numpy.memmap(self.idsfile % id, mode='r', offset=28, dtype='u8')
        return rt
    
class SubHaloDepot(Depot):
    def __init__(self, tabfile):
        self.tabfile = tabfile
        first = Snapshot(tabfile %0, 'cmugadget.SubHaloTab', idtype='u8',
            floattype='f4')
        self.Nhalo = int(first.C['NtotSubgroups'])
        Depot.__init__(self, subdtype, first.C['Nfiles'])

    def fetch(self, id):
        f = Snapshot(self.tabfile % id , 'cmugadget.SubHaloTab', 
                 idtype='u8', floattype='f4')
        rt = numpy.empty(f.C['N'][1], dtype=subdtype)
        rt['pos'] = f[1, 'halopos']
        rt['vel'] = f[1, 'halovel']
        rt['len'] = f[1, 'halolen']
        rt['mass'] = f[1, 'halomass']
        rt['vdisp'] = f[1, 'haloveldisp']
        rt['vcirc'] = f[1, 'halovmax']
        rt['rcirc'] = f[1, 'halovmaxrad']
        rt['grnr'] = f[1,'haloparent']
        rt['lenbytype'] = 0
        rt['massbytype'] = 0
        return rt

class GroupDepot(Depot):
    def __init__(self, tabfile):
        self.tabfile = tabfile
        first = Snapshot(tabfile %0, 'cmugadget.GroupTab', idtype='u8',
            floattype='f4')
        Depot.__init__(self, groupdtype, first.C['Nfiles'])
        self.Ngroups = first.C['Ntot'][0]
    def fetch(self, id):
        f = Snapshot(self.tabfile % id , 'cmugadget.GroupTab', 
                 idtype='u8', floattype='f4')
        rt = numpy.empty(f.C['N'][0], dtype=groupdtype)
        rt['pos'] = f[0, 'pos']
        rt['vel'] = f[0, 'vel']
        rt['mass'] = f[0, 'mass']
        rt['len'] = f[0, 'length']
        rt['massbytype'] = f[0, 'massbytype']
        rt['lenbytype'] = f[0, 'lenbytype']
        rt['nhalo'] = 0
        return rt

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="commands")
    split_parser = subparsers.add_parser("split", help="split particles to six files")
    split_parser.add_argument("snap", type=SnapDir)
    split_parser.set_defaults(func=splitMain)
    lookup_parser = subparsers.add_parser("lookup", help="lookup properties")
    lookup_parser.add_argument("snap", type=SnapDir)
    lookup_parser.add_argument("type", type=int)
    lookup_parser.add_argument("field")
    lookup_parser.set_defaults(func=lookupMain)
    args = parser.parse_args()

    args.func(args)

def splitMain(args):
    snap = args.snap

    for type in range(6):
        try:
            os.makedirs(snap.subhalodir + '/%d' % type)
        except OSError:
            pass
    
    file(snap.subhalodir + '/header.txt', mode='w').write(str(snap.C))
    depot_par = SubHaloParDepot(snap.tabfile)
    depot_sub = SubHaloDepot(snap.tabfile)
    depot_group = GroupDepot(snap.grouptabfile)
    
    if os.path.isfile(snap.groupfile):
        size = os.path.getsize(snap.groupfile)
        N = depot_group.Ngroups
        if size == groupdtype.itemsize  * N:
            raise Exception("File seems to be right. quit!")
        else:
            print 'ovewrite', size, groupdtype.itemsize, N

    F = {}
    for type in range(6):
        F[type] = {}
        for field in pdtype.fields:
            F[type][field] = snap.open(type, field, mode='w')
        

    halodump = file(snap.subhalofile, 'w')
    groupdump = file(snap.groupfile, 'w')
    
    realgroupid = 0
    for i in range(len(depot_sub)):
        tab = Snapshot(snap.tabfile % i, 'cmugadget.SubHaloTab', idtype='u8',
                floattype='f4')
        group = depot_group.consume(tab.C['N'][0])
        for groupid in range(tab.C['N'][0]):
            if snap.hascollective(realgroupid, group[groupid]):
                fakenhalo = tab[0, 'nhalo'][groupid]
                depot_subcol = SubHaloDepot(snap.collectivetabfile(realgroupid))
                depot_parcol = SubHaloParDepot(snap.collectivetabfile(realgroupid))
                pars = depot_parcol.consume()
                realhalo = depot_subcol.consume()
                nhalo = len(realhalo)
                # skip the pars from ids file (it's of wrong order)
                junk = depot_par.consume(group['len'][groupid])
                assert len(pars) == group['len'][groupid]
            else:
                nhalo = tab[0, 'nhalo'][groupid]
                pars = depot_par.consume(group['len'][groupid])
                realhalo = depot_sub.consume(nhalo)

            halo = numpy.empty(nhalo + 1, dtype=realhalo.dtype)

            halo[:-1] = realhalo
            # last halo is contamination
            halo[-1]['mass'] = numpy.nan
            halo[-1]['pos'] = numpy.nan
            halo[-1]['grnr'] = halo[0]['grnr']
            halo[-1]['len'] = group['len'][groupid] - halo['len'][:-1].sum()
            assert (halo['len'] >= 0).all()
            parstart = numpy.concatenate(([0], halo['len'].cumsum()))
            parend = parstart[1:]
            #make sure group and tab are of the same file
            assert group[groupid]['len'] == tab[0, 'length'][groupid]
    
            for haloid in range(nhalo + 1):
                parinhalo = pars[parstart[haloid]:parend[haloid]]
                parinhalo.sort(order='type')
                start = parinhalo['type'].searchsorted(range(6), side='left')
                end = parinhalo['type'].searchsorted(range(6), side='right')
                for type in range(6):
                    thistype = parinhalo[start[type]:end[type]]
                    halo[haloid]['massbytype'][type] = thistype['mass'].sum(dtype='f8')
                    halo[haloid]['lenbytype'][type] = len(thistype)
                    for field in pdtype.fields:
                        thistype[field].tofile(F[type][field])
            group[groupid]['nhalo'] = nhalo
            print i, 'group', groupid, 'nhalo', \
                    group['nhalo'][groupid], \
                    'npar', group['len'][groupid]
            halo.tofile(halodump)        
            halodump.flush()
            if not numpy.allclose(group[groupid]['massbytype'],
                    halo['massbytype'].sum(axis=0)):
                print group[groupid]['massbytype'], halo['massbytype'].sum(axis=0)
            assert numpy.all(group[groupid]['lenbytype'] ==
                   halo['lenbytype'].sum(axis=0))
            realgroupid += 1
        group.tofile(groupdump)
        groupdump.flush()
        for type in range(6):
            for field in pdtype.fields:
                F[type][field].flush()

def lookupMain(args):
    type = args.type
    field = args.field
    snap = args.snap

    Nfiles = snap.C['Nfiles']
    id = snap.load(type, 'id')
    dtype = snap.schema[field].dtype

    if os.path.isfile(snap.filename(type, field)):
        size = os.path.getsize(snap.filename(type, field))
        N = len(id)
        if size == dtype.itemsize  * N:
            raise Exception("File seems to be right. quit!")
        else:
            print 'ovewrite', size, dtype.itemsize, N

    arg = id.argsort()
    idsorted = id[arg]

    val = numpy.empty(len(id), dtype)
#    outfile = numpy.memmap(outputdir + '/%d/%s.raw' %(type, field), mode='w',
#            dtype=first.schema[field].dtype, shape=len(id))

    Nfound = 0
    for fid in range(Nfiles):
        f = snap.opensnap(fid)
        idlookup = f[type, 'id']
        valuelookup = f[type, field]
        ind = idsorted.searchsorted(idlookup)
        found = idsorted.take(ind, mode='wrap') == idlookup
        assert (id[arg[ind[found]]] == idlookup[found]).all()
        val.put(arg[ind[found]], valuelookup[found])
        Nfound += found.sum()
        print fid, found.sum()
    assert Nfound == len(id)
    del idsorted

    val.tofile(snap.open(type, field, 'w'))

if __name__ == '__main__':
    main()
