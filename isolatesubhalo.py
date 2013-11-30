from gaepsi.snapshot import Snapshot
from gaepsi.cosmology import Cosmology
from lib import * 

import numpy
import argparse
from sys import argv
import os.path
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
        # i don't understand this; we do not want halo parent
        # it's always zero
        rt['parent'] = f[1,'haloparent']
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

    sp = subparsers.add_parser("split", help="split particles to six files")
    sp.add_argument("snap", type=SnapDir)
    sp.add_argument("--check", action='store_true', default=False)
    sp.set_defaults(func=splitMain)

    sp = subparsers.add_parser("lookup", help="lookup properties")
    sp.add_argument("snap", type=SnapDir)
    sp.add_argument("type", type=int)
    sp.add_argument("field")
    sp.set_defaults(func=lookupMain)

    sp = subparsers.add_parser("fixgroupid", help="fix group id")
    sp.add_argument("snap", type=SnapDir)
    sp.set_defaults(func=fixgroupidMain)

    sp = subparsers.add_parser("subhalotype", help="add type property to subhalo")
    sp.add_argument("snap", type=SnapDir)
    sp.set_defaults(func=subhalotypeMain)

    args = parser.parse_args()

    args.func(args)

def fixgroupidMain(args):
    """ no need for newer isolatehalo """
    snap = args.snap
    g = snap.readgroup()
    sub = numpy.memmap(snap.subhalofile, mode='r+', dtype=subdtype)
    sub['groupid'] = numpy.repeat(numpy.arange(len(g)), g['nhalo'] + 1)
    sub.flush()
    print 'fixed', snap.snapid

def subhalotypeMain(args):
    snap = args.snap
    g = snap.readsubhalo()
    cmask = numpy.memmap(snap.filename('subhalo', 'type'), 
            mode='w+', dtype='u1', shape=len(g))
    cmask[:] = numpy.isnan(g['mass'])
    cmask.flush()

def splitMain(args):
    snap = args.snap

    for type in range(6):
        try:
            os.makedirs(snap.subhalodir + '/%d' % type)
        except OSError:
            pass
    
    file(snap.subhalodir + '/header.txt', mode='w').write(str(snap.C))
    print snap.tabfile
    depot_sub = SubHaloDepot(snap.tabfile)
    depot_group = GroupDepot(snap.grouptabfile)
    
    if not args.check:
        wrong_file_or_die(snap.groupfile, depot_group.Ngroups *
                groupdtype.itemsize)

        depot_par = SubHaloParDepot(snap.tabfile)
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
                realhalo = depot_subcol.consume()
                nhalo = len(realhalo)
                # skip the pars from ids file (it's of wrong order)
                if not args.check:
                    depot_parcol = SubHaloParDepot(snap.collectivetabfile(realgroupid))
                    pars = depot_parcol.consume()
                    junk = depot_par.consume(group['len'][groupid])
                    assert len(pars) == group['len'][groupid]
            else:
                nhalo = tab[0, 'nhalo'][groupid]
                if not args.check:
                    pars = depot_par.consume(group['len'][groupid])
                realhalo = depot_sub.consume(nhalo)

            halo = numpy.empty(nhalo + 1, dtype=realhalo.dtype)

            halo[:-1] = realhalo
            # last halo is contamination
            halo[-1]['mass'] = numpy.nan
            halo[-1]['pos'] = numpy.nan
            halo[:]['groupid'] = realgroupid
            halo[-1]['len'] = group['len'][groupid] - halo['len'][:-1].sum()
            assert (halo['len'] >= 0).all()
            assert (halo['len'][0:-2] >= halo['len'][1:-1]).all()
            parstart = numpy.concatenate(([0], halo['len'].cumsum()))
            parend = parstart[1:]
            #make sure group and tab are of the same file
            assert group[groupid]['len'] == tab[0, 'length'][groupid]
    
            if not args.check:
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
                    'npar', group['len'][groupid], halo[-1]['len'], \
                    'collective', snap.hascollective(realgroupid, group[groupid])

            if not args.check:
                halo.tofile(halodump)        
                halodump.flush()
                assert numpy.all(group[groupid]['lenbytype'] ==
                       halo['lenbytype'].sum(axis=0))
                if not numpy.allclose(group[groupid]['massbytype'],
                        halo['massbytype'].sum(axis=0)):
                    print group[groupid]['massbytype'], halo['massbytype'].sum(axis=0)
            realgroupid += 1

        if not args.check:
            group.tofile(groupdump)
            groupdump.flush()
            for type in range(6):
                for field in pdtype.fields:
                    F[type][field].flush()
    assert depot_sub.Nhalo == depot_sub.total_read

def lookupMain(args):
    type = args.type
    field = args.field
    snap = args.snap

    Nfiles = snap.C['Nfiles']
    id = snap.load(type, 'id')
    dtype = snap.schema[field].dtype

    wrong_file_or_die(snap.filename(type, field),
            len(id) * dtype.itemsize)

    arg = id.argsort()
    idsorted = id[arg]

    val = numpy.empty(len(id), dtype)
    if len(id) > 0:
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
