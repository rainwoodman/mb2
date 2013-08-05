from gaepsi.snapshot import Snapshot
from gaepsi.cosmology import Cosmology
import numpy
import argparse
from sys import argv
import os.path
ROOT = '/physics/yfeng1/mb2'


def makeparser():
    parser = argparse.ArgumentParser()
    class action(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, self.dest, values)
            snapid = values
            namespace.snapname = ROOT + '/snapdir/snapdir_%03d/snapshot_%03d.%%d' %(snapid, snapid)
            namespace.outputdir = ROOT + '/subhalos/%03d' % snapid
            first = Snapshot(namespace.snapname % 0, 'cmugadget')
            namespace.first = first
            namespace.C = first.C
            C = first.C
            namespace.cosmology = Cosmology(h=C['h'], M=C['OmegaM'], L=C['OmegaL'])
            namespace.U = namespace.cosmology.U
            namespace.schema = first.schema
            namespace.open = lambda type, comp, schema=first.schema: \
                    numpy.fromfile(namespace.outputdir + '/%d/%s.raw' % (type, comp),
                            dtype=schema[comp].dtype)

    parser.add_argument('snapid', action=action, type=int)
    return parser

def main():
    parser = makeparser()
    subparsers = parser.add_subparsers(help="commands")
    split_parser = subparsers.add_parser("split", help="split particles to six files")
    split_parser.set_defaults(func=splitMain)
    lookup_parser = subparsers.add_parser("lookup", help="lookup properties")
    lookup_parser.add_argument("type", type=int)
    lookup_parser.add_argument("field")
    lookup_parser.set_defaults(func=lookupMain)
    args = parser.parse_args()

    args.func(args)

subdtype = numpy.dtype([
    ('mass', 'f4'), 
    ('len', 'i4'), 
    ('pos', ('f4', 3)), 
    ('vel', ('f4', 3)),   
    ('vdisp', 'f4'),  
    ('vcirc', 'f4'),  
    ('rcirc', 'f4'),  
    ('grnr', 'i4'),  
    ('massbytype', ('f4', 6)),  
    ('lenbytype', ('u4',6)), 
    ('unused', 'f4'),    
    ('unused2', 'u4'), 
   ])
groupdtype = numpy.dtype([
    ('mass', 'f4'), 
    ('len', 'i4'), 
    ('pos', ('f4', 3)), 
    ('vel', ('f4', 3)), 
    ('nhalo', 'i4'),   
    ('massbytype', ('f4', 6)),  
    ('lenbytype', ('u4',6)), 
   ])

pdtype = numpy.dtype([('pos', ('f4', 3)),
    ('vel', ('f4', 3)),
    ('mass', 'f4'),
    ('id', 'u8'),
    ('type', 'u1')
    ])

class Depot:
    def __init__(self, dtype, nshipment, fetch=lambda x: []):
        self.next_shipment = 0
        self.max_shipment = nshipment
        self.pool = numpy.empty(0, dtype)
        self.fetch = fetch

    def consume(self, n):
        while len(self.pool) < n:
            try: 
                self.supply()
            except StopIteration:
                break
        if len(self.pool) < n:
            raise ValueError("insufficient supply limit reached")
        rt = self.pool[:n].copy()
        self.pool = self.pool[n:].copy()
        return rt

    def supply(self):
        if self.next_shipment == self.max_shipment:
            raise StopIteration
        self.pool = numpy.append(self.pool, 
                self.fetch(self.next_shipment), axis=0)
        self.next_shipment += 1

def splitMain(args):
    snapid = args.snapid
    outputdir = args.outputdir

    tabfile = ROOT + \
        '/groups_and_subgroups/groups_%03d/subhalo_tab_%03d.%%d' %(snapid, snapid)
    posvelfile = tabfile.replace('_tab_', '_posvel_')
    idsfile = tabfile.replace('_tab_', '_ids_')
    groupfile = tabfile.replace('subhalo', 'group')
    
    for type in range(6):
        try:
            os.makedirs(outputdir + '/%d' % type)
        except OSError:
            pass
    
    def readsubhalo_par(id):
        f = file(posvelfile % id , mode='r')
        header = numpy.fromfile(f, dtype='i4', count=9)
        size = header[2]
        rt = numpy.empty(shape=size, dtype=pdtype)
    
        rt['pos'] = numpy.fromfile(f, dtype=('f4', 3), count=size)
        rt['vel'] = numpy.fromfile(f, dtype=('f4', 3), count=size)
        rt['type'] = numpy.fromfile(f, dtype='u1', count=size)
        rt['mass'] = numpy.fromfile(f, dtype='f4', count=size)
        # do not use it
        age = numpy.fromfile(f, dtype='f4', count=size)
        rt['id'] = numpy.memmap(idsfile % id, mode='r', offset=28, dtype='u8')
        return rt
    
    def readsubhalo_tab(id):
        f = Snapshot(tabfile % id , 'cmugadget.SubHaloTab', 
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

    def readgroup_tab(id):
        f = Snapshot(groupfile % id , 'cmugadget.GroupTab', 
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
    
    first = Snapshot(tabfile %0, 'cmugadget.SubHaloTab', idtype='u8',
            floattype='f4')
    Nfiles = first.C['Nfiles'] 
    firstgrp = Snapshot(groupfile %0, 'cmugadget.GroupTab', idtype='u8',
            floattype='f4')
    Nfiles2 = firstgrp.C['Nfiles'] 
    depot_par = Depot(dtype=pdtype, nshipment=Nfiles, fetch=readsubhalo_par)
    depot_sub = Depot(dtype=subdtype, nshipment=Nfiles, fetch=readsubhalo_tab)
    depot_group = Depot(dtype=groupdtype, nshipment=Nfiles2, fetch=readgroup_tab)
    
    F = {}
    for type in range(6):
        F[type] = {}
        for field in pdtype.fields:
            F[type][field] = file(outputdir + '/%d/%s.raw' %(type, field), mode='w')
        
    halodump = file(outputdir + '/subhalotab.raw', 'w')
    groupdump = file(outputdir + '/grouphalotab.raw', 'w')
    
    for i in range(Nfiles):
        tab = Snapshot(tabfile % i, 'cmugadget.SubHaloTab', idtype='u8',
                floattype='f4')
        group = depot_group.consume(tab.C['N'][0])
        group['nhalo'] = tab[0, 'nhalo']
        for groupid in range(tab.C['N'][0]):
            nhalo = tab[0, 'nhalo'][groupid]
            realhalo = depot_sub.consume(nhalo)
            halo = numpy.empty(nhalo + 1, dtype=realhalo.dtype)
            halo[:-1] = realhalo
            # last halo is contamination
            halo[-1]['mass'] = numpy.nan
            halo[-1]['pos'] = numpy.nan
            halo[-1]['grnr'] = halo[0]['grnr']
            halo[-1]['len'] = group['len'][groupid] - halo['len'][:-1].sum()
            parstart = numpy.concatenate(([0], halo['len'].cumsum()))
            parend = parstart[1:]
            assert group[groupid]['len'] == tab[0, 'length'][groupid]
            pars = depot_par.consume(group['len'][groupid])
    
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
        group.tofile(groupdump)
        groupdump.flush()
        for type in range(6):
            for field in pdtype.fields:
                F[type][field].flush()

def lookupMain(args):
    type = args.type
    field = args.field
    outputdir = args.outputdir

    Nfiles = args.C['Nfiles']
    id = numpy.fromfile(outputdir + '/%d/id.raw' %(type), dtype='u8')
    arg = id.argsort()
    idsorted = id[arg]
    dtype = args.schema[field].dtype
    valsorted = numpy.empty(len(id), dtype)
#    outfile = numpy.memmap(outputdir + '/%d/%s.raw' %(type, field), mode='w',
#            dtype=first.schema[field].dtype, shape=len(id))

    Nfound = 0
    for fid in range(Nfiles):
        f = Snapshot(args.snapname % fid, 'cmugadget', template=args.first)
        idlookup = f[type, 'id']
        valuelookup = f[type, field]
        ind = idsorted.searchsorted(idlookup)
        found = idsorted.take(ind, mode='wrap') == idlookup
        valsorted.put(ind[found], valuelookup[found])
        Nfound += found.sum()
        print fid, found.sum()
    assert Nfound == len(id)
    del idsorted

    val = numpy.empty_like(valsorted)
    val.put(arg, valsorted)
    val.tofile(outputdir + '/%d/%s.raw' % (type, field))

if __name__ == '__main__':
    main()
