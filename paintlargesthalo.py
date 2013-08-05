import numpy
from gaepsi.gaplot import GaplotContext
from gaepsi.gaplot import nl_, n_
from gaepsi.cosmology import WMAP7
from sys import argv
from nishitab import subhalodtype as nishitabdtype
from matplotlib.figure import Figure
from matplotlib import cm
snapid = int(argv[1])
if len(argv) == 3:
   dostar = True
else:
   dostar = False
snappath = '/physics/yfeng1/mb2/snapdir/snapdir_%03d/snapshot_%03d.%%d' %(snapid,
        snapid)
mappath = '/physics/yfeng1/mb2/snapmesh/snapmesh_%03d.map' % (snapid)
nishitabpath = '/physics/yfeng1/mb2/galaxy_props/%03d/subhalo_tab_%03d.bin' \
          % (snapid, snapid)
N = numpy.fromfile(nishitabpath, 'i4', 1)[0]
nishitab = \
        numpy.memmap('/physics/yfeng1/mb2/galaxy_props/%03d/subhalo_tab_%03d.bin' \
          % (snapid, snapid), mode='r', offset=4, dtype=nishitabdtype, shape=N)

ct = GaplotContext()
ct.use(snappath, 'cmugadget', periodic=True)
print ct.C['Ntot']
print ct.C['redshift']

halomass = numpy.array([100, 200, 500, 800, 1000, 2000, 5000, 10000, 20000, 50000], dtype='f8')
haloarg = [numpy.abs(nishitab['mass'] / mass - 1.0).argmin()\
         for mass in halomass]

center = nishitab['pos'][haloarg]
vel = nishitab['vel'][haloarg]
mass = nishitab['mass'][haloarg]
#radius = WMAP7.Rvir(nishitab['mass'][haloarg], z=5.0) * 1.5
radius = WMAP7.Rvir(nishitab['mass'][haloarg], z=ct.C['redshift'])
print center, mass, radius
raise
def paintgas(center, size):
    ct.schema('gas', 0, ['mass', 'sml', 'ie', 'ye'])
    C = numpy.zeros((len(center), ct.shape[0], ct.shape[1]))
    L = numpy.zeros((len(center), ct.shape[0], ct.shape[1]))
    
    print len(ct.fids())
    
    fids = ct.fids()
    for i in range(0, len(fids), 16):
        print 'doing files', i
        gas = ct.read('gas', fids=fids[i:i+16], np=8)
        ct.makeT('gas')
        for c, s, cc, ll in zip(center, size, C, L): 
            ct.view(center=c, size=s)
            dC, dL = ct.paint('gas', 'T', 'mass', 'sml')
            dC[dL == 0] = 0
            cc += dC * dL
            ll += dL
    return C/L, L

def paintstar(center, size, vel):
    ct.schema('star', 4, ['mass', 'sft'])
    C = numpy.zeros((len(center), ct.shape[0], ct.shape[1]))
    L = numpy.zeros((len(center), ct.shape[0], ct.shape[1]))
    
    star = ct.read('star', np=8)
    star.smooth(ct.T['star'])
    star['T'] = ct.cosmology.a2t(star['sft'])
    for c, s, cc, ll, v in zip(center, size, C, L, vel): 
        ct.view(center=c, size=s)
        dC, dL = ct.paint('star', 'T', 'mass', 'sml')
        dC[dL == 0] = 0
        cc += dC * dL
        ll += dL
    return C/L, L


if dostar:
  Cs, Ls = paintstar(center, radius * 0.2, vel)
  for i, c, s, cs, ls in zip(range(len(Cs)), center, radius * 0.2, Cs, Ls):
    figure = Figure(dpi=200)
    ax = ct.newaxes(figure)
    ct.view(center=c, size=s)
    now = ct.cosmology.z2t(ct.C['redshift'])
    past = ct.cosmology.z2t(11.0)
    ct.imshow(n_(cs, vmin=past, vmax=now), nl_(ls, vmin='25db'),
            cmap=cm.coolwarm_r)
    ct.frame(axis=False, scale={'color':'w'})
    ct.print_svg('%03d/star-%05d.svg' % (snapid, int(mass[i])))
else:  
  Cg, Lg = paintgas(center, radius * 10)
  for i, c, s, cg, lg in zip(range(len(Cg)), center, radius * 10, Cg, Lg):
    figure = Figure(dpi=200)
    ax = ct.newaxes(figure)
    ct.view(center=c, size=s)
    ct.imshow(nl_(cg, vmin=3.5, vmax=7.5), nl_(lg, vmin='40db'))
    ct.frame(axis=False, scale={'color':'w'})
    ct.print_svg('%03d/gas-%05d.svg' % (snapid, int(mass[i])))

