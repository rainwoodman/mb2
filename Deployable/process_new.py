




import numpy as numpy
import sys
import os
from scipy import integrate
import math


import parameters_new as p

sys.path.append('subroutines/')

sys.path.append(p.snapshot_dir + '/scripts')
import readsubhalo
from nebulaemission import *
#import readMB
import cosmocalc
import population_synthesis_mod as population_synthesis
import utilities as u
import continuum
from sys import argv

snap = readsubhalo.SnapDir(argv[1], p.snapshot_dir)


#------------------------------------------------
# make output directories
output_dir = snap.subhalodir + '/SED'
os.system('mkdir '+ output_dir)


#------------------------------------------------
#------------------------------------------------
#determine relations between redshift and age
aou_cz=cosmocalc.cosmocalc(snap.redshift)['zage_Gyr']
redshifts=numpy.arange(snap.redshift,20,0.1)
ages=[]
for z in redshifts:
    aou_z=cosmocalc.cosmocalc(z)['zage_Gyr']
    ages.append(aou_cz-aou_z)
ages=numpy.array(ages)


#------------------------------------------------
#------------------------------------------------
# read in SSPs
mets=numpy.array([0.05,0.02,0.008,0.004,0.0004])
metsl={0.05:'005',0.02:'002',0.008:'0008',0.004:'0004',0.0004:'00004'}

ssp={}
for Z in mets:ssp[Z]=getattr(population_synthesis,'pegase')('Salpeter/i.z'+metsl[Z])

# now we convert the ssps to an array

sspdtype = numpy.dtype([
    ('Z', 'f8'),
    ('frac', ('f8', 52)),
    ('lam', ('f8', 1220)),
    ('sed', ('f8', (52, 1220))),
    ('bolometric_luminosity', ('f8', 52)),
    ('age', ('f8', 52)),
    ('Hbeta_flux', ('f8', 52)),
    ('nebula', ('f8', len(nebula_lines))),
    ])

nebula_lam = numpy.empty(len(nebula_lines), 'f8')
nebula_name = numpy.empty(len(nebula_lines), 'S20')
sspa = numpy.empty(len(mets), dtype=sspdtype)
lines = list(nebula_lines.keys())

for j, line in enumerate(lines):
    nebula_lam[j] = nebula_lines[line]['l']
    nebula_name[j] = line

for i, Z in enumerate(mets):
    sspitem = sspa[i]
    sspitem['Z'] = Z
    sspi = getattr(population_synthesis,'pegase')('Salpeter/i.z'+metsl[Z])
    sspitem['frac'][:] = sspi['frac']
    sspitem['lam'][:] = sspi['lam']
    sspitem['sed'][:] = sspi['sed']
    sspitem['bolometric_luminosity'][:] = sspi['bolometric_luminosity']
    sspitem['age'][:] = sspi['age']
    sspitem['Hbeta_flux'][:] = sspi['nebula_lines']['fluxes']['Hbeta']
    for j, line in enumerate(lines):
        sspitem['nebula'][j] = nebula_lines[line][Z]

lam=sspa[0]['lam']

lamz=lam*(1.+snap.redshift)

#------------------------------------------------
# save wavelength grid
numpy.save('lam.npy',lam)

#------------------------------------------------
# build IGM array
igm=continuum.expteff(lam*(1.+snap.redshift),snap.redshift)


subtab = snap.readsubhalo()

#------------------------------------------------
#------------------------------------------------
# Define observed Filters

#------------------------------------------------
#Read in filter transmission curves, and interpolate onto wavelength grid
obs_fT={}
obs_filters=[]
obs_outputs = {}
for filter_set_key in p.obs_filter_sets.keys():
    obs_outputs[filter_set_key] = {}
    os.system('mkdir -p '+ snap.subhalodir + '/subhalo/ObsFilter/' + filter_set_key)
    for f in p.obs_filter_sets[filter_set_key]:     
        obs_outputs[filter_set_key][f] = numpy.memmap(snap.filename('subhalo',
            'ObsFilter/' + filter_set_key + '/' + f), shape=len(subtab), dtype='f4', mode='w+')
    
        obs_filters.append((filter_set_key, f))
    
        fname='filters/'+filter_set_key.split('.')[0]+'/'+filter_set_key.split('.')[1]+'/'+f+'.txt'
        d=u.readc(fname,5)
        filter_trans_l=numpy.array(map(float,d[0]))
        filter_trans_T=numpy.array(map(float,d[1]))            
        maxT=numpy.max(filter_trans_T)
        filter_trans_T[numpy.where(filter_trans_T<0.05*maxT)]=0.0                                
        obs_fT[f]=numpy.interp(lamz,filter_trans_l,filter_trans_T)

#------------------------------------------------
#------------------------------------------------
# Define rest_frame filters - suffixed '_r'
  
#------------------------------------------------
#Read in filter transmission curves, and interpolate onto wavelength grid
rest_fT={}
rf_filters=[]
rf_outputs = {}
for filter_set_key in p.rf_filter_sets.keys():
    rf_outputs[filter_set_key] = {}
    os.system('mkdir -p '+ snap.subhalodir + '/subhalo/RfFilter/' + filter_set_key)
    for f in p.rf_filter_sets[filter_set_key]:   
        rf_outputs[filter_set_key][f] = numpy.memmap(snap.filename('subhalo',
            'RfFilter/' + filter_set_key + '/' + f), shape=len(subtab), dtype='f4', mode='w+')
    
        rf_filters.append((filter_set_key, f))
       
        fname='filters/'+filter_set_key.split('.')[0]+'/'+filter_set_key.split('.')[1]+'/'+f+'.txt'
        d=u.readc(fname,5)
        filter_trans_l=numpy.array(map(float,d[0]))
        filter_trans_T=numpy.array(map(float,d[1]))            
        maxT=numpy.max(filter_trans_T)
        filter_trans_T[numpy.where(filter_trans_T<0.05*maxT)]=0.0                                
        rest_fT[f]=numpy.interp(lam,filter_trans_l,filter_trans_T)





#------------------------------------------------
#------------------------------------------------
# Main Code

starmass = snap.load(4, 'mass', subtab)
starmet = snap.load(4, 'met', subtab)
starsft = snap.load(4, 'sft', subtab)

print 'Total Number of Halos:',len(subtab)

print '-------------------------------'
print '-------------------------------'


selectedhalos = ((subtab['lenbytype'][:, 4] > 0) & ~numpy.isnan(subtab['mass'])).nonzero()[0]
SEDindex = numpy.memmap(snap.subhalodir + '/SED/index.raw',
        shape=len(subtab), dtype=('i8'), mode='w+')
SEDindex[:] = -1
SEDindex[selectedhalos] = numpy.arange(len(selectedhalos))
SEDindex.flush()

StellarSEDs = numpy.memmap(snap.subhalodir + '/SED/stellar.raw',
        shape=len(selectedhalos), dtype=('f4', len(lam)), mode='w+')
FullSEDs = numpy.memmap(snap.subhalodir + '/SED/full.raw',
        shape=len(selectedhalos), dtype=('f4', len(lam)), mode='w+')

for i, id in enumerate(selectedhalos):
        
    halo=subtab[id]
    stellar_mass=halo['massbytype'][4]*10**10
    n_DM=halo['lenbytype'][1]

    #------------------------------------------------
    # construct stellar SED and nebula lines
       
    #print 'Halo ID:',id, halo['lenbytype'][4]
    if i % 1000 == 0:
        print 'Halo', i, id, halo['lenbytype'][4]
    stellar_mass_wrec=0.0 #stellar masses corrected for recycling
    
    halostarmass = starmass[id]
    halostarmet = starmet[id]
    halostarsft = starsft[id]

    z =(1./halostarsft)-1
    age = numpy.interp(z, redshifts, ages) * 1000.

    Zi = numpy.abs(sspa['Z'][:, None] - halostarmet[None, :]).argmin(axis=0)
    agei = numpy.empty(len(Zi), dtype='intp')
    for k in range(len(sspa)):
        mask = k == Zi
        agei[mask] = \
            numpy.abs(sspa['age'][k][None, :] - age[mask][:, None]).argmin(axis=-1)
    # use rind to quickly access the Zi, agei from the raveled arrays)
    rind = numpy.ravel_multi_index((Zi, agei), sspa['age'].shape)

    stellar_mass_wrec = 1e10 * (halostarmass * sspa['frac'].take(rind)).sum(dtype='f8')

    # use bincount to combine sum on the same sed
    wt = numpy.bincount(rind, weights=1e10 * halostarmass, 
            minlength=numpy.prod(sspa['age'].shape))
    ste_sed = (sspa['sed'].reshape(-1, len(lam)) * wt[:, None]).sum(axis=0)
    Hbeta_flux=halostarmass*(1e10)*sspa['Hbeta_flux'].take(rind)

    # use bincount to combine sum on the same nebula line set
    # bin on metalicity
    wt = numpy.bincount(Zi, weights=Hbeta_flux, 
            minlength=len(sspa))
    nbl = (wt[:, None] * sspa['nebula']).sum(axis=0)

    stellar_sed = ste_sed.copy()
    #------------------------------------------------
    # save stellar mass with recycling
    
    #print 'log10(M_*_rec):',numpy.round(numpy.log10(stellar_mass_wrec),3)
    #print 'fraction of initial mass remaining:',numpy.round(stellar_mass_wrec/stellar_mass,3)
    
    #------------------------------------------------
    # add nebula lines 
            
    if p.include_nebular_emission==True:
    
        stellar_sed_lam=stellar_sed*(3.*10**8)*(10**10)/(lam**2)

        FWHM=velocity*nebula_lam/(299792.) #l in \AA, velocity in kms, c in kms
        sigma=FWHM/2.3548                        
        sed_lam =  stellar_sed_lam + \
                (1. / (sigma[None, :]* numpy.sqrt(2.*math.pi)) * \
                numpy.exp(-((lam[:, None]-nebula_lam[None,
                    :])**2)/(2.*sigma[None, :]**2))*nbl[None, :]).sum(axis=-1)
   
        sed=sed_lam/((3.*10**8)*(10**10)/(lam**2)) #convert back to f_nu
    
    else:
        
        sed=stellar_sed
    
     
            
    #-------
    # apply IGM absorption
    
    if p.apply_IGM_absorption==True:       
        sed=sed*igm
    
    
    # save seds
    StellarSEDs[i][:] = stellar_sed
    FullSEDs[i][:] = sed
     
    #--------
    # determine fluxes in each band 
    for fs, f in obs_filters:
        flux = integrate.trapz(1 / lamz * sed*obs_fT[f],x=lamz)/integrate.trapz(1 / lamz * obs_fT[f],x=lamz)
        obs_outputs[fs][f][id] = flux / 1e28
    for fs, f in rf_filters:
        flux = integrate.trapz(1 / lam * sed*rest_fT[f],x=lam)/integrate.trapz(1 / lam * rest_fT[f],x=lam)
        rf_outputs[fs][f][id] = flux / 1e28

StellarSEDs.flush()
FullSEDs.flush()
for fs in obs_outputs:
    for f in obs_outputs[fs]:
        obs_outputs[fs][f].flush()
for fs in rf_outputs:
    for f in rf_outputs[fs]:
        rf_outputs[fs][f].flush()
        
