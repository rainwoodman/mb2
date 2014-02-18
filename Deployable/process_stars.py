




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
    # for each age
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

Nbins = numpy.prod(sspa['age'].shape)


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
    os.system('mkdir -p '+ snap.subhalodir + '/4/ObsFilter/' + filter_set_key)
    for f in p.obs_filter_sets[filter_set_key]:     
        obs_outputs[filter_set_key][f] = numpy.memmap(snap.filename('4/ObsFilter',
            filter_set_key + '/' + f), shape=Nbins, dtype='f4', mode='w+')
    
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
    os.system('mkdir -p '+ snap.subhalodir + '/4/RfFilter/' + filter_set_key)
    for f in p.rf_filter_sets[filter_set_key]:   
    
        rf_outputs[filter_set_key][f] = numpy.memmap(snap.filename('4/RfFilter',
            filter_set_key + '/' + f), shape=Nbins, dtype='f4', mode='w+')
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

starmet = snap.load(4, 'met')
starsft = snap.load(4, 'sft')


print 'Total Number of Stars:',len(starmet)

print '-------------------------------'
print '-------------------------------'

z =(1./starsft)-1
age = numpy.interp(z, redshifts, ages) * 1000.

Zi = numpy.abs(sspa['Z'][:, None] - starmet[None, :]).argmin(axis=0)
agei = numpy.empty(len(Zi), dtype='intp')
for k in range(len(sspa)):
    mask = k == Zi
    agei[mask] = \
        numpy.abs(sspa['age'][k][None, :] - age[mask][:, None]).argmin(axis=-1)

# use rind to quickly access the Zi, agei from the raveled arrays)
rind = numpy.int64(numpy.ravel_multi_index((Zi, agei), sspa['age'].shape))
rind.tofile(snap.filename(4, 'SEDindex'))
stellar_mass_rec_frac = numpy.float32(sspa['frac']).take(rind)
stellar_mass_rec_frac.tofile(snap.filename(4, 'recfrac'))
print stellar_mass_rec_frac.dtype
# stars in the same rind have identical sed.
for rind0 in numpy.arange(Nbins):

    print 'doing', rind0, '/', Nbins
    Zi0, agei0 = numpy.unravel_index(rind0, sspa['age'].shape)

    # get the idential sed template
    ste_sed = sspa['sed'].reshape(-1, len(lam))[rind0]

    # use bincount to combine sum on the same nebula line set
    Hbeta_flux=sspa['Hbeta_flux'].ravel()[rind0]
    nebula = sspa['nebula'][Zi0]

    #------------------------------------------------
    # add nebula lines 
            
    if p.include_nebular_emission==True:
        stellar_sed_lam=ste_sed*(3.*10**8)*(10**10)/(lam**2)

        FWHM=velocity*nebula_lam/(299792.) #l in \AA, velocity in kms, c in kms
        sigma=FWHM/2.3548                        
        sed_lam =  stellar_sed_lam + \
                (1. / (sigma[None, :]* numpy.sqrt(2.*math.pi)) * \
                numpy.exp(-((lam[:, None]-nebula_lam[None,
                    :])**2)/(2.*sigma[None, :]**2))* Hbeta_flux * nebula[None, :]).sum(axis=-1)
   
        sed=sed_lam/((3.*10**8)*(10**10)/(lam**2)) #convert back to f_nu
    else:
        sed=ste_sed

            
    #-------
    # apply IGM absorption
    
    if p.apply_IGM_absorption==True:       
        sed=sed*igm
    
    #--------
    # determine fluxes in each band 
    for fs, f in obs_filters:
        flux = integrate.trapz(1 / lamz * sed*obs_fT[f],x=lamz)/integrate.trapz(1 / lamz * obs_fT[f],x=lamz)
        # flux is a scalar!
        # multiply by starmass to get real filtervalue
        obs_outputs[fs][f][rind0] = 1e10 * flux / 1e28
    for fs, f in rf_filters:
        flux = integrate.trapz(1 / lam * sed*rest_fT[f],x=lam)/integrate.trapz(1 / lam * rest_fT[f],x=lam)
        # flux is a scalar!
        rf_outputs[fs][f][rind0] = 1e10 * flux / 1e28

for fs in obs_outputs:
    for f in obs_outputs[fs]:
        obs_outputs[fs][f].flush()
        
for fs in rf_outputs:
    for f in rf_outputs[fs]:
        rf_outputs[fs][f].flush()
