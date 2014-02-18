
import numpy as numpy
import sys
import os
from scipy import integrate
import os.path

sys.path.append(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'subroutines/'))

from nebulaemission import *
#import readMB
import cosmocalc
import population_synthesis_mod as population_synthesis
import utilities as u
import continuum
from sys import argv

def readsspa():
    #------------------------------------------------
    #------------------------------------------------
    # read in SSPs
    mets=numpy.array([0.05,0.02,0.008,0.004,0.0004])
    metsl={0.05:'005',0.02:'002',0.008:'0008',0.004:'0004',0.0004:'00004'}

    dir = os.path.dirname(os.path.abspath(__file__))

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
        sspi = population_synthesis.pegase(os.path.join('Salpeter',
            'i.z'+metsl[Z]), dir=dir)
        sspitem['frac'][:] = sspi['frac']
        sspitem['lam'][:] = sspi['lam']
        sspitem['sed'][:] = sspi['sed']
        sspitem['bolometric_luminosity'][:] = sspi['bolometric_luminosity']
        sspitem['age'][:] = sspi['age']
        # for each age
        sspitem['Hbeta_flux'][:] = sspi['nebula_lines']['fluxes']['Hbeta']
        for j, line in enumerate(lines):
            sspitem['nebula'][j] = nebula_lines[line][Z]
    return sspa

def measurefilters(sspa, redshift, obs_filter_sets={}, rf_filter_sets={}, 
        nebula_sets={}, include_nebula_emission=True, apply_IGM_absorption=True):
    dir = os.path.dirname(__file__)
    Nbins = numpy.prod(sspa['age'].shape)
    lam=sspa[0]['lam']
    nebula_lam = numpy.empty(len(nebula_lines), 'f8')
    nebula_name = numpy.empty(len(nebula_lines), 'S20')
    nebula_use = numpy.empty(len(nebula_lines), '?')

    lines = list(nebula_lines.keys())

    for j, line in enumerate(lines):
        nebula_lam[j] = nebula_lines[line]['l']
        nebula_name[j] = line
        if line in nebula_sets:
            nebula_use[j] = True
        else:
            nebula_use[j] = False
    lamz=lam*(1.+redshift)
    #------------------------------------------------
    # build IGM array
    igm=continuum.expteff(lam*(1.+redshift), redshift)

    #------------------------------------------------
    #------------------------------------------------
    # Define observed Filters

    #------------------------------------------------
    #Read in filter transmission curves, and interpolate onto wavelength grid
    obs_fT={}
    obs_filters=[]
    obs_outputs = {}
    for filter_set_key in obs_filter_sets.keys():
        obs_outputs[filter_set_key] = {}
        for f in obs_filter_sets[filter_set_key]:     
            obs_outputs[filter_set_key][f] = numpy.empty(shape=Nbins, dtype='f4')
        
            obs_filters.append((filter_set_key, f))
        
            n1, n2 = filter_set_key.split('.')
            fname=os.path.join(dir, 'filters', n1, n2, f+'.txt')
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
    for filter_set_key in rf_filter_sets.keys():
        rf_outputs[filter_set_key] = {}
        for f in rf_filter_sets[filter_set_key]:   
        
            rf_outputs[filter_set_key][f] = numpy.empty(shape=Nbins, dtype='f4')
            rf_filters.append((filter_set_key, f))
           
            n1, n2 = filter_set_key.split('.')
            fname=os.path.join(dir, 'filters', n1, n2, f+'.txt')
            d=u.readc(fname,5)
            filter_trans_l=numpy.array(map(float,d[0]))
            filter_trans_T=numpy.array(map(float,d[1]))            
            maxT=numpy.max(filter_trans_T)
            filter_trans_T[numpy.where(filter_trans_T<0.05*maxT)]=0.0                                
            rest_fT[f]=numpy.interp(lam,filter_trans_l,filter_trans_T)

    nebula_outputs = {}
    for key in nebula_sets:
        nebula_outputs[key] = numpy.zeros(shape=Nbins, dtype='f4')

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
                
        if include_nebula_emission==True:
            stellar_sed_lam=ste_sed*(3.*10**8)*(10**10)/(lam**2)

            FWHM=velocity*nebula_lam/(299792.) #l in \AA, velocity in kms, c in kms
            sigma=FWHM/2.3548                        
            sed_lam =  stellar_sed_lam + \
                    (1. / (sigma[None, :]* numpy.sqrt(2.*numpy.pi)) * \
                    numpy.exp(-((lam[:, None]-nebula_lam[None,
                        :])**2)/(2.*sigma[None, :]**2))* Hbeta_flux * nebula[None, :]).sum(axis=-1)
       
            sed=sed_lam/((3.*10**8)*(10**10)/(lam**2)) #convert back to f_nu
        else:
            sed=ste_sed

                
        #-------
        # apply IGM absorption
        
        if apply_IGM_absorption==True:       
            sed=sed*igm
            nebula_igm_factor = continuum.expteff(nebula_lam*(1. + redshift), redshift)
        else:
            neubla_igm_factor = 1
        for j in range(len(nebula)):
            if nebula_use[j]:
                nebula_outputs[nebula_name[j]][rind0] = \
                        Hbeta_flux * nebula[j] * 1e10 / 1e28 * \
                        nebula_igm_factor[j]

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
    return obs_outputs, rf_outputs, nebula_outputs

def process_stars(sspa, redshift, starmet, starsft):
    """ sft is in scale factor a """

    #------------------------------------------------
    #------------------------------------------------
    #determine relations between redshift and age
    aou_cz=cosmocalc.cosmocalc(redshift)['zage_Gyr']
    redshifts=numpy.arange(redshift,20,0.1)
    ages=[]
    for z in redshifts:
        aou_z=cosmocalc.cosmocalc(z)['zage_Gyr']
        ages.append(aou_cz-aou_z)
    ages=numpy.array(ages)

#------------------------------------------------
#------------------------------------------------
# Main Code

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
    stellar_mass_rec_frac = numpy.float32(sspa['frac']).take(rind)

    print stellar_mass_rec_frac.dtype
    return rind, stellar_mass_rec_frac
