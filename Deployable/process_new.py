




import numpy as numpy
import sys
import os
from scipy import integrate
import math


import parameters as p

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
output_dir = p.output_dir + '/%03d/' % snap.snapid
os.system('mkdir '+ output_dir)

if p.output_SEDs==True:
    os.system('mkdir '+ output_dir+'/neb')
    os.system('mkdir '+ output_dir+'/SEDs')
    os.system('mkdir '+ output_dir+'/StellarSEDs')


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
mets=[0.05,0.02,0.008,0.004,0.0004]
metsl={0.05:'005',0.02:'002',0.008:'0008',0.004:'0004',0.0004:'00004'}

ssp={}
for Z in mets:ssp[Z]=getattr(population_synthesis,'pegase')('Salpeter/i.z'+metsl[Z])

lam=ssp[0.02]['lam']
lamz=lam*(1.+snap.redshift)

#------------------------------------------------
# save wavelength grid
numpy.save('lam.npy',lam)

#------------------------------------------------
# build IGM array
igm=continuum.expteff(lam*(1.+snap.redshift),snap.redshift)


#------------------------------------------------
#------------------------------------------------
# Define observed Filters

#------------------------------------------------
#Read in filter transmission curves, and interpolate onto wavelength grid
obs_fT={}
obs_filters=[]
for filter_set_key in p.obs_filter_sets.keys():
    for f in p.obs_filter_sets[filter_set_key]:     
    
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
for filter_set_key in p.rf_filter_sets.keys():
    for f in p.rf_filter_sets[filter_set_key]:   
    
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


#------------------------------------------------
#------------------------------------------------
# setup outputs
out_cols=['index','gas_mass','BH_mass','DM_mass','sfr','stellar_mass']+['stellar_mass_wrec']
properties={}


#MBHS = readMB.MBHaloStars(p.halo_file,p.stars_file)
subtab = snap.readsubhalo()
subtabsfr = snap.load('subhalo', 'sfr')
outputs = {}
for filter_set_key in p.obs_filter_sets.keys():
    outputs[filter_set_key] = numpy.memmap(snap.filename('subhalo',
        filter_set_key), shape=len(subtab), dtype=readsubhalo.bandfilterdtype[filter_set_key], mode='w+')

for filter_set_key in p.rf_filter_sets.keys():
    outputs[filter_set_key] = numpy.memmap(snap.filename('subhalo',
        filter_set_key), shape=len(subtab), dtype=readsubhalo.bandfilterdtype[filter_set_key], mode='w+')

starmass = snap.load(4, 'mass', subtab)
starmet = snap.load(4, 'met', subtab)
starsft = snap.load(4, 'sft', subtab)

print 'Total Number of Halos:',len(subtab)

print '-------------------------------'
print '-------------------------------'


if p.halo_ID==False:
    do_halos=len(subtab) #do all halos
else:
    do_halos=1


for i in range(do_halos):


    if p.halo_ID==False:
        id=i    
    else:
        id=p.halo_ID
        
    halo=subtab[id]
    sfr=subtabsfr[id]
    stellar_mass=halo['massbytype'][4]*10**10
    gas_mass=halo['massbytype'][0]*10**10
    BH_mass=halo['massbytype'][5]*10**10
    DM_mass=halo['massbytype'][1]*10**10
    n_DM=halo['lenbytype'][1]

    
    #------------------------------------------------
    # check the halo meets the stellar mass and number of DM particles limits
    
    halo_selected=True
    
    if p.stellar_mass_limit != False:
        if numpy.log10(stellar_mass)<p.stellar_mass_limit:
            halo_selected=False
            
    if p.dark_matter_particle_limit != False:
        if n_DM<p.dark_matter_particle_limit:
            halo_selected=False    
    
    
    #------------------------------------------------
    # construct stellar SED and nebula lines
       
    if halo_selected==True:
        
        properties['index']=id
        properties['gas_mass']=gas_mass
        properties['DM_mass']=DM_mass
        properties['stellar_mass']=stellar_mass
        properties['BH_mass']=BH_mass
        properties['sfr']=sfr
        
        
        print '-------------------------------'
        print '-------------------------------'
        
        print 'Halo ID:',id
        print 'log10(M_DM):',numpy.round(numpy.log10(DM_mass),3)
        print 'log10(M_*):',numpy.round(numpy.log10(stellar_mass),3)
        
        
        stellar_mass_wrec=0.0 #stellar masses corrected for recycling
        
        
        halostarmass = starmass[id]
        halostarmet = starmet[id]
        halostarsft = starsft[id]
        for si in range(halo['lenbytype'][4]):
            mass=halostarmass[si]
            z=(1./halostarsft[si])-1
            
            Z=halostarmet[si]
            Zi=u.find_nearest(mets,Z)
            Zn=mets[Zi]
            
            age=numpy.interp(z,redshifts,ages)*1000.
            agei=u.find_nearest(ssp[Zn]['age'],age)
            agen=ssp[Zn]['age'][agei]
            
            stellar_mass_wrec+=mass*ssp[Zn]['frac'][agei]
            
                        
            ste_sed=(10**10)*mass*ssp[Zn]['sed'][agei] # stellar SED for just this star!
            
            Hbeta_flux=mass*(10**10)*ssp[Zn]['nebula_lines']['fluxes']['Hbeta'][agei]
                                                              
            if si==0:
                stellar_sed=ste_sed
                lam=ssp[Zn]['lam']                
                nbl={}
                for line in nebula_lines.keys():nbl[line]=Hbeta_flux*nebula_lines[line][Zn]                                                       
            else:
                stellar_sed=stellar_sed+ste_sed
                for line in nebula_lines.keys():nbl[line]=nbl[line]+Hbeta_flux*nebula_lines[line][Zn]  

                    

        #------------------------------------------------
        # save stellar mass with recycling
        
        stellar_mass_wrec=stellar_mass_wrec*10**10
        properties['stellar_mass_wrec']=stellar_mass_wrec
        
        print 'log10(M_*_rec):',numpy.round(numpy.log10(stellar_mass_wrec),3)
        print 'fraction of initial mass remaining:',numpy.round(stellar_mass_wrec/stellar_mass,3)
        
        
                
        #------------------------------------------------
        # add nebula lines 
                
        if p.include_nebular_emission==True:
        
            stellar_sed_lam=stellar_sed*(3.*10**8)*(10**10)/(lam**2)

            sed_lam=stellar_sed_lam
            for line in nebula_lines.keys():
                l=nebula_lines[line]['l']
                FWHM=velocity*l/(299792.) #l in \AA, velocity in kms, c in kms
                sigma=FWHM/2.3548                        
                sed_lam=sed_lam+nbl[line]*(1./(sigma*numpy.sqrt(2.*math.pi)))*(numpy.exp(-((lam-l)**2)/(2.*sigma**2)))
       
            sed=sed_lam/((3.*10**8)*(10**10)/(lam**2)) #convert back to f_nu
        
        else:
            
            sed=stellar_sed
        
         
                
        #-------
        # apply IGM absorption
        
        if p.apply_IGM_absorption==True:       
            sed=sed*igm
        
        
        if p.output_SEDs==True:
         
             #-------
             # save full SED
             numpy.save(output_dir+ '/SEDs/%09d' % id +'.npy',sed) 
             numpy.save(output_dir+ '/StellarSEDs/%09d' % id +'.npy',stellar_sed) 
         
        #--------
        # determine fluxes in each band 
        for fs, f in obs_filters:
            properties[f]=integrate.trapz(sed*obs_fT[f],x=lamz)/integrate.trapz(obs_fT[f],x=lamz)
            outputs[fs][f][id] = properties[f] / 1e28
        for fs, f in rf_filters:
            properties[f+'_r']=integrate.trapz(sed*rest_fT[f],x=lam)/integrate.trapz(rest_fT[f],x=lam)
            outputs[fs][f][id] = properties[f + '_r'] / 1e28
        
        print '---- --------------------- ----'
        print '---- rest-frame photometry ----'
        
        for fs, f in rf_filters:
            print f,numpy.round(properties[f+'_r']/10**28,3)      

        print '-- observed-frame photometry --'
        
        for fs, f in obs_filters:
            print f,numpy.round(properties[f]/10**28,3)    

