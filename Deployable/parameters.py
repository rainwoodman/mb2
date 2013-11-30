


# modify and then do python process.py



snapshot_id=18 # z=10
snap_z=10.0

halo_ID=0 #choose one specific halo. Set to False for all.

snapshot_dir='/home/yfeng1/physics/mb2/'   #<--------------------------------- YOU WILL NEED TO CHANGE THIS
halo_file=snapshot_dir + '/galaxy_props/%03d/' % snapshot_id + 'subhalo_tab_%03d' %snapshot_id+'.bin'
stars_file=snapshot_dir+ '/galaxy_props/%03d/' %snapshot_id + 'star_props_%03d' %snapshot_id+'.bin'



# -------------------------------------------
# include nebular emission (See Wilkins et al. (2013d)

include_nebular_emission=True



# -------------------------------------------
# apply Madau et al. (1998) IGM absorption correction

apply_IGM_absorption=True



# -------------------------------------------
# creates full SED for each halo in output_dir/SEDs

output_SEDs=False


# -------------------------------------------
# define output director

output_dir='output'


# -------------------------------------------
# Halo selection criteria

stellar_mass_limit=False #only calculate SED for galaxies with log10(stellar_mass)>stellar_mass_limit. Set to False if no limit.

dark_matter_particle_limit=False  #only calculate SED for galaxies with n_DM>dark_matter_particle_limit. Set to False if no limit.




# -------------------------------------------
# define observed-frame filters - see filters/ for choices (HST.WFC3 is resolved into filters/HST/WFC3/)


obs_filter_sets={}
obs_filter_sets['HST.WFC3']=['f105w','f125w','f160w']
obs_filter_sets['HST.ACS']=['f435w','f606w','f775w','f850lp']
obs_filter_sets['Spitzer.Irac']=['ch1','ch2','ch3','ch4']

# -------------------------------------------
# define rest-frame filters

rf_filter_sets={}
rf_filter_sets['FAKE.FAKE']=['1500','2000','2500','V','Un','Vn','Vw']
rf_filter_sets['SDSS.SDSS']=['u','g','r','i','z']
rf_filter_sets['UKIRT.WFCAM']=['Y','J','H','K']
rf_filter_sets['GALEX.GALEX']=['FUV','NUV']


