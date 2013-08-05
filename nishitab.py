import numpy
subhalodtype = numpy.dtype([
    ('mass', 'f4'),      #mass of the halo in units of 10^10Msun/h
    ('len', 'i4'),        #number of particles in this halo
    ('pos', ('f4', 3)),   #position x,y,z of the potential minimum of the halo in units of kpc/h
    ('vel', ('f4', 3)),   #velocity of the halo in km/s
    ('vdisp', 'f4'),    #velocity dispersion of the halo in km/s
    ('vcirc', 'f4'),    #maximum circular velocity of the halo in km/s
    ('rcirc', 'f4'),    #radius at the maximum circular velocity of the halo in km/s
    ('grnr', 'i4'),       #parent fof group (You may not need this information)
    ('massbytype', ('f4', 6)),  #mass of each  component (0-gas, 1-dark matter, 4-stars, 5-BH) in 10^10Msun/h
    ('lenbytype', ('u4',6)),  #number of particles of each type (type same as above)
    ('sfr', 'f4'),      #SFR in Msun/yr 
    ('staroffset', 'u4'), #the starting index for reading stellar properties in star_props_0xx.bin
    ])

