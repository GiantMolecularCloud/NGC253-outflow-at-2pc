from __future__ import print_function

########################
# CASA ANALYSIS SCRIPT #
########################

# Get estimates for various quantities for disk (inside disk mask) and the rest (mainly outflows).

###################################################################################################

# import required modules
execfile('scripts/casa_imports.py')
execfile('NGC253/project_info.py')


###################################################################################################

# load sub project info
execfile(os.path.join(projectdir, 'info.py'))
ratedir = ratedir+'.improved3'
os.system('mkdir -p '+resultdir)


###################################################################################################

# regrid and rotate model
#########################

# The model needs to be regridded to the various datasets.

def regrid_model(dataset, ax_type='major'):

    # replace any potential NaNs
    im = fits.open(os.path.join(ratedir, 'diskfit.total_model.fits'))[0]
    im.data[np.isnan(im.data)] = 0.0
    im.writeto(os.path.join(ratedir, 'diskfit.total_model.zeros.fits'), overwrite=True)

    imregrid(imagename = os.path.join(ratedir, 'diskfit.total_model.fits'),
        template = os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.3.0s.mom1'),
        output   = os.path.join(ratedir, 'diskfit.total_model.regrid_'+dataset['line']+'.'+ax_type+'.image'),
        overwrite = True
        )

    rotate_by = str( 90*u.degree - disk_PA )

    os.system('rm -rf '+os.path.join(ratedir, 'diskfit.total_model.regrid_'+dataset['line']+'.'+ax_type+'.rotate.image'))
    ia.open(os.path.join(ratedir, 'diskfit.total_model.regrid_'+dataset['line']+'.'+ax_type+'.image'))
    temp = ia.rotate(pa = rotate_by,
        outfile = os.path.join(ratedir, 'diskfit.total_model.regrid_'+dataset['line']+'.'+ax_type+'.rotate.image'),
        overwrite = True
        )
    temp.close()
    ia.done()
    exportfits(imagename = os.path.join(ratedir, 'diskfit.total_model.regrid_'+dataset['line']+'.'+ax_type+'.rotate.image'),
        fitsimage    = os.path.join(ratedir, 'diskfit.total_model.regrid_'+dataset['line']+'.'+ax_type+'.rotate.fits'),
        dropstokes   = True,
        dropdeg      = True,
        overwrite    = True
        )


###################################################################################################

# rotate data
#############

# Aligning major axis and y-axis simplifies the task of getting distances significantly.

def rotate_image(dataset, ax_type='major', SNR=5.0):

    # replace any potential NaNs
    im = fits.open(os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.fits'))[0]
    im.data[np.isnan(im.data)] = 0.0
    im.writeto(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.zeros.fits'), overwrite=True)

    rotate_by = str( 90*u.degree - disk_PA )

    os.system('rm -rf '+os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.rotate.image'))
    print("rotating "+os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.fits'))
    ia.open(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.zeros.fits'))
    temp = ia.rotate(pa = rotate_by,
        outfile  = os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.rotate.image')
        )
    temp.close()
    ia.done()
    exportfits(imagename = os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.rotate.image'),
        fitsimage    = os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.rotate.fits'),
        velocity     = True,
        optical      = True,
        dropstokes   = True,
        dropdeg      = True,
        overwrite    = True
        )
    os.system('rm -rf '+os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.rotate.image'))


###################################################################################################

# get mass cube
###############

def get_mass_cube(dataset, ax_type='major', SNR=5.0):

    # load data
    nondisk = os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.rotate.fits')
    im = fits.open(nondisk)[0]
    mass_header = copy.deepcopy(im.header)

    # brightness temperature to mass conversion
    pix_scale  = u.Quantity(str(np.abs(im.header['cdelt1']))+im.header['cunit1']).to(u.arcsec)
    chan_width = u.Quantity(str(np.abs(im.header['cdelt3']))+im.header['cunit3']).to(u.km/u.s)
    area_cm2 = ( distance.to(u.cm)*np.sin(pix_scale) )**2
    Tb_to_coldens = 0.5e20*chan_width.value/dataset['Xco_corr']/u.cm**2
    coldens_to_mass  = ((area_cm2 * (1.36*2.0*u.u).to(u.kg) )/(1*u.Msun)).to(1.0*u.cm**2)       # atomic weight of H2: 2    # 1.36 to account for helium

    Tb_to_mass = (Tb_to_coldens*coldens_to_mass).value
    print("1K in 1 pixel "+str(chan_width)+" channels corresponds to "+str(Tb_to_mass)+"Msun")

    # convert to mass
    mass_data = im.data*Tb_to_mass

    # fix header
    mass_header['bunit'] = 'Msun'

    # save mass cube
    fits.writeto(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.mass_cube.fits'), data=mass_data, header=mass_header, overwrite=True)


###################################################################################################

# find major axis
#################

# Get the y coordinate that corresponds to the major axis.
def find_major_axis(dataset):

    im = fits.open(os.path.join(ratedir, 'diskfit.total_model.regrid_'+dataset['line']+'.'+ax_type+'.rotate.fits'))[0]
    im_wcs = WCS(im.header)
    x,y = im_wcs.all_world2pix(kin_center.ra.value, kin_center.dec.value,1)
    major_axis = int(np.round(y))
    return major_axis


###################################################################################################

# get distance map
##################

# create a map of projected distance from the major axis
def get_distance_map(dataset, ax_type='major'):

    # map of correct dimensions
    template = fits.open(os.path.join(ratedir, 'diskfit.total_model.regrid_'+dataset['line']+'.'+ax_type+'.rotate.fits'))[0]

    # three maps for the three launching site estimates
    dmap   = np.full_like(template.data, np.nan)
    dmap_p = np.full_like(template.data, np.nan)
    dmap_m = np.full_like(template.data, np.nan)

    ma_y = find_major_axis(dataset)

    # pixel to km conversion
    pix_scale = u.Quantity(str(np.abs(template.header['cdelt1']))+template.header['cunit1'])
    pix_to_km = ((distance*np.sin(pix_scale)).to(u.km)).value

    # uncertainty in origin/major axis: 1.25"
    delta_pix = int(1.25*u.arcsec/pix_scale)

    ny,nx = dmap.shape
    for y in np.arange(ny):

        # distance in pixel and physical
        dist   = pix_to_km * np.abs(y-ma_y)                   # km
        dist_p = pix_to_km * np.abs(y-ma_y+delta_pix)       # +1.25"
        dist_m = pix_to_km * np.abs(y-ma_y-delta_pix)       # -1.25"

        # write to map
        dmap[y]   = np.full_like(dmap[y], dist)
        dmap_p[y] = np.full_like(dmap[y], dist_p)
        dmap_m[y] = np.full_like(dmap[y], dist_m)

    # correct header
    template.header['bunit'] = 'km'

    # write to disk
    fits.writeto(os.path.join(ratedir, dataset['cube']+'.distance.'+ax_type+'.fits'), data=dmap, header=template.header, overwrite=True)
    fits.writeto(os.path.join(ratedir, dataset['cube']+'.distance_p.'+ax_type+'.fits'), data=dmap_p, header=template.header, overwrite=True)
    fits.writeto(os.path.join(ratedir, dataset['cube']+'.distance_m.'+ax_type+'.fits'), data=dmap_m, header=template.header, overwrite=True)


###################################################################################################

# rotate mask
#############

def rotate_mask(dataset, ax_type, mask_file):
    """
    Rotate the mask to the major axis frame. Preparation for ouflow rate calculation.
    """

    # replace any potential NaNs
    im = fits.open(mask_file)[0]
    im.data[np.isnan(im.data)] = 0.0
    im.writeto(mask_file+'.zeros', overwrite=True)

    # regrid mask to match the cube
    imregrid(imagename = mask_file,
        template = os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.3.0s.mom1'),
        output   = os.path.join(ratedir, mask_file+'.mask.regrid'),
        overwrite = True
        )

    # rotate mask to match the rotated cube
    rotate_by = str( 90*u.degree - disk_PA )
    os.system('rm -rf '+os.path.join(ratedir, mask_file+'.mask.rotate'))
    ia.open(os.path.join(ratedir, mask_file+'.mask.regrid'))
    temp = ia.rotate(pa = rotate_by,
        outfile  = os.path.join(ratedir, mask_file+'.mask.rotate'),
        overwrite = True
        )
    temp.close()
    ia.done()
    exportfits(imagename = os.path.join(ratedir, mask_file+'.mask.rotate'),
        fitsimage    = os.path.join(ratedir, mask_file+'.mask.rotate.fits'),
        velocity     = True,
        optical      = True,
        dropstokes   = True,
        dropdeg      = True,
        overwrite    = True
        )


###################################################################################################

# get relative velocity cube
############################

def get_velo(dataset, ax_type='major'):
    """
    A cube containing the relative velocity between each pixel and the corresponding model velocity
    """

    print("Rate cube: "+dataset['line']+" "+ax_type)

    # load data
    mass  = fits.open(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.5.0s.mass_cube.fits'))[0]
    model = fits.open(os.path.join(ratedir, 'diskfit.total_model.regrid_'+dataset['line']+'.'+ax_type+'.rotate.fits'))[0]

    # cube to store mass outflow rates
    velo   = np.full_like(mass.data, np.nan)
    velo_p = np.full_like(mass.data, np.nan)
    velo_m = np.full_like(mass.data, np.nan)
    velo_header = copy.deepcopy(mass.header)
    del velo_header['history']
    velo_header['bunit'] = 'km/s'

    ma_y = find_major_axis(dataset)

    # uncertainty in origin/major axis: 1.25", get corresponding velocity shift
    pix_scale = u.Quantity(str(np.abs(mass.header['cdelt3']))+mass.header['cunit3'])
    delta_pix = int(50*u.km/u.s/pix_scale)

    # cube extent
    nv, ny, nx = velo.shape
    for v in np.arange(nv):
        vel = (u.Quantity(str((v-mass.header['crpix3'])*mass.header['cdelt3']+mass.header['crval3'])+mass.header['cunit3']).to(u.km/u.s)).value     # km/s

        for y in np.arange(ny):
            for x in np.arange(nx):
                print(dataset['line']+" "+ax_type+": pixel ["+str(v)+","+str(x)+","+str(y)+"] of ["+str(nv)+","+str(nx)+","+str(ny)+"]", end='\r')

                velo[v,y,x]   = np.abs(vel - model.data[ma_y,x])
                velo_p[v,y,x] = np.abs(vel - model.data[ma_y+delta_pix,x])
                velo_m[v,y,x] = np.abs(vel - model.data[ma_y-delta_pix,x])
    print("\n")

    fits.writeto(os.path.join(ratedir, dataset['cube']+'.velo.fits'), data=velo, header=velo_header, overwrite=True)
    fits.writeto(os.path.join(ratedir, dataset['cube']+'.velo_p.fits'), data=velo_p, header=velo_header, overwrite=True)
    fits.writeto(os.path.join(ratedir, dataset['cube']+'.velo_m.fits'), data=velo_m, header=velo_header, overwrite=True)


###################################################################################################

# re-order error cubes to make computation easier
################################################

def restructure_errors(dataset, ax_type='major'):
    """
    Re-order the error maps/cubes such that the positive error cube contains all the upper errors
    and the negative error cube contains all the lower errors.
    """

    for a1,a2 in [[os.path.join(ratedir, dataset['cube']+'.distance_p.'+ax_type+'.fits'),
                   os.path.join(ratedir, dataset['cube']+'.distance_m.'+ax_type+'.fits')],
                  [os.path.join(ratedir, dataset['cube']+'.velo_p.fits'),
                   os.path.join(ratedir, dataset['cube']+'.velo_m.fits')]]:

        # load cubes
        x1 = fits.open(a1)[0]
        x2 = fits.open(a2)[0]

        # create new cubes that store the reordered values, copy/deepcopy is NOT possible ("I/O error on closed file")
        p = fits.open(a1)[0]
        m = fits.open(a2)[0]

        # reorder values
        p.data[ (x1.data > x2.data) ] = x1.data[ (x1.data > x2.data) ]
        p.data[ (x1.data < x2.data) ] = x2.data[ (x1.data < x2.data) ]
        m.data[ (x1.data < x2.data) ] = x1.data[ (x1.data < x2.data) ]
        m.data[ (x1.data > x2.data) ] = x2.data[ (x1.data > x2.data) ]

        # save to disk
        p.writeto(a1+'.test', overwrite=True)
        m.writeto(a2+'.test', overwrite=True)


###################################################################################################
# get deprojected outflow rate, energy and momentum
###################################################################################################

# bootstrap errors by assuming a flat dstribution of inclinations

def time_evolution(dataset, ax_type='major', SNR=5.0, binwidth=0.1, maskfile='', sampling=1.0):

    print("Dataset: "+dataset['line']+" "+ax_type+" "+str(SNR)+"sigma "+maskfile+" sampling: "+str(sampling)+"degrees")

    rates_j    = []
    energies_j = []
    momenta_j  = []

    # assume flat distribution between 48degrees and 108 degrees
    inclinations = np.arange(48,108.01,sampling)
    n_incl       = len(inclinations)
    for j in inclinations:
        print("\tinclination ["+str(j)+"/48-108]: "+dataset['line']+" "+ax_type+" "+str(SNR)+"sigma "+maskfile+" sampling: "+str(sampling)+"degrees")

        # load data
        mass   = fits.open(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.mass_cube.fits'))[0]
        dmap   = fits.open(os.path.join(ratedir, dataset['cube']+'.distance.'+ax_type+'.fits'))[0]
        dmap_p = fits.open(os.path.join(ratedir, dataset['cube']+'.distance_p.'+ax_type+'.fits'))[0]
        dmap_m = fits.open(os.path.join(ratedir, dataset['cube']+'.distance_m.'+ax_type+'.fits'))[0]
        velo   = fits.open(os.path.join(ratedir, dataset['cube']+'.velo.fits'))[0]
        velo_p = fits.open(os.path.join(ratedir, dataset['cube']+'.velo_p.fits'))[0]
        velo_m = fits.open(os.path.join(ratedir, dataset['cube']+'.velo_m.fits'))[0]

        # correct inclinations
        dmap.data   = dmap.data   /np.abs(np.sin(j/180.*np.pi))
        dmap_p.data = dmap_p.data /np.abs(np.sin(j/180.*np.pi))
        dmap_m.data = dmap_m.data /np.abs(np.sin(j/180.*np.pi))
        velo.data   = velo.data   /np.abs(np.cos(j/180.*np.pi))
        velo_p.data = velo_p.data /np.abs(np.cos(j/180.*np.pi))
        velo_m.data = velo_m.data /np.abs(np.cos(j/180.*np.pi))

        # mask the maps/cubes if mask is set
        if not ( maskfile == '' ):
            mask = fits.open(os.path.join(ratedir, 'nondisk.'+dataset['line']+'.'+maskfile+'.mask.fits.mask.rotate.fits'))[0]

            # mask distance map
            dmap.data   *= mask.data
            dmap_p.data *= mask.data
            dmap_m.data *= mask.data

            # mask cubes
            nv,ny,nx = mass.data.shape
            for v in np.arange(nv):
                mass.data[v]   *= mask.data
                velo.data[v]   *= mask.data
                velo_p.data[v] *= mask.data
                velo_m.data[v] *= mask.data

        # convert distance map to cube to get matching array shapes
        nv, ny, nx = mass.shape
        dcube   = np.array([dmap.data   for v in np.arange(nv)])
        dcube_p = np.array([dmap_p.data for v in np.arange(nv)])
        dcube_m = np.array([dmap_m.data for v in np.arange(nv)])

        # pixel crossing time scale, the time over which the gas was ejected
        pix_scale = u.Quantity(str(np.abs(mass.header['cdelt1']))+mass.header['cunit1'])
        pix_km    = (distance*np.sin(pix_scale)).to(u.km).value
        t_pix     = (pix_km / velo.data)   /60/60/24/365                # sec -> yr
        t_pix_p   = (pix_km / velo_p.data) /60/60/24/365                # sec -> yr
        t_pix_m   = (pix_km / velo_m.data) /60/60/24/365                # sec -> yr

        # ejection time scale
        t_eject   = dcube/velo.data     /60/60/24/365/1e6          # sec -> Myr
        t_eject_p = dcube_p/velo_p.data /60/60/24/365/1e6          # sec -> Myr
        t_eject_m = dcube_m/velo_m.data /60/60/24/365/1e6          # sec -> Myr

        # mass outflow rate at pixel x at time t_eject
        mdot   = mass.data/t_pix
        mdot_p = mass.data/t_pix_p
        mdot_m = mass.data/t_pix_m

        # energy at pixel x at time t_eject
        # the conversion factor from Msun km^2 s^-2 to erg is 1.988475e+43
        energy   = 0.5*mass.data*velo.data**2     *1.988475e+43       # Msun km^2 s^-2 -> erg
        energy_p = 0.5*mass.data*velo_p.data**2 *1.988475e+43       # Msun km^2 s^-2 -> erg
        energy_m = 0.5*mass.data*velo_m.data**2 *1.988475e+43       # Msun km^2 s^-2 -> erg

        # momentum at pixel x at time t_eject
        momentum   = mass.data*velo.data                        # Msun km/s
        momentum_p = mass.data*velo_p.data                    # Msun km/s
        momentum_m = mass.data*velo_m.data                    # Msun km/s

        # bin outflow rates/energies/momenta by ejection time in bins of binwidth Myr
        hbw = binwidth/2.
        for b in np.arange(0,10,binwidth):
            selection   = (t_eject  >b-hbw) & (t_eject  <b+hbw)
            selection_p = (t_eject_p>b-hbw) & (t_eject_p<b+hbw)
            selection_m = (t_eject_m>b-hbw) & (t_eject_m<b+hbw)
            rate   = np.abs(np.nansum( mdot[selection]     *t_pix[selection]     ) /(binwidth*1e6))
            rate_1 = np.abs(np.nansum( mdot_p[selection_p] *t_pix_p[selection_p] ) /(binwidth*1e6))
            rate_2 = np.abs(np.nansum( mdot_m[selection_m] *t_pix_m[selection_m] ) /(binwidth*1e6))
            ene    = np.abs(np.nansum( energy[selection]     ))
            ene_1  = np.abs(np.nansum( energy_p[selection_p] ))
            ene_2  = np.abs(np.nansum( energy_m[selection_m] ))
            mom    = np.abs(np.nansum( momentum[selection]     ))
            mom_1  = np.abs(np.nansum( momentum_p[selection_p] ))
            mom_2  = np.abs(np.nansum( momentum_m[selection_m] ))

            rates_j.append([b,j,rate])
            rates_j.append([b,j,rate_1])
            rates_j.append([b,j,rate_2])
            energies_j.append([b,j,ene])
            energies_j.append([b,j,ene_1])
            energies_j.append([b,j,ene_2])
            momenta_j.append([b,j,mom])
            momenta_j.append([b,j,mom_1])
            momenta_j.append([b,j,mom_2])

    # convert to array for easier handling
    rates_j    = np.array(rates_j)
    energies_j = np.array(energies_j)
    momenta_j  = np.array(momenta_j)

    fnpickle(rates_j,    os.path.join(resultdir, dataset['cube']+'.'+ax_type+'.'+str(SNR)+'s.rate_time_bootstrap.'+maskfile+'.pickle'))
    fnpickle(energies_j, os.path.join(resultdir, dataset['cube']+'.'+ax_type+'.'+str(SNR)+'s.energy_time_bootstrap.'+maskfile+'.pickle'))
    fnpickle(momenta_j,  os.path.join(resultdir, dataset['cube']+'.'+ax_type+'.'+str(SNR)+'s.momentum_time_bootstrap.'+maskfile+'.pickle'))

    # get percentiles in each time bin as error estimation
    rates    = []
    energies = []
    momenta  = []
    hbw = binwidth/2.
    for b in np.arange(0,10,binwidth):
        bin_rates = np.abs(rates_j[:,2][(rates_j[:,0] > b-hbw) & (rates_j[:,0] < b+hbw)])
        rate      = np.percentile(bin_rates, 50)
        rate_m    = np.percentile(bin_rates, 16)
        rate_p    = np.percentile(bin_rates, 84)
        rates.append([b,rate,rate_p,rate_m])
    for b in np.arange(0,10,binwidth):
        bin_energies = np.abs(energies_j[:,2][(energies_j[:,0] > b-hbw) & (energies_j[:,0] < b+hbw)])
        rate      = np.percentile(bin_energies, 50)
        rate_m    = np.percentile(bin_energies, 16)
        rate_p    = np.percentile(bin_energies, 84)
        energies.append([b,rate,rate_p,rate_m])
    for b in np.arange(0,10,binwidth):
        bin_momenta = np.abs(momenta_j[:,2][(momenta_j[:,0] > b-hbw) & (momenta_j[:,0] < b+hbw)])
        rate      = np.percentile(bin_momenta, 50)
        rate_m    = np.percentile(bin_momenta, 16)
        rate_p    = np.percentile(bin_momenta, 84)
        momenta.append([b,rate,rate_p,rate_m])

    # write out values
    np.savetxt(os.path.join(resultdir, dataset['cube']+'.'+ax_type+'.'+str(SNR)+'s.rate_time_bootstrap.'+maskfile+'.txt'),
               rates,
               fmt = ['%7.1f','%9.3f','%9.3f','%9.3f'],
               header = " time      rate    rate_p    rate_m\n[Myr] [Msun/yr] [Msun/yr] [Msun/yr]"
              )
    np.savetxt(os.path.join(resultdir, dataset['cube']+'.'+ax_type+'.'+str(SNR)+'s.energy_time_bootstrap.'+maskfile+'.txt'),
               energies,
               fmt = ['%7.1f','%9.3f','%9.3f','%9.3f'],
               header = " time      energy energy_p energy_m\n[Myr]     [erg]     [erg]     [erg]"
              )
    np.savetxt(os.path.join(resultdir, dataset['cube']+'.'+ax_type+'.'+str(SNR)+'s.momentum_time_bootstrap.'+maskfile+'.txt'),
               momenta,
               fmt = ['%7.1f','%9.3f','%9.3f','%9.3f'],
               header = " time     momentum  mom_p     mom_m\n[Myr][Msunkm/s][Msunkm/s][Msunkm/s]"
              )

    print("Finished: "+dataset['line']+" "+ax_type+" "+str(SNR)+"sigma "+maskfile+" sampling: "+str(sampling)+"degrees")


###################################################################################################

# get distance evolution with errors
####################################

# bootstrap errors by assuming a flat dstribution of inclinations

def distance_evolution(dataset, ax_type='major', SNR=5.0, binwidth=1.0, maskfile='', sampling=1.0):

    print("Dataset: "+dataset['line']+" "+ax_type+" "+str(SNR)+"sigma "+maskfile+" sampling: "+str(sampling)+"degrees")

    rates_j    = []
    energies_j = []
    momenta_j  = []

    # assume flat distribution between 48degrees and 108 degrees
    inclinations = np.arange(48,108.01,sampling)
    n_incl       = len(inclinations)
    for j in inclinations:
        print("inclination ["+str(j)+"/48-108]: "+dataset['line']+" "+ax_type+" "+str(SNR)+"sigma "+maskfile+" sampling: "+str(sampling)+"degrees")

        # load data
        mass   = fits.open(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.mass_cube.fits'))[0]
        dmap   = fits.open(os.path.join(ratedir, dataset['cube']+'.distance.'+ax_type+'.fits'))[0]
        dmap_p = fits.open(os.path.join(ratedir, dataset['cube']+'.distance_p.'+ax_type+'.fits'))[0]
        dmap_m = fits.open(os.path.join(ratedir, dataset['cube']+'.distance_m.'+ax_type+'.fits'))[0]
        velo   = fits.open(os.path.join(ratedir, dataset['cube']+'.velo.fits'))[0]
        velo_p = fits.open(os.path.join(ratedir, dataset['cube']+'.velo_p.fits'))[0]
        velo_m = fits.open(os.path.join(ratedir, dataset['cube']+'.velo_m.fits'))[0]

        # correct inclinations
        dmap.data   = dmap.data   /np.abs(np.sin(j/180.*np.pi))
        dmap_p.data = dmap_p.data /np.abs(np.sin(j/180.*np.pi))
        dmap_m.data = dmap_m.data /np.abs(np.sin(j/180.*np.pi))
        velo.data   = velo.data   /np.abs(np.cos(j/180.*np.pi))
        velo_p.data = velo_p.data /np.abs(np.cos(j/180.*np.pi))
        velo_m.data = velo_m.data /np.abs(np.cos(j/180.*np.pi))

        # mask the maps/cubes if mask is set
        if not ( maskfile == '' ):
            mask = fits.open(os.path.join(ratedir, 'nondisk.'+dataset['line']+'.'+maskfile+'.mask.fits.mask.rotate.fits'))[0]

            # mask distance map
            dmap.data   *= mask.data
            dmap_p.data *= mask.data
            dmap_m.data *= mask.data

            # mask cubes
            nv,ny,nx = mass.data.shape
            for v in np.arange(nv):
                mass.data[v]   *= mask.data
                velo.data[v]   *= mask.data
                velo_p.data[v] *= mask.data
                velo_m.data[v] *= mask.data

        # convert distance map to cube to get matching array shapes
        nv, ny, nx = mass.shape
        dcube   = np.array([dmap.data   for v in np.arange(nv)])
        dcube_p = np.array([dmap_p.data for v in np.arange(nv)])
        dcube_m = np.array([dmap_m.data for v in np.arange(nv)])

        # pixel crossing time scale, the time over which the gas was ejected
        pix_scale = u.Quantity(str(np.abs(mass.header['cdelt1']))+mass.header['cunit1'])
        pix_km    = (distance*np.sin(pix_scale)).to(u.km).value
        t_pix     = np.abs( (pix_km / velo.data)   /60/60/24/365 )                # sec -> yr
        t_pix_p   = np.abs( (pix_km / velo_p.data) /60/60/24/365 )                # sec -> yr
        t_pix_m   = np.abs( (pix_km / velo_m.data) /60/60/24/365 )                # sec -> yr

        # mass outflow rate at pixel x at time t_eject
        mdot   = mass.data/t_pix
        mdot_p = mass.data/t_pix_p
        mdot_m = mass.data/t_pix_m

        # energy at pixel x at time t_eject
        # the conversion factor from Msun km^2 s^-2 to erg is 1.988475e+43
        energy   = 0.5*mass.data*velo.data**2   *1.988475e+43       # Msun km^2 s^-2 -> erg
        energy_p = 0.5*mass.data*velo_p.data**2 *1.988475e+43       # Msun km^2 s^-2 -> erg
        energy_m = 0.5*mass.data*velo_m.data**2 *1.988475e+43       # Msun km^2 s^-2 -> erg

        # momentum at pixel x at time t_eject
        momentum   = mass.data*velo.data                        # Msun km/s
        momentum_p = mass.data*velo_p.data                    # Msun km/s
        momentum_m = mass.data*velo_m.data                    # Msun km/s

        # to bin over distance, it needs to be in units of arcsec
        pix_km_arcsec = pix_scale.to(u.arcsec).value/pix_km
        pix_scale_as  = pix_scale.to(u.arcsec).value
        dcube   *= pix_km_arcsec
        dcube_p *= pix_km_arcsec
        dcube_m *= pix_km_arcsec

        # bin 0.1Myr averaged outflow rates by distance in bins of 1.0arcsec
        hbw = binwidth/2.
        rates    = []
        energies = []
        momenta  = []
        for b in np.arange(0,40,binwidth):
            print("\tbin ["+str(b-hbw)+","+str(b+hbw)+"]", end='\r')
            selection   = (dcube  >b-hbw) & (dcube  <b+hbw)
            selection_p = (dcube_p>b-hbw) & (dcube_p<b+hbw)
            selection_m = (dcube_m>b-hbw) & (dcube_m<b+hbw)
            rate   = np.abs(np.nansum( mdot[selection]     ) *pix_scale_as/binwidth)
            rate_1 = np.abs(np.nansum( mdot_p[selection_p] ) *pix_scale_as/binwidth)
            rate_2 = np.abs(np.nansum( mdot_m[selection_m] ) *pix_scale_as/binwidth)
            ene    = np.abs(np.nansum( energy[selection]     ))
            ene_1  = np.abs(np.nansum( energy_p[selection_p] ))
            ene_2  = np.abs(np.nansum( energy_m[selection_m] ))
            mom    = np.abs(np.nansum( momentum[selection]     ))
            mom_1  = np.abs(np.nansum( momentum_p[selection_p] ))
            mom_2  = np.abs(np.nansum( momentum_m[selection_m] ))

            rates_j.append([b,j,rate])
            rates_j.append([b,j,rate_1])
            rates_j.append([b,j,rate_2])
            energies_j.append([b,j,ene])
            energies_j.append([b,j,ene_1])
            energies_j.append([b,j,ene_2])
            momenta_j.append([b,j,mom])
            momenta_j.append([b,j,mom_1])
            momenta_j.append([b,j,mom_2])

    # convert to array for easier handling
    rates_j    = np.array(rates_j)
    energies_j = np.array(energies_j)
    momenta_j  = np.array(momenta_j)

    fnpickle(rates_j,    os.path.join(resultdir, dataset['cube']+'.'+ax_type+'.'+str(SNR)+'s.rate_dist_bootstrap.'+maskfile+'.pickle'))
    fnpickle(energies_j, os.path.join(resultdir, dataset['cube']+'.'+ax_type+'.'+str(SNR)+'s.energy_dist_bootstrap.'+maskfile+'.pickle'))
    fnpickle(momenta_j,  os.path.join(resultdir, dataset['cube']+'.'+ax_type+'.'+str(SNR)+'s.mometum_dist_bootstrap.'+maskfile+'.pickle'))

    # get percentiles in each time bin as error estimation
    rates    = []
    energies = []
    momenta  = []
    hbw = binwidth/2.
    for b in np.arange(0,40,binwidth):
        bin_rates = np.abs(rates_j[:,2][(rates_j[:,0] > b-hbw) & (rates_j[:,0] < b+hbw)])
        rate      = np.percentile(bin_rates, 50)
        rate_m    = np.percentile(bin_rates, 16)
        rate_p    = np.percentile(bin_rates, 84)
        rates.append([b,rate,rate_p,rate_m])
    for b in np.arange(0,40,binwidth):
        bin_energies = np.abs(energies_j[:,2][(energies_j[:,0] > b-hbw) & (energies_j[:,0] < b+hbw)])
        rate      = np.percentile(bin_energies, 50)
        rate_m    = np.percentile(bin_energies, 16)
        rate_p    = np.percentile(bin_energies, 84)
        energies.append([b,rate,rate_p,rate_m])
    for b in np.arange(0,40,binwidth):
        bin_momenta = np.abs(momenta_j[:,2][(momenta_j[:,0] > b-hbw) & (momenta_j[:,0] < b+hbw)])
        rate      = np.percentile(bin_momenta, 50)
        rate_m    = np.percentile(bin_momenta, 16)
        rate_p    = np.percentile(bin_momenta, 84)
        momenta.append([b,rate,rate_p,rate_m])

    # write out values
    np.savetxt(os.path.join(resultdir, dataset['cube']+'.'+ax_type+'.'+str(SNR)+'s.rate_dist_bootstrap.'+maskfile+'.txt'),
               rates,
               fmt = ['%7.1f','%9.3f','%9.3f','%9.3f'],
               header = " dist      rate    rate_p    rate_m\n [''] [Msun/yr] [Msun/yr] [Msun/yr]"
              )
    np.savetxt(os.path.join(resultdir, dataset['cube']+'.'+ax_type+'.'+str(SNR)+'s.energy_dist_bootstrap.'+maskfile+'.txt'),
               energies,
               fmt = ['%7.1f','%9.3f','%9.3f','%9.3f'],
               header = " dist      energy energy_p energy_m\n ['']     [erg]     [erg]     [erg]"
              )
    np.savetxt(os.path.join(resultdir, dataset['cube']+'.'+ax_type+'.'+str(SNR)+'s.momentum_dist_bootstrap.'+maskfile+'.txt'),
               momenta,
               fmt = ['%7.1f','%9.3f','%9.3f','%9.3f'],
               header = " dist     momentum  mom_p     mom_m\n''[][Msunkm/s][Msunkm/s][Msunkm/s]"
              )

    print("Finished: "+dataset['line']+" "+ax_type+" "+str(SNR)+"sigma "+maskfile+" sampling: "+str(sampling)+"degrees")


###################################################################################################

# recent outflow rate
#####################

def averaged_rate(qtype, ax_type='major', SNR=5.0, window=[0.,1.], err_type='value'):
    """
    Display the outflow rate averaged over the given window (in Myr for time or arcsec for distance).
    """

    def load_quantity(qtype):
        q = {}
        for dataset in datasets:
            q[dataset['line']] = {}
            for ax_type in ['major']: #,'minor']:
                q[dataset['line']][ax_type] = {}
                for SNR in [5.0]: #SNRs:
                    q[dataset['line']][ax_type][SNR] = {}
                    for mask in ['','real_outflow','superbubble','cospatial_disk']:

                        try:
                            q[dataset['line']][ax_type][SNR][mask] = np.genfromtxt(os.path.join(resultdir, dataset['cube']+'.'+ax_type+'.'+str(SNR)+'s.'+qtype+'.'+mask+'.txt'),
                                              names = ('bin','q','q_p','q_m')
                                             )
                        except:
                            q[dataset['line']][ax_type][SNR][mask] = []
        return q

    def avg_q(q,dataset,ax_type,SNR,mask,window):

        bins   = q[dataset['line']][ax_type][SNR][mask]['bin']
        rate   = q[dataset['line']][ax_type][SNR][mask]['q']
        rate_p = q[dataset['line']][ax_type][SNR][mask]['q_p']
        rate_m = q[dataset['line']][ax_type][SNR][mask]['q_m']

        binwidth = bins[1] - bins[0]

        selection  = (bins>window[0]) & (bins<window[1])
        avg   = np.nansum( rate[selection] )   *binwidth / (window[1]-window[0])
        avg_p = np.nansum( rate_p[selection] ) *binwidth / (window[1]-window[0])
        avg_m = np.nansum( rate_m[selection] ) *binwidth / (window[1]-window[0])

        return avg, avg_p, avg_m

    def dex_err_q(avg,avg_p,avg_m):
        err_p_dex = np.log10(avg_p/avg)
        err_m_dex = np.log10(avg_m/avg)
        return avg,err_p_dex,err_m_dex

    def lin_err_q(avg,avg_p,avg_m):
        err_p = avg_p-avg
        err_m = avg-avg_m
        return avg,err_p,err_m

    def print_q(q):
        print("\nDe-projected "+qtype+" rate averaged over ["+str(window[0])+","+str(window[1])+"]")
        print(' '*19+'   {:<24}   {:<24}   {:<24}'.format(*[d['line'] for d in datasets]))
        if err_type=='value':
            print(' '*16+('   '+'{:^8}{:^8}{:^8}'.format('value','value+','value-'))*3)
        elif err_type=='linear':
            print(' '*16+('   '+'{:^8}{:^8}{:^8}'.format('value','err+','err-'))*3)
        elif err_type=='dex':
            print(' '*16+('   '+'{:^8}{:^8}{:^8}'.format('value','dex+','dex-'))*3)

        for maskfile in ['', 'real_outflow', 'superbubble', 'cospatial_disk']:
            print('{:>15}'.format(maskfile), end='')
            for dataset in datasets:
                avg, avg_p, avg_m = avg_q(q,dataset,ax_type,SNR,maskfile,window)
                if err_type=='value':
                    avg, err_p, err_m = (avg, avg_p, avg_m)
                    fmt = '{:8.1e}{:8.1e}{:8.1e}'
                elif err_type=='linear':
                    avg, err_p, err_m = lin_err_q(avg, avg_p, avg_m)
                    fmt = '{:8.1e}{:8.1e}{:8.1e}'
                elif err_type=='dex':
                    avg, err_p, err_m = dex_err_q(avg, avg_p, avg_m)
                    fmt = '{:8.1e}{:8.2f}{:8.2f}'

                print('   ', end='')
                if 'rate' in qtype:
                    print('{:8.1f}{:8.2f}{:8.2f}'.format( avg, err_p, err_m ), end='')
                elif 'energy' in qtype:
                    print(fmt.format( avg, err_p, err_m ), end='')
                elif 'momentum' in qtype:
                    print(fmt.format( avg, err_p, err_m ), end='')
            print('')

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        print_q(load_quantity(qtype))


###################################################################################################
# get time evolution
###################################################################################################

def parallel_times(inp):
    dataset  = inp[0]
    ax_type  = inp[1]
    SNR      = inp[2]
    maskfile = inp[3]
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            time_evolution(dataset, ax_type, SNR, binwidth=0.1, maskfile=maskfile, sampling=1.0)
    except:
        print("Failed: "+dataset['line']+" "+ax_type+" "+str(SNR)+" "+maskfile+" deproject: "+str(deproject))

inps = [[d, a, SNR, maskfile] for d in datasets for a in ['major'] for SNR in [5.0] for maskfile in ['', 'real_outflow', 'superbubble', 'cospatial_disk']]
pool = Pool(20)
pool.map(parallel_times, inps)
pool.close()


###################################################################################################
# get distance evolution
###################################################################################################

def parallel_distances(inp):
    dataset  = inp[0]
    ax_type  = inp[1]
    SNR      = inp[2]
    maskfile = inp[3]
    try:
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            distance_evolution(dataset, ax_type, SNR, binwidth=1.0, maskfile=maskfile, sampling=1.0)
    except:
        print("Failed: "+dataset['line']+" "+ax_type+" "+str(SNR)+" "+maskfile+" deproject: "+str(deproject))

inps = [[d, a, SNR, maskfile] for d in datasets for a in ['major'] for SNR in [5.0] for maskfile in ['', 'real_outflow', 'superbubble', 'cospatial_disk']]
pool = Pool(20)
pool.map(parallel_distances, inps)
pool.close()


###################################################################################################
# calculate average quantities
###################################################################################################

averaged_rate('rate_time_bootstrap',     window=[0.,1.],  err_type='value')
averaged_rate('energy_time_bootstrap',   window=[0.,1.],  err_type='value')
averaged_rate('momentum_time_bootstrap', window=[0.,1.],  err_type='value')
averaged_rate('rate_dist_bootstrap',     window=[0.,20.], err_type='value')
averaged_rate('energy_dist_bootstrap',   window=[0.,20.], err_type='value')
averaged_rate('momentum_dist_bootstrap', window=[0.,20.], err_type='value')

averaged_rate('rate_time_bootstrap',     window=[0.,1.],  err_type='linear')
averaged_rate('energy_time_bootstrap',   window=[0.,1.],  err_type='linear')
averaged_rate('momentum_time_bootstrap', window=[0.,1.],  err_type='linear')
averaged_rate('rate_dist_bootstrap',     window=[0.,20.], err_type='linear')
averaged_rate('energy_dist_bootstrap',   window=[0.,20.], err_type='linear')
averaged_rate('momentum_dist_bootstrap', window=[0.,20.], err_type='linear')

averaged_rate('rate_time_bootstrap',     window=[0.,1.],  err_type='dex')
averaged_rate('energy_time_bootstrap',   window=[0.,1.],  err_type='dex')
averaged_rate('momentum_time_bootstrap', window=[0.,1.],  err_type='dex')
averaged_rate('rate_dist_bootstrap',     window=[0.,20.], err_type='dex')
averaged_rate('energy_dist_bootstrap',   window=[0.,20.], err_type='dex')
averaged_rate('momentum_dist_bootstrap', window=[0.,20.], err_type='dex')


###################################################################################################
###################################################################################################
