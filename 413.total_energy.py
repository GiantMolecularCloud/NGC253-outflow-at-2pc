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
os.system('mkdir -p '+ratedir)


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

# get energy cube
#################

def get_energy_cube(dataset, ax_type='major', SNR=3.0):

    # load data
    mass_cube = fits.open(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.mass_cube.fits'))[0]
    model     = fits.open(os.path.join(ratedir, 'diskfit.total_model.regrid_'+dataset['line']+'.'+ax_type+'.rotate.fits'))[0]

    # get pixel position of major axis
    ma_y = find_major_axis(dataset)

    # pixel to km conversion
    pix_scale = u.Quantity(str(np.abs(mass_cube.header['cdelt1']))+mass_cube.header['cunit1'])
    pix_to_km = ((distance*np.sin(pix_scale)).to(u.km)).value

    # cube to store mass outflow rates
    energy_cube   = np.full_like(mass_cube.data, np.nan)
    energy_cube_p = np.full_like(mass_cube.data, np.nan)
    energy_cube_m = np.full_like(mass_cube.data, np.nan)
    energy_header = copy.deepcopy(mass_cube.header)
    del energy_header['history']
    energy_header['bunit'] = 'erg'

    # cube extent
    nv, ny, nx = energy_cube.shape

    for v in np.arange(nv):
        velo = (u.Quantity(str((v-mass_cube.header['crpix3'])*mass_cube.header['cdelt3']+mass_cube.header['crval3'])+mass_cube.header['cunit3']).to(u.km/u.s)).value     # km/s

        for y in np.arange(ny):
            for x in np.arange(nx):
                print(dataset['line']+" "+ax_type+" "+str(SNR)+"sigma : pixel ["+str(v)+","+str(x)+","+str(y)+"] of ["+str(nv)+","+str(nx)+","+str(ny)+"]", end='\r')

                mass   = mass_cube.data[v,y,x]*1.98847542e+33   # Msun -> g
                modv   = model.data[ma_y,x]*1e5                 # km/s -> cm/s
                modv_p = model.data[ma_y+25,x]*1e5              # km/s -> cm/s
                modv_m = model.data[ma_y-25,x]*1e5              # km/s -> cm/s

                # energy needs to be saved modulo 1e51 to not overflow a 32bit fits image
                energy   = 0.5*mass*np.abs(velo-modv)**2 /1e51
                energy_p = 0.5*mass*np.abs(velo-modv_p)**2 /1e51
                energy_m = 0.5*mass*np.abs(velo-modv_m)**2 /1e51

                energy_cube[v,y,x]   = energy
                energy_cube_p[v,y,x] = energy_p
                energy_cube_m[v,y,x] = energy_m
    print("\n")

    # mask inf values if they happen somehow
    energy_cube[np.isinf(energy_cube)]     = np.nan
    energy_cube_p[np.isinf(energy_cube_p)] = np.nan
    energy_cube_m[np.isinf(energy_cube_m)] = np.nan

    fits.writeto(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.energy_cube.fits'), data=energy_cube, header=energy_header, overwrite=True)
    fits.writeto(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.energy_cube_p.fits'), data=energy_cube_p, header=energy_header, overwrite=True)
    fits.writeto(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.energy_cube_m.fits'), data=energy_cube_m, header=energy_header, overwrite=True)


###################################################################################################

# get integrated outflow rate
#############################

def integrate_energy(dataset, ax_type='major', SNR=3.0):

    energy_cube   = fits.open(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.energy_cube.fits'))[0]
    energy_cube_p = fits.open(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.energy_cube_p.fits'))[0]
    energy_cube_m = fits.open(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.energy_cube_m.fits'))[0]
    energy   = np.nansum(energy_cube.data)
    energy_p = np.nansum(energy_cube_p.data)
    energy_m = np.nansum(energy_cube_m.data)
    return [energy,energy_p,energy_m]


###################################################################################################

# get integrated outflow rate with defined region
####################################################

def integrate_energy_region(dataset, mask_file, ax_type='major', SNR=3.0):

    randname = 'temp.'+str(int(np.random.rand()*1e4))

    # regrid mask to match the cube
    imregrid(imagename = mask_file,
        template = os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.3.0s.mom1'),
        output   = os.path.join(ratedir, randname+'.mask.regrid'),
        overwrite = True
        )

    # rotate mask to match the rotated cube
    rotate_by = str( 90*u.degree - disk_PA )
    os.system('rm -rf '+os.path.join(ratedir, randname+'.mask.rotate'))
    ia.open(os.path.join(ratedir, randname+'.mask.regrid'))
    temp = ia.rotate(pa = rotate_by,
        outfile  = os.path.join(ratedir, randname+'.mask.rotate'),
        overwrite = True
        )
    temp.close()
    ia.done()
    exportfits(imagename = os.path.join(ratedir, randname+'.mask.rotate'),
        fitsimage    = os.path.join(ratedir, randname+'.mask.rotate.fits'),
        velocity     = True,
        optical      = True,
        dropstokes   = True,
        dropdeg      = True,
        overwrite    = True
        )

    # multiply each channel with the mask
    energy_cube = fits.open(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.energy_cube.fits'))[0]
    mask      = fits.open(os.path.join(ratedir, randname+'.mask.rotate.fits'))[0]
    masked_rate = copy.deepcopy(energy_cube.data)
    nv,ny,nx = energy_cube.data.shape
    for v in np.arange(nv):
        masked_rate[v] *= mask.data

    os.system('rm -rf '+os.path.join(ratedir, randname+'.mask.rotate'))
    os.system('rm -rf '+os.path.join(ratedir, randname+'.mask.rotate.fits'))
    os.system('rm -rf '+os.path.join(ratedir, randname+'.mask.regrid'))

    energy = np.nansum(masked_rate)
    return energy


###################################################################################################

# get outflow energy
####################

# At 100sigma there are only a few very small regions of non-disk emission left in the CO(3-2) cube
# which cannot be rotated. Rotation loses emission by resampling the data leading to an empty image
# in that case. Thus try/except statements are needed to keep the code running.


# generate energy cubes in parallel
# This is pure python and can be done in parallel
def parallel_energy_cubes(inp):
    try:
        print("Executing: "+inp[0]['line']+" "+inp[1]+" "+str(inp[2]))
        get_energy_cube(inp[0], ax_type=inp[1], SNR=inp[2])
    except:
        print("Failed: "+inp[0]['line']+" "+inp[1]+" "+str(inp[2]))

for dataset in datasets:
#    for ax_type in ['major','minor']:
    for ax_type in ['major']:
        plist = [[dataset, ax_type, SNR] for SNR in SNRs]
        pool = Pool(len(SNRs))
        pool.map(parallel_energy_cubes, plist)
        pool.close()


# sequentially read in energies
energies = {}
for dataset in datasets:
    energies[dataset['line']] = {}
    for ax_type in ['major','minor']:
        energies[dataset['line']][ax_type] = {}
        for SNR in SNRs:

            try:
                E, E_p, E_m           = integrate_energy(dataset, ax_type=ax_type, SNR=SNR)
                E_real_outflow        = integrate_energy_region(dataset, os.path.join(ratedir, 'nondisk.'+dataset['line']+'.real_outflow_new.mask'), ax_type=ax_type, SNR=SNR)
                E_real_superbubble    = integrate_energy_region(dataset, os.path.join(ratedir, 'nondisk.'+dataset['line']+'.superbubble_new.mask'), ax_type=ax_type, SNR=SNR)
                E_real_cospatial_disk = integrate_energy_region(dataset, os.path.join(ratedir, 'nondisk.'+dataset['line']+'.cospatial_disk_new.mask'), ax_type=ax_type, SNR=SNR)
            except:
                E, E_p, E_m           = np.nan, np.nan, np.nan
                E_real_outflow        = np.nan
                E_real_superbubble    = np.nan
                E_real_cospatial_disk = np.nan

            energies[dataset['line']][ax_type][SNR] = [E,E_p,E_m,E_real_outflow,E_real_superbubble,E_real_cospatial_disk]

fnpickle(energies, 'energies.pickle')


###################################################################################################

# deprojected energies
######################

dep_energies = {}

for dataset in datasets:
    dep_energies[dataset['line']] = {}
    for ax_type in ['major','minor']:
        dep_energies[dataset['line']][ax_type] = {}
        for SNR in SNRs:
            dep_energies[dataset['line']][ax_type][SNR] = {}

            E, E_p, E_m, E_outflow, E_bubble, E_codisk = energies[dataset['line']][ax_type][SNR]

            dep_E         = [ E/(np.sin(j*u.degree).value)**2 for j in np.arange(48,108) ]
            dep_E_outflow = [ E_outflow/(np.sin(j*u.degree).value)**2 for j in np.arange(48,108) ]
            dep_E_bubble  = [ E_bubble/(np.sin(j*u.degree).value)**2 for j in np.arange(48,108) ]
            dep_E_codisk  = [ E_codisk/(np.sin(j*u.degree).value)**2 for j in np.arange(48,108) ]

            dep_energies[dataset['line']][ax_type][SNR]['total']   = np.nanpercentile(dep_E, [16,50,81])
            dep_energies[dataset['line']][ax_type][SNR]['outflow'] = np.nanpercentile(dep_E_outflow, [16,50,81])
            dep_energies[dataset['line']][ax_type][SNR]['bubble']  = np.nanpercentile(dep_E_bubble, [16,50,81])
            dep_energies[dataset['line']][ax_type][SNR]['codisk']  = np.nanpercentile(dep_E_codisk, [16,50,81])

fnpickle(dep_energies, 'energies_deprojected.pickle')

###################################################################################################
