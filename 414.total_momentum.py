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

# get momentum cube
#################

# Create a cube with the mass outflow rate of each pixel as an intermediate step.

def get_momentum_cube(dataset, ax_type='major', SNR=3.0):

    # load data
    mass_cube = fits.open(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.mass_cube.fits'))[0]
    model     = fits.open(os.path.join(ratedir, 'diskfit.total_model.regrid_'+dataset['line']+'.'+ax_type+'.rotate.fits'))[0]

    # get pixel position of major axis
    ma_y = find_major_axis(dataset)

    # pixel to km conversion
    pix_scale = u.Quantity(str(np.abs(mass_cube.header['cdelt1']))+mass_cube.header['cunit1'])
    pix_to_km = ((distance*np.sin(pix_scale)).to(u.km)).value

    # cube to store mass outflow rates
    momentum_cube   = np.full_like(mass_cube.data, np.nan)
    momentum_cube_p = np.full_like(mass_cube.data, np.nan)
    momentum_cube_m = np.full_like(mass_cube.data, np.nan)
    momentum_header = copy.deepcopy(mass_cube.header)
    del momentum_header['history']
    momentum_header['bunit'] = 'erg'

    # cube extent
    nv, ny, nx = momentum_cube.shape

    for v in np.arange(nv):
        velo = (u.Quantity(str((v-mass_cube.header['crpix3'])*mass_cube.header['cdelt3']+mass_cube.header['crval3'])+mass_cube.header['cunit3']).to(u.km/u.s)).value     # km/s

        for y in np.arange(ny):
            for x in np.arange(nx):
                print(dataset['line']+" "+ax_type+" "+str(SNR)+"sigma : pixel ["+str(v)+","+str(x)+","+str(y)+"] of ["+str(nv)+","+str(nx)+","+str(ny)+"]", end='\r')

                mass   = mass_cube.data[v,y,x]      # Msun -> g
                modv   = model.data[ma_y,x]         # km/s -> cm/s
                modv_p = model.data[ma_y+25,x]      # km/s -> cm/s
                modv_m = model.data[ma_y-25,x]      # km/s -> cm/s

                momentum   = mass*np.abs(velo-modv)
                momentum_p = mass*np.abs(velo-modv_p)
                momentum_m = mass*np.abs(velo-modv_m)

                momentum_cube[v,y,x]   = momentum
                momentum_cube_p[v,y,x] = momentum_p
                momentum_cube_m[v,y,x] = momentum_m
    print("\n")

    # mask inf values if they happen somehow
    momentum_cube[np.isinf(momentum_cube)]     = np.nan
    momentum_cube_p[np.isinf(momentum_cube_p)] = np.nan
    momentum_cube_m[np.isinf(momentum_cube_m)] = np.nan

    fits.writeto(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.momentum_cube.fits'), data=momentum_cube, header=momentum_header, overwrite=True)
    fits.writeto(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.momentum_cube_p.fits'), data=momentum_cube_p, header=momentum_header, overwrite=True)
    fits.writeto(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.momentum_cube_m.fits'), data=momentum_cube_m, header=momentum_header, overwrite=True)


###################################################################################################

# get integrated outflow rate
#############################

def integrate_momentum(dataset, ax_type='major', SNR=3.0):

    momentum_cube   = fits.open(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.momentum_cube.fits'))[0]
    momentum_cube_p = fits.open(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.momentum_cube_p.fits'))[0]
    momentum_cube_m = fits.open(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.momentum_cube_m.fits'))[0]
    momentum   = np.nansum(momentum_cube.data)
    momentum_p = np.nansum(momentum_cube_p.data)
    momentum_m = np.nansum(momentum_cube_m.data)
    return [momentum,momentum_p,momentum_m]


###################################################################################################

# get integrated outflow rate without defined region
####################################################

def integrate_momentum_region(dataset, mask_file, ax_type='major', SNR=3.0):

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
    momentum_cube = fits.open(os.path.join(ratedir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s.momentum_cube.fits'))[0]
    mask      = fits.open(os.path.join(ratedir, randname+'.mask.rotate.fits'))[0]
    masked_rate = copy.deepcopy(momentum_cube.data)
    nv,ny,nx = momentum_cube.data.shape
    for v in np.arange(nv):
        masked_rate[v] *= mask.data

    os.system('rm -rf '+os.path.join(ratedir, randname+'.mask.rotate'))
    os.system('rm -rf '+os.path.join(ratedir, randname+'.mask.rotate.fits'))
    os.system('rm -rf '+os.path.join(ratedir, randname+'.mask.regrid'))

    momentum = np.nansum(masked_rate)
    return momentum


###################################################################################################

# get outflow momentum
####################

# At 100sigma there are only a few very small regions of non-disk emission left in the CO(3-2) cube
# which cannot be rotated. Rotation loses emission by resampling the data leading to an empty image
# in that case. Thus try/except statements are needed to keep the code running.


# generate momentum cubes in parallel
# This is pure python and can be done in parallel
def parallel_momentum_cubes(inp):
    try:
        print("Executing: "+inp[0]['line']+" "+inp[1]+" "+str(inp[2]))
        get_momentum_cube(inp[0], ax_type=inp[1], SNR=inp[2])
    except:
        print("Failed: "+inp[0]['line']+" "+inp[1]+" "+str(inp[2]))

for dataset in datasets:
#    for ax_type in ['major','minor']:
    for ax_type in ['major']:
        plist = [[dataset, ax_type, SNR] for SNR in SNRs]
        pool = Pool(len(SNRs))
        pool.map(parallel_momentum_cubes, plist)
        pool.close()


# sequentially read in momenta
momenta = {}
for dataset in datasets:
    momenta[dataset['line']] = {}
    for ax_type in ['major','minor']:
        momenta[dataset['line']][ax_type] = {}
        for SNR in SNRs:

            try:
                P, P_p, P_m           = integrate_momentum(dataset, ax_type=ax_type, SNR=SNR)
                P_real_outflow        = integrate_momentum_region(dataset, os.path.join(ratedir, 'nondisk.'+dataset['line']+'.real_outflow_new.mask'), ax_type=ax_type, SNR=SNR)
                P_real_superbubble    = integrate_momentum_region(dataset, os.path.join(ratedir, 'nondisk.'+dataset['line']+'.superbubble_new.mask'), ax_type=ax_type, SNR=SNR)
                P_real_cospatial_disk = integrate_momentum_region(dataset, os.path.join(ratedir, 'nondisk.'+dataset['line']+'.cospatial_disk_new.mask'), ax_type=ax_type, SNR=SNR)
            except:
                P, P_p, P_m           = np.nan, np.nan, np.nan
                P_real_outflow        = np.nan
                P_real_superbubble    = np.nan
                P_real_cospatial_disk = np.nan

            momenta[dataset['line']][ax_type][SNR] = [P,P_p,P_m,P_real_outflow,P_real_superbubble,P_real_cospatial_disk]

fnpickle(momenta, 'momenta.pickle')


###################################################################################################

# deprojected momenta
######################

dep_momenta = {}

for dataset in datasets:
    dep_momenta[dataset['line']] = {}
    for ax_type in ['major','minor']:
        dep_momenta[dataset['line']][ax_type] = {}
        for SNR in SNRs:
            dep_momenta[dataset['line']][ax_type][SNR] = {}

            P, P_p, P_m, P_outflow, P_bubble, P_codisk = momenta[dataset['line']][ax_type][SNR]

            dep_P         = [ P/(np.sin(j*u.degree).value)**2 for j in np.arange(48,108) ]
            dep_P_outflow = [ P_outflow/(np.sin(j*u.degree).value)**2 for j in np.arange(48,108) ]
            dep_P_bubble  = [ P_bubble/(np.sin(j*u.degree).value)**2 for j in np.arange(48,108) ]
            dep_P_codisk  = [ P_codisk/(np.sin(j*u.degree).value)**2 for j in np.arange(48,108) ]

            dep_momenta[dataset['line']][ax_type][SNR]['total']   = np.nanpercentile(dep_P, [16,50,81])
            dep_momenta[dataset['line']][ax_type][SNR]['outflow'] = np.nanpercentile(dep_P_outflow, [16,50,81])
            dep_momenta[dataset['line']][ax_type][SNR]['bubble']  = np.nanpercentile(dep_P_bubble, [16,50,81])
            dep_momenta[dataset['line']][ax_type][SNR]['codisk']  = np.nanpercentile(dep_P_codisk, [16,50,81])

fnpickle(dep_momenta, 'momenta_deprojected.pickle')

###################################################################################################
