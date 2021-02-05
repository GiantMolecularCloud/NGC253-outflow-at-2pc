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

os.system('mkdir -p '+dnddir)


###################################################################################################

# separate disk and other
#########################

def separate_disk_other(inps):
    dataset = inps[0]
    SNR     = inps[1]
    ax_type = inps[2]
    cube_file = os.path.join(datdir, dataset['cube']+'.mask_'+str(SNR)+'s')
    mask_file = os.path.join(sepdir, dataset['cube']+'.ppV_mask_'+ax_type)
    disk_file = os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.disk.'+str(SNR)+'s')
    nondisk_file = os.path.join(dnddir, dataset['cube']+'.ppV_mask_'+ax_type+'.non-disk.'+str(SNR)+'s')

    if os.path.exists(mask_file):

        # disk
        ######

        casalog.post("*"*80)
        casalog.post("getting disk for "+cube_file)
        casalog.post("*"*80)

        os.system('rm -rf '+disk_file+'.image')
        immath(imagename = [cube_file, mask_file],
            mode    = 'evalexpr',
            expr    = 'IM0*IM1',
            outfile = disk_file+'.image'
            )
        exportfits(imagename = disk_file+'.image',
            fitsimage = disk_file+'.fits',
            dropdeg   = True,
            velocity  = True,
            overwrite = True
            )
        # mask 0 values
        cube = fits.open(disk_file+'.fits')[0]
        cube.data[cube.data==0.0] = np.nan
        fits.writeto(disk_file+'.fits', data=cube.data, header=cube.header, overwrite=True)


        # other (outflows)
        ##################

        # CASA does not accept the perfectly well defined mask even after regridding to exactly the image
        # use astropy instead because it just works, unlike casa

        casalog.post("*"*80)
        casalog.post("getting non-disk for "+cube_file)
        casalog.post("*"*80)

        # load data
        cube           = fits.open(cube_file+'.fits')[0]
        mask           = fits.open(mask_file+'.fits')[0]
        disk_outline   = fits.open(os.path.join(sepdir, dataset['cube']+'.disk_outline.mask.fits'))[0]
        overall_region = fits.open(os.path.join(sepdir, dataset['cube']+'.overall_region.mask.fits'))[0]

        # replace nan with 0 for multiplication
        disk_outline.data[np.isnan(disk_outline.data)] = 0
        overall_region.data[np.isnan(overall_region.data)] = 0

        # mask data
        inverted_mask = (mask.data-1.)*-1.
        masked_data   = cube.data*inverted_mask
        masked_data   = np.array([masked_data[i]*disk_outline.data*overall_region.data for i in np.arange(len(masked_data))])

        masked_data[masked_data == 0.0] = np.nan

        fits.writeto(nondisk_file+'.fits', data=masked_data, header=cube.header, overwrite=True)

        importfits(fitsimage = nondisk_file+'.fits',
            imagename = nondisk_file+'.image',
            overwrite = True)


###################################################################################################

# separate disk from non-disk
#############################

# do not consider emission outside the sliced up area
sidecenter1 = SkyCoord(center.ra-minor_length/2.*np.cos(disk_PA)/np.cos(center.dec), center.dec+minor_length/2.*np.sin(disk_PA), frame='icrs')
sidecenter2 = SkyCoord(center.ra+minor_length/2.*np.cos(disk_PA)/np.cos(center.dec), center.dec-minor_length/2.*np.sin(disk_PA), frame='icrs')
edge1 = SkyCoord(sidecenter1.ra+major_length/2.*np.sin(disk_PA)/np.cos(sidecenter1.dec),sidecenter1.dec+major_length/2.*np.cos(disk_PA), frame='icrs')
edge2 = SkyCoord(sidecenter1.ra-major_length/2.*np.sin(disk_PA)/np.cos(sidecenter1.dec),sidecenter1.dec-major_length/2.*np.cos(disk_PA), frame='icrs')
edge3 = SkyCoord(sidecenter2.ra-major_length/2.*np.sin(disk_PA)/np.cos(sidecenter2.dec),sidecenter2.dec-major_length/2.*np.cos(disk_PA), frame='icrs')
edge4 = SkyCoord(sidecenter2.ra+major_length/2.*np.sin(disk_PA)/np.cos(sidecenter2.dec),sidecenter2.dec+major_length/2.*np.cos(disk_PA), frame='icrs')

with open(os.path.join(sepdir, 'overall_region.region'), 'w') as reg_file:
    reg_file.write("#CRTFv0 CASA Region Text Format version 0\n")
    reg_file.write("poly [["+str(edge1.ra)+","+str(edge1.dec)+"],["+str(edge2.ra)+","+str(edge2.dec)+"],["+str(edge3.ra)+","+str(edge3.dec)+"],["+str(edge4.ra)+","+str(edge4.dec)+"]] coord=J2000")

for dataset in datasets:
    os.system('cp '+os.path.join(sepdir, 'overall_region.region')+' '+os.path.join(sepdir, dataset['cube']+'.overall_region.region'))
    region_to_mask(region=os.path.join(sepdir, dataset['cube']+'.overall_region.region'), template=os.path.join(momdir, dataset['cube']+'.mask_3.0s.mom0'), overwrite=True)
    exportfits(imagename = os.path.join(sepdir, dataset['cube']+'.overall_region.mask'),
            fitsimage = os.path.join(sepdir, dataset['cube']+'.overall_region.mask.fits'),
            dropdeg   = True,
            velocity  = True,
            overwrite = True
            )
    os.system('cp '+os.path.join(sepdir, 'disk_outline.region')+' '+os.path.join(sepdir, dataset['cube']+'.disk_outline.region'))
    region_to_mask(region=os.path.join(sepdir, dataset['cube']+'.disk_outline.region'), template=os.path.join(momdir, dataset['cube']+'.mask_3.0s.mom0'), overwrite=True)
    exportfits(imagename = os.path.join(sepdir, dataset['cube']+'.disk_outline.mask'),
            fitsimage = os.path.join(sepdir, dataset['cube']+'.disk_outline.mask.fits'),
            dropdeg   = True,
            velocity  = True,
            overwrite = True
            )


# major axis first
inps = []
for dataset in datasets:
    for SNR in SNRs:
        inps.append([dataset, SNR, 'major'])
pool = Pool(10)                             # cannot run all in parallel because they access the same files
pool.map(separate_disk_other, inps)
pool.close()

# minor axis when major is finished
inps = []
for dataset in datasets:
    for SNR in SNRs:
        inps.append([dataset, SNR, 'minor'])
pool = Pool(10)                             # cannot run all in parallel because they access the same files
pool.map(separate_disk_other, inps)
pool.close()


###################################################################################################
