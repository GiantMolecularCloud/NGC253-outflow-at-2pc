########################
# CASA ANALYSIS SCRIPT #
########################

# Get pV diagram along the previously defined slits.

###################################################################################################

# import required modules
execfile('scripts/casa_imports.py')
execfile('NGC253/project_info.py')


###################################################################################################

# load sub project info
execfile(os.path.join(projectdir, 'info.py'))

major_slices = fnunpickle('slices_major.pickle')
minor_slices = fnunpickle('slices_minor.pickle')


###################################################################################################

# use the cube masked at 3.0 sigma
SNR = 3.0

os.system('mkdir -p '+pvdir)


# regrid cubes
##############

for dataset in datasets:
    masked_cube = os.path.join(datdir, dataset['cube']+'.mask_'+str(SNR)+'s')
    regrid_cube = os.path.join(datdir, dataset['cube']+'.mask_'+str(SNR)+'s.regrid')
    imregrid(imagename = masked_cube,
        template   = os.path.join(datdir, 'template_cube.image'),
        output     = regrid_cube,
        asvelocity = True,
        overwrite  = True
        )


# get pV diagram
################

for dataset in datasets:
    regrid_cube = os.path.join(datdir, dataset['cube']+'.mask_'+str(SNR)+'s.regrid')
    for a_slice in major_slices:
        pV_file = os.path.join(pvdir, dataset['cube']+'.pV_major.slice_major_'+str(a_slice['slice']))

        os.system('rm -rf '+pV_file)
        os.system('rm -rf '+pV_file+'.fits')
        impv(imagename = regrid_cube,
            outfile    = pV_file,
            overwrite  = True,
            mode       = 'length',
            center     = str(a_slice['center'].to_string('hmsdms')).split(' '),
            length     = "{0.value}{0.unit}".format(major_length),
            pa         = "{0.value}{0.unit}".format(a_slice['PA']),
            width      = "{0.value}{0.unit}".format(a_slice['slicewidth']),
            unit       = 'arcsec',
            chans      = '',
            mask       = ''
            )
        exportfits(imagename = pV_file,
            fitsimage    = pV_file+'.fits',
            velocity     = True,
            dropstokes   = True,
            dropdeg      = True,
            overwrite    = True
            )


###################################################################################################

# stitch pVs
############

# FAILS IN CASA FOR UNKNOWN REASON BUT WORKS IN PYTHON3

concatlist = []
empty  = np.full_like(fits.open(os.path.join(pvdir, 'NGC253.CO_3-2.pV_major.slice_major_10.fits'))[0].data, np.nan)
header = fits.open(os.path.join(pvdir, 'NGC253.CO_3-2.pV_major.slice_major_10.fits'))[0].header

for a_slice in major_slices:
    for dataset in datasets:
        pV_file = os.path.join(pvdir, dataset['cube']+'.pV_major.slice_major_'+str(a_slice['slice'])+'.fits')
        if os.path.exists(pV_file):
            concatlist.append(fits.open(pV_file)[0].data)
        else:
            concatlist.append(empty)

header['ctype3'] = 'SLICE'
header['cunit3'] = 'NONE'
header['cdelt3'] = 'NONE'
header['crval3'] = 'NONE'
header['crpix3'] = 'NONE'
fits.writeto('NGC253.CO.pV_major.pVs.fits', data=np.array(concatlist), header=header, overwrite=True)


###################################################################################################
