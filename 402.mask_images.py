########################
# CASA ANALYSIS SCRIPT #
########################

# Calculate masks at certain SNR thresholds.

###################################################################################################

# import required modules
execfile('scripts/casa_imports.py')
execfile('NGC253/project_info.py')


###################################################################################################

# load sub project info
execfile(os.path.join(projectdir, 'info.py'))


###################################################################################################

# SNR masking
#############

for dataset in datasets:
    for SNR in SNRs:
        cube        = os.path.join(datdir, ''+dataset['cube'])
        masked_cube = os.path.join(datdir, ''+dataset['cube']+'.mask_'+str(SNR)+'s')
        rms         = dataset['rms'].value

        os.system('cp -r '+cube+' '+masked_cube)
        ia.open(masked_cube)
        ia.calcmask(mask='"'+masked_cube+'" > '+str(SNR)+'*'+str(rms), name=str(SNR)+'sigma')
        ia.done()

        exportfits(imagename = masked_cube,
            fitsimage    = masked_cube+'.fits',
            velocity     = True,
            dropstokes   = True,
            dropdeg      = True,
            overwrite    = True
            )

        # export mask to mask image
        makemask(mode = 'copy',
            inpimage  = masked_cube,
            inpmask   = masked_cube+':'+str(SNR)+'sigma',
            output    = masked_cube+'.mask',
            overwrite = True
            )


###################################################################################################
