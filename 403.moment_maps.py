########################
# CASA ANALYSIS SCRIPT #
########################

# Make moment maps using the mask.

###################################################################################################

# import required modules
execfile('scripts/casa_imports.py')
execfile('NGC253/project_info.py')


###################################################################################################

# load sub project info
execfile(os.path.join(projectdir, 'info.py'))


###################################################################################################

# make moment maps
##################

os.system('mkdir -p '+momdir)

for dataset in datasets:
    for SNR in SNRs:
        masked_cube = os.path.join(datdir, ''+dataset['cube']+'.mask_'+str(SNR)+'s')

        for mom in [0,1,2,4]:
            mom_map = os.path.join(momdir, dataset['cube']+'.mask_'+str(SNR)+'s.mom'+str(mom))
            os.system('rm -rf '+mom_map)
            os.system('rm -rf '+mom_map+'.fits')
            immoments(imagename = masked_cube,
                outfile    = mom_map,
                moments    = [mom],
                axis       = 'spectral',
                region     = '',
                chans      = '',
                includepix = [0,1000]       # is needed for moment 4 despite the map being positive already
                )
            exportfits(imagename = mom_map,
                fitsimage    = mom_map+'.fits',
                dropstokes   = True,
                dropdeg      = True,
                overwrite    = True
                )


###################################################################################################
