########################
# CASA ANALYSIS SCRIPT #
########################

# Create ppv mask files from the pV masks.

###################################################################################################

# import required modules
execfile('scripts/casa_imports.py')
execfile('NGC253/project_info.py')


###################################################################################################

# load sub project info
execfile(os.path.join(projectdir, 'info.py'))


###################################################################################################

# Regrid mask
#############

os.system('mkdir -p '+sepdir)

# smooth the mask to cover single pixels that are missed by the pV-to-ppV conversion
for dataset in datasets:
    print("\nRegrid masks ("+dataset['line']+")\n")

    for mtype in ['major','minor']:
        imregrid(imagename = os.path.join(sepdir, 'NGC253.CO.ppV_mask_'+mtype+'.mask_model_variable'),
            template = os.path.join(datdir, dataset['cube']),
            output   = os.path.join(sepdir, 'temp.'+dataset['cube']+'.ppV_mask_'+mtype),
            asvelocity = True,
            axes       = [-1],
            overwrite  = True
            )
        os.system('rm -rf '+os.path.join(sepdir, dataset['cube']+'.ppV_mask_'+mtype))
        immath(imagename = os.path.join(sepdir, 'temp.'+dataset['cube']+'.ppV_mask_'+mtype),
            outfile = os.path.join(sepdir, dataset['cube']+'.ppV_mask_'+mtype),
            mode    = 'evalexpr',
            expr    = 'iif(IM0>0.1,1,0)'
            )
        exportfits(imagename = os.path.join(sepdir, dataset['cube']+'.ppV_mask_'+mtype),
            fitsimage    = os.path.join(sepdir, dataset['cube']+'.ppV_mask_'+mtype+'.fits'),
            dropstokes   = True,
            dropdeg      = True,
            overwrite    = True
            )


###################################################################################################
