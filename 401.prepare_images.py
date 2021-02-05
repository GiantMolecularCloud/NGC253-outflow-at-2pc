########################
# CASA ANALYSIS SCRIPT #
########################

# Reframe and regrid the original datasets to a common structure.

###################################################################################################

# import required modules
execfile('scripts/casa_imports.py')
execfile('NGC253/project_info.py')


###################################################################################################

# load sub project info
execfile(os.path.join(projectdir, 'info.py'))


###################################################################################################

# Import
########

# Import the fits files into casa image format
os.system('mkdir -p '+datdir)
os.system('rm -rf '+os.path.join(datdir, 'temp.CO_1-0.import.image'))
os.system('rm -rf '+os.path.join(datdir, 'temp.CO_2-1.import.image'))
os.system('rm -rf '+os.path.join(datdir, 'temp.CO_3-2.import.image'))
importfits(fitsimage = os.path.join(rawdir, 'NGC253.CO_1-0.fits'),
    imagename        = os.path.join(datdir, 'temp.CO_1-0.import.image'),
    defaultaxes      = True,
    defaultaxesvalues = ['','','','I']
    )
importfits(fitsimage = os.path.join(rawdir, 'NGC253.CO_2-1.fits'),
    imagename        = os.path.join(datdir, 'temp.CO_2-1.import.image'),
    defaultaxes      = True,
    defaultaxesvalues = ['','','','I']
    )
importfits(fitsimage = os.path.join(rawdir, 'NGC253.CO_3-2.fits'),
    imagename        = os.path.join(datdir, 'temp.CO_3-2.import.image'),
    defaultaxes      = True,
    defaultaxesvalues = ['','','','I']
    )


###################################################################################################

# Reframe
#########

# Match the spectral coordinate systems
os.system('rm -rf '+os.path.join(datdir, 'temp.CO_1-0.reframe.image'))
os.system('rm -rf '+os.path.join(datdir, 'temp.CO_2-1.reframe.image'))
os.system('rm -rf '+os.path.join(datdir, 'temp.CO_3-2.reframe.image'))
imreframe(imagename = os.path.join(datdir, 'temp.CO_1-0.import.image'),
    output          = os.path.join(datdir, 'temp.CO_1-0.reframe.image'),
    outframe        = 'LSRK',
    restfreq        = '115.27120180GHz'
    )
imreframe(imagename = os.path.join(datdir, 'temp.CO_2-1.import.image'),
    output          = os.path.join(datdir, 'temp.CO_2-1.reframe.image'),
    outframe        = 'LSRK',
    restfreq        = '230.53800000GHz'
    )
imreframe(imagename = os.path.join(datdir, 'temp.CO_3-2.import.image'),
    output          = os.path.join(datdir, 'temp.CO_3-2.reframe.image'),
    outframe        = 'LSRK',
    restfreq        = '345.79598990GHz'
    )

###################################################################################################

# To Kelvin
###########

# convert to Kelvin if not already
for line in ['CO_1-0', 'CO_2-1', 'CO_3-2']:
    bunit = imhead(imagename=os.path.join(datdir, 'temp.'+line+'.reframe.image'), mode='get', hdkey='bunit')
    os.system('rm -rf '+os.path.join(datdir, 'temp.'+line+'.K.image'))
    if (bunit == 'Jy/beam'):
        restfreq = ((imhead(imagename=os.path.join(datdir, 'temp.'+line+'.reframe.image'), mode='get', hdkey='restfreq')['value'])*u.Hz).to(u.GHz)
        bmin     = (imhead(imagename=os.path.join(datdir, 'temp.'+line+'.reframe.image'), mode='get', hdkey='bmin')['value'])*u.arcsec
        bmaj     = (imhead(imagename=os.path.join(datdir, 'temp.'+line+'.reframe.image'), mode='get', hdkey='bmaj')['value'])*u.arcsec
        immath(imagename = os.path.join(datdir, 'temp.'+line+'.reframe.image'),
            outfile = os.path.join(datdir, 'temp.'+line+'.K.image'),
            mode    = 'evalexpr',
        	expr    = '1.226e6*IM0/('+str(restfreq.value)+'^2*'+str(bmin.value)+'*'+str(bmaj.value)+')'
            )
        imhead(imagename = os.path.join(datdir, 'temp.'+line+'.K.image'),
            mode   = 'put',
            hdkey  = 'bunit',
            hdvalue = 'K'
            )
    elif (bunit == 'K'):
        print(line+" is Kelvin already")
        os.system('cp -r '+os.path.join(datdir, 'temp.'+line+'.reframe.image')+' '+os.path.join(datdir, 'temp.'+line+'.K.image'))
    else:
        raise Exception("Unknown bunit in file "+os.path.join(datdir, 'temp.'+line+'.reframe.image'))



###################################################################################################

# Fix conversions
#################

# some conversions are written into the header instead of being applied to the data
# apply these changes by exporting to fits and reimporting

for line in ['CO_1-0', 'CO_2-1', 'CO_3-2']:
    exportfits(imagename = os.path.join(datdir, 'temp.'+line+'.K.image'),
        fitsimage  = os.path.join(datdir, 'temp.'+line+'.K.fits'),
        velocity   = True,
        optical    = True,
        overwrite  = True,
        dropstokes = False,
        stokeslast = True,
        history    = True,
        dropdeg    = False
        )
    importfits(fitsimage = os.path.join(datdir, 'temp.'+line+'.K.fits'),
        imagename = os.path.join(datdir, 'NGC253.'+line),
        overwrite = True
        )
    exportfits(imagename = os.path.join(datdir, 'NGC253.'+line),
        fitsimage        = os.path.join(datdir, 'NGC253.'+line+'.fits'),
        overwrite = True,
        dropdeg   = True
        )

print("\n\nUpdate the noise values in the line definition file ('outflow_comparison_info.py')!\n\n")


###################################################################################################

# create empty cube
###################

for line in ['CO_1-0', 'CO_2-1', 'CO_3-2']:
    os.system('rm -rf '+os.path.join(datdir, 'NGC253.'+line+'.empty.image'))
    immath(imagename = os.path.join(datdir, 'NGC253.'+line+'.fits'),
        outfile      = os.path.join(datdir, 'NGC253.'+line+'.empty.image'),
        mode         = 'evalexpr',
        expr         = 'IM0*0.0'
        )
    exportfits(imagename = os.path.join(datdir, 'NGC253.'+line+'.empty.image'),
        fitsimage        = os.path.join(datdir, 'NGC253.'+line+'.empty.fits'),
        velocity         = True,
        dropdeg          = True,
        dropstokes       = True,
        overwrite        = True
        )


###################################################################################################

# clean up
##########

if not (debug == True):
    os.system('rm -rf '+os.path.join(datdir, 'temp.CO_1-0.import.image'))
    os.system('rm -rf '+os.path.join(datdir, 'temp.CO_2-1.import.image'))
    os.system('rm -rf '+os.path.join(datdir, 'temp.CO_3-2.import.image'))
    os.system('rm -rf '+os.path.join(datdir, 'temp.CO_1-0.reframe.image'))
    os.system('rm -rf '+os.path.join(datdir, 'temp.CO_2-1.reframe.image'))
    os.system('rm -rf '+os.path.join(datdir, 'temp.CO_3-2.reframe.image'))
    os.system('rm -rf '+os.path.join(datdir, 'temp.CO_1-0.K.image'))
    os.system('rm -rf '+os.path.join(datdir, 'temp.CO_2-1.K.image'))
    os.system('rm -rf '+os.path.join(datdir, 'temp.CO_3-2.K.image'))
    os.system('rm -rf '+os.path.join(datdir, 'temp.CO_1-0.K.fits'))
    os.system('rm -rf '+os.path.join(datdir, 'temp.CO_2-1.K.fits'))
    os.system('rm -rf '+os.path.join(datdir, 'temp.CO_3-2.K.fits'))


###################################################################################################
