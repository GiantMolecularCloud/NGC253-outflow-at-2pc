########################
# CASA ANALYSIS SCRIPT #
########################

# Calculate wider or more narrow masks to quantify the error associated with the pV mask.

###################################################################################################

# import required modules
execfile('NGC253/project_info.py')


###################################################################################################

# load sub project info
execfile(os.path.join(projectdir, 'info.py'))
major_slices = fnunpickle(os.path.join(subprojectdir,'pickle_jar','slices_major.py3.pickle'))

tmodel_major = os.path.join(varmdir,'diskfit.total_model.major_slice_')

def gauss(xvalues, amp, sigma, shift, const):
    return [(amp*np.exp(-1.*((x+shift)/sigma)**2))+const for x in xvalues]


slice_list_major = [s['slice'] for s in major_slices]
mask_width_major = gauss(xvalues = slice_list_major,
                         amp   = 30.,
                         sigma = 2.5,
                         shift = -1.*np.mean(slice_list_major),
                         const = 25.
                         )

# generate masks that are -50% to +50% relative to the actual mask
mask_factor = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]
mask_widths = [list(np.array(mask_width_major)*x) for x in mask_factor]


###################################################################################################
# variable mask model velocity
###################################################################################################

os.system('mkdir -p '+varmdir)

def pixels_to_offsets(fitsimage, pixels):
    """
    Convert list of pixels to offset values using the header information of a given image.
    """
    header = fits.open(fitsimage)[0].header
    crval = header['CRVAL1']
    crpix = header['CRPIX1']
    cdelt = header['CDELT1']
    return [(i-crpix)*cdelt+crval if not np.isnan(i) else np.nan for i in pixels]

def offsets_to_pixels(fitsimage, offsets):
    """
    Convert list of offset values to pixels using the header information of a given image.
    """
    header = fits.open(fitsimage)[0].header
    crval = header['CRVAL1']
    crpix = header['CRPIX1']
    cdelt = header['CDELT1']
    return [int(np.around((i-crval)/cdelt+crpix)) if not np.isnan(i) else np.nan for i in offsets]


def velocities_to_pixels(fitsimage, velocities):
    """
    Convert list of velocity values to pixels using the header information of a given image.
    """
    header = fits.open(fitsimage)[0].header
    crval = header['CRVAL2']
    crpix = header['CRPIX2']
    cdelt = header['CDELT2']
    return [int(np.around((i-crval)/cdelt+crpix)) if not np.isnan(i) else np.nan for i in velocities]

def grow_median():
    """
    grow the median velocity by mask_radius pixels in each direction.
    NUMPY DOES BANKING ROUNDING WHICH OMITS EVERY OTHER NUMBER FOR LISTS OF n.5 AND THUS BREAKS THE CODE.
    A SIMPLE WORKAROUND IS np.floor(x+0.5) TO GET TRUE ROUNDING FOR VALUES EXACLTY IN BETWEEN TWO
    WHOLE NUMBERS.
    """
    pV_togrow_data   = fits.open('temp2.fits')[0].data
    pV_grown_header  = fits.open('temp2.fits')[0].header
    pV_grown_data    = copy.deepcopy(pV_togrow_data)
    for x in np.arange(len(pV_togrow_data)):                                    # loop over x axis (velocity)
        for y in np.arange(len(pV_togrow_data[0])):                             # loop over y axis (offset)
            if (pV_togrow_data[x,y] == 1.0):                                    # grow mask around median velocity
                #print("pixel ["+str(x)+','+str(y)+"] is unity.")               # debug output
                for i in [int(np.floor(ii+0.5)) for ii in np.arange(x-mask_radius,x+mask_radius,1)]:              # loop over +- mask_radius pixels in x
                    for j in [int(np.floor(jj+0.5)) for jj in np.arange(y-mask_radius,y+mask_radius,1)]:          # loop over +- mask_radius pixels in y
                        if (i >= 0) and (i < len(pV_togrow_data)):              # do not go beyond image edges
                            if (j >= 0) and (j < len(pV_togrow_data[0])):       # do not go beyond image edges
                                distance = np.sqrt( np.abs(x-i)**2 + np.abs(y-j)**2 )
                                if (distance < mask_radius):                             # grow mask within a radius of 25 pixels
                                    pV_grown_data[i,j] = 1.0
    return pV_grown_data, pV_grown_header

# start loop
for idx,mask_width in enumerate(mask_widths):
    mfac = mask_factor[idx]
    print("Mask factor "+str(mfac))
    for a_slice in tqdm(major_slices):

        modelV = tmodel_major+str(a_slice['slice'])+'.fits'
        pV_file = os.path.join(projectdir,'410.disk_non-disk','04.pVs','NGC253.CO_1-0.pV_minor.slice_minor_'+str(a_slice['slice'])+'.fits')      # Need to use old files to get correct header
        mask_radius = mask_width[a_slice['slice']]                                                                                               # and later regrid the ppV mask to the actual images
                                                                                                                                                 # like I did in all other 420 scripts.
        if os.path.exists(modelV) and os.path.exists(pV_file):

            os.system('rm -rf temp.fits')
            os.system('rm -rf temp2.fits')

            # set the pixels at the median velocity to unity
            modelV_data    = fits.open(modelV)[0].data[0]
            offsets        = pixels_to_offsets(modelV, np.arange(1,len(modelV_data)+1))
            offsets_pix    = offsets_to_pixels(modelV, offsets)
            velocities_pix = velocities_to_pixels(pV_file, [i*1e3 for i in modelV_data])     # needs correction for 1e3 wrong velocity infos

            # copy the empty pV for each median pV_file
            os.system('cp -r '+os.path.join(varmdir,'pV_empty.fits temp.fits'))

            # coordinates are now offsets_pix and velocities_pix
            pV_tofill_data   = fits.open('temp.fits')[0].data
            pV_tofill_header = fits.open(pV_file)[0].header
            for p,V in zip(offsets_pix, velocities_pix):
                if not (np.isnan(V)):
                    pV_tofill_data[V-1,p-1] = 1.0
            fits.writeto('temp2.fits', data=pV_tofill_data, header=pV_tofill_header)

            pV_grown_data,pV_grown_header = grow_median()
            fits.writeto(os.path.join(varmdir,'NGC253.CO.pV_major.slice_major_'+str(a_slice['slice'])+'.mask_factor_'+str(mfac)+'.fits'), data=pV_grown_data, header=pV_grown_header, overwrite=True)

            os.system('rm -rf temp.fits')
            os.system('rm -rf temp2.fits')


###################################################################################################
# plot alternativ masks
###################################################################################################

paperdir = os.path.join(projectdir,'paper_K18')
plotdir  = os.path.join('plots/NGC253/420.disk_non-disk/')

easy_aplpy.settings.velo_fontsize           = 20      # unit: point
easy_aplpy.settings.colorbar_label_fontsize = 20      # unit: point
easy_aplpy.settings.colorbar_ticks_fontsize = 20      # unit: point
easy_aplpy.settings.grid_label_fontsize     = 20      # unit: point
easy_aplpy.settings.colorbar_width          = 0.08    # relative to panel size
easy_aplpy.settings.scalebar_frame          = False
easy_aplpy.settings.scalebar_linestyle      = 'solid' # or any other plt.plot linestyle
easy_aplpy.settings.scalebar_linewidth      = 4       # unit: points
easy_aplpy.settings.scalebar_color          = 'red'   # any named color or mpl.color instance
easy_aplpy.settings.ticks_color             = 'black' # this setting overrules the matplotlibrc defaults
easy_aplpy.settings.frame_color             = 'black'
easy_aplpy.settings.tick_label_fontsize     = 20      # unit: point
easy_aplpy.settings.axis_label_fontsize     = 20      # unit: point
easy_aplpy.settings.props                   = {'boxstyle': "round", 'facecolor': "w", 'edgecolor': "black", 'linewidth': 1.0, 'alpha': 1.0}
easy_aplpy.settings.ticks_xspacing          = Angle('0 0 1.0', unit='hourangle')
easy_aplpy.settings.ticks_yspacing          = 10.0*u.arcsec
easy_aplpy.settings.ticks_minor_frequency   = 5

# plot alternativ masks in the style of appendix figure B
for idx,mask_width in enumerate(mask_widths):
    mfac = mask_factor[idx]

    contours = []
    labels   = []
    lines    = []
    for a_slice in major_slices:
        for CO in ['CO_1-0','CO_2-1','CO_3-2']:
            c = []
            pV       = os.path.join(pvdir, 'NGC253.'+CO+'.pV_major.slice_major_'+str(a_slice['slice'])+'.fits')
            var_mask = os.path.join(varmdir,'NGC253.CO.pV_major.slice_major_'+str(a_slice['slice'])+'.mask_factor_'+str(mfac)+'.fits')
            model    = tmodel_major+str(a_slice['slice'])+'.fits'

            if os.path.exists(pV):
                c.append([pV, [1,4,16], 'black'])
            if os.path.exists(var_mask):
                c.append([var_mask, [0.5,1], 'gold', {'filled': True, 'alpha': 0.25, 'smooth': 1}])
                c.append([var_mask, [0.5], 'gold', {'filled': False, 'alpha': 0.4, 'smooth': 1}])
            contours.append(c)

            if os.path.exists(model):
                model_data = fits.open(model)[0].data[0]
                lines.append([[[i*u.arcsec for i in np.arange(-25.,25.1,0.1)], [(j*u.km/u.s).to(u.m/u.s) for j in model_data], {'color': 'red'}]])

        labels.append([[[0.025, 0.95], str((a_slice['slice']-int(len(major_slices)/2.))*(a_slice['slicewidth']/2.).to(u.arcsec).value)+"$^{\prime\prime}$", {'size': 20, 'bbox': easy_aplpy.settings.props, 'ha': 'left', 'va': 'top'}]])
        labels.append([[]])
        labels.append([[]])
    labels[0] = [[[0.5, 0.95], 'CO (1-0) '+str(mfac*100)+'\%', {'size': 24, 'bbox': easy_aplpy.settings.props, 'ha': 'center', 'va': 'top'}]]
    labels[1] = [[[0.5, 0.95], 'CO (2-1) '+str(mfac*100)+'\%', {'size': 24, 'bbox': easy_aplpy.settings.props, 'ha': 'center', 'va': 'top'}]]
    labels[2] = [[[0.5, 0.95], 'CO (3-2) '+str(mfac*100)+'\%', {'size': 24, 'bbox': easy_aplpy.settings.props, 'ha': 'center', 'va': 'top'}]]

    easy_aplpy.plot.grid(os.path.join(paperdir,'NGC253.CO.pV_major.pVs.fits'),
        [3,9],
        np.arange(0,27),
        figsize  = (15,18),
        cmap     = 'binary',
        colorbar = ['right', 'K'],
        vmin     = 0.1,
        vmax     = 30,
        stretch  = 'log',
        recenter = [0*u.arcsec, 260*u.km/u.s, minor_length, 500*u.km/u.s],
        labels   = ['offset [$^{\prime\prime}$]','velocity [km\,s$^{-1}$]'],
        channel_label = None,
        texts    = labels[0:27],
        contours = contours[0:27],
        lines    = lines[0:27],
        out      = os.path.join(plotdir,'mask_factor_'+str(mfac)+'_a.pdf')
        )
    easy_aplpy.plot.grid(os.path.join(paperdir,'NGC253.CO.pV_major.pVs.fits'),
        [3,10],
        np.arange(27,57),
        figsize  = (15,20),
        cmap     = 'binary',
        colorbar = ['right', 'K'],
        vmin     = 0.1,
        vmax     = 30,
        stretch  = 'log',
        recenter = [0*u.arcsec, 260*u.km/u.s, minor_length, 500*u.km/u.s],
        labels   = ['offset ["]','velocity [km\,s$^{-1}$]'],
        channel_label = None,
        texts    = labels[27:57],
        contours = contours[27:57],
        lines    = lines[27:57],
        out      = os.path.join(plotdir,'mask_factor_'+str(mfac)+'_b.pdf')
        )


###################################################################################################
# plot alternativ masks comparison
###################################################################################################

a_slice = major_slices[9]

# generate new pV file and mask file with the correct order
shape     = fits.open(os.path.join(pvdir, 'NGC253.CO_1-0.pV_major.slice_major_9.fits'))[0].data.shape
concat_pv = np.zeros((3*len(mask_widths),shape[0],shape[1]))
concat_hd = fits.getheader(os.path.join(pvdir, 'NGC253.CO_1-0.pV_major.slice_major_9.fits'))
concat_hd['CTYPE3'] = 'slice'
concat_hd['CUNIT3']  = None
concat_hd['CDELT3'] = None
concat_hd['CRVAL3'] = None
concat_hd['CRPIX3'] = None
for idx,mask_width in enumerate(mask_widths):
    for idx2,CO in enumerate(['CO_1-0','CO_2-1','CO_3-2']):
        pv_data = fits.getdata(os.path.join(pvdir, 'NGC253.'+CO+'.pV_major.slice_major_'+str(a_slice['slice'])+'.fits'))
        concat_pv[idx*3+idx2] = pv_data
fits.writeto(os.path.join(varmdir,'NGC253.CO.pV_major.pVs.fits'), data=concat_pv, header=concat_hd, overwrite=True)

# overlays
contours = []
labels   = []
lines    = []
for idx,mask_width in enumerate(mask_widths):
    mfac = mask_factor[idx]
    for CO in ['CO_1-0','CO_2-1','CO_3-2']:
        c = []
        pV       = os.path.join(pvdir, 'NGC253.'+CO+'.pV_major.slice_major_'+str(a_slice['slice'])+'.fits')
        var_mask = os.path.join(varmdir,'NGC253.CO.pV_major.slice_major_'+str(a_slice['slice'])+'.mask_factor_'+str(mfac)+'.fits')
        model    = tmodel_major+str(a_slice['slice'])+'.fits'

        if os.path.exists(pV):
            c.append([pV, [1,4,16], 'black'])
        if os.path.exists(var_mask):
            c.append([var_mask, [0.5,1], 'gold', {'filled': True, 'alpha': 0.25, 'smooth': 1}])
            c.append([var_mask, [0.5], 'gold', {'filled': False, 'alpha': 0.4, 'smooth': 1}])
        contours.append(c)

        if os.path.exists(model):
            model_data = fits.open(model)[0].data[0]
            lines.append([[[i*u.arcsec for i in np.arange(-25.,25.1,0.1)], [(j*u.km/u.s).to(u.m/u.s) for j in model_data], {'color': 'red'}]])

    labels.append([[[0.025, 0.95], str(np.round(mfac*100,1))+'\%', {'size': 20, 'bbox': easy_aplpy.settings.props, 'ha': 'left', 'va': 'top'}]])
    labels.append([[]])
    labels.append([[]])
labels[0] += [[[0.5, 0.95], 'CO (1-0) ', {'size': 24, 'bbox': easy_aplpy.settings.props, 'ha': 'center', 'va': 'top'}]]
labels[1] += [[[0.5, 0.95], 'CO (2-1) ', {'size': 24, 'bbox': easy_aplpy.settings.props, 'ha': 'center', 'va': 'top'}]]
labels[2] += [[[0.5, 0.95], 'CO (3-2) ', {'size': 24, 'bbox': easy_aplpy.settings.props, 'ha': 'center', 'va': 'top'}]]

# plot
easy_aplpy.plot.grid(os.path.join(varmdir,'NGC253.CO.pV_major.pVs.fits'),
    [3,11],
    np.arange(0,33),
    figsize  = (15,23),
    cmap     = 'binary',
    colorbar = ['right', 'K'],
    vmin     = 0.1,
    vmax     = 30,
    stretch  = 'log',
    recenter = [0*u.arcsec, 260*u.km/u.s, minor_length, 500*u.km/u.s],
    labels   = ['offset [$^{\prime\prime}$]','velocity [km\,s$^{-1}$]'],
    channel_label = None,
    texts    = labels,
    contours = contours,
    lines    = lines,
    out      = os.path.join(plotdir,'mask_factor_comparison.pdf')
    )


###################################################################################################
# plot alternativ masks comparison
###################################################################################################

paperdir = os.path.join(projectdir,'paper_K18')
plotdir  = os.path.join('plots/NGC253/paper_K+19/v0.16')

easy_aplpy.settings.velo_fontsize           = 16      # unit: point
easy_aplpy.settings.colorbar_label_fontsize = 16      # unit: point
easy_aplpy.settings.colorbar_ticks_fontsize = 16      # unit: point
easy_aplpy.settings.grid_label_fontsize     = 16      # unit: point
easy_aplpy.settings.colorbar_width          = 0.08    # relative to panel size
easy_aplpy.settings.scalebar_frame          = False
easy_aplpy.settings.scalebar_linestyle      = 'solid' # or any other plt.plot linestyle
easy_aplpy.settings.scalebar_linewidth      = 4       # unit: points
easy_aplpy.settings.scalebar_color          = 'red'   # any named color or mpl.color instance
easy_aplpy.settings.ticks_color             = 'black' # this setting overrules the matplotlibrc defaults
easy_aplpy.settings.frame_color             = 'black'
easy_aplpy.settings.tick_label_fontsize     = 20      # unit: point
easy_aplpy.settings.axis_label_fontsize     = 20      # unit: point
easy_aplpy.settings.props                   = {'boxstyle': "round", 'facecolor': "w", 'edgecolor': "black", 'linewidth': 1.0, 'alpha': 1.0}
easy_aplpy.settings.ticks_xspacing          = Angle('0 0 1.0', unit='hourangle')
easy_aplpy.settings.ticks_yspacing          = 10.0*u.arcsec
easy_aplpy.settings.ticks_minor_frequency   = 5

a_slice = major_slices[9]
few_mask_factor = [0.8,0.9,1.0,1.1,1.2]

# generate new pV file and mask file with the correct order
shape     = fits.open(os.path.join(pvdir, 'NGC253.CO_1-0.pV_major.slice_major_9.fits'))[0].data.shape
concat_pv = np.zeros((3*5,shape[0],shape[1]))
concat_hd = fits.getheader(os.path.join(pvdir, 'NGC253.CO_1-0.pV_major.slice_major_9.fits'))
concat_hd['CTYPE3'] = 'slice'
concat_hd['CUNIT3']  = None
concat_hd['CDELT3'] = None
concat_hd['CRVAL3'] = None
concat_hd['CRPIX3'] = None
for idx,mask_width in enumerate(few_mask_factor):
    for idx2,CO in enumerate(['CO_1-0','CO_2-1','CO_3-2']):
        pv_data = fits.getdata(os.path.join(pvdir, 'NGC253.'+CO+'.pV_major.slice_major_'+str(a_slice['slice'])+'.fits'))
        concat_pv[idx*3+idx2] = pv_data
fits.writeto(os.path.join(varmdir,'NGC253.CO.pV_major.pVs2.fits'), data=concat_pv, header=concat_hd, overwrite=True)

# overlays
contours = []
labels   = []
lines    = []
for idx,mask_width in enumerate(few_mask_factor):
    mfac = few_mask_factor[idx]
    for CO in ['CO_1-0','CO_2-1','CO_3-2']:
        c = []
        pV       = os.path.join(pvdir, 'NGC253.'+CO+'.pV_major.slice_major_'+str(a_slice['slice'])+'.fits')
        var_mask = os.path.join(varmdir,'NGC253.CO.pV_major.slice_major_'+str(a_slice['slice'])+'.mask_factor_'+str(mfac)+'.fits')
        model    = tmodel_major+str(a_slice['slice'])+'.fits'

        if os.path.exists(pV):
            c.append([pV, [1,4,16], 'black'])
        if os.path.exists(var_mask):
            c.append([var_mask, [0.5,1], 'gold', {'filled': True, 'alpha': 0.25, 'smooth': 1}])
            c.append([var_mask, [0.5], 'gold', {'filled': False, 'alpha': 0.4, 'smooth': 1}])
        contours.append(c)

        if os.path.exists(model):
            model_data = fits.open(model)[0].data[0]
            lines.append([[[i*u.arcsec for i in np.arange(-25.,25.1,0.1)], [(j*u.km/u.s).to(u.m/u.s) for j in model_data], {'color': 'red'}]])

    labels.append([[[0.025, 0.95], str(int(mfac*100))+'\%', {'size': 20, 'bbox': easy_aplpy.settings.props, 'ha': 'left', 'va': 'top'}]])
    labels.append([[]])
    labels.append([[]])
labels[0] += [[[0.5, 0.95], 'CO (1-0) ', {'size': 20, 'bbox': easy_aplpy.settings.props, 'ha': 'center', 'va': 'top'}]]
labels[1] += [[[0.5, 0.95], 'CO (2-1) ', {'size': 20, 'bbox': easy_aplpy.settings.props, 'ha': 'center', 'va': 'top'}]]
labels[2] += [[[0.5, 0.95], 'CO (3-2) ', {'size': 20, 'bbox': easy_aplpy.settings.props, 'ha': 'center', 'va': 'top'}]]

# plot
easy_aplpy.plot.grid(os.path.join(varmdir,'NGC253.CO.pV_major.pVs2.fits'),
    [3,5],
    np.arange(0,3*5),
    figsize  = (15,10),
    cmap     = 'binary',
    colorbar = ['right', 'K'],
    vmin     = 0.1,
    vmax     = 30,
    stretch  = 'log',
    recenter = [0*u.arcsec, 260*u.km/u.s, minor_length, 500*u.km/u.s],
    labels   = ['offset [$^{\prime\prime}$]','velocity [km\,s$^{-1}$]'],
    channel_label = None,
    texts    = labels,
    contours = contours,
    lines    = lines,
    out      = os.path.join(plotdir,'figB1.pdf')
    )


###################################################################################################
# pV to ppV mask
###################################################################################################

def pv_slice_point_to_ppv_coord(pv_points,a_slice,pix_scale):
    """
    Convert a point in pV space to a list of ppV (RA,DEC, V) points using
    the information about the orientation and size of the slice.

    Parameters
    ----------
    pv_points : list
        List of pV points as astropy.units objects.
    a_slice : dictionary
        Dictionary describing the slice. in particular: center, width and PA.
    pix_scale : astropy.units angular
        Scale on which the ppV line that corresponds to a pV point will be sampled.

    Returns
    -------
    list
        List of ppV points as astropy.SkyCoord and astropy.units objects.
    """
    ppv_points = []
    for pv_point in pv_points:
        center     = a_slice['center']
        slicewidth = a_slice['slicewidth']
        PA = a_slice['PA']
        p  = pv_point[0]
        v  = pv_point[1]

        # offset position p to pp position in RA/DEC
        pp = SkyCoord( center.ra - p*np.sin(PA)/np.cos(center.dec),         # right ascension
                       center.dec - p*np.cos(PA),                           # declination
                       frame = 'icrs')

        # a single p point corresponds to a line in pp space perpendicular to the slice
        # sample this line by a half the pixel scale to catch each pixel
        off = np.arange(slicewidth.to(u.arcsec).value/-2., slicewidth.to(u.arcsec).value/2.+pix_scale.to(u.arcsec).value, pix_scale.to(u.arcsec).value)*u.arcsec
        for x in off:
            ppv_points.append([SkyCoord( pp.ra + x*np.cos(PA)/np.cos(pp.dec),
                                         pp.dec - x*np.sin(PA),
                                         frame = 'icrs'),
                               v]
                             )
    return ppv_points

def pv_mask_to_ppv_mask(parameters):
    """
    Convert a pV mask fits file to a ppV mask fits file.

    Parameters
    ----------
    a_slice : dictionary
        Slice object as defined by script 405.
    ax_type : string
        Specifying if a_slice is a major axis or minor axis slice.
    mfac : mask factor
    """

    a_slice = parameters[0]
    ax_type = parameters[1]
    mfac    = parameters[2]
    pV = os.path.join(varmdir,'NGC253.CO.pV_major.slice_major_'+str(a_slice['slice'])+'.mask_factor_'+str(mfac)+'.fits')
    print("pV -> ppV for slice "+str(a_slice['slice'])+", mfac="+str(mfac))

    if not os.path.exists(pV):
        print("does not exist: "+pV)
    else:
        pv_mask    = fits.open(pV)[0].data
        pv_header  = fits.open(pV)[0].header
        pv_wcs     = WCS(pv_header)
        ppv_mask   = fits.open(empty_pV)[0].data
        ppv_header = fits.open(empty_pV)[0].header
        ppv_wcs    = WCS(ppv_header)
        all_ppv_points = []
        crpix,crval,cdelt = get_axis_info(empty_pV,axis=1)
        pv_mask_true = np.transpose(np.where(pv_mask == 1.0))           # get indices where pv_mask is True
        for (v,p) in pv_mask_true:

            # get pV coordinate corresponding to pixel
            pv = pix_to_coord([p,v], pv_header)

            # get ppV coordinate of pV coordinate
            ppv_points = pv_slice_point_to_ppv_coord([pv],a_slice,np.abs(cdelt))
            all_ppv_points.append(ppv_points)

        # flatten list
        all_ppv_points = [item for sublist in all_ppv_points for item in sublist]

        # convert ppv points to ppv coordinates
        ppv_pix = []
        for ppv in all_ppv_points:
            ppv_pix.append(ppv_wcs.all_world2pix(ppv[0].ra.value, ppv[0].dec.value, ppv[1].to(u.m/u.s).value,1))
        ppv_pix = [ [int(round(i.tolist())) for i in j]  for j in ppv_pix]

        # set corresponding pixel in empty cube to 1
        for (x,y,z) in ppv_pix:
            ppv_mask[z,y,x] = 1.0

        print("writing "+'temp.ppV_mask.slice_major_'+str(a_slice['slice'])+'.mask_factor_'+str(mfac)+'.fits')
        fits.writeto(os.path.join(varmdir,'temp.ppV_mask.slice_'+ax_type+'_'+str(a_slice['slice'])+'.mask_factor_'+str(mfac)+'.fits'), data=ppv_mask, header=ppv_header, overwrite=True)


###################################################################################################
# generate the masks

empty_pV = os.path.join(projectdir,'410.disk_non-disk','02.data','NGC253.empty.fits')
execfile('scripts/NGC253/pix_to_coord.py')
execfile('scripts/casa_scripts/get_axis_info.py')

# run in parallel
parameters = [[a_slice, 'major', mfac] for a_slice in major_slices for mfac in mask_factor]
p = Pool(25)
p.map(pv_mask_to_ppv_mask, parameters)


####################################################################################################
# smooth mask

# smooth the mask to cover single pixels that are missed by the pV-to-ppV conversion
for mfac in mask_factor:
    for a_slice in major_slices:
        temp_mask_file = os.path.join(varmdir,'temp.ppV_mask.slice_major_'+str(a_slice['slice'])+'.mask_factor_'+str(mfac)+'.fits')
        ppV_mask_file  = os.path.join(varmdir,'NGC253.CO.ppV_mask.slice_major_'+str(a_slice['slice'])+'.mask_factor_'+str(mfac))
        if os.path.exists(temp_mask_file):
            crpix,crval,cdelt = get_axis_info(temp_mask_file,1)
            smooth = str(abs(cdelt))
            os.system('rm -rf '+temp_mask_file+'.smooth')
            imsmooth(imagename = temp_mask_file,
                outfile = temp_mask_file+'.smooth',
                major   = smooth,
                minor   = smooth,
                pa      = '0deg'
                )
            os.system('rm -rf '+ppV_mask_file)
            immath(imagename = temp_mask_file+'.smooth',
                outfile = ppV_mask_file,
                mode    = 'evalexpr',
                expr    = 'iif(IM0>0.1,1,0)'
                )
            exportfits(imagename = ppV_mask_file,
                fitsimage = ppV_mask_file+'.fits',
                velocity  = True,
                optical   = True,
                overwrite = True,
                dropdeg   = True
                )


###################################################################################################
# merge ppV masks

for mfac in mask_factor:
    ppv_slice_mask_data    = []
    ppv_slice_mask_headers = []
    for a_slice in major_slices:
        ppV_mask_file = os.path.join(varmdir,'NGC253.CO.ppV_mask.slice_major_'+str(a_slice['slice'])+'.mask_factor_'+str(mfac)+'.fits')
        if os.path.exists(ppV_mask_file):
            ppv_slice_mask_data.append( fits.open(ppV_mask_file)[0].data )
            ppv_slice_mask_headers.append( fits.open(ppV_mask_file)[0].header )
    collapsed_mask_data = copy.deepcopy(ppv_slice_mask_data[0])
    for ppv_slice_mask in tqdm(ppv_slice_mask_data):
        collapsed_mask_data += ppv_slice_mask
    collapsed_mask_data[collapsed_mask_data >= 1.0] = 1.0
    fits.writeto(os.path.join(varmdir,'NGC253.CO.ppV_mask_major.mask_factor_'+str(mfac)+'.fits'), data=collapsed_mask_data, header=ppv_slice_mask_headers[0], overwrite=True)


###################################################################################################
# regrid masks
###################################################################################################

for dataset in datasets:
    print("\nRegrid masks ("+dataset['line']+")\n")

    for mfac in mask_factor:
        imregrid(imagename = os.path.join(varmdir,'NGC253.CO.ppV_mask_major.mask_factor_'+str(mfac)+'.fits'),
            template = os.path.join(datdir, dataset['cube']),
            output   = os.path.join(varmdir,'temp.'+dataset['cube']+'.mask_factor_'+str(mfac)),
            asvelocity = True,
            axes       = [-1],
            overwrite  = True
            )
        os.system('rm -rf '+os.path.join(varmdir, dataset['cube']+'.ppV_mask_major'))
        immath(imagename = os.path.join(varmdir,'temp.'+dataset['cube']+'.mask_factor_'+str(mfac)),
            outfile = os.path.join(varmdir, dataset['cube']+'.ppV_mask_major.mask_factor_'+str(mfac)),
            mode    = 'evalexpr',
            expr    = 'iif(IM0>0.1,1,0)'
            )
        exportfits(imagename = os.path.join(varmdir, dataset['cube']+'.ppV_mask_major.mask_factor_'+str(mfac)),
            fitsimage    = os.path.join(varmdir, dataset['cube']+'.ppV_mask_major.mask_factor_'+str(mfac)+'.fits'),
            dropstokes   = True,
            dropdeg      = True,
            overwrite    = True
            )



###################################################################################################
# separate disk from non-disk
###################################################################################################

# taken from 406
# code that is not needed here is removed

def separate_disk_other(dataset, mfac):
    cube_file = os.path.join(datdir, dataset['cube']+'.mask_5.0s')
    mask_file = os.path.join(varmdir, dataset['cube']+'.ppV_mask_major.mask_factor_'+str(mfac))
    disk_file = os.path.join(varmdir, dataset['cube']+'.ppV_mask_major.mask_factor_'+str(mfac)+'.disk.5.0s')
    nondisk_file = os.path.join(varmdir, dataset['cube']+'.ppV_mask_major.mask_factor_'+str(mfac)+'.non-disk.5.0s')

    if os.path.exists(mask_file):
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

        casalog.post("*"*80)
        casalog.post("getting non-disk for "+cube_file)
        casalog.post("*"*80)
        cube           = fits.open(cube_file+'.fits')[0]
        mask           = fits.open(mask_file+'.fits')[0]
        disk_outline   = fits.open(os.path.join(sepdir, dataset['cube']+'.disk_outline.mask.fits'))[0]
        overall_region = fits.open(os.path.join(sepdir, dataset['cube']+'.overall_region.mask.fits'))[0]
        disk_outline.data[np.isnan(disk_outline.data)] = 0
        overall_region.data[np.isnan(overall_region.data)] = 0
        inverted_mask = (mask.data-1.)*-1.
        masked_data   = cube.data*inverted_mask
        masked_data   = np.array([masked_data[i]*disk_outline.data*overall_region.data for i in np.arange(len(masked_data))])
        masked_data[masked_data == 0.0] = np.nan
        fits.writeto(nondisk_file+'.fits', data=masked_data, header=cube.header, overwrite=True)
        importfits(fitsimage = nondisk_file+'.fits',
            imagename = nondisk_file+'.image',
            overwrite = True)

# separate disk from non-disk
for dataset in datasets:
    for mfac in mask_factor:
        separate_disk_other(dataset, mfac)


###################################################################################################
# get flux in disk and non-disk components
###################################################################################################

flux = {}
for dataset in datasets:
    flux[dataset['line']] = {}
    for mfac in mask_factor:
        disk = fits.getdata(os.path.join(varmdir, dataset['cube']+'.ppV_mask_major.mask_factor_'+str(mfac)+'.disk.5.0s.fits'))
        nond = fits.getdata(os.path.join(varmdir, dataset['cube']+'.ppV_mask_major.mask_factor_'+str(mfac)+'.non-disk.5.0s.fits'))
        f_disk = np.nansum(disk)
        f_nond = np.nansum(nond)
        flux[dataset['line']][str(mfac)] = [f_disk, f_nond]

# absolute change due to mask
print('{:>6}'.format('mask')+ '{:^20}'.format('CO (1,0)')+ '{:^20}'.format('CO (2,1)')+ '{:^20}'.format('CO (3,2)'))
print('{:>7}'.format('')+ 3*('{:^10}'.format('disk')+'{:^10}'.format('non-disk')))
for mfac in mask_factor:
    print('{:3d}'.format(int(mfac*100))+'%  '+
          '{:10.2e}'.format(flux['CO_1-0'][str(mfac)][0])+ '{:10.2e}'.format(flux['CO_1-0'][str(mfac)][1])+
          '{:10.2e}'.format(flux['CO_2-1'][str(mfac)][0])+ '{:10.2e}'.format(flux['CO_2-1'][str(mfac)][1])+
          '{:10.2e}'.format(flux['CO_3-2'][str(mfac)][0])+ '{:10.2e}'.format(flux['CO_3-2'][str(mfac)][1]))

# relative change due to mask
print('{:>6}'.format('mask')+ '{:^20}'.format('CO (1,0)')+ '{:^20}'.format('CO (2,1)')+ '{:^20}'.format('CO (3,2)'))
print('{:>7}'.format('')+ 3*('{:^10}'.format('disk')+'{:^10}'.format('non-disk')))
for mfac in mask_factor:
    print('{:3d}'.format(int(mfac*100))+'%  '+
          '{:10.1f}'.format((flux['CO_1-0'][str(mfac)][0]/flux['CO_1-0']['1.0'][0])*100)+ '{:10.1f}'.format((flux['CO_1-0'][str(mfac)][1]/flux['CO_1-0']['1.0'][1])*100)+
          '{:10.1f}'.format((flux['CO_2-1'][str(mfac)][0]/flux['CO_2-1']['1.0'][0])*100)+ '{:10.1f}'.format((flux['CO_2-1'][str(mfac)][1]/flux['CO_2-1']['1.0'][1])*100)+
          '{:10.1f}'.format((flux['CO_3-2'][str(mfac)][0]/flux['CO_3-2']['1.0'][0])*100)+ '{:10.1f}'.format((flux['CO_3-2'][str(mfac)][1]/flux['CO_3-2']['1.0'][1])*100))

# relative change in dex
print('{:>6}'.format('mask')+ '{:^20}'.format('CO (1,0)')+ '{:^20}'.format('CO (2,1)')+ '{:^20}'.format('CO (3,2)'))
print('{:>7}'.format('')+ 3*('{:^10}'.format('disk')+'{:^10}'.format('non-disk')))
for mfac in mask_factor:
    print('{:3d}'.format(int(mfac*100))+'%  '+
          '{:10.2f}'.format(np.log10(flux['CO_1-0'][str(mfac)][0]/flux['CO_1-0']['1.0'][0]))+ '{:10.2f}'.format(np.log10(flux['CO_1-0'][str(mfac)][1]/flux['CO_1-0']['1.0'][1]))+
          '{:10.2f}'.format(np.log10(flux['CO_2-1'][str(mfac)][0]/flux['CO_2-1']['1.0'][0]))+ '{:10.2f}'.format(np.log10(flux['CO_2-1'][str(mfac)][1]/flux['CO_2-1']['1.0'][1]))+
          '{:10.2f}'.format(np.log10(flux['CO_3-2'][str(mfac)][0]/flux['CO_3-2']['1.0'][0]))+ '{:10.2f}'.format(np.log10(flux['CO_3-2'][str(mfac)][1]/flux['CO_3-2']['1.0'][1])))



###################################################################################################
#
###################################################################################################
