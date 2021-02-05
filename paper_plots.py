##############################
# EASY APLPY PLOTTING SCRIPT #
##############################

# Generate the plots of the paper Krieger+19

###################################################################################################

# import required modules
execfile('NGC253/project_info.py')
import matplotlib.ticker as ticker

###################################################################################################

# directories
paperdir = os.path.join(projectdir,'paper_K19')
plotdir  = os.path.join('XXXXX','paper_K19')
os.system('mkdir -p '+plotdir)

# change plot preferences
# This script uses super sampling with a factor of two to get thin contour lines!
easy_aplpy.settings.velo_fontsize           = 36      # unit: point
easy_aplpy.settings.colorbar_label_fontsize = 36      # unit: point
easy_aplpy.settings.colorbar_ticks_fontsize = 36      # unit: point
easy_aplpy.settings.grid_label_fontsize     = 36      # unit: point
easy_aplpy.settings.colorbar_width          = 0.15    # relative to panel size
easy_aplpy.settings.scalebar_frame          = False
easy_aplpy.settings.scalebar_linestyle      = 'solid' # or any other plt.plot linestyle
easy_aplpy.settings.scalebar_linewidth      = 4       # unit: points
easy_aplpy.settings.scalebar_color          = 'red'   # any named color or mpl.color instance
easy_aplpy.settings.scalebar_fontsize       = 36      # only used in channel map to prevent bar sliding over map
easy_aplpy.settings.beam_frame              = False
easy_aplpy.settings.beam_color              = 'black'
easy_aplpy.settings.ticks_color             = 'black' # this setting overrules the matplotlibrc defaults
easy_aplpy.settings.frame_color             = 'black'
easy_aplpy.settings.tick_label_fontsize     = 36      # unit: point
easy_aplpy.settings.axis_label_fontsize     = 36      # unit: point
easy_aplpy.settings.props                 = {'boxstyle': "round", 'facecolor': "w", 'edgecolor': "black", 'linewidth': 1.0, 'alpha': 1.0}

easy_aplpy.settings.ticks_xspacing = Angle('0 0 1.0', unit='hourangle')
easy_aplpy.settings.ticks_yspacing = 10.0*u.arcsec
easy_aplpy.settings.ticks_minor_frequency = 5

default_recenter = [SkyCoord('00h47m33.01s -25d17m20.0s'), 46.0*u.arcsec, 35.0*u.arcsec]
default_scalebar = [14.73*u.arcsec, '250\,pc', 'bottom left']


###################################################################################################

# figure 1: CO channel map
##########################

easy_aplpy.settings.colorbar_width          = 0.05

#COcube = os.path.join(paperdir,'paper.CO.cube.fits')
COcube = os.path.join(paperdir,'paper.CO32.cube.mask_3s.fits')
easy_aplpy.plot.grid(COcube,
    [3,4],
    [20*u.km/u.s,60*u.km/u.s,100*u.km/u.s,140*u.km/u.s,180*u.km/u.s,220*u.km/u.s,260*u.km/u.s,300*u.km/u.s,340*u.km/u.s,380*u.km/u.s,420*u.km/u.s,460*u.km/u.s],
    figsize  = (30,30),
    cmap     = 'jet', #easy_aplpy.custom_colormaps.viridis_cropped,
    colorbar = ['right', 'T$_{\mathrm{b}}$ [K]'],
    vmin     = 1,
    vmax     = 75,
    stretch  = 'log',
    recenter = default_recenter,
    channel_label = 'physical',
    contours = [[[COcube, ch, [x*0.37 for x in [10,20,40,80]], 'black']] for ch in np.arange(1,10+3*4*16,16)],
    scalebar = default_scalebar,
    out      = os.path.join(plotdir,'fig1.pdf')
    )

easy_aplpy.settings.colorbar_width          = 0.15


###################################################################################################

# figure 2: CO moment maps
##########################

mommap = os.path.join(paperdir,'paper.CO.mom0.fits')
easy_aplpy.plot.map(mommap,
    figsize   = (16,16),
    cmap      = easy_aplpy.custom_colormaps.viridis_cropped,
    stretch   = 'linear',
    vmin      = -800,
    vmax      = 8000,
    recenter  = [SkyCoord('00h47m33.01s -25d17m17.0s'), 46.0*u.arcsec, 36.0*u.arcsec],
    contours  = [[mommap, [250, 500, 1000, 2000, 4000, 8000], 'black']],
    colorbar  = ['right', 'F [K\,km\,s$^{-1}$]'],
    scalebar  = default_scalebar,
    beam      = 'bottom right',
    texts     = [[[0.0,1.0], 'moment 0', {'size':60, 'ha':'left', 'va':'top', 'style':'oblique', 'weight':'black', 'color':'black', 'bbox': easy_aplpy.settings.props}]],
    out       = os.path.join(plotdir,'fig2a.pdf')
    )

mommap = os.path.join(paperdir,'paper.CO.mom1.fits')
easy_aplpy.plot.map(mommap,
    figsize   = (16,16),
    cmap      = 'RdBu',
    stretch   = 'linear',
    vmin      = 50,
    vmax      = 450,
    recenter  = [SkyCoord('00h47m33.01s -25d17m17.0s'), 46.0*u.arcsec, 36.0*u.arcsec],
    contours  = [[mommap, [100,150,200,250,300,350,400], 'black']],
    colorbar  = ['right', '$v$ [km\,s$^{-1}$]'],
    scalebar  = default_scalebar,
    beam      = 'bottom right',
    texts     = [[[0.0,1.0], 'moment 1', {'size':60, 'ha':'left', 'va':'top', 'style':'oblique', 'weight':'black', 'color':'black', 'bbox': easy_aplpy.settings.props}]],
    out       = os.path.join(plotdir,'fig2b.pdf')
    )

mommap = os.path.join(paperdir,'paper.CO.mom2.fits')
easy_aplpy.plot.map(mommap,
    figsize   = (16,16),
    cmap      = 'magma',
    stretch   = 'linear',
    vmin      = -10,
    vmax      = 100,
    recenter  = [SkyCoord('00h47m33.01s -25d17m17.0s'), 46.0*u.arcsec, 36.0*u.arcsec],
    contours  = [[mommap, [20,40,60,80,100], 'black']],
    colorbar  = ['right', '$\Delta v$ [km\,s$^{-1}$]'],
    scalebar  = default_scalebar,
    beam      = 'bottom right',
    texts     = [[[0.0,1.0], 'moment 2', {'size':60, 'ha':'left', 'va':'top', 'style':'oblique', 'weight':'black', 'color':'black', 'bbox': easy_aplpy.settings.props}]],
    out       = os.path.join(plotdir,'fig2c.pdf')
    )


###################################################################################################

# figure 3: position-velocity diagrams along major axis
#######################################################

pvmap = os.path.join(paperdir, 'paper.CO.pv_major.fits')
#easy_aplpy.helpers.m_to_km(pvmap, overwrite=True)
easy_aplpy.plot.map(pvmap,
    figsize   = (16,16),
    cmap      = easy_aplpy.custom_colormaps.viridis_cropped,
    stretch   = 'linear',
    vmin      = -3,
    vmax      = 30,
    recenter  = [0.0*u.arcsec, 260*u.km/u.s, 50*u.arcsec, 500*u.km/u.s],
    contours  = [[pvmap, [x*0.37 for x in [10,20,40,80]], 'black']],
    colorbar  = ['right', 'T$_{\mathrm{b}}$ [K]'],
    labels    = [r'offset along major axis [$^{\prime\prime}$]', r'velocity [km\,s$^{-1}$]'],
    out       = os.path.join(plotdir, 'fig3.pdf')
    )


###################################################################################################

# figure 4: disk/non-disk slice overview
########################################

slicewidth   = 5.0*u.arcsec
major_length = 50.*u.arcsec
minor_length = 50.*u.arcsec
nslices_major = 2*int(minor_length.value/slicewidth.value)-1
major_slices = fnunpickle(os.path.join(paperdir,'slices_major.pickle'))

polygons = []
labels = []
for a_slice in major_slices:
    offset    = (a_slice['slice']+1)*(major_length-3*u.arcsec)/(nslices_major+1)-(major_length-3*u.arcsec)/2.
    label_pos = SkyCoord(a_slice['center'].ra+offset*np.sin(a_slice['PA'])/np.cos(a_slice['center'].dec),a_slice['center'].dec+offset*np.cos(a_slice['PA']), frame='icrs')
    polygons.append( [ [a_slice['edge1'], a_slice['edge2'], a_slice['edge3'], a_slice['edge4']], {'facecolor':"black", 'alpha':0.25} ] )
    polygons.append( [ [a_slice['edge1'], a_slice['edge2'], a_slice['edge3'], a_slice['edge4']], {'edgecolor':"black"} ] )
    labels.append( [ label_pos, str((a_slice['slice']-int(len(major_slices)/2.))*(a_slice['slicewidth']/2.).to(u.arcsec).value)+'$^{\prime\prime}$', {'rotation':(a_slice['PA']-90*u.degree).value, 'ha':"center", 'va':"center", 'color':"white", 'size':40} ] )

for panel,line,name in [['a','CO_1-0','CO (1-0)'], ['b','CO_2-1','CO (2-1)'], ['c','CO_3-2','CO (3-2)']]:
    mommap = os.path.join(paperdir,'NGC253.'+line+'.mask_3.0s.mom0.fits')
    labels.append([[0.05,0.95], name, {'size':48, 'ha':'left', 'va':'top', 'style':'oblique', 'weight':'black', 'color':'black', 'bbox':easy_aplpy.settings.props}]),
    easy_aplpy.plot.map(mommap,
        figsize   = (16,16),
        cmap      = easy_aplpy.custom_colormaps.viridis_cropped,
        stretch   = 'log',
        vmin      = 1,
        vmax      = 8000,
        recenter  = [SkyCoord('00h47m33.18s -25d17m20.5s'), 1.2*u.arcmin, 1.2*u.arcmin],
        contours  = [[mommap, [250, 500, 1000, 2000, 4000, 8000], 'black']],
        colorbar  = ['right', 'F [K\,km\,s$^{-1}$]'],
        scalebar  = default_scalebar,
        beam      = 'bottom right',
        texts     = labels,
        polygons  = polygons,
        out       = os.path.join(plotdir,'fig4'+panel+'.pdf')
        )




###################################################################################################

# figure 5: example pV plots for unbiased masking
#################################################

# labels need to be corrected by hand!

execfile('scripts/NGC253/410.disk_non-disk/400.disk_non-disk.helper_functions.py')
pV_path      = 'NGC253/410.disk_non-disk/04.pVs/'
kin_path     = 'NGC253/410.disk_non-disk/05.kinematics/'
total_model  = 'NGC253/410.disk_non-disk/07.model/diskfit.total_model.fits'
tmodel_major = 'NGC253/410.disk_non-disk/07.model/diskfit.total_model.major_slice_'
tmodel_minor = 'NGC253/410.disk_non-disk/07.model/diskfit.total_model.minor_slice_'

for line,name,panel in [['CO_1-0', 'CO (1-0)', 'a'],['CO_2-1', 'CO (2-1)', 'b'],['CO_3-2', 'CO (3-2)', 'c']]:
    a_slice = major_slices[9]

    contours = []
    lines    = []
    texts    = []
    pV    = pV_path+'NGC253.'+line+'.pV_major.slice_major_'+str(a_slice['slice'])+'.fits'
    model = tmodel_major+str(a_slice['slice'])+'.fits'
    var_mask   = pV_path+'NGC253.CO.pV_major.slice_major_'+str(a_slice['slice'])+'.mask_model_variable.fits'

    if not os.path.exists(pV):
        pV = pV_path+'pV_empty.fits'

    if os.path.exists(pV):
        contours.append([pV, [1,2,4,8,16,32], 'black'])
    if os.path.exists(var_mask):
        contours.append([var_mask, [0.5,1], 'gold', {'filled': True, 'alpha': 0.25, 'smooth': 1}])
        contours.append([var_mask, [0.5], 'gold', {'filled': False, 'alpha': 0.4, 'smooth': 1}])

    if os.path.exists(model):
        model_data = fits.open(model)[0].data[0]
        lines.append([[i*u.arcsec for i in np.arange(-25.,25.1,0.1)], [(j*u.km/u.s).to(u.m/u.s) for j in model_data], {'color': 'red'}])

    texts.append([[0.05,0.95], name, {'size': 28, 'ha':'left', 'va':'top', 'color':'black', 'bbox': easy_aplpy.settings.props}])

    fig = easy_aplpy.plot.map(pV,
        figsize  = (10,8),
        recenter = [0*u.arcsec, 260*u.km/u.s, major_length, 500*u.km/u.s],
        labels   = ['offset along major axis [$^{\prime\prime}$]','velocity [km\,s$^{-1}$]'],
        contours = contours,
        lines    = lines,
        cmap     = 'binary',
        vmin     = -5,
        vmax     = 30,
        stretch  = 'linear',
        texts    = texts,
        out      = os.path.join(plotdir,'fig5'+panel+'.pdf')
        )


###################################################################################################

# figure 6: separated moment maps
#################################

# uses old aplpy_plotting code
ap._velo_fontsize           = 36      # unit: point
ap._colorbar_fontsize       = 36      # unit: point
ap._colorbar_ticks_fontsize = 36      # unit: point
ap._grid_label_fontsize     = 36      # unit: point
ap._colorbar_width          = 0.15    # relative to panel size
ap._scalebar_frame          = False
ap._scalebar_linestyle      = 'solid' # or any other plt.plot linestyle
ap._scalebar_linewidth      = 4       # unit: points
ap._scalebar_color          = 'red'   # any named color or mpl.color instance
ap._scalebar_fontsize       = 36      # only used in channel map to prevent bar sliding over map
ap._beam_frame              = False
ap._beam_color              = 'black'
ap._ticks_color             = 'black' # this setting overrules the matplotlibrc defaults
ap._frame_color             = 'black'
ap._tick_label_fontsize     = 36      # unit: point
ap._axis_label_fontsize     = 36      # unit: point
ap._props                 = {'boxstyle': "round", 'facecolor': "w", 'edgecolor': "black", 'linewidth': 1.0, 'alpha': 1.0}

ap.ticks_xspacing = Angle('0 0 1.0', unit='hourangle')
ap.ticks_yspacing = 10.0*u.arcsec
ap.ticks_minor_frequency = 5



# moment 0
infiles = ['NGC253/420.disk_non-disk/03.moments/NGC253.CO_'+CO+'.mask_5.0s.mom0.log.fits' for CO in ['1-0','2-1','3-2']]
infiles += ['NGC253/420.disk_non-disk/06.disk_non-disk/NGC253.CO_'+CO+'.ppV_mask_major.'+kin_type+'.5.0s.mom0.log.fits' for kin_type in ['disk','non-disk'] for CO in ['1-0','2-1','3-2']]
contours = []
for idx,f in enumerate(infiles):
    if idx%3==2:
        contours.append([f, [2.0,2.5,3.0,3.5], 'black'])
    else:
        contours.append([f, [1.7,2.0,2.3,2.6,2.9,3.2,3.5], 'black'])


ap.aplpy_moments_grid(infiles,
    scaling  = [[1, 3.8, r'log (F [K\,km\,s$^{-1}$])',  'linear', ap.viridis_cropped],
                [1, 3.8, r'log (F [K\,km\,s$^{-1}$])',  'linear', ap.viridis_cropped],
                [1, 3.8, r'log (F [K\,km\,s$^{-1}$])',  'linear', ap.viridis_cropped]],
    contours = contours,
    #labels   = [None, 'disk', None, None, 'non-disk', None],
    labels   = ['CO ('+CO+') ' for CO in ['1-0','2-1','3-2']] + ['','',''] + ['','',''],
    label_kwargs = {'bbox': ap._props},
    figsize  = (30.0,32.0),
    recenter = [SkyCoord('00h47m33.18s -25d17m20.5s'), 0.8*u.arcmin, 0.8*u.arcmin],
    execute_code = [["fig.add_label(0.05,0.05, relative=True, text='disk + non-disk', size=36, rotation='vertical', ha='left', va='bottom', style='oblique', weight='black', color='black', bbox={'boxstyle': 'round', 'facecolor': 'w', 'edgecolor': 'black', 'linewidth': 2.0, 'alpha': 1.0})"],[],[],
                    ["fig.add_label(0.05,0.05, relative=True, text='disk', size=36, rotation='vertical', ha='left', va='bottom', style='oblique', weight='black', color='black', bbox={'boxstyle': 'round', 'facecolor': 'w', 'edgecolor': 'black', 'linewidth': 2.0, 'alpha': 1.0})"],[],[],
                    ["fig.add_label(0.05,0.05, relative=True, text='non-disk', size=36, rotation='vertical', ha='left', va='bottom', style='oblique', weight='black', color='black', bbox={'boxstyle': 'round', 'facecolor': 'w', 'edgecolor': 'black', 'linewidth': 2.0, 'alpha': 1.0})"],[],[]],
    out      = os.path.join(plotdir,'fig6a.pdf')
    )

# moment 1
infiles = ['NGC253/420.disk_non-disk/03.moments/NGC253.CO_'+CO+'.mask_10.0s.mom1.fits' for CO in ['1-0','2-1','3-2']]
infiles += ['NGC253/420.disk_non-disk/06.disk_non-disk/NGC253.CO_'+CO+'.ppV_mask_major.'+kin_type+'.10.0s.mom1.fits' for kin_type in ['disk','non-disk'] for CO in ['1-0','2-1','3-2']]
contours = [[f, [150,200,250,300,350], 'black'] for f in infiles]
ap.aplpy_moments_grid(infiles,
    scaling  = [[50, 450,  r'v$_{rad}$ [km\,s$^{-1}$]', 'linear', 'RdBu'],
                [50, 450,  r'v$_{rad}$ [km\,s$^{-1}$]', 'linear', 'RdBu'],
                [50, 450,  r'v$_{rad}$ [km\,s$^{-1}$]', 'linear', 'RdBu']],
    contours = contours,
    #labels   = [None, 'disk', None, None, 'non-disk', None],
    labels   = ['CO ('+CO+') ' for CO in ['1-0','2-1','3-2']] + ['','',''] + ['','',''],
    label_kwargs = {'bbox': ap._props},
    figsize  = (30,32),
    recenter = [SkyCoord('00h47m33.18s -25d17m20.5s'), 0.8*u.arcmin, 0.8*u.arcmin],
    execute_code = [["fig.add_label(0.05,0.05, relative=True, text='disk + non-disk', size=36, rotation='vertical', ha='left', va='bottom', style='oblique', weight='black', color='black', bbox={'boxstyle': 'round', 'facecolor': 'w', 'edgecolor': 'black', 'linewidth': 2.0, 'alpha': 1.0})"],[],[],
                    ["fig.add_label(0.05,0.05, relative=True, text='disk', size=36, rotation='vertical', ha='left', va='bottom', style='oblique', weight='black', color='black', bbox={'boxstyle': 'round', 'facecolor': 'w', 'edgecolor': 'black', 'linewidth': 2.0, 'alpha': 1.0})"],[],[],
                    ["fig.add_label(0.05,0.05, relative=True, text='non-disk', size=36, rotation='vertical', ha='left', va='bottom', style='oblique', weight='black', color='black', bbox={'boxstyle': 'round', 'facecolor': 'w', 'edgecolor': 'black', 'linewidth': 2.0, 'alpha': 1.0})"],[],[]],
    out      = os.path.join(plotdir,'fig6b.pdf')
    )

# moment 2
infiles = ['NGC253/420.disk_non-disk/03.moments/NGC253.CO_'+CO+'.mask_5.0s.mom2.fits' for CO in ['1-0','2-1','3-2']]
infiles += ['NGC253/420.disk_non-disk/06.disk_non-disk/NGC253.CO_'+CO+'.ppV_mask_major.'+kin_type+'.5.0s.mom2.fits' for kin_type in ['disk'] for CO in ['1-0','2-1','3-2']]
contours = [[f, [20,40,60,80,100], 'black'] for f in infiles]
ap.aplpy_moments_grid(infiles,
    scaling  = [[0,   80,   r'$\sigma$ [km\,s$^{-1}$]', 'linear', 'magma'],
                [0,   80,   r'$\sigma$ [km\,s$^{-1}$]', 'linear', 'magma'],
                [0,   80,   r'$\sigma$ [km\,s$^{-1}$]', 'linear', 'magma']],
    contours = contours,
    #labels   = [None, 'disk', None, None, 'non-disk', None],
    labels   = ['CO ('+CO+') ' for CO in ['1-0','2-1','3-2']] + ['','',''] + ['','',''],
    label_kwargs = {'bbox': ap._props},
    figsize  = (30,22),
    recenter = [SkyCoord('00h47m33.18s -25d17m20.5s'), 0.8*u.arcmin, 0.8*u.arcmin],
    execute_code = [["fig.add_label(0.05,0.05, relative=True, text='disk + non-disk', size=36, rotation='vertical', ha='left', va='bottom', style='oblique', weight='black', color='black', bbox={'boxstyle': 'round', 'facecolor': 'w', 'edgecolor': 'black', 'linewidth': 2.0, 'alpha': 1.0})"],[],[],
                    ["fig.add_label(0.05,0.05, relative=True, text='disk', size=36, rotation='vertical', ha='left', va='bottom', style='oblique', weight='black', color='black', bbox={'boxstyle': 'round', 'facecolor': 'w', 'edgecolor': 'black', 'linewidth': 2.0, 'alpha': 1.0})"],[],[]],
    out      = os.path.join(plotdir,'fig6c.pdf')
    )


###################################################################################################

# figure 7: zoom-in on CO (3-2) non-disk
########################################

sepdir  = os.path.join(projectdir, '06.disk_non-disk')
maskdir = os.path.join(projectdir, '07.mass_outflow_rate')
bubmap = os.path.join(maskdir, 'nondisk.CO_3-2.superbubble_new.mask.fits')
codmap = os.path.join(maskdir, 'nondisk.CO_3-2.cospatial_disk_new.mask.fits')

mommap = os.path.join(sepdir,'NGC253.CO_3-2.ppV_mask_major.non-disk.5.0s.mom0.log.fits')
easy_aplpy.plot.map(mommap,
    figsize   = (24,24),
    cmap      = easy_aplpy.custom_colormaps.viridis_cropped,
    stretch   = 'linear',
    vmin      = 1,
    vmax      = 3.0,
    recenter  = [SkyCoord('00h47m33.01s -25d17m18.5s'), 48.0*u.arcsec, 40.0*u.arcsec],
    contours  = [[mommap, [2.0,2.5,3.0], 'black'],
                 [bubmap, [0.5], 'red', {'linewidths': 5}],
                 [codmap, [0.5], 'red', {'linewidths': 5}]],
    colorbar  = ['right', r'log (F [K\,km\,s$^{-1}$])'],
    scalebar  = default_scalebar,
    beam      = 'bottom right',
    texts     = [[[0.0,1.0], 'moment 0', {'size':48, 'ha':'left', 'va':'top', 'style':'oblique', 'weight':'black', 'color':'black', 'bbox': easy_aplpy.settings.props}]],
    out       = os.path.join(plotdir,'fig7a.pdf')
    )

mommap = os.path.join(sepdir,'NGC253.CO_3-2.ppV_mask_major.non-disk.5.0s.mom1.fits')
easy_aplpy.plot.map(mommap,
    figsize   = (24,24),
    cmap      = 'RdBu',
    stretch   = 'linear',
    vmin      = 50,
    vmax      = 450,
    recenter  = [SkyCoord('00h47m33.01s -25d17m18.5s'), 46.0*u.arcsec, 40.0*u.arcsec],
    contours  = [#[mommap, [150,200,250,300,350], 'black'],
                 [bubmap, [0.5], 'black', {'linewidths': 5}],
                 [codmap, [0.5], 'black', {'linewidths': 5}]],
    colorbar  = ['right', '$v$ [km\,s$^{-1}$]'],
    scalebar  = default_scalebar,
    beam      = 'bottom right',
    texts     = [[[0.0,1.0], 'moment 1', {'size':48, 'ha':'left', 'va':'top', 'style':'oblique', 'weight':'black', 'color':'black', 'bbox': easy_aplpy.settings.props}]],
    out       = os.path.join(plotdir,'fig7b.pdf')
    )


###################################################################################################
# figure 8: outflow rates
###################################################################################################

datasets = [{'line': 'CO_1-0', 'cube': 'NGC253.CO_1-0', 'plot_color': 'aqua'},
            {'line': 'CO_2-1', 'cube': 'NGC253.CO_2-1', 'plot_color': 'dodgerblue'},
            {'line': 'CO_3-2', 'cube': 'NGC253.CO_3-2', 'plot_color': 'midnightblue'}]

# plot rate evolution
rates_td = {}
for dataset in datasets:
    rates_td[dataset['line']] = np.genfromtxt(os.path.join(projectdir, '07.mass_outflow_rate.improved3', dataset['cube']+'.major.5.0s.rate_bootstrap.real_outflow.txt'),
                                              names = ('bin','rate','rate_p','rate_m')
                                             )
fig = plt.figure(figsize=(8,4))
ax  = fig.subplots(1)
ylim = 0.0
lowlim = 0.2
for dataset,lim in zip(datasets,[5,4,3]):
    bin    = rates_td[dataset['line']]['bin']
    rate   = rates_td[dataset['line']]['rate']
    rate_p = rates_td[dataset['line']]['rate_p']
    rate_m = rates_td[dataset['line']]['rate_m']
    ax.plot(bin[bin<=lowlim], rate[bin<=lowlim],
            lw=3, ls=':',
            color=dataset['plot_color']
           )
    ax.plot(bin[(bin<=lim) & (bin>=lowlim)], rate[(bin<=lim) & (bin>=lowlim)],
            lw=3, ls='-',
            color=dataset['plot_color'],
            label=dataset['line'].replace('_','(').replace('-','--')+')'
           )
    ax.plot(bin[bin>=lim], rate[bin>=lim],
            lw=3, ls=':',
            color=dataset['plot_color']
           )
    ax.fill_between(bin, rate_p, rate_m,
                           color=dataset['plot_color'],
                           alpha=0.25
                          )
    if (ylim<np.nanmax([rate,rate_p,rate_m])): ylim = np.nanmax([rate,rate_p,rate_m])
ax.grid(True)
ax.set_axisbelow(True)
ax.set_xlabel(r't$_\mathrm{eject}$ [Myr]')
ax.set_ylabel(r'$\dot{\mathrm{M}}_{out}$ [M$_\odot$\,yr$^{-1}$]')
ax.set_xlim(-0.1,5.1)
ax.set_ylim(0.1, 2*ylim)
ax.set_yscale('log')
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
ax.legend(loc=3, bbox_to_anchor=(0.0,1.15,1.0,0.1), borderaxespad=0., ncol=3, numpoints=1, fancybox=True, mode='expand', fontsize=12)
fig.savefig(os.path.join(plotdir,'fig8a.pdf'), bbox_inches='tight')


# plot distance evolution
rates_dd = {}
for dataset in datasets:
    rates_dd[dataset['line']] = np.genfromtxt(os.path.join(projectdir, '07.mass_outflow_rate.improved3', dataset['cube']+'.major.5.0s.dist_bootstrap.real_outflow.txt'),
                                              names = ('bin','rate','rate_p','rate_m')
                                             )

fig = plt.figure(figsize=(8,4))
ax = fig.subplots(1)
ylim = 0.0
lowlim = 33.9369576
for dataset,lim in zip(datasets,[509.054364,322.4010972,203.6217456]):
    bin    = rates_dd[dataset['line']]['bin']*16.9684788
    rate   = rates_dd[dataset['line']]['rate']
    rate_p = rates_dd[dataset['line']]['rate_p']
    rate_m = rates_dd[dataset['line']]['rate_m']
    rate_m = [r if r>0. else 0. for r in rate_m]
    ax.plot(bin[bin<=lowlim], rate[bin<=lowlim],
            lw=3, ls=':',
            color=dataset['plot_color']
           )
    ax.plot(bin[(bin<=lim) & (bin>=lowlim)], rate[(bin<=lim) & (bin>=lowlim)],
            lw=3, ls='-',
            color=dataset['plot_color'],
            label=dataset['line'].replace('_','(').replace('-','--')+')'
           )
    ax.plot(bin[bin>=lim], rate[bin>=lim],
            lw=3, ls=':',
            color=dataset['plot_color']
           )
    ax.fill_between(bin, rate_p, rate_m,
                    color=dataset['plot_color'],
                    alpha=0.25
                   )
    if (ylim<np.nanmax([rate,rate_p,rate_m])): ylim = np.nanmax([rate,rate_p,rate_m])
ax.grid(True)
ax.set_axisbelow(True)
ax.set_xlabel(r'd [pc]')
ax.set_ylabel(r'$\dot{\mathrm{M}}_{out}$ [M$_\odot$\,yr$^{-1}$]')
ax.set_xlim(-10,510)
ax.set_ylim(0.1, 2*ylim)
ax.set_yscale('log')
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
fig.savefig(os.path.join(plotdir,'fig8b.pdf'), bbox_inches='tight')


###################################################################################################
# figure 9: energy and momentum evolution
###################################################################################################

datasets = [{'line': 'CO_1-0', 'cube': 'NGC253.CO_1-0', 'plot_color': 'aqua'},
            {'line': 'CO_2-1', 'cube': 'NGC253.CO_2-1', 'plot_color': 'dodgerblue'},
            {'line': 'CO_3-2', 'cube': 'NGC253.CO_3-2', 'plot_color': 'midnightblue'}]

# plot rate evolution
energies_td = {}
for dataset in datasets:
    energies_td[dataset['line']] = np.genfromtxt(os.path.join(projectdir, '08.results', dataset['cube']+'.major.5.0s.energy_time_bootstrap.real_outflow.txt'),
                                              names = ('bin','q','q_p','q_m')
                                             )
fig = plt.figure(figsize=(8,4))
ax  = fig.subplots(1)
ylim = 0.0
lowlim = 0.2
for dataset,lim in zip(datasets,[5,4,3]):
    bin    = energies_td[dataset['line']]['bin']
    energy   = energies_td[dataset['line']]['q']
    energy_p = energies_td[dataset['line']]['q_p']
    energy_m = energies_td[dataset['line']]['q_m']
    ax.plot(bin[bin<=lowlim], energy[bin<=lowlim],
            lw=3, ls=':',
            color=dataset['plot_color']
           )
    ax.plot(bin[(bin<=lim) & (bin>=lowlim)], energy[(bin<=lim) & (bin>=lowlim)],
            lw=3, ls='-',
            color=dataset['plot_color'],
            label=dataset['line'].replace('_','(').replace('-','--')+')'
           )
    ax.plot(bin[bin>=lim], energy[bin>=lim],
            lw=3, ls=':',
            color=dataset['plot_color']
           )
    ax.fill_between(bin, energy_p, energy_m,
                           color=dataset['plot_color'],
                           alpha=0.25
                          )
    if (ylim<np.nanmax([energy,energy_p,energy_m])): ylim = np.nanmax([energy,energy_p,energy_m])
ax.set_axisbelow(True)
ax.grid(True)
ax.set_xlabel(r't$_\mathrm{eject}$ [Myr]')
ax.set_ylabel(r'E$_\mathrm{kin}$ [erg]')
ax.set_xlim(-0.1,5.1)
ax.set_ylim(1e48, 2*ylim)
ax.set_yscale('log')
ax.legend(loc=3, bbox_to_anchor=(0.0,1.15,1.0,0.1), borderaxespad=0., ncol=3, numpoints=1, fancybox=True, mode='expand', fontsize=12)
fig.savefig(os.path.join(plotdir,'fig9a.pdf'), bbox_inches='tight')


# plot distance evolution
energies_dd = {}
for dataset in datasets:
    energies_dd[dataset['line']] = np.genfromtxt(os.path.join(projectdir, '08.results', dataset['cube']+'.major.5.0s.energy_dist_bootstrap.real_outflow.txt'),
                                              names = ('bin','q','q_p','q_m')
                                             )
fig = plt.figure(figsize=(8,4))
ax  = fig.subplots(1)
ylim = 0.0
lowlim = 33.9369576
for dataset,lim in zip(datasets,[509.054364,322.4010972,203.6217456]):
    bin    = energies_dd[dataset['line']]['bin']*16.9684788
    energy   = energies_dd[dataset['line']]['q']
    energy_p = energies_dd[dataset['line']]['q_p']
    energy_m = energies_dd[dataset['line']]['q_m']
    energy_m = [r if r>0. else 0. for r in energy_m]
    ax.plot(bin[bin<=lowlim], energy[bin<=lowlim],
            lw=3, ls=':',
            color=dataset['plot_color']
           )
    ax.plot(bin[(bin<=lim) & (bin>=lowlim)], energy[(bin<=lim) & (bin>=lowlim)],
            lw=3, ls='-',
            color=dataset['plot_color'],
            label=dataset['line'].replace('_','(').replace('-','--')+')'
           )
    ax.plot(bin[bin>=lim], energy[bin>=lim],
            lw=3, ls=':',
            color=dataset['plot_color']
           )
    ax.fill_between(bin, energy_p, energy_m,
                           color=dataset['plot_color'],
                           alpha=0.25
                          )
    if (ylim<np.nanmax([energy,energy_p,energy_m])): ylim = np.nanmax([energy,energy_p,energy_m])
ax.set_axisbelow(True)
ax.grid(True)
ax.set_xlabel(r'd [pc]')
ax.set_ylabel(r'E$_\mathrm{kin}$ [erg]')
ax.set_xlim(-10,510)
ax.set_ylim(1e51, 2*ylim)
ax.set_yscale('log')
fig.savefig(os.path.join(plotdir,'fig9c.pdf'), bbox_inches='tight')


# plot momentum evolution
momenta_td = {}
for dataset in datasets:
    momenta_td[dataset['line']] = np.genfromtxt(os.path.join(projectdir, '08.results', dataset['cube']+'.major.5.0s.momentum_time_bootstrap.real_outflow.txt'),
                                              names = ('bin','q','q_p','q_m')
                                             )
fig = plt.figure(figsize=(8,4))
ax  = fig.subplots(1)
ylim = 0.0
lowlim = 0.2
for dataset,lim in zip(datasets,[5,4,3]):
    bin    = momenta_td[dataset['line']]['bin']
    momentum   = momenta_td[dataset['line']]['q']
    momentum_p = momenta_td[dataset['line']]['q_p']
    momentum_m = momenta_td[dataset['line']]['q_m']
    ax.plot(bin[bin<=lowlim], momentum[bin<=lowlim],
            lw=3, ls=':',
            color=dataset['plot_color']
           )
    ax.plot(bin[(bin<=lim) & (bin>=lowlim)], momentum[(bin<=lim) & (bin>=lowlim)],
            lw=3, ls='-',
            color=dataset['plot_color'],
            label=dataset['line'].replace('_','(').replace('-','--')+')'
           )
    ax.plot(bin[bin>=lim], momentum[bin>=lim],
            lw=3, ls=':',
            color=dataset['plot_color']
           )
    ax.fill_between(bin, momentum_p, momentum_m,
                           color=dataset['plot_color'],
                           alpha=0.25
                          )
    if (ylim<np.nanmax([momentum,momentum_p,momentum_m])): ylim = np.nanmax([momentum,momentum_p,momentum_m])
ax.set_axisbelow(True)
ax.grid(True)
ax.set_xlabel(r't$_\mathrm{eject}$ [Myr]')
ax.set_ylabel(r'P [M$_\odot$\,km\,s$^{-1}$]')
ax.yaxis.set_label_position('right')
ax.tick_params('y', left=False, labelleft=False, right=True, labelright=True)
ax.set_xlim(-0.1,5.1)
ax.set_ylim(1e4, 2*ylim)
ax.set_yscale('log')
ax.legend(loc=3, bbox_to_anchor=(0.0,1.15,1.0,0.1), borderaxespad=0., ncol=3, numpoints=1, fancybox=True, mode='expand', fontsize=12)
fig.savefig(os.path.join(plotdir,'fig9b.pdf'), bbox_inches='tight')


# plot distance evolution
momenta_dd = {}
for dataset in datasets:
    momenta_dd[dataset['line']] = np.genfromtxt(os.path.join(projectdir, '08.results', dataset['cube']+'.major.5.0s.momentum_dist_bootstrap.real_outflow.txt'),
                                              names = ('bin','q','q_p','q_m')
                                             )
fig = plt.figure(figsize=(8,4))
ax  = fig.subplots(1)
ylim = 0.0
lowlim = 33.9369576
for dataset,lim in zip(datasets,[509.054364,322.4010972,203.6217456]):
    bin    = momenta_dd[dataset['line']]['bin']*16.9684788
    momentum   = momenta_dd[dataset['line']]['q']
    momentum_p = momenta_dd[dataset['line']]['q_p']
    momentum_m = momenta_dd[dataset['line']]['q_m']
    momentum_m = [r if r>0. else 0. for r in momentum_m]
    ax.plot(bin[bin<=lowlim], momentum[bin<=lowlim],
            lw=3, ls=':',
            color=dataset['plot_color']
           )
    ax.plot(bin[(bin<=lim) & (bin>=lowlim)], momentum[(bin<=lim) & (bin>=lowlim)],
            lw=3, ls='-',
            color=dataset['plot_color'],
            label=dataset['line'].replace('_','(').replace('-','--')+')'
           )
    ax.plot(bin[bin>=lim], momentum[bin>=lim],
            lw=3, ls=':',
            color=dataset['plot_color']
           )
    ax.fill_between(bin, momentum_p, momentum_m,
                           color=dataset['plot_color'],
                           alpha=0.25
                          )
    if (ylim<np.nanmax([momentum,momentum_p,momentum_m])): ylim = np.nanmax([momentum,momentum_p,momentum_m])
ax.set_axisbelow(True)
ax.grid(True)
ax.set_xlabel(r'd [pc]')
ax.set_ylabel(r'P [M$_\odot$\,km\,s$^{-1}$]')
ax.yaxis.set_label_position('right')
ax.tick_params('y', left=False, labelleft=False, right=True, labelright=True)
ax.set_xlim(-10,510)
ax.set_ylim(1e6, 2*ylim)
ax.set_yscale('log')
fig.savefig(os.path.join(plotdir,'fig9d.pdf'), bbox_inches='tight')


###################################################################################################

# Appendix A: model
###################

easy_aplpy.plot.map(os.path.join(paperdir,'NGC253.CO_1-0.mom1.fits'),
    figsize   = (16,16),
    vmin    = 100,
    vmax    = 400,
    stretch = 'linear',
    cmap    = 'RdBu',
    beam    = None,
    recenter = [SkyCoord('00h47m33.01s -25d17m17.0s'), 80.0*u.arcsec, 80.0*u.arcsec],
    contours = [[os.path.join(paperdir,'NGC253.CO_1-0.mom1.fits'), np.arange(100,450,50), 'black']],
    colorbar = ['right', 'km\,s$^{-1}$'],
    texts     = [[[0.0,1.0], 'data', {'size':60, 'ha':'left', 'va':'top', 'style':'oblique', 'weight':'black', 'color':'black', 'bbox': easy_aplpy.settings.props}]],
    out     = os.path.join(plotdir,'figA1a.pdf')
    )
easy_aplpy.plot.map(os.path.join(paperdir,'diskfit.total_model.fits'),
    figsize   = (16,16),
    vmin    = 100,
    vmax    = 400,
    stretch = 'linear',
    cmap    = 'RdBu',
    beam    = None,
    recenter = [SkyCoord('00h47m33.01s -25d17m17.0s'), 80.0*u.arcsec, 80.0*u.arcsec],
    contours = [[os.path.join(paperdir,'diskfit.total_model.fits'), np.arange(100,450,50), 'black']],
    colorbar = ['right', 'km\,s$^{-1}$'],
    texts     = [[[0.0,1.0], 'model', {'size':60, 'ha':'left', 'va':'top', 'style':'oblique', 'weight':'black', 'color':'black', 'bbox': easy_aplpy.settings.props}]],
    out     = os.path.join(plotdir,'figA1b.pdf')
    )
easy_aplpy.plot.map(os.path.join(paperdir,'NGC253.CO_1-0.mom1.fits'),
    figsize   = (16,16),
    vmin    = 100,
    vmax    = 400,
    stretch = 'linear',
    cmap    = 'RdBu',
    beam    = None,
    recenter = [SkyCoord('00h47m33.01s -25d17m17.0s'), 80.0*u.arcsec, 80.0*u.arcsec],
    contours = [[os.path.join(paperdir,'diskfit.total_model.fits'), np.arange(100,450,50), 'black']],
    colorbar = ['right', 'km\,s$^{-1}$'],
    texts     = [[[0.0,1.0], 'data + model contours', {'size':60, 'ha':'left', 'va':'top', 'style':'oblique', 'weight':'black', 'color':'black', 'bbox': easy_aplpy.settings.props}]],
    out     = os.path.join(plotdir,'figA1c.pdf')
    )
easy_aplpy.plot.map(os.path.join(paperdir,'diskfit.total_residual.fits'),
    figsize   = (16,16),
    vmin    = -150,
    vmax    = 150,
    stretch = 'linear',
    cmap    = 'RdBu',
    beam    = None,
    recenter = [SkyCoord('00h47m33.01s -25d17m17.0s'), 80.0*u.arcsec, 80.0*u.arcsec],
    contours = [[os.path.join(paperdir,'diskfit.total_residual.fits'), np.arange(-50,75,25), 'black']],
    colorbar = ['right', 'km\,s$^{-1}$'],
    texts     = [[[0.0,1.0], 'residual', {'size':60, 'ha':'left', 'va':'top', 'style':'oblique', 'weight':'black', 'color':'black', 'bbox': easy_aplpy.settings.props}]],
    out     = os.path.join(plotdir,'figA1d.pdf')
    )


###################################################################################################

# Appendix B: all pV slices
###########################

easy_aplpy.settings.colorbar_label_fontsize = 20      # unit: point
easy_aplpy.settings.colorbar_ticks_fontsize = 20      # unit: point
easy_aplpy.settings.grid_label_fontsize     = 20      # unit: point
easy_aplpy.settings.tick_label_fontsize     = 20      # unit: point
easy_aplpy.settings.axis_label_fontsize     = 20      # unit: point
easy_aplpy.settings.colorbar_width          = 0.08    # relative to panel size

execfile('scripts/NGC253/410.disk_non-disk/400.disk_non-disk.helper_functions.py')
pvdir        = 'NGC253/420.disk_non-disk/04.pVs/'
pV_path      = 'NGC253/410.disk_non-disk/04.pVs/'
kin_path     = 'NGC253/410.disk_non-disk/05.kinematics/'
total_model  = 'NGC253/410.disk_non-disk/07.model/diskfit.total_model.fits'
tmodel_major = 'NGC253/410.disk_non-disk/07.model/diskfit.total_model.major_slice_'
tmodel_minor = 'NGC253/410.disk_non-disk/07.model/diskfit.total_model.minor_slice_'

contours = []
labels   = []
lines    = []
for a_slice in major_slices:
    for CO in ['CO_1-0','CO_2-1','CO_3-2']:
        c = []
        pV       = os.path.join(pvdir, 'NGC253.'+CO+'.pV_major.slice_major_'+str(a_slice['slice'])+'.fits')
        var_mask = os.path.join(pV_path,'NGC253.CO.pV_major.slice_major_'+str(a_slice['slice'])+'.mask_model_variable.fits')
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
labels[0] = [[[0.5, 0.95], 'CO (1-0)', {'size': 24, 'bbox': easy_aplpy.settings.props, 'ha': 'center', 'va': 'top'}]]
labels[1] = [[[0.5, 0.95], 'CO (2-1)', {'size': 24, 'bbox': easy_aplpy.settings.props, 'ha': 'center', 'va': 'top'}]]
labels[2] = [[[0.5, 0.95], 'CO (3-2)', {'size': 24, 'bbox': easy_aplpy.settings.props, 'ha': 'center', 'va': 'top'}]]

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
    out      = os.path.join(plotdir,'figB1a.pdf')
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
    out      = os.path.join(plotdir,'figB1b.pdf')
    )

easy_aplpy.settings.colorbar_label_fontsize = 36      # unit: point
easy_aplpy.settings.colorbar_ticks_fontsize = 36      # unit: point
easy_aplpy.settings.grid_label_fontsize     = 36      # unit: point
easy_aplpy.settings.tick_label_fontsize     = 36      # unit: point
easy_aplpy.settings.axis_label_fontsize     = 36      # unit: point
easy_aplpy.settings.colorbar_width          = 0.15    # relative to panel size



###################################################################################################

# crop images
#############

# for all images in plotdir ...
for image in os.listdir(plotdir):

    # ... trim to cut off unnecessary white space
    os.system('pdfcrop '+os.path.join(plotdir,image))
    os.system('rm -rf '+os.path.join(plotdir,image))
    os.system('mv '+' '+os.path.join(plotdir,image[:-4]+'-crop.pdf')+' '+os.path.join(plotdir,image))

    # ... make rastered png image for manual adjustment of labels
    #os.system('/usr/bin/convert -density 300 '+plotdir+image+' '+plotdir+image[:-3]+'png')

    # compress images to save space
    os.system('convert '+os.path.join(plotdir,image)+' -compress Zip '+os.path.join(plotdir,image))


###################################################################################################
