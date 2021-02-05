##########################
# PYTHON PLOTTING SCRIPT #
##########################

# Plot the outflow rate evolution in time and distance.

###################################################################################################

execfile('scripts/plotting_imports.py')
execfile('NGC253/project_info.py')
execfile(os.path.join(projectdir, 'info.py'))

ratedir = ratedir+'.improved3'

###################################################################################################

# load data
###########

def load_rates(rate_type):
    rates = {}
    for dataset in datasets:
        rates[dataset['line']] = {}
        for ax_type in ['major']: #,'minor']:
            rates[dataset['line']][ax_type] = {}
            for SNR in [5.0]: #SNRs:
                rates[dataset['line']][ax_type][SNR] = {}
                for mask in ['','real_outflow','superbubble','cospatial_disk']:

                    try:
                        rates[dataset['line']][ax_type][SNR][mask] = np.genfromtxt(os.path.join(ratedir, dataset['cube']+'.'+ax_type+'.'+str(SNR)+'s.'+rate_type+'.'+mask+'.txt'),
                                          names = ('bin','rate','rate_p','rate_m')
                                         )
                    except:
                        rates[dataset['line']][ax_type][SNR][mask] = []
    return rates


###################################################################################################

# plot rate evolution
#####################

def plot_rate_time(rates, savepath):

    # set up figure
    fig = plt.figure(figsize=(8,16))
    gs  = GridSpec(4,1)
    gs.update(left=0.1, right=0.95, top=0.85, bottom=0.1, hspace=0.0, wspace=0.0)
    axes    = [[],[],[],[]]
    axes[0] = plt.subplot(gs[0, 0:])
    axes[1] = plt.subplot(gs[1, 0:], sharex=axes[0])
    axes[2] = plt.subplot(gs[2, 0:], sharex=axes[0])
    axes[3] = plt.subplot(gs[3, 0:], sharex=axes[0])


    # plot data
    ax_type = 'major'
    SNR     = 5.0
    ylims   = {mask: 0.0 for mask in rates[datasets[0]['line']][ax_type][SNR].keys()}

    for idx,mask in enumerate(rates[datasets[0]['line']][ax_type][SNR].keys()):
        axes[idx].text(0.95,0.9,
                       mask.replace('_',' '),
                       fontsize=16,
                       va='top', ha='right', transform=axes[idx].transAxes, bbox=props
                      )

        text = r'$\int \dot{\mathrm{M}} \, \mathrm{d}t$'

        for dataset in datasets:
            try:
                bin    = rates[dataset['line']][ax_type][SNR][mask]['bin']
                rate   = rates[dataset['line']][ax_type][SNR][mask]['rate']
                rate_p = rates[dataset['line']][ax_type][SNR][mask]['rate_p']
                rate_m = rates[dataset['line']][ax_type][SNR][mask]['rate_m']

                ymax = np.nanmax([rate,rate_p,rate_m])
                if ( ylims[mask] < ymax ):
                    ylims[mask] = ymax

                # plot solid inside trusted range and dashed outside the trusted range
                axes[idx].plot(bin[bin<=2.0], rate[bin<=2.0],
                               lw=3, ls='-',
                               color=dataset['plot_color'],
                               label=dataset['line'].replace('_','(').replace('-','--')+')'
                               )
                axes[idx].plot(bin[bin>=2.0], rate[bin>=2.0],
                               lw=3, ls='--',
                               color=dataset['plot_color'],
                               label=dataset['line'].replace('_','(').replace('-','--')+')'
                               )
                axes[idx].fill_between(bin, rate_p, rate_m,
                                       color=dataset['plot_color'],
                                       alpha=0.25
                                      )

                text += '\n$'+('{:.2e}'.format(np.nansum(rate)*(bin[1]-bin[0])*1e6)).replace('e+0',r'\times 10^')+'\,\mathrm{M}_\odot$'
            except:
                pass
        axes[idx].text(0.95,0.75,
                       text,
                       fontsize=16,
                       va='top', ha='right', transform=axes[idx].transAxes, bbox=props
                      )

    # format figure
    axes[0].xaxis.set_label_position('top')
    axes[0].set_xlabel(r't$_\mathrm{eject}$ [Myr]')
    axes[3].set_xlabel(r't$_\mathrm{eject}$ [Myr]')
    axes[0].tick_params(labelbottom=False, labeltop=True)
    axes[1].tick_params(labelbottom=False)
    axes[2].tick_params(labelbottom=False)
    for ax,ylim in zip(axes,ylims.values()):
        ax.xaxis.set_ticks_position('both')
        ax.set_ylabel(r'$\dot{\mathrm{M}}_{out}$ [M$_\odot$\,yr$^{-1}$]')
        ax.grid(True)
        ax.set_xlim(-0.1,3.1)
        ax.set_ylim(-0.1*ylim,1.1*ylim)

    axes[0].legend(loc=3, bbox_to_anchor=(0.0,1.15,1.0,0.1), borderaxespad=0., ncol=3, numpoints=1, fancybox=True, mode='expand', fontsize=16)

    fig.savefig(savepath, bbox_inches='tight')


###################################################################################################

# plot rate evolution
#####################

def plot_rate_distance(rates, savepath):

    # set up figure
    fig = plt.figure(figsize=(8,16))
    gs  = GridSpec(4,1)
    gs.update(left=0.1, right=0.95, top=0.85, bottom=0.1, hspace=0.0, wspace=0.0)
    axes    = [[],[],[],[]]
    axes[0] = plt.subplot(gs[0, 0:])
    axes[1] = plt.subplot(gs[1, 0:], sharex=axes[0])
    axes[2] = plt.subplot(gs[2, 0:], sharex=axes[0])
    axes[3] = plt.subplot(gs[3, 0:], sharex=axes[0])


    # plot data
    ax_type = 'major'
    SNR     = 5.0
    ylims   = {mask: 0.0 for mask in rates[datasets[0]['line']][ax_type][SNR].keys()}

    for idx,mask in enumerate(rates[datasets[0]['line']][ax_type][SNR].keys()):
        axes[idx].text(0.95,0.9,
                       mask.replace('_',' '),
                       fontsize=16,
                       va='top', ha='right', transform=axes[idx].transAxes, bbox=props
                      )

        for dataset in datasets:
            try:
                bin    = rates[dataset['line']][ax_type][SNR][mask]['bin']
                rate   = rates[dataset['line']][ax_type][SNR][mask]['rate']
                rate_p = rates[dataset['line']][ax_type][SNR][mask]['rate_p']
                rate_m = rates[dataset['line']][ax_type][SNR][mask]['rate_m']

                ymax = np.nanmax([rate,rate_p,rate_m])
                if ( ylims[mask] < ymax ):
                    ylims[mask] = ymax

                # plot solid inside trusted range and dashed outside the trusted range
                axes[idx].plot(bin, rate,
                               lw=3, ls='-',
                               color=dataset['plot_color'],
                               label=dataset['line'].replace('_','(').replace('-','--')+')'
                               )
                axes[idx].fill_between(bin, rate_p, rate_m,
                                       color=dataset['plot_color'],
                                       alpha=0.25
                                      )
            except:
                pass

    # format figure
    axes[0].xaxis.set_label_position('top')
    axes[0].set_xlabel(r'd [$^{\prime\prime}$]')
    axes[3].set_xlabel(r'd [$^{\prime\prime}$]')
    axes[0].tick_params(labelbottom=False, labeltop=True)
    axes[1].tick_params(labelbottom=False)
    axes[2].tick_params(labelbottom=False)
    for ax,ylim in zip(axes,ylims.values()):
        ax.xaxis.set_ticks_position('both')
        ax.set_ylabel(r'$\dot{\mathrm{M}}_{out}$ [M$_\odot$\,yr$^{-1}$]')
        ax.grid(True)
        ax.set_xlim(-0.6,30.6)
        ax.set_ylim(-0.1*ylim,1.1*ylim)

    axes[0].legend(loc=3, bbox_to_anchor=(0.0,1.15,1.0,0.1), borderaxespad=0., ncol=3, numpoints=1, fancybox=True, mode='expand', fontsize=16)

    fig.savefig(savepath, bbox_inches='tight')


###################################################################################################

# plot bootstrap rates
######################

def plot_bootstrap_time(maskfile, savepath):

    # load data
    rates = []
    for dataset in datasets:
        with open(os.path.join(ratedir, dataset['cube']+'.major.5.0s.rate_bootstrap.'+maskfile+'.pickle'), "rb") as f:
            rates.append( pickle.load(f, encoding='latin1') )

    # set up figure
    fig = plt.figure(figsize=(8,12))
    gs  = GridSpec(3,1)
    gs.update(left=0.1, right=0.95, top=0.85, bottom=0.1, hspace=0.0, wspace=0.0)
    axes    = [[],[],[]]
    axes[0] = plt.subplot(gs[0, 0:])
    axes[1] = plt.subplot(gs[1, 0:], sharex=axes[0])
    axes[2] = plt.subplot(gs[2, 0:], sharex=axes[0])

    # plot data
    ax_type  = 'major'
    SNR      = 5.0
    binwidth = 0.1

    # plot all inclinations colorcoded
    for idx,dataset in enumerate(datasets):
        inclinations = np.arange(48,108.01,1.0)

        axes[idx].text(0.95,0.9,
                       dataset['line'].replace('_','(').replace('-','--')+')',
                       fontsize=16,
                       va='top', ha='right', transform=axes[idx].transAxes, bbox=props
                      )

        for j in inclinations:

                j_bin  =        rates[idx][:,0][(rates[idx][:,1] == j) & (rates[idx][:,1] == j)]
                j_rate = np.abs(rates[idx][:,2][(rates[idx][:,1] == j) & (rates[idx][:,1] == j)])

                axes[idx].plot(j_bin, j_rate,
                               lw=1, ls='-',
                               color = mpl.cm.viridis(np.linspace(0,1,len(inclinations)))[list(inclinations).index(j)],
                               label = str(j)
                               )

        # plot 16, 50, 84 percentiles
        perc = []
        hbw = binwidth/2.
        for b in np.arange(0,10,binwidth):
            bin_rates = rates[idx][:,2][(rates[idx][:,0] > b-hbw) & (rates[idx][:,0] < b+hbw) & (rates[idx][:,2] > 0.)]
            perc.append([b,np.percentile(bin_rates,16),np.percentile(bin_rates, 50),np.percentile(bin_rates, 84)])
        perc = np.array(perc)
        for i in [1,2,3]:
            axes[idx].plot(perc[:,0], perc[:,i],
                           lw=3, ls='-',
                           color = 'r'
                           )

    # format figure
    axes[0].xaxis.set_label_position('top')
    axes[0].set_xlabel(r't$_\mathrm{eject}$ [Myr]')
    axes[2].set_xlabel(r't$_\mathrm{eject}$ [Myr]')
    axes[0].tick_params(labelbottom=False, labeltop=True)
    axes[1].tick_params(labelbottom=False)
    for ax in axes:
        ax.xaxis.set_ticks_position('both')
        ax.set_ylabel(r'$\dot{\mathrm{M}}_{out}$ [M$_\odot$\,yr$^{-1}$]')
        ax.grid(True)
        ax.set_xlim(-0.1,3.1)
        ax.set_ylim(-5,100)
    axes[0].legend(loc=3, bbox_to_anchor=(0.0,1.15,1.0,0.1), borderaxespad=0., ncol=10, numpoints=1, fancybox=True, mode='expand', fontsize=16)

    fig.savefig(savepath, bbox_inches='tight')


###################################################################################################

# plot bootstrap rates
######################

def plot_bootstrap_dist(maskfile, savepath):

    # load data
    rates = []
    for dataset in datasets:
        with open(os.path.join(ratedir, dataset['cube']+'.major.5.0s.dist_bootstrap.'+maskfile+'.pickle'), "rb") as f:
            rates.append( pickle.load(f, encoding='latin1') )

    # set up figure
    fig = plt.figure(figsize=(8,12))
    gs  = GridSpec(3,1)
    gs.update(left=0.1, right=0.95, top=0.85, bottom=0.1, hspace=0.0, wspace=0.0)
    axes    = [[],[],[]]
    axes[0] = plt.subplot(gs[0, 0:])
    axes[1] = plt.subplot(gs[1, 0:], sharex=axes[0])
    axes[2] = plt.subplot(gs[2, 0:], sharex=axes[0])

    # plot data
    ax_type  = 'major'
    SNR      = 5.0
    binwidth = 1.0

    # plot all inclinations colorcoded
    for idx,dataset in enumerate(datasets):
        inclinations = np.arange(48,108.01,1.0)

        axes[idx].text(0.95,0.9,
                       dataset['line'].replace('_','(').replace('-','--')+')',
                       fontsize=16,
                       va='top', ha='right', transform=axes[idx].transAxes, bbox=props
                      )

        for j in inclinations:

                j_bin  =        rates[idx][:,0][(rates[idx][:,1] == j) & (rates[idx][:,1] == j)]
                j_rate = np.abs(rates[idx][:,2][(rates[idx][:,1] == j) & (rates[idx][:,1] == j)])

                axes[idx].plot(j_bin, j_rate,
                               lw=1, ls='-',
                               color = mpl.cm.viridis(np.linspace(0,1,len(inclinations)))[list(inclinations).index(j)],
                               label = str(j)
                               )

        # plot 16, 50, 84 percentiles
        perc = []
        hbw = binwidth/2.
        for b in np.arange(0,40,binwidth):
            bin_rates = rates[idx][:,2][(rates[idx][:,0] > b-hbw) & (rates[idx][:,0] < b+hbw) & (rates[idx][:,2] > 0.)]
            if not ( bin_rates.size == 0 ):
                perc.append([b,np.percentile(bin_rates,16),np.percentile(bin_rates, 50),np.percentile(bin_rates, 84)])
        perc = np.array(perc)
        for i in [1,2,3]:
            axes[idx].plot(perc[:,0], perc[:,i],
                           lw=3, ls='-',
                           color = 'r'
                           )

    # format figure
    axes[0].xaxis.set_label_position('top')
    axes[0].set_xlabel(r'd [$^{\prime\prime}$]')
    axes[2].set_xlabel(r'd [$^{\prime\prime}$]')
    axes[0].tick_params(labelbottom=False, labeltop=True)
    axes[1].tick_params(labelbottom=False)
    for ax in axes:
        ax.xaxis.set_ticks_position('both')
        ax.set_ylabel(r'$\dot{\mathrm{M}}_{out}$ [M$_\odot$\,yr$^{-1}$]')
        ax.grid(True)
        ax.set_xlim(-0.6,30.6)
        ax.set_ylim(-5,100)
    axes[0].legend(loc=3, bbox_to_anchor=(0.0,1.15,1.0,0.1), borderaxespad=0., ncol=10, numpoints=1, fancybox=True, mode='expand', fontsize=16)

    fig.savefig(savepath, bbox_inches='tight')


###################################################################################################
###################################################################################################

# carry out plotting
####################

rates_tp = load_rates('rate')
rates_td = load_rates('rate.deprojected')
rates_tboot = load_rates('rate_bootstrap')
rates_dp = load_rates('distance_evolution')
rates_dd = load_rates('distance_evolution.deprojected')
rates_dboot = load_rates('dist_bootstrap')

plot_rate_time(    rates_tp,    os.path.join(plotdir,'rate_function.pdf'))
plot_rate_time(    rates_td,    os.path.join(plotdir,'rate_function.deprojected.pdf'))
plot_rate_time( rates_tboot,    os.path.join(plotdir,'rate_function_bootstrap.pdf'))
plot_rate_distance(rates_dp,    os.path.join(plotdir,'distance_function.pdf'))
plot_rate_distance(rates_dd,    os.path.join(plotdir,'distance_function.deprojected.pdf'))
plot_rate_distance(rates_dboot, os.path.join(plotdir,'distance_function_bootstrap.pdf'))
plot_bootstrap_time('real_outflow', os.path.join(plotdir,'bootstrap_distribution_time.pdf'))
plot_bootstrap_dist('real_outflow', os.path.join(plotdir,'bootstrap_distribution_dist.pdf'))


###################################################################################################
###################################################################################################
