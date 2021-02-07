'''
Some Plots for the SNOWPACK data...
Based on *.nc files generated with pro2netcdf.py
Author: @Hannes Hohenwarter
'''

import warnings
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
import xarray as xr
import pandas as pd
import numpy as np
import os
import SnowTools as st
warnings.simplefilter(action='ignore', category=FutureWarning)


def vert_profile(dds, var, profile, type, mrun='multi', time=None, freq='6H'):
    ''' Plot vert. Profiles of Snowpack at time= 'time' or at intervals with freq= 'freq'
    :param dds: data (xarray)
    :param var: Variable to Print
    :param profile: name of profile (str)
    :param type: name of type (str)
    :param mrun: choose specific model run by number or multible runs 'multi'
    :param time: print only one timestemp e.g. 'YYYY-MM-DDThh:mm:ss' (string) or tuple (start, end)
    :param freq: print profile every e.g. 6 hours '6H' (string)
    '''

    if mrun == 'multi':
        ds = dds
    elif dds.run.isin([mrun]).any():  # check if mrun is in the range of ds.run
        ds = dds.sel(run=mrun)
    else:
        print('mrun must be "multi" or within the range of' + str(dds.run.values))

    if time is None:  # take whole periode of *.smet file
        start = '2016-01-31T22:00:00'
        end = '2016-02-02T22:00:00'
        #  create time vector with 'freq' intervals
        time = pd.date_range(start=start, end=end, freq=freq)
    elif len(time) == 2:  # use defined start and end time
        start = time[0]
        end = time[1]
        time = pd.date_range(start=start, end=end, freq=freq)

    #  calculate var if necessary and add to Dataset
    if var == 'Sn_wet':
        ds[var] = st.Sn_wet(ds)
        ds.attrs[var] = ''
    elif var == 'Sn38_hand':
        ds[var] = st.Sn38(ds)
        ds.attrs[var] = ''
    elif var == 'Sn_wet_manu':
        ds[var] = st.Sn_wet(ds, kind='manu')
        ds.attrs[var] = ''
    elif var == 'sigma_wet':
        ds[var] = ds['shear_strength'] * np.exp(-0.235 * ds.lwc)
        ds.attrs[var] = 'kPa'

    (row, col) = sub_shape(len(time))  # get row, col for subplot
    max_height = float(ds.height.max())  # maximum snow height
    wl_bottom = 0
    wl_top = 0

    fig, axs = plt.subplots(row, col)
    fig.set_size_inches(8, 10)
    fig.set_tight_layout(True)
    colors = plt.cm.winter(np.linspace(0, 1, len(ds.run)))
    p = 0
    # plot..
    for i in np.arange(row):
        for j in np.arange(col):
            if p < len(time):
                df = ds.sel(date=time[p])
                #max = ds.lwc.sel(date=time).max()
                for r in np.arange(len(ds.run)):
                    ddf = df.sel(run=r+1).to_dataframe()
                    ddf.plot.line(ax=axs[i, j], x=var, y='height', c=colors[r], linewidth=0.6)

                    #  plot mean_LWC
                    if var == 'lwc':
                        axs[i, j].vlines(ddf.lwc.mean(), 0, max_height, colors=colors[r],
                                         linestyles='dashed', linewidth=0.5, alpha=0.4)

                        if r == 0:
                            ax = axs[i, j].inset_axes([0.5, 0.6, 0.45, 0.3])
                            max_bar(df, ax, colors, max=35)

                #  Plot critical value and ax Limits for variables
                if var == 'lwc':
                    axs[i, j].vlines(6, 0, max_height)
                    axs[i, j].vlines(3, 0, max_height, linewidth=0.5)
                    plt.setp(axs, xlim=(-0.1, 30))
                    plt.setp(axs, xticks=[0, 10, 20, 30])
                elif var in ['Sn38', 'Sk38']:
                    axs[i, j].vlines(1, 0, max_height)
                    plt.setp(axs, xlim=(0, 6.5))
                elif var == 'shear_strength':
                    df_orig = df.sel(run=st.orig_run(profile, type)).to_dataframe()
                    df_orig['sigma'] = st.shear_str(df_orig)
                    axs[i, j].vlines(1, 0, max_height)
                    df_orig.plot(ax=axs[i, j], x='sigma', y='height', c='black',
                                 linewidth=0.6, linestyle='dashed')
                    plt.setp(axs, xlim=(0, 10))
                elif var == 'sigma_wet':
                    df_orig = df.sel(run=st.orig_run(profile, type)).to_dataframe()
                    axs[i, j].vlines(1, 0, max_height)
                    df_orig.plot(ax=axs[i, j], x='shear_strength', y='height', c='black',
                                 linewidth=0.6, linestyle='dashed')
                    plt.setp(axs, xlim=(0, 10))
                elif var in ['Sn_wet', 'Sn_wet_manu']:
                    axs[i, j].vlines(1, 0, max_height)
                    plt.setp(axs, xlim=(0, 10))
                    df_orig = df.sel(run=st.orig_run(profile, type)).to_dataframe()
                    if var == 'Sn_wet_manu':
                        df_orig['Sn_dry'] = st.Sn_dry(df_orig, kind='manu')
                    else:
                        df_orig['Sn_dry'] = st.Sn_dry(df_orig)
                    df_orig.plot(ax=axs[i, j], x='Sn_dry', y='height', c='black',
                                 linewidth=0.6, linestyle='dashed')

                #  Plot height of WL
                if profile == 'SchichtW_R':
                    wl_bottom = 0
                    wl_top = 20
                elif profile in ['MeltCrust_R', 'SurfHoar_R', 'SurfHoaR_hr_R']:
                    wl_bottom = 20
                    wl_top = 20.5
                    if type == 'dWL':
                        left, right = plt.xlim()
                        axs[i, j].hlines(17.5, left, right, colors='gray', linestyles='dashed',
                                         linewidth=0.5, alpha=0.4)
                elif profile == 'DepthHoar_R':
                    wl_bottom = 0
                    wl_top = 3
                axs[i, j].axhspan(wl_bottom, wl_top, facecolor='gray', alpha=0.4)

                #  plot time and remove legend and xlabels
                axs[i, j].legend().remove()
                axs[i, j].text(.98, .02, time[p].strftime('%H:%M'), horizontalalignment='right',
                               transform=axs[i, j].transAxes, fontsize=10)
                axs[i, j].set_xlabel('', visible=False)

                #  get rid of tick labels at middle plots
                if i != row-1:
                    plt.setp(axs[i, j].get_xticklabels(), visible=False)
                if j != 0:
                    plt.setp(axs[i, j].get_yticklabels(), visible=False)

            else:
                # remove empty subplots
                axs[i, j].set_axis_off()
            p += 1


    plt.setp(axs, ylim=(0, np.ceil(max_height)), yticks=np.arange(0, max_height, step=10))
    y_lab = fig.text(-0.01, .5, 'Snow Height (cm)', horizontalalignment='left', verticalalignment='center',
                     rotation=90, fontsize=12)
    x_lab = fig.text(0.5, -0.01, var_label[var], horizontalalignment='center', fontsize=12)

    if type == 'roh':
        fig.text(1, 0.60, r'$\frac{\rho_{WL}}{\rho_{Slab}} = $' + str(st.density_rel[profile]),
                 horizontalalignment='left', fontsize=10)
    handles, labels = axs[0, 0].get_legend_handles_labels()
    labels = getattr(st, type)[profile[:-2]]
    lgd = fig.legend(handles, labels.round(decimals=2).astype(str), loc='center left', title=legend_title[type],
                     bbox_to_anchor=(1, 0.5))

    tit = fig.suptitle(var_title[var] + '\n vertical profiles from ' + profile[:-2] + ' changing ' +
                       type_title[type], y=1.05, fontsize=15)

    return fig, lgd, tit, x_lab, y_lab


def snowpack(ds, var):
    #  querschniit der schneedecke
    #  var als color über dz als stacked bar plot (width=1) !!!

    if var == 'Sn_wet':
        import SnowTools as st
        ds[var] = st.Sn_wet(ds)
    elif var == 'Sn38_hand':
        import SnowTools as st
        ds[var] = st.Sn38(ds)

    data = ds[var].values.transpose()
    norm = plt.Normalize(vmin=np.nanmin(data), vmax=np.nanmax(data)) 
    my_cmap = my_colormap(var)
    colors = my_cmap(norm(data))  
    dz = ds.dz.to_pandas()
    fig, ax = plt.subplots()
    dz.plot.bar(ax=ax, stacked=True, width=1, color=colors, legend=None)
    sm = ScalarMappable(cmap=my_cmap, norm=norm)
    sm.set_array(np.array([]))
    cbar = plt.colorbar(sm)

    newday = dz.index.strftime('%H:%M') == '00:00'
    noon = dz.index.strftime('%H:%M') == '12:00'
    ticklabels = dz.index.strftime('%b %d %Y \n%H:%M')
    minorlabels = dz.index.strftime('%H:%M')
    ticklabels[noon] = minorlabels[noon]
    tickloc = np.arange(len(ds.date)) + 1
    ax.xaxis.set_major_locator(ticker.FixedLocator(tickloc[newday]))
    ax.xaxis.set_minor_locator(ticker.FixedLocator(tickloc[noon]))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter(ticklabels[newday]))
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(ticklabels[noon]))
    for label in ax.get_xticklabels():
        label.set_rotation(0)
    plt.show()


def max_timeline(ds, var, profile, type, value='max'):

    if var == 'Sn_wet':
        ds[var] = st.Sn_wet(ds)
        ds.attrs[var] = ''

    if value == 'max':
        dmax = ds.max(dim='elem')
    elif value == 'min':
        dmax = ds.min(dim='elem')
    elif value == 'mean':
        dmax = ds.mean(dim='elem')

    colors = plt.cm.winter(np.linspace(0, 1, len(dmax.run)))
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 5)
    fig.set_tight_layout(True)
    lines = [None] * len(dmax.run)
    for r in dmax.run:
        lines[int(r)-1] = dmax[var].sel(run=int(r)).plot.line(x='date', c=colors[int(r)-1], linewidth=0.6)

    if value == 'max':
        plt.ylim(0, 35)
        plt.axhline(y=6, color='k', linestyle='dashed')  # hline at lwc=6
        tit = plt.title(max_title[var] + ' maximum of ' + profile[:-2] + ' changing ' + type_title[type], fontsize=15)
    elif value == 'min':
        plt.ylim(0, 15)
        plt.axhline(y=1, color='k', linestyle='dashed')  # hline at Sn=1
        tit = plt.title('Minimum of ' + max_title[var] + ' from ' + profile[:-2] + ' changing ' + type_title[type],
                        fontsize=15)
    elif value == 'mean':
        plt.ylim(0, 15)
        plt.axhline(y=3, color='k', linestyle='dashed')  # hline at lwc=3
        tit = plt.title('Mean of ' + max_title[var] + ' from ' + profile[:-2] + ' changing ' + type_title[type],
                        fontsize=15)

    axs = ax.inset_axes([0.75, 0.65, 0.2, 0.2])
    max_bar(dmax, axs, colors, value, dim='date', fontsize=6)

    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d \n%H:%M'))
    #ax.xaxis.set_minor_formatter(mdates.DateFormatter('%H:%M'))
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    _ = plt.xticks(rotation=0, horizontalalignment='center')

    plt.autoscale(enable=True, axis='x', tight=True)
    plt.ylabel(var_label[var], fontsize=12)
    plt.xlabel('Date', fontsize=12)
    labels = getattr(st, type)[profile[:-2]]
    lgd = plt.legend(labels.round(decimals=2).astype(str), loc=2, title=legend_title[type])
    if type == 'roh':
        fig.text(0.1, 0.45, r'$\frac{\rho_{WL}}{\rho_{Slab}} = $' + str(st.density_rel[profile]),
                     horizontalalignment='left', fontsize=10)

    return fig, lgd, tit


def sub_shape(n):
    '''
    :param n: length of array
    :return: #row and #col of squared Matrix to fit all n elements
    e.g: n=7 -> row=3, col=3
    '''
    # reshape to closest to square matrix
    col = np.ceil(np.sqrt(n)).astype(int)
    row = np.ceil(n/col).astype(int)

    return row, col


def my_colormap(var):
    #define colormaps for all variables
    if var in ['lwc', 'Sn_wet']:
        return plt.cm.get_cmap('Blues')
    else:
        return plt.cm.get_cmap('Jet')
        # return defalt colormap for everything else


def max_bar(df, ax, colors, value='max', max=35, dim='elem', fontsize=5):

    if value in ['max', 'mean']:
        maxima = df.max(dim=dim).to_dataframe()
        maxima.plot.bar(ax=ax, y='lwc', width=0.9, color=colors, legend=None)
        plt.setp(ax, ylim=(0, max + 1))
        [ax.text(i, v, ' {:.1f}'.format(v), fontsize=fontsize, color='black', horizontalalignment='center',
                 rotation=90) for i, v in enumerate(maxima.lwc)]
    elif value == 'min':
        minima = df.min(dim=dim).to_dataframe()
        minima.plot.bar(ax=ax, y='Sn_wet', width=0.9, color=colors, legend=None)
        plt.setp(ax, ylim=(0, 6))
        [ax.text(i, v, ' {:.1f}'.format(v), fontsize=fontsize, color='black', horizontalalignment='center',
                 rotation=90) for i, v in enumerate(minima.Sn_wet)]

    ax.text(.5, -0.1, value + ' values', fontsize=7, horizontalalignment='center', transform=ax.transAxes)
    ax.patch.set_alpha(0)
    ax.axis('off')


def wever_alpha_n():

    profiles = ['LayerCh', 'SurfHoar']
    marker = {'LayerCh': 'o', 'SurfHoar': 'x'}

    fig, axs = plt.subplots()
    for profile in profiles:
        alpha = 7.3 * st.kWL[profile] + 1.9
        n = 15.68 * np.exp(-0.46 * st.kWL[profile]) + 1
        axs.scatter(st.kWL[profile], alpha, marker=marker[profile], color='tab:blue')
        axs.scatter(st.kWL[profile], n, marker=marker[profile], color='tab:orange')
    axs.set_xlabel('grain diameter')
    axs.set_title('Parameterizations')
    axs.legend([r'$\alpha$ (LayerCh)', 'n (LayerCh)', r'$\alpha$ (SurfHoar)', 'n (SurfHoar)'])
    fig.savefig('parametrizations' + '.pdf')


def plot_extr_time(var='lwc', extr='max'):
    profiles = ['MeltCrust_R', 'SurfHoar_R', 'DepthHoar_R', 'LayerCh_R']
    types = ['kWL', 'hWL', 'dWL', 'roh']
    mittel = ['ArithMean', 'GeoMean']
    path = '/Users/fuzi/Documents/Simulations/Snowpack-3.5.0/data/output/'

    start = '2016-01-31T22:00:00'
    end = '2016-02-02T22:00:00'

    for profile in profiles:
        outpath = '/Users/fuzi/Documents/Simulations/Snowpack-3.5.0/data/Plots/extr_times/'
        if not os.path.exists(outpath):
            os.makedirs(outpath)

        y_ticklab = [None] * 4
        for i, typ in zip(np.arange(4), types):
            y_ticklab[i] = profile[:-2] + '\n' + typ
        if profile == 'SchichtW':
            y_ticklab.pop(2)

        fig, axs = plt.subplots()
        fig.set_size_inches(9, 4)

        i = 0
        j = np.array(0)
        for typ in types:
            for mit in mittel:
                netfile = path + mit + '/' + 'NetCDF/' + profile + '_' + typ + '.nc'
                if not os.path.exists(netfile):
                    print(netfile + ' does not exist')
                    continue
                else:
                    ds = st.get_ext_time(xr.open_dataset(netfile), var, extr)

                colors = plt.cm.winter(np.linspace(0, 1, len(ds.run)))
                for r in np.arange(len(ds.run)):
                    df = ds.sel(run=r+1) + i/2  # separate lines for different runs
                    df.plot.line(ax=axs, c=colors[r])
                    i += 1
                i += 1
                j = np.append(j, i/2)
            i += 1
        axs.set_title('')
        axs.set_xlim(pd.Timestamp(start), pd.Timestamp(end))
        axs.xaxis.set_major_formatter(mdates.DateFormatter('%b %d \n%H:%M'))
        #axs.xaxis.set_minor_formatter(mdates.DateFormatter('%H:%M'))
        axs.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
        _ = plt.xticks(rotation=0, horizontalalignment='center')
        axs.set_yticks(j[1::2]+0.5)
        axs.set_ylabel('')
        axs.set_xlabel('Date', fontsize=12)
        axs.set_yticklabels(y_ticklab)
        axs.set_ylim(0.5, i / 2)
        plt.hlines(j[2::2]+0.75, pd.Timestamp(start), pd.Timestamp(end), alpha=0.2, linestyles='--')
        plt.hlines(j[1::2]+0.5, pd.Timestamp(start), pd.Timestamp(end), alpha=0.1, linestyles='dotted')
        #tit = plt.title('Timespan when ' + extr_tit[extr] + ' reaches critical value', fontsize=15)

        for y, mit in zip((j[:-1] + j[1:])/2, mittel*4):
            axs.text(pd.Timestamp('2016-01-31T23:00:00'), y, mit, fontsize=8)

        #fig.savefig(outpath + profile + '_' + var + '_' + extr + '.pdf',  # safe with title
         #           bbox_extra_artists=(tit,), bbox_inches='tight')  # Plot with Title
        fig.savefig(outpath + profile + '_' + var + '_' + extr + '.pdf', bbox_inches='tight')


def main(var, mean_typ, plot_typ='vert'):
    profiles = ['MeltCrust_R', 'SurfHoar_R', 'DepthHoar_R', 'LayerCh_R']
    types = ['kWL', 'hWL', 'dWL', 'roh']
    path = '/Users/fuzi/Documents/Simulations/Snowpack-3.5.0/data/output/' + mean_typ + '/NetCDF/'
    outpath = '/Users/fuzi/Documents/Simulations/Snowpack-3.5.0/data/Plots/' + mean_typ + \
              '/' + var + '_' + plot_typ + '/'
    freq = '30T'

    if not os.path.exists(outpath):
        os.makedirs(outpath)

    for profile in profiles:

        if profile == 'DepthHoar_R':
            start = {'dWL': '2016-02-01T11:00:00', 'hWL': '2016-02-01T10:00:00', 'roh': '2016-02-01T11:00:00',
                     'kWL': '2016-02-01T11:00:00'}
            end = {'dWL': '2016-02-01T18:30:00', 'hWL': '2016-02-01T17:30:00', 'roh': '2016-02-01T18:30:00',
                   'kWL': '2016-02-01T18:30:00'}
        else:
            start = {'dWL': '2016-02-01T09:00:00', 'hWL': '2016-02-01T09:00:00', 'roh': '2016-02-01T09:00:00',
                     'kWL': '2016-02-01T09:00:00'}
            end = {'dWL': '2016-02-01T16:30:00', 'hWL': '2016-02-01T16:30:00', 'roh': '2016-02-01T16:30:00',
                   'kWL': '2016-02-01T16:30:00'}

        for type in types:
            netfile = path + profile + '_' + type + '.nc'
            if not os.path.exists(netfile):
                print(netfile + ' does not exist')
                continue
            else:
                ds = xr.open_dataset(netfile)

                if plot_typ == 'vert':
                    fig, lgd, tit, x_lab, y_lab = vert_profile(ds, var, profile, type,
                                                               time=(start[type], end[type]), freq=freq)
                    # add "tit" to bbox_extra_artists to include title
                    fig.savefig(outpath + profile + '_' + type + '_' + var + '.pdf',
                                bbox_extra_artists=(lgd, x_lab, y_lab,), bbox_inches='tight')

                    plt.close(fig)
                elif plot_typ == 'max':
                    fig, lgd, tit = max_timeline(ds, var, profile, type)
                    # add "tit" to bbox_extra_artists to include title
                    fig.savefig(outpath + profile + '_' + type + '_' + var + '_max.pdf',
                                bbox_extra_artists=(lgd,), bbox_inches='tight')
                    plt.close(fig)
                elif plot_typ == 'min':
                    fig, lgd, tit = max_timeline(ds, var, profile, type, value='min')
                    # add "tit" to bbox_extra_artists to include title
                    fig.savefig(outpath + profile + '_' + type + '_' + var + '_min.pdf',
                                bbox_extra_artists=(lgd,), bbox_inches='tight')
                    plt.close(fig)
                elif plot_typ == 'mean':
                    fig, lgd, tit = max_timeline(ds, var, profile, type, value='mean')
                    # add "tit" to bbox_extra_artists to include title
                    fig.savefig(outpath + profile + '_' + type + '_' + var + '_mean.pdf',
                                bbox_extra_artists=(lgd,), bbox_inches='tight')
                    plt.close(fig)


# Some dictionaries for the plot labels

legend_title = {'dWL': '$h_{WL}$ (cm)', 'hWL': '$h_{slab}$ (cm)', 'roh': '$density$ (%)',
                'kWL': '$d_{grain}$ (mm)'}
type_title = {'dWL': 'height of WL', 'hWL': 'height of Slab',
              'roh': 'density of WL', 'kWL': 'grain diameter'}
var_title = {'Sn_wet': 'Natural Stability $(Sn_{wet})$', 'Sn_wet_manu': '$(Sn_{wet,manu})$',
             'lwc': 'Liquid Water Content $(\\theta_{w,v})$', 'shear_strength': 'Shear Strength $(\sigma_d)$',
             'sigma_wet': 'Shear Strength $(\sigma_{wet})$', 'stress': 'Stress $(\\tau)$'}
max_title = {'Sn_wet': '$Sn_{wet}$', 'Sn_wet_manu': '$Sn_{wet,manu}$', 'lwc': '$\\theta_{w,v}$',
             'shear_strength': '$\sigma_{dry})$', 'sigma_wet': '$\sigma_{wet}$',
             'stress': '$(\\tau)$'}
var_label = {'lwc': '$\\theta_{w,v}$ (%$_{vol}$)', 'Sn_wet': '$Sn_{wet}$', 'Sn_wet_manu': '$Sn_{wet,manu}$',
             'shear_strength': '$\sigma_d$ (kPa)', 'stress': '$\\tau$ (kPa)', 'sigma_wet': '$\sigma_{wet}$ (kPa)'}
max_label = {'lwc': '$\\theta_{max}$ (%$_{vol}$)', 'Sn_wet': '$Sn_{wet, max}$', 'Sn_wet_manu': '$Sn_{wet,manu}$',
             'shear_strength': '$\sigma_{max}$ (kPa)', 'stress': '$\\tau_{max}$ (kPa)'}
extr_tit = {'max': '$MAX_{LWC}$', 'mean': 'LWC$_{index}$', 'min': '$Sn_{wet}$'}
extr_val = {'max': '6%$_{vol}$', 'mean': '3%$_{vol}$', 'min': '1'}

if __name__ == "__main__":
    for mean_typ in ['ArithMean', 'GeoMean', 'hRes/arith', 'hRes/geo']:
        # 'ArithMean', 'GeoMean', 'hRes/arith', 'hRes/geo'
        #main('lwc', mean_typ)
        main('lwc', mean_typ, 'max')
        #main('Sn_wet', mean_typ)
        main('Sn_wet', mean_typ, 'min')
        #main('Sn_wet_manu', mean_typ)
        #main('shear_strength', mean_typ)
        #main('sigma_wet', mean_typ)
        #main('stress', mean_typ)
    #plot_extr_time()
    #plot_extr_time('lwc', 'mean')
    #plot_extr_time('Sn_wet', 'min')
