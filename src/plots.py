from astropy.cosmology import LambdaCDM
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from src.lensed_qso import RADIO_CUTOFF, XRAY_CUTOFF, PROPERTIES as LQSO_PROPERTIES, AGNFITTER_FIELDS

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os.path
from src.lensed_qso import XRAY_CUTOFF, RADIO_CUTOFF
import src.ms_data_reader


LCDM = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)  # Cosmological constants as Speagle uses them

log_m_stars = np.linspace(7, 13, num=10000)

# Speagle GMS constants
a = 0.84
sa = 0.02
b = -0.026
sb = 0.003
c = 6.51
sc = 0.24
d = -0.11
sd = 0.03


def speagle_gms(log_m_star, t, log_m_star_err=None, t_err=None):
    if log_m_star_err is None:
        log_m_star_err = 0
    if t_err is None:
        t_err = 0

    log_sfr = (a + b * t) * log_m_star - (c + d * t)
    log_sfr_err = np.sqrt(
        np.power(sa * log_m_star, 2.) +
        np.power(sb * t *log_m_star, 2.) +
        np.power(sc, 2.) +
        np.power(t * sd, 2.) +
        np.power((a + b * t) * log_m_star_err, 2.) +
        np.power((b * log_m_star + d) * t_err, 2.)
    )

    return log_sfr, log_sfr_err


def speagle_gms_residual(df):
    df['univ_age'] = LCDM.age(df['redshift']).value

    lsfr, lsfr_err = speagle_gms(df['logMstar'], df['univ_age'])

    return df['logSFR'] - lsfr, lsfr_err


def plot_speagle_residual(df, fig=None, ax=None, label=None, save_name='speagle_res', x_field='univ_age', x_label='Age of universe [Gyr]', errorbar_kwargs=None):
    if errorbar_kwargs is None:
        errorbar_kwargs = {}

    if fig == ax == None:
        fig, ax = plt.subplots()

        ax.set_xlabel(x_label)
        ax.set_ylabel('Residual with Speagle, logSFR')

    speagle_res = speagle_gms_residual(df)

    ax.errorbar(df[x_field], speagle_res[0], yerr=None if np.sum([df['logSFR_me'],\
                   df[f'logSFR_pe']]) == 0 else np.reshape([df['logSFR_me'], \
                    df[f'logSFR_pe']], (2, len(df[f'logSFR_pe']))), label=label, fmt='o', **errorbar_kwargs)

    lgd = ax.legend(loc='center right', bbox_to_anchor=(1.7, 0.5),
              ncol=1, fancybox=True, shadow=True)

    ax.axhline(0, xmin=0, xmax=1, color='grey', ls='--')
    ax.axhline(0.65, xmin=0, xmax=1, color='grey', ls='--')
    ax.axhline(-0.65, xmin=0, xmax=1, color='grey', ls='--')

    fig.savefig(os.path.join('plots', f'{save_name}.pdf'), bbox_extra_artists=(lgd,), bbox_inches='tight')

    return fig, ax


def plot_lqsos_in_speagle(df, fig=None, ax=None, label=None, save_name='speagle', sfr_type='logSFR', group=True, errorbar_kwargs=None):
    """
    label must be a list when group=False
    """
    if errorbar_kwargs is None:
        errorbar_kwargs = {}
    # Get age
    #t = np.power(10., lqso.agn_fitter_output()['age'].iloc[2])
    #t_pe = (lqso.agn_fitter_output()['age'].iloc[3] - t) * np.log(10.) * t
    #t_me = (t - lqso.agn_fitter_output()['age'].iloc[1]) * np.log(10.) * t
    #t = t * np.power(10., -9)
    #t_pe = t_pe * np.power(10., -9)
    #t_me = t_me * np.power(10., -9)

    z_min = 1.019
    z_max = 1.589

    if ax is None:
        # These things are only done if no ax is given
        fig, ax = plt.subplots(figsize=(10,8))

        sp_ms_max, sp_ms_err_max = speagle_gms(log_m_stars, LCDM.age(z_max).value)

        # ax.fill_between(log_m_stars, sp_ms_max - sp_ms_err_max, sp_ms_max + sp_ms_err_max, alpha=.6, label=f'Speagle+2014, z={z_max}')

        sp_ms, sp_ms_err = speagle_gms(log_m_stars, LCDM.age(z_min).value)

        # ax.fill_between(log_m_stars, sp_ms - sp_ms_err, sp_ms + sp_ms_err, alpha=.6, label=f'Speagle+2014, z={z_min}', color='fuchsia')

        ax.fill_between(log_m_stars, sp_ms - sp_ms_err, sp_ms_max + sp_ms_err_max, alpha=.4, color='grey', label=f'Speagle+2014, z=[{z_min},{z_max}]')

        ax.set_ylabel(sfr_type)
        ax.set_xlabel('log$M_*$')

    # Plot galaxy
    if group:
        ax.errorbar(df['logMstar'], df[sfr_type],
                 xerr=None if np.sum([df['logMstar_me'], df['logMstar_pe']]) == 0 else
                 np.reshape([df['logMstar_me'], df['logMstar_pe']], (2, len(df['logMstar_pe']))),
                 yerr=None if np.sum([df[f'{sfr_type}_me'], df[f'{sfr_type}_pe']]) == 0 else
                 np.reshape([df[f'{sfr_type}_me'], df[f'{sfr_type}_pe']], (2, len(df[f'{sfr_type}_pe']))), label=label, fmt='o',
                 **errorbar_kwargs)
    else:
        for lab, r in zip(label, df.iterrows()):
            _, row = r
            ax.errorbar(row['logMstar'], row[sfr_type], xerr=[[row['logMstar_me']], [row['logMstar_pe']]],
                        yerr=[[row[f'{sfr_type}_me']], [row[f'{sfr_type}_pe']]], label=lab, fmt='o', **errorbar_kwargs)

    ax.set_ylim(ymax=4.1)
    lgd = ax.legend(loc='center right', bbox_to_anchor=(1.65, 0.5),
              ncol=1, fancybox=True, shadow=True)

    fig.savefig(os.path.join('plots', f'{save_name}.pdf'), bbox_extra_artists=(lgd,), bbox_inches='tight')
    # fig.savefig(os.path.join('plots', 'speagle.svg'))

    return fig, ax

def plot_lqsos_in_speagle_z_scaled(df, fig=None, ax=None, label=None, save_name='speagle', group=True, errorbar_kwargs=None):
    """
    Plots all entries in df on Speagle ms, as if they are at redshift 1.5 using their residual to the Speagle main sequence.
    label must be a list when group=False
    """
    if errorbar_kwargs is None:
        errorbar_kwargs = {}

    z = 1.5

    gms, gms_err = speagle_gms(log_m_stars, LCDM.age(z).value)

    df['gms_logsfr'] = np.interp(df['logMstar'], log_m_stars, gms)
    df['gms_res'] = speagle_gms_residual(df)[0]

    if ax is None:
        # These things are only done if no ax is given
        fig, ax = plt.subplots(figsize=(10,8))

        ax.plot(log_m_stars, gms, color='black', label='Speagle MS $\pm 1\sigma$')
        ax.plot(log_m_stars, gms + gms_err, color='black', linestyle='--')
        ax.plot(log_m_stars, gms - gms_err, color='black', linestyle='--')

        ax.set_ylabel('logSFR')
        ax.set_xlabel('log$M_*$')

    # Plot galaxy
    if group:
        ax.errorbar(df['logMstar'], df['gms_logsfr'] + df['gms_res'],
                 xerr=None if np.sum([df['logMstar_me'], df['logMstar_pe']]) == 0 else
                 np.reshape([df['logMstar_me'], df['logMstar_pe']], (2, len(df['logMstar_pe']))),
                 yerr=None if np.sum([df[f'logSFR_me'], df[f'logSFR_pe']]) == 0 else
                 np.reshape([df[f'logSFR_me'], df[f'logSFR_pe']], (2, len(df[f'logSFR_pe']))), label=label, fmt='o',
                 **errorbar_kwargs)
    else:
        for lab, r in zip(label, df.iterrows()):
            _, row = r
            ax.errorbar(row['logMstar'], row['logSFR'], xerr=[[row['logMstar_me']], [row['logMstar_pe']]],
                        yerr=[[row[f'logSFR_me']], [row[f'logSFR_pe']]], label=lab, fmt='o', **errorbar_kwargs)

    ax.set_ylim(ymax=4.1)
    lgd = ax.legend(loc='center right', bbox_to_anchor=(1.65, 0.5),
              ncol=1, fancybox=True, shadow=True)

    fig.savefig(os.path.join('plots', f'{save_name}.pdf'), bbox_extra_artists=(lgd,), bbox_inches='tight')
    # fig.savefig(os.path.join('plots', 'speagle.svg'))

    return fig, ax


def hist_stellarmass(df, fig, ax,label, zorder=1, binwidth=0.20, alpha=0.5, density=True):
    if fig is None:
        fig, ax = plt.subplots()

    ax.hist(df['logMstar'], zorder=zorder, bins=np.arange(int(min(df['logMstar'])), 12.5, binwidth), alpha=alpha,density=density, label=label, edgecolor='black')
    ax.set_xlabel('Log M_star', fontsize=14)
    ax.set_ylabel('normalised number density', fontsize=14)
    lgd = ax.legend(loc='center right', bbox_to_anchor=(1.65, 0.5),
              ncol=1, fancybox=True, shadow=True)

    fig.savefig(os.path.join('plots', 'hist_stellarmass.pdf'), bbox_extra_artists=(lgd,), bbox_inches='tight')
    return fig, ax


def plot_agnf_output(lqsos, field_1, field_2, color_scale_field=None, component='_sub', equals_line=False, logx=False, logy=False, unique_markers=True):
    cm = plt.cm.get_cmap('RdYlBu')
    fig, ax = plt.subplots()
    markers=['o', 'v', '8', 's','h', 'p', 'D', 'X', '>', '<']

    if not unique_markers:
        ax_scatter = ax.scatter(lqsos[field_1], lqsos[field_2], c=lqsos[color_scale_field] if color_scale_field is not None else None, cmap=cm, zorder=100)#, vmin=0, vmax=10000)
        ax.errorbar(lqsos[field_1], lqsos[field_2],
                    xerr=None if np.sum([lqsos[f'{field_1}_me'], lqsos[f'{field_1}_pe']]) == 0 else
                    np.reshape([lqsos[f'{field_1}_me'], lqsos[f'{field_1}_pe']], (2, len(lqsos[f'{field_1}_pe']))),
                    yerr=None if np.sum([lqsos[f'{field_2}_me'], lqsos[f'{field_2}_pe']]) == 0 else
                    np.reshape([lqsos[f'{field_2}_me'], lqsos[f'{field_2}_pe']], (2, len(lqsos[f'{field_2}_pe']))),
                    zorder=0, fmt='o', marker=None)

        if color_scale_field is not None:
            cbar = fig.colorbar(ax_scatter, cmap=cm)
            cbar.set_label(color_scale_field)
    else:
        vmin, vmax = None, None
        if color_scale_field is not None:
            vmin, vmax = np.min(lqsos[color_scale_field]), np.max(lqsos[color_scale_field])

        for i, row in lqsos.iterrows():
            ax_scatter = ax.scatter(row[field_1], row[field_2], c=row[color_scale_field] if color_scale_field is not None else None, cmap=cm, zorder=100, vmin=vmin, vmax=vmax, marker=markers[i] if color_scale_field else None, label=row['name'], edgecolors='black' if color_scale_field else None)
            ax.errorbar(row[field_1], row[field_2],
                        xerr=None if np.sum([row[f'{field_1}_me'], row[f'{field_1}_pe']]) == 0 else
                        np.reshape([row[f'{field_1}_me'], row[f'{field_1}_pe']], (2, 1)),
                        yerr=None if np.sum([row[f'{field_2}_me'], row[f'{field_2}_pe']]) == 0 else
                        np.reshape([row[f'{field_2}_me'], row[f'{field_2}_pe']], (2, 1)),
                        zorder=0, fmt='', marker=None, color=ax_scatter.get_fc())

        if color_scale_field is not None:
            cbar = fig.colorbar(ScalarMappable(Normalize(vmin=vmin, vmax=vmax), cmap=cm))
            cbar.set_label(color_scale_field)

    if unique_markers:
        if color_scale_field:
            bbox_h = 1.65
        else:
            bbox_h = 1.35
        lgd = ax.legend(loc='center right', bbox_to_anchor=(bbox_h, 0.5),
                  ncol=1, fancybox=True, shadow=True)

    if equals_line:
        x = np.linspace(np.min(lqsos[field_1] - lqsos[f'{field_1}_me']), np.max(lqsos[field_1] + lqsos[f'{field_1}_pe']), 10000)
        ax.plot(x, x, linestyle='--', color='black')

    ax.set_xlabel(field_1)
    ax.set_ylabel(field_2)

    if field_1 == 'EBVbbb' and field_2 == 'Nh':
        # Add Type1/2 AGN separation line as found in AGNfitter paper
        ax.vlines(0.2, ymin=21.5, ymax=25, color='black', ls='--')
        ax.hlines(21.5, xmin=0.2, xmax=1, color='black', ls='--')

    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')

    fig.savefig(os.path.join('plots', f'{field_1}_{field_2}.pdf'), bbox_extra_artists=(lgd,), bbox_inches='tight')
    # fig.savefig(os.path.join('plots', f'{field_1}_{field_2}.svg'))

    return fig, ax


def plot_n_runs_pars(lqso, n=10, nrows=4, sub_folder=None):
    pars = AGNFITTER_FIELDS

    ncols = len(pars) // nrows
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(16, 10))
    fig.suptitle(lqso.name)
    for i in range(n):
        lqso.agn_fitter_output(run_time=i, sub_folder=sub_folder)
        for j, p in enumerate(pars):
            r = int(j / ncols)
            c = j % ncols

            val, pe, me = lqso.get_agnf_output_field(p, demag=True)
            axs[r,c].errorbar([val], [1], xerr=[[me], [pe]], fmt='o')
            axs[r,c].set_title(p)
            axs[r,c].tick_params(axis='y', left=False, labelleft=False)

    fig.tight_layout()
    fig.savefig(os.path.join('plots', f'{lqso.name}_pars.pdf'))


def residual_plot(lqso, errors=False):
    #readin in the data
    realizations = lqso.agnf_sed_realizations
    data = lqso.agnf_sed_data
    rest_nu_rea = realizations['nu_rest']
    rest_nu_data = data['nu_rest']

    #setting the general figure
    fig, (ax1,ax2) = plt.subplots(nrows=2,ncols=1, figsize=(10,10), gridspec_kw={'hspace':0, 'height_ratios': [2, 1]}, sharex=True)
    ax1.set_xscale ('log')
    ax1.set_yscale('log')
    ax2.set_xscale ('log')
    # ax2.set_yscale('log')

    # ax1.set_title(f'{lqso.name}', fontsize=15)
    ax1.set_ylabel('$\\nu \\rm{L}(\\nu)[erg \ s^{-1}]$')
    ax2.set_ylabel('$\sigma$')
    # ax2.ticklabel_format(style='scientific', axis='y')
    ax2.set_xlabel('rest frame $\\nu$[Hz]')

    #For getting the right colors
    terms=['SBnuLnu','BBnuLnu','GAnuLnu', 'TOnuLnu','TOTALnuLnu']
    color= ['forestgreen','darkblue', 'gold', 'purple', 'red']

    #first plot
    for i in range(10):
        for term in range(5) :
            #print(realizations[f'{term}{i}'])
            rea = realizations[f'{terms[term]}{i}']
            ax1.plot(rest_nu_rea, rea, color=color[term], alpha=0.5)

    upper = data['upp']==True
    nupper= data['upp']== False
    ax1.errorbar(rest_nu_data[upper], data['nuLnu'][upper],yerr=data['nuLnu_err'][upper],color='black',zorder=2, fmt='v')
    ax1.errorbar(rest_nu_data[nupper], data['nuLnu'][nupper],yerr=data['nuLnu_err'][nupper],color='black',zorder=2, fmt='o')

    ax1.set_ylim(ymin=1e44)
    # ax1.set_xlim(xmin=5e11, xmax=1.2e16)
    ax1.legend(['Starburst','Accretion disk','Stars', 'Torus','Total'])


    #second plot
    for i in range(10):
        rea_interp = np.interp(x = rest_nu_data, xp=rest_nu_rea, fp=realizations[f'TOTALnuLnu{i}'])
        ax2.errorbar(rest_nu_data[nupper], (data['nuLnu'][nupper] - rea_interp[nupper]) / data['nuLnu_err'][nupper], fmt='o' , color='black')
        ax2.errorbar(rest_nu_data[upper], (data['nuLnu'][upper] - rea_interp[upper]) / data['nuLnu_err'][upper], fmt='v' , color='black')
    ax2.set_xlim(xmin=3e18 / XRAY_CUTOFF, xmax=3e18 / RADIO_CUTOFF)
    ax2.axhline(0, xmin=0, xmax=1, color='black')

    fig.tight_layout()
    fig.savefig(os.path.join('plots', f'{lqso.name}_agnf_residuals.pdf'))


def plot_lqsos_vs_stacey(lqsos):
    fig, ax = plt.subplots()
    ax.errorbar(range(len(lqsos['name'])), np.zeros(len(lqsos['name'])),
                yerr=np.reshape([lqsos['logmu_SFR_IR_me'], lqsos['logmu_SFR_IR_pe']], (2, len(lqsos['name']))),
                label='This work', color='blue', fmt='o', zorder=50, alpha=.6)
    ax.errorbar(range(len(lqsos['name'])), lqsos['stacey_sfr'] - lqsos['logmu_SFR_IR'],
                yerr=np.reshape([lqsos['stacey_sfr_me'], lqsos['stacey_sfr_pe']], (2, len(lqsos['name']))),
                label='Stacey+2018', color='green', fmt='o', zorder=0, alpha=.6)

    ax.axhline(y=0, linestyle='dashed', color='grey')

    ax.set_xticks(range(len(lqsos['name'])), lqsos['name'], rotation=90)
    ax.set_ylabel('$\Delta\log\mu\mathrm{SFR_{IR}}$')
    ax.legend()

    fig.tight_layout()
    fig.savefig(os.path.join('plots', 'SFR_IR_stacey_res.pdf'))







def plot_evolution_df(df, fig=None, ax=None, context=True):
    """
    This function plots the stellar mass evolution assuming constant sfr
    under the formula M_star (t) = SFR * t + (M_star(age) - SFR * age)

    This becomes 0 at t = SFR * age - M_star(age) / SFR

    The goal of this is to give a rough indication of whether your sfr's make sense
    """
     #redshift of the galaxy
    z = df['redshift']


    #stellar mass of the galaxy
    logM , pe_logM, me_logM = df['logMstar'],  df['logMstar_pe'], df['logMstar_me']
    M = 10 ** logM


    #total SFR of the galaxy
    sfr_tot = df['SFR']
    sfr_tot_pe = df['SFR_pe']
    sfr_tot_me = df['SFR_me']

    #age of the individual galaxy, plotting this on the plot
    age = LCDM.age(z).value * 1e9 #in yrs

    #setting up the plot
    if ax is None:
        fig,ax= plt.subplots(figsize=(10,8))
    ax.set_xlabel('age of universe [yr]')
    ax.set_ylabel('$M_\star /M_{\odot}$')


    #TODO: add error prop
    M_gas= df['Mgas']
    M_gas_err=df['Mgas_err']

    #make a range of ages in order to make the evolution line
    #lower limit = where no solar mass had been formed
    #upper limit = where all the gas mass has depleted
    if context == True:
        axins = ax.inset_axes([0.9, 0.39, 0.8, 0.6])
    for i in range(len(z)):
        #the constant in the formula
        b = M[i] - (sfr_tot[i] * age[i])
        name = df['name'].values[i]
        ax.scatter(age[i], M[i], label = f'{name}', zorder=100, s=49) #placing the galaxy


        age_range = np.linspace((-(b)/sfr_tot[i]), age[i], 100)
        age_range2 = np.linspace( age[i],(M_gas[i] + M[i] - b)/sfr_tot[i], 100)

        #formula of the M_star assuming constant sfr

        M_range = sfr_tot[i] * age_range + b
        M_range2 = sfr_tot[i] * age_range2 + b

        if i==0:
            ax.plot(age_range, M_range, color='fuchsia', label='time until galaxy formed' )
            ax.plot(age_range2, M_range2, color='blue', label='time until gas depletes' )
        else:
            ax.plot(age_range, M_range, color='fuchsia')
            ax.plot(age_range2, M_range2, color='blue')

        if context == True:
            axins.plot(age_range, M_range, color='fuchsia')
            axins.plot(age_range2, M_range2, color='blue')
            axins.scatter(age[i], M[i], label = f'{name}', zorder=100, s=49) #placing the galaxy

    if context == True:

        ax.set_xlim(8e8,6e9)

        for label, df in src.ms_data_reader.FILES.items():
            if 'Birkin' in label:
                continue
            if 'Sun' in label:
                continue
            ax.scatter(LCDM.age(df[df['redshift'] > 0]['redshift']) * 1e9, np.power(10., df[df['redshift'] > 0]['logMstar']), label=label,
                zorder=50, s=50, alpha=.7)

            axins.scatter(LCDM.age(df[df['redshift'] > 0]['redshift']) * 1e9, np.power(10., df[df['redshift'] > 0]['logMstar']), label=label, zorder=50, s=50, alpha=.7)

        x1, x2, y1, y2 = 3.3e9, 5.3e9, 0, 1.1e11
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        ax.indicate_inset_zoom(axins, edgecolor="black")



    if context==True:
        lgd = ax.legend(loc='center right', bbox_to_anchor=(1.7, 0.19),
                  ncol=2, fancybox=True, shadow=True)
        fig.savefig(os.path.join('plots', 'total_evolution_withdata.pdf'), bbox_extra_artists=(lgd,), bbox_inches='tight')
<<<<<<< HEAD

=======


>>>>>>> ec27e85d231b654e21e880146378f7fbb660d625
    if context==False:
        lgd = ax.legend(loc='center right', bbox_to_anchor=(1.3, 0.5),
                  ncol=1, fancybox=True, shadow=True)
        fig.savefig(os.path.join('plots', 'evolution.pdf'), bbox_extra_artists=(lgd,), bbox_inches='tight')

    return fig, ax


