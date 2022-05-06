from astropy.cosmology import LambdaCDM

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os.path


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


def plot_speagle_residual(df, fig=None, ax=None, label=None, save_name='speagle_res', x_field='univ_age', x_label='Age of universe [Gyr]', errorbar_kwargs=None):
    if errorbar_kwargs is None:
        errorbar_kwargs = {}

    if fig == ax == None:
        fig, ax = plt.subplots()

        ax.set_xlabel(x_label)
        ax.set_ylabel('Residual with Speagle, logSFR')

    df['univ_age'] = LCDM.age(df['redshift']).value

    speagle_log_sfr = speagle_gms(df['logMstar'], df['univ_age'])

    ax.errorbar(df[x_field], df['logSFR'] - speagle_log_sfr[0], yerr=None if np.sum([df['logSFR_me'], df[f'logSFR_pe']]) == 0 else np.reshape([df['logSFR_me'], df[f'logSFR_pe']], (2, len(df[f'logSFR_pe']))), label=label, fmt='o', **errorbar_kwargs)
    ax.legend()

    ax.axhline(0, xmin=0, xmax=1, color='grey', ls='--')

    fig.savefig(os.path.join('plots', f'{save_name}.pdf'))

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

        ax.set_ylabel('log SFR')
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

    ax.legend(ncol=2)

    fig.savefig(os.path.join('plots', f'{save_name}.pdf'))
    # fig.savefig(os.path.join('plots', 'speagle.svg'))

    return fig, ax


def plot_agnf_output(lqsos, field_1, field_2, color_scale_field=None, component='_sub', equals_line=False, logx=False, logy=False, unique_markers=True):
    f1vs = []
    f1es = [[],[]]

    f2vs = []
    f2es = [[],[]]

    fcs = []
    

    for lqso in lqsos:
        f1v, f1pe, f1me = lqso.get_agnf_output_field(field_1, component=component, demag=True)

        f1vs.append(f1v)
        f1es[1].append(f1pe)
        f1es[0].append(f1me)

        f2v, f2pe, f2me = lqso.get_agnf_output_field(field_2, component=component, demag=True)

        f2vs.append(f2v)
        f2es[1].append(f2pe)
        f2es[0].append(f2me)

        if color_scale_field is not None:
            fcs.append(lqso.get_agnf_output_field(color_scale_field, component=component, demag=True)[0])

    cm = plt.cm.get_cmap('RdYlBu')
    fig, ax = plt.subplots()
    markers=['o', 'v', '8', 's','h', 'p', 'D', 'X', '>', '<']
    
    if color_scale_field is not None:
        ax_scatter = ax.scatter(f1vs, f2vs, c=fcs, cmap=cm, zorder=100)#, vmin=0, vmax=10000)
        ax.errorbar(f1vs, f2vs, xerr=f1es, yerr=f2es, zorder=0, fmt='o')#, vmin=0, vmax=10000)
        cbar = fig.colorbar(ax_scatter, cmap=cm)

        cbar.set_label(color_scale_field)
    else:
        ax.errorbar(f1vs, f2vs, xerr=f1es, yerr=f2es, zorder=0, fmt='o')

    if equals_line:
        x = np.linspace(min(np.array(f1vs) - f1es[0]), max(np.array(f1vs) + f1es[1]), 10000)
        ax.plot(x, x, linestyle='--', color='black')

    ax.set_xlabel(field_1)
    ax.set_ylabel(field_2)

    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')

    fig.tight_layout()
    fig.savefig(os.path.join('plots', f'{field_1}_{field_2}.pdf'))
    # fig.savefig(os.path.join('plots', f'{field_1}_{field_2}.svg'))

    return fig, ax


def plot_n_runs_pars(lqso, n=10, nrows=4, sub_folder=None):
    pars = ['tau', 'age', 'Nh', 'irlum', 'SB', 'BB', 'GA', 'TO', 'EBVbbb', 'EBVgal', 'logMstar', 'SFR_opt', 'LIR(8-1000)', 'Lbb(0.1-1)', 'Lbbdered(0.1-1)', 'Lga(01-1)', 'Ltor(1-30)', 'Lsb(1-30)', 'SFR_IR', '-ln_like']

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

    ax1.set_title('plot of the residuals', fontsize=15)
    ax1.set_ylabel('$\\nu \\rm{L}(\\nu)[erg \ s^{-1}]$', fontsize=14)
    ax2.set_ylabel('$\\nu \\rm{L}(\\nu)[erg \ s^{-1}]$', fontsize=14)
    # ax2.ticklabel_format(style='scientific', axis='y')
    ax2.set_xlabel('rest frame $\\nu$[Hz]', fontsize=14)

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
    ax1.errorbar(rest_nu_data[upper], data['nuLnu'][upper],yerr=data['nuLnu_err'][upper],color='black',zorder=2, fmt='o', marker='v')
    ax1.errorbar(rest_nu_data[nupper], data['nuLnu'][nupper],yerr=data['nuLnu_err'][nupper],color='black',zorder=2, fmt='o')

    ax1.set_ylim(ymin=1e44)
    # ax1.set_xlim(xmin=5e11, xmax=1.2e16)
    ax1.legend(['Starburst','Accretion disk','Stars', 'Torus','Total'])


    #second plot
    for i in range(10):
        rea_interp = np.interp(x = rest_nu_data, xp=rest_nu_rea, fp=realizations[f'TOTALnuLnu{i}'])
        if errors:
            ax2.errorbar(rest_nu_data, data['nuLnu'] - rea_interp, yerr=data['nuLnu_err'], fmt='o' , color='black')
        else:
            ax2.scatter(rest_nu_data, data['nuLnu'] - rea_interp, color='black')
    ax2.set_xlim(xmin=3e18 / XRAY_CUTOFF, xmax=3e18 / RADIO_CUTOFF)
    ax2.axhline(0, xmin=0, xmax=1, color='black')

    fig.savefig(os.path.join('plots', f'{lqso.name}_agnf_residuals.pdf'))


def plot_lqsos_vs_stacey(lqsos):
    x_labels = []
    agnf_sfr_ir = []
    agnf_sfr_ir_err = [[], []]
    stacey_sfr_ir = []
    stacey_sfr_ir_err = [[], []]

    for lqso in lqsos:
        if not pd.isnull(lqso.props['stacey_sfr'].values[0]):
            # Append Stacey data
            stacey_sfr_ir.append(lqso.props['stacey_sfr'].values[0])
            stacey_sfr_ir_err[0].append(lqso.props['stacey_sfr_me'].values[0])
            stacey_sfr_ir_err[1].append(lqso.props['stacey_sfr_pe'].values[0])

            # Append our data
            v, pe, me = lqso.get_agnf_output_field('SFR_IR', demag=False)  # Don't demag since Stacey doesn't either

            # Take log
            agnf_sfr_ir.append(np.log10(v))
            agnf_sfr_ir_err[0].append(me / (v * np.log(10)))
            agnf_sfr_ir_err[1].append(pe / (v * np.log(10)))

            x_labels.append(lqso.name)

    fig, ax = plt.subplots()
    ax.errorbar(range(len(x_labels)), np.zeros(len(x_labels)), yerr=agnf_sfr_ir_err, label='This work', color='blue', fmt='o', zorder=50, alpha=.6)
    ax.errorbar(range(len(x_labels)), np.array(stacey_sfr_ir) - np.array(agnf_sfr_ir), yerr=stacey_sfr_ir_err, label='Stacey+2018', color='green', fmt='o', zorder=0, alpha=.6)
    ax.axhline(y=0, linestyle='dashed', color='grey')

    ax.set_xticks(range(len(x_labels)), x_labels, rotation=90)
    ax.set_ylabel('$\Delta\log\mu\mathrm{SFR_{IR}}$')
    ax.legend()

    fig.tight_layout()
    fig.savefig(os.path.join('plots', 'SFR_IR_stacey_res.pdf'))



def plot_evolution(lqso, fig=None, ax=None, single=False):
    """
    This function plots the stellar mass evolution assuming constant sfr
    under the formula M_star (t) = SFR * t + (M_star(age) - SFR * age)

    This becomes 0 at t = SFR * age - M_star(age) / SFR

    The goal of this is to give a rough indication of whether your sfr's make sense
    """
     #redshift of the galaxy
    z = lqso.props.z_qso.values[0]

    #stellar mass of the galaxy
    logM , pe_logM, me_logM = lqso.get_agnf_output_field('logMstar', demag=True)
    M = 10 ** logM

    #TODO: calculate the non-log Mstar errors (if we want them)

    #total SFR of the galaxy
    sfr_ir, sfr_ir_pe, sfr_ir_me = lqso.get_agnf_output_field('SFR_IR', demag=True)
    sfr_opt, sfr_opt_pe, sfr_opt_me = lqso.get_agnf_output_field('SFR_opt', demag=True)

    sfr_tot = sfr_opt + sfr_ir
    sfr_tot_pe = np.sqrt(sfr_ir_pe ** 2. + sfr_opt_pe ** 2.)
    sfr_tot_me = np.sqrt(sfr_ir_me ** 2. + sfr_opt_me ** 2.)

    #age of the individual galaxy, plotting this on the plot
    age = LCDM.age(z).value * 1e9 #in yrs

    #setting up the plot
    if ax is None:
        fig,ax= plt.subplots(figsize=(10,8))
    ax.set_xlabel('age of universe [yr]')
    ax.set_ylabel('Stellar mass [solar mass]')

    ax.scatter(age, M, label = f'{lqso.name}', zorder=100, s=49) #placing the galaxy

    #the constant in the formula
    b = M - (sfr_tot * age)

    #TODO: add real gas masses
    M_gas= 0.4*M #gas mass in solar masses

    #make a range of ages in order to make the evolution line
    #lower limit = where no solar mass had been formed
    #upper limit = where all the gas mass has depleted
    age_range = np.linspace((-(M - sfr_tot * age)/sfr_tot), age, 100)
    age_range2 = np.linspace( age,(M_gas + M - b)/sfr_tot, 100)

    #formula of the M_star assuming constant sfr
    M_range = sfr_tot * age_range + b
    M_range2 = sfr_tot * age_range2 + b

    if single:
        ax.plot(age_range, M_range, color='fuchsia', label='time until galaxy formed' )
        ax.plot(age_range2, M_range2, color='blue', label='time until gas depletes' )
        ax.set_title(f'Stellar mass evolution of {lqso.name}')

        fig.savefig(os.path.join('plots', f'{lqso.name}_evolution.pdf'))


    elif lqso.name == 'J0806+2006':
        ax.plot(age_range, M_range, color='fuchsia', label='time until galaxy formed' )
        ax.plot(age_range2, M_range2, color='blue', label='time until gas depletes' )
        ax.set_title('Stellar mass evolution')
    else:
        ax.plot(age_range, M_range, color='fuchsia')
        ax.plot(age_range2, M_range2, color='blue')
    ax.legend()

    if single == False:
        fig.savefig(os.path.join('plots', 'total_evolution.pdf'))
        
   

    return fig, ax




