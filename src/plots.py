from astropy.cosmology import LambdaCDM

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os.path
import pandas as pd

from src.lensed_qso import LensedQSO

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


def plot_lqso_in_speagle(lqso, fig=None, ax=None, label=None):
    # Get age
    #t = np.power(10., lqso.agn_fitter_output()['age'].iloc[2])
    #t_pe = (lqso.agn_fitter_output()['age'].iloc[3] - t) * np.log(10.) * t
    #t_me = (t - lqso.agn_fitter_output()['age'].iloc[1]) * np.log(10.) * t
    #t = t * np.power(10., -9)
    #t_pe = t_pe * np.power(10., -9)
    #t_me = t_me * np.power(10., -9)

    z_min = 1.019
    z_max = 1.589

    # Star formation rate
    sfr_ir, sfr_ir_pe, sfr_ir_me = lqso.get_agnf_output_field('SFR_IR', demag=True)

    sfr_opt, sfr_opt_pe, sfr_opt_me = lqso.get_agnf_output_field('SFR_opt', demag=True)

    sfr_tot = sfr_opt + sfr_ir
    sfr_tot_pe = np.sqrt(sfr_ir_pe ** 2. + sfr_opt_pe ** 2.)
    sfr_tot_me = np.sqrt(sfr_ir_me ** 2. + sfr_opt_me ** 2.)

    # logM_star
    log_m_star, log_m_star_pe,log_m_star_me = lqso.get_agnf_output_field('logMstar', demag=True)

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

    # Format errors properly
    xerr = np.array([log_m_star_me, log_m_star_pe]).reshape((2, 1))

    yerr_tot = np.array([sfr_tot_me, sfr_tot_pe]).reshape((2, 1))
    yerr_tot = yerr_tot / (sfr_tot * np.log(10.))  # Calculate error in total SFR from S

    yerr_ir = np.array([sfr_ir_me, sfr_ir_pe]).reshape((2, 1))
    yerr_ir = yerr_ir / (sfr_ir * np.log(10.))  # Calculate error in SFR_IR from S

    yerr_opt = np.array([sfr_opt_me, sfr_opt_pe]).reshape((2, 1))
    yerr_opt = yerr_opt / (sfr_opt * np.log(10.))  # Calculate error in SFR_opt from S

    # Plot galaxy
    ax.errorbar(log_m_star, np.log10(sfr_tot), xerr=xerr, yerr=yerr_tot, label=f'{lqso.name} SFR_tot, z={lqso.props.z_qso.values[0]:.3f}' if label is None else label, fmt='o', capsize=4)
    # ax.errorbar(log_m_star, np.log10(sfr_ir), xerr=xerr, yerr=yerr_ir, label=f'{lqso.name} SFR_ir, z={lqso.props.z_qso.values[0]:.3f}', fmt='o', capsize=4)
    # ax.errorbar(log_m_star, np.log10(sfr_opt), xerr=xerr, yerr=yerr_opt, label=f'{lqso.name} SFR_opt, z={lqso.props.z_qso.values[0]:.3f}', fmt='o', capsize=4)

    ax.legend(ncol=2)

    fig.savefig(os.path.join('plots', 'speagle.pdf'))
    # fig.savefig(os.path.join('plots', 'speagle.svg'))

    return fig, ax


def plot_agnf_output(lqsos, field_1, field_2, color_scale_field=None, component='_sub', equals_line=False):
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


def residual_plot(lqso):
    realizations = pd.read_csv( os.path.join('data',lqso.name,'agnfitter', 'sed_realizations.csv'))
    data = pd.read_csv( os.path.join('data',lqso.name,'agnfitter', 'sed_data.csv'))
    rest_nu_rea = realizations['nu_rest']
    rest_nu_data = data['nu_rest']
    
    fig, ax = plt.subplots(nrows=1,ncols=1, figsize=(10,8))
    ax.set_xscale ('log')
    ax.set_yscale('log')
    
    terms=['SBnuLnu','BBnuLnu','GAnuLnu', 'TOnuLnu','TOTALnuLnu']
    color= ['forestgreen','darkblue', 'gold', 'purple', 'red']
    
    for i in range(10):
        for term in range(5) :
            #print(realizations[f'{term}{i}'])
            rea = realizations[f'{terms[term]}{i}']
            ax.plot(rest_nu_rea, rea, color=color[term], alpha=0.5)
            
    ax.errorbar(rest_nu_data, data['nuLnu'],yerr=data['nuLnu_err'],color='black',zorder=2, fmt='o')
            
    plt.legend(terms)
    ax.set_ylim(ymin=1e44)
    ax.set_title('plot of the residuals', fontsize=15)
    ax.set_ylabel('$\\nu \\rm{L}(\\nu)[erg \ s^{-1}]$', fontsize=14)
    
    #TODO upperlimits, residuals, x label, pas als allerlaatste lars zn layout doen, opslaan


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
    # ax.errorbar(range(len(x_labels)), agnf_sfr_ir, yerr=agnf_sfr_ir_err, label='This work', color='black', fmt='o', zorder=50)
    ax.errorbar(range(len(x_labels)), np.array(stacey_sfr_ir) - np.array(agnf_sfr_ir), yerr=stacey_sfr_ir_err, label='Stacey+2018 - This work', color='black', fmt='o', zorder=0)
    ax.axhline(y=0, linestyle='dashed', color='grey')

    ax.set_xticks(range(len(x_labels)), x_labels, rotation=90)
    ax.set_ylabel('$\log\mu\mathrm{SFR_{IR}}$')
    ax.legend()

    fig.savefig(os.path.join('plots', 'SFR_IR_stacey_res.pdf'))
