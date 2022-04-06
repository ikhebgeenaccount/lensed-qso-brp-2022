from astropy.cosmology import LambdaCDM

import numpy as np
import matplotlib.pyplot as plt

import os.path

LCDM = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)  # Cosmological constants as Speagle uses them

log_m_stars = np.linspace(8, 14, num=10000)

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


def plot_lqso_in_speagle(lqso, fig=None, ax=None):
    if lqso.agn_fitter_output() is None:
        print('No AGNfitter output found for', lqso.name)
        return
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
    sfr_ir = lqso.agn_fitter_output()['SFR_IR'].iloc[2]
    sfr_ir_pe = lqso.agn_fitter_output()['SFR_IR'].iloc[3] - sfr_ir
    sfr_ir_me = sfr_ir - lqso.agn_fitter_output()['SFR_IR'].iloc[1]

    sfr_opt = lqso.agn_fitter_output()['SFR_opt'].iloc[2]
    sfr_opt_pe = lqso.agn_fitter_output()['SFR_opt'].iloc[3] - sfr_opt
    sfr_opt_me = sfr_opt - lqso.agn_fitter_output()['SFR_opt'].iloc[1]
    
    sfr_tot = sfr_opt + sfr_ir 
    sfr_tot_pe = np.sqrt(sfr_ir_pe **2 + sfr_opt_pe**2)
    sfr_tot_me = np.sqrt(sfr_ir_me **2 + sfr_opt_me**2)

    # logM_star
    log_m_star = lqso.agn_fitter_output()['logMstar'].iloc[2]
    log_m_star_pe = lqso.agn_fitter_output()['logMstar'].iloc[3] - log_m_star
    log_m_star_me = log_m_star - lqso.agn_fitter_output()['logMstar'].iloc[1]

    # Magnification
    mu = lqso.props.magnification.values[0]

    if ax is None:
        # These things are only done if no ax is given
        fig, ax = plt.subplots(figsize=(10,8))

        sp_ms, sp_ms_err = speagle_gms(log_m_stars, LCDM.age(z_max).value)

        ax.fill_between(log_m_stars, sp_ms - sp_ms_err, sp_ms + sp_ms_err, alpha=.6, label=f'Speagle+2014, z={z_max}')

        sp_ms, sp_ms_err = speagle_gms(log_m_stars, LCDM.age(z_min).value)

        ax.fill_between(log_m_stars, sp_ms - sp_ms_err, sp_ms + sp_ms_err, alpha=.6, label=f'Speagle+2014, z={z_min}', color='fuchsia')

        ax.set_ylabel('log SFR')
        ax.set_xlabel('log$M_*$')


    # Format errors properly
    xerr = np.array([log_m_star_me, log_m_star_pe]).reshape((2, 1))
    
    yerr_ir = np.array([sfr_ir_me, sfr_ir_pe]).reshape((2, 1))
    yerr_ir = yerr_ir / (sfr_ir * np.log(10.))  # Calculate error in log SFR from SFR

    yerr_opt = np.array([sfr_opt_me, sfr_opt_pe]).reshape((2, 1))
    yerr_opt = yerr_opt / (sfr_opt * np.log(10.))  # Calculate error in log SFR from SFR
    
    yerr_tot = np.array([sfr_tot_me, sfr_tot_pe]).reshape((2, 1))
    yerr_tot = yerr_tot / (sfr_tot * np.log(10.))  # Calculate error in total SFR from S    

    # Plot galaxy
    #ax.errorbar(log_m_star - np.log10(mu), np.log10(sfr_ir/mu), xerr=xerr, yerr=yerr_ir, label=f'{lqso.name} SFR_IR, z={lqso.props.z_qso.values[0]:.3f}', fmt='o', capsize=4)
    #ax.errorbar(log_m_star - np.log10(mu), np.log10(sfr_opt/mu), xerr=xerr, yerr=yerr_opt, label=f'{lqso.name} SFR_opt, z={lqso.props.z_qso.values[0]:.3f}', fmt='o', capsize=4)
    ax.errorbar(log_m_star - np.log10(mu), np.log10(sfr_tot/mu), xerr=xerr, yerr=yerr_tot, label=f'{lqso.name} SFR_tot, z={lqso.props.z_qso.values[0]:.3f}', fmt='o', capsize=4)

    ax.legend()

    fig.savefig(os.path.join('plots', f'speagle.pdf'))
    fig.savefig(os.path.join('plots', f'speagle.svg'))

    return fig, ax
