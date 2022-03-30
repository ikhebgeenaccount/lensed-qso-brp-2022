#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 11:15:07 2022

@author: abbo
"""
import numpy as np
import matplotlib.pyplot as plt


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


def plot_lqso_in_speagle(lqso, ax=None):
    # Get age
    t = lqso.agn_fitter_output()['age'].iloc[2]
    t_pe = lqso.agn_fitter_output()['age'].iloc[3] - t
    t_me = t - lqso.agn_fitter_output()['age'].iloc[1]
    # Star formation rate
    sfr = lqso.agn_fitter_output()['SFR_IR'].iloc[2]
    sfr_pe = lqso.agn_fitter_output()['SFR_IR'].iloc[3] - sfr
    sfr_me = sfr - lqso.agn_fitter_output()['SFR_IR'].iloc[1]
    # logM_star
    log_m_star = lqso.agn_fitter_output()['logMstar'].iloc[2]
    log_m_star_pe = lqso.agn_fitter_output()['logMstar'].iloc[3] - log_m_star
    log_m_star_me = log_m_star - lqso.agn_fitter_output()['logMstar'].iloc[1]

    if ax is None:
        # These things are only done if no ax is given
        fig, ax = plt.subplots()

        sp_ms, sp_ms_err = speagle_gms(log_m_stars, t, t_err=t_me)

        ax.fill_between(log_m_stars, sp_ms - sp_ms_err, sp_ms + sp_ms_err, alpha=.6, label='Speagle+2014 GMS')

        ax.set_ylabel('log SFR')
        ax.set_xlabel('log$M_*$')

    # Format errors properly
    xerr = np.array([log_m_star_me, log_m_star_pe]).reshape((2, 1))
    yerr = np.array([sfr_me, sfr_pe]).reshape((2, 1))
    yerr = yerr / (sfr * np.log(10.))  # Calculate error in log SFR from SFR

    # Plot galaxy
    ax.errorbar(log_m_star, np.log10(sfr), xerr=xerr, yerr=yerr, label=lqso.name, fmt='o')

    ax.legend()
