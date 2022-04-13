from src.lensed_qso import LensedQSO
from src.model_subtraction import model_subtraction
from src.agn_fitter_automated import run_agn_fitter
from src.plots import plot_lqso_in_speagle

import argparse
import os
import json


GALAXIES = ['J0806+2006', 'J0924+0219', 'B1152+200', 'J1330+1810', 'J1455+1447', 'J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']


def agnfitter(galaxies, run_ten=False, git_push=False, rX=False, copy=False, model_sub=False, settings=None, component='_sub_demag', speagle=False):

    ax = None

    for g in galaxies:
        lqso = LensedQSO(g)

        if model_sub:
            model_subtraction(lqso)
        run_agn_fitter([g], run_ten=run_ten, rX=rX, settings=settings, component=component)
        lqso.agn_fitter_output(copy=copy)

        if speagle:
            _, ax = plot_lqso_in_speagle(lqso, ax=ax)

        if git_push:
            os.system('git add data/*')
            os.system(f'git commit -m "AGNfitter output automated {lqso.name}"')
            os.system('git push')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # args:
    #   --push: flag whether to git push or not
    parser.add_argument('--push', help='Push to GitHub repo', action='store_true')
    parser.add_argument('--rX', help='Run rX version of AGNfitter', action='store_true')
    parser.add_argument('--tenmode', help='Run all galaxies 10 times', action='store_true')
    parser.add_argument('--copy', help='Copy files from AGNfitter OUTPUT to git repo', action='store_true')
    parser.add_argument('--modelsub', help='Run model subtraction for galaxies', action='store_true')
    parser.add_argument('--component', help='Flux component to use, default "_sub_demag"', type=str, default='_sub_demag')
    parser.add_argument('--speagle', help='Plot run galaxies in a Speagle MS', action='store_true')
    #   --single
    parser.add_argument('--single', type=str, help='Run a single galaxy, give name, if not given, runs all galaxies')
    parser.add_argument('--settings', type=str, help='Settings file to use, if not given uses default settings')

    args = parser.parse_args()

    gals = GALAXIES
    if args.single:
        gals = [args.single]

    settings = None
    if args.settings:
        with open(args.settings, 'r') as f:
            settings = json.load(f)

    agnfitter(gals, run_ten=args.tenmode, git_push=args.push,
              rX=args.rX, copy=args.copy, model_sub=args.modelsub,
              settings=settings, component=args.component, speagle=args.speagle
    )
