from src.lensed_qso import LensedQSO
from src.lens_subtraction.model_subtraction import model_subtraction
from src.agnfitter.agn_fitter_automated import run_agn_fitter

import argparse
import os
import json


GALAXIES = ['J0806+2006', 'J0924+0219', 'B1152+200', 'J1330+1810', 'J1455+1447', 'J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']


def agnfitter(galaxies, run_times=1, git_push=False, rX=False, copy=False, model_sub=False, settings=None, component='_sub_demag'):

    for g in galaxies:
        lqso = LensedQSO(g)

        if model_sub:
            model_subtraction(lqso)
        run_agn_fitter([g], run_times=run_times, rX=rX, component=component)

        if copy:
            lqso.agn_fitter_output(copy=copy)

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
    parser.add_argument('--times', help='Run all galaxies n times', type=int, default=1)
    parser.add_argument('--copy', help='Copy files from AGNfitter OUTPUT to git repo', action='store_true')
    parser.add_argument('--modelsub', help='Run model subtraction for galaxies', action='store_true')
    parser.add_argument('--component', help='Flux component to use, default _sub', type=str, default='_sub')
    #   --single
    parser.add_argument('--selection', type=str, help='Run a selection of galaxies, separated by commas')
    # parser.add_argument('--settings', type=str, help='Settings file to use, if not given uses default settings')

    args = parser.parse_args()

    gals = GALAXIES
    if args.selection:
        gals = args.selection.split(',')

    # settings = None
    # if args.settings:
    #     with open(args.settings, 'r') as f:
    #         settings = json.load(f)

    agnfitter(gals, run_times=args.times, git_push=args.push,
              rX=args.rX, copy=args.copy, model_sub=args.modelsub,
              # settings=settings,
              component=args.component
    )
