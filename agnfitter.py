from src.lensed_qso import LensedQSO
from src.model_subtraction import model_subtraction
from src.agn_fitter_automated import run_agn_fitter

import argparse
import os


GALAXIES = ['J0806+2006', 'J0924+0219', 'B1152+200', 'J1330+1810', 'J1455+1447', 'J1524+4409', 'B1600+434', 'B1608+656', 'J1633+3134', 'J1650+4251']


def agnfitter(galaxies, run_ten=False, git_push=False, rX=False, copy=False):
    for g in galaxies:
        lqso = LensedQSO(g)
        model_subtraction(lqso)
        run_agn_fitter([g], run_ten=run_ten, rX=rX)
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
    parser.add_argument('--tenmode', help='Run all galaxies 10 times', action='store_true')
    parser.add_argument('--copy', help='Copy files from AGNfitter OUTPUT to git repo', action='store_true')
    #   --single
    parser.add_argument('--single', type=str, help='Run a single galaxy, give name, if not given, runs all galaxies')

    args = parser.parse_args()

    gals = GALAXIES
    if args.single:
        gals = [args.single]

    agnfitter(gals, run_ten=args.tenmode, git_push=args.push, rX=args.rX, copy=args.copy)
