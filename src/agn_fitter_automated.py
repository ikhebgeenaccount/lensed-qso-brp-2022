from src.lensed_qso import LensedQSO
import os.path
import os

RUN_TEN_MODE = 'python ../AGNfitter/RUN_AGNfitter_multi.py ../AGNfitter/example/SETTINGS_AGNfitter_{name}.py --overwrite'
RUN_SINGLE_MODE = RUN_TEN_MODE + ' -n 0'

AGN_FITTER_PATH = os.path.join(os.pardir, 'AGNfitter')
AGN_FITTER_RX_PATH = os.path.join(os.pardir, 'AGNfitter-rX_v0.1', 'AGNfitter')


def run_agn_fitter(galaxies, rX=False, run_ten=False, settings=None, demag=False):
    if rX:
        path = AGN_FITTER_RX_PATH
    else:
        path = AGN_FITTER_PATH

    for g in galaxies:
        lqso = LensedQSO(g)

        # Update catalog
        print(f'Updating catalog for {lqso.name}')
        cat, l = lqso.sed_to_agn_fitter(rX=rX, component='_sub_demag' if demag else '_sub')

        with open(os.path.join(path, 'data', f'{lqso.name}.txt'), 'w') as cf:
            cf.write(cat)

        # Update settings
        print(f'Updating settings for {lqso.name}')
        settings = lqso.agn_settings(rX=rX, settings=settings)
        with open(os.path.join(path, 'example', f'SETTINGS_AGNfitter_{lqso.name}.py'), 'w') as sf:
            sf.write(settings)

        command = RUN_SINGLE_MODE.format(**{'name': lqso.name})
        if run_ten:
            print('Selected run_ten mode')
            command = RUN_TEN_MODE.format(**{'name': lqso.name})

        print(f'Running AGNfitter for {lqso.name}')
        os.system(command)
