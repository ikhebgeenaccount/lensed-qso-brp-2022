import astropy.io.fits

import numpy as np
import pandas as pd


def sdss_fits_to_csv(galaxy, file, url='', units='Ang,e-17 erg/s/cm2/Ang,e-17 erg/s/cm2/Ang', hdu_index=1):
    """
    Extracts the wavelength, flux and flux error from a SDSS fits file obtained with the BOSS spectograph.
    Fits file should be stored in data/[galaxy]/[file], where file is the fits file.
    Csv file is saved to data/[galaxy]/SDSS_[file].csv where file is the fits file without the .fits extension.
    :param galaxy: Name of the galaxy
    :param file: File name of fits file
    :param url: Optional, url from where the fits file was obtained
    :param units: Units of the fields, default: 'Ang,e-17 erg/s/cm2/Ang,e-17 erg/s/cm2/Ang'
    :param hdu_index: the index of the hdu containing the wavelength and flux data in the fits file, default: 1
    :return: A Pandas DataFrame with the fields wavelength, flux and fluxerr
    """
    f = astropy.io.fits.open(f'data/{galaxy}/{file}')

    df = pd.DataFrame()

    # Source for SDSS fields: http://www.sdss3.org/dr9/spectro/spectro_basics.php

    # ivar is inverse variance, defined as 1/sigma^2
    # ivar = 0 indicated bad data in that row
    # Select only that data with ivar > 0
    selection = f[hdu_index].data['ivar'] > 0

    # loglam is 1og_10(lambda), convert
    df['wavelength'] = np.power(10., f[hdu_index].data['loglam'])[selection]
    df['flux'] = f[hdu_index].data['flux'][selection]
    # Convert ivar to flux error
    df['fluxerr'] = 1. / np.sqrt(f[hdu_index].data['ivar'][selection])

    info_string = f"# Obtained from SDSS, fits file {file}, url={url}\n#{units}\n"

    csv_file = f'data/{galaxy}/SDSS_{file[0:-5]}.csv'
    with open(csv_file, 'w') as file:
        file.write(info_string)

    df.to_csv(csv_file, mode='a', index=False)

    return df


