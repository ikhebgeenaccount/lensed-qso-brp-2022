import time
import warnings

import pandas as pd
import requests

import astropy.units as u
from bs4 import BeautifulSoup


BASE_URL = 'http://vo.ned.ipac.caltech.edu/services/{service}'

last_request_time = 0


def _generate_url(service):
    return BASE_URL.format(service=service)


def _get_page(url, params):
    # Make sure we do not request data more than once per second
    global last_request_time
    time_since_last_request = time.time() - last_request_time
    if time_since_last_request < 1:
        time.sleep((1 - time_since_last_request))
    last_request_time = time.time()

    return requests.get(url, params=params)


def _response_to_dataframe(res):
    xml_str = res.text
    soup = BeautifulSoup(xml_str, features='xml')

    # Check if object found
    msg = soup.find('PARAM', {'name': 'Message'})  # PARAM message is only there when no object found
    if msg is not None:
        warnings.warn(msg['value'] + ', ' + res.url)

    data = {}

    # Find column titles
    cols = []
    for f in soup.find_all('FIELD'):
        data[f['ID']] = []
        cols.append(f['ID'])

    # Loop through all table rows and add them to lists
    for tr in soup.find_all('TR'):
        for i, td in enumerate(tr.find_all('TD')):
            data[cols[i]].append(td.text.strip())

    out_df = pd.DataFrame(data)
    out_df.rename(columns={
        'DataSpectralPassBand': 'filter',
        'DataSpectralValue': 'frequency',
        'DataFrequencyMode': 'fremode',
        'DataFluxValue': 'flux_total',
        'DataFluxStatErr': 'flux_err',
        'DataFluxUnit': 'unit',
        'DataRefcode': 'bibcode',
        'DataQualifiers': 'qualifiers',
    }, inplace=True)

    out_df['wavelength'] = 3e18 / out_df['frequency'].astype(float)
    out_df['flux_total'] = _convert_array_units(out_df['flux_total'], out_df['unit'], 'mJy')
    out_df['flux_err'] = _convert_array_units(out_df['flux_err'], out_df['unit'], 'mJy')

    return out_df[['filter', 'wavelength', 'flux_total', 'flux_err', 'bibcode', 'qualifiers']]


def _convert_array_units(arr, units, to_unit):
    """
    Converts all values in array arr from units in unit to units in to_unit
    :param arr: Array of values
    :param units: Array of strings, those being units of values
    :param to_unit: Single string of unit
    :return: Array containing unitless values converted to to_unit
    """
    return [(v * getattr(u, unit)).to(getattr(u, to_unit)) for v, unit in zip(arr, units)]


def query_sed(target_name):
    return _get_page(_generate_url('querySED'), params={'TARGETNAME': target_name, 'REQUEST': 'queryData'})


def access_sed(target_name):
    """
    Retrieves all SED data from NED and returns them as a pandas DataFrame.
    :param target_name:
    :return:
    """
    res = _get_page(_generate_url('accessSED'), params={'TARGETNAME': target_name, 'REQUEST': 'getData'})

    return _response_to_dataframe(res)
