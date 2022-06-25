import time

import pandas as pd
import requests

from bs4 import BeautifulSoup


BASE_URL = 'http://vo.ned.ipac.caltech.edu/services/{service}'

last_request_time = 0


def _generate_url(service):
    return BASE_URL.format(service=service)


def _get_page(url, params):
    # Make sure we do not request data more than once per second
    time_since_last_request = time.time() * 1000 - last_request_time
    if time_since_last_request < 1000:
        time.sleep((1000 - time_since_last_request) / 1000)

    return requests.get(url, params=params)


def _xml_to_dataframe(xml_str):
    soup = BeautifulSoup(xml_str, features='xml')

    data = {}

    # Find column titles
    cols = []
    for f in soup.find_all('FIELD'):
        data[f['ID']] = []
        cols.append(f['ID'])

    for tr in soup.find_all('TR'):
        for i, td in enumerate(tr.find_all('TD')):
            data[cols[i]].append(td.text.strip())

    return pd.DataFrame(data)


def query_sed(target_name):
    return _get_page(_generate_url('querySED'), params={'TARGETNAME': target_name, 'REQUEST': 'queryData'})


def access_sed(target_name):
    req = _get_page(_generate_url('accessSED'), params={'TARGETNAME': target_name, 'REQUEST': 'getData'})
    return _xml_to_dataframe(req.text)


pd.set_option('max_columns', None)
print(access_sed('SDSS0806+2006')[['DataSpectralPassBand', 'DataFluxValue', 'DataFluxStatErr', 'DataFluxUnit', 'DataSpectralPublishedValue']])