import time
import warnings

import pandas as pd
import requests

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

    return pd.DataFrame(data)


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


# pd.set_option('max_columns', None)
# ned_sed = access_sed('J0806+2006')
# print(ned_sed.columns)
# print(ned_sed)
