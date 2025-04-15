import astropy.table as atpy
import copy
import httpio
import numpy as np
import astropy.io.fits as pyfits
import urllib3
import os
from pylru import lrudecorator, lrucache
import pyarrow.parquet as pq
import fsspec
import pickle
import requests
import traceback
from ..utils import read_spectra, read_models

urllib3.disable_warnings()


class si:
    DESI_USER = None
    DESI_PASSWD = None
    cache = lrucache(100)
    gaiaIndex = None


class GaiaIndex:

    def __init__(self, key_arr, pos1, pos2, bin_url=None, columns=None):
        self.key_arr = key_arr
        self.pos1 = pos1
        self.pos2 = pos2
        self.bin_url = bin_url
        self.columns = columns

    def search_id(self, key_val):
        xid = np.searchsorted(self.key_arr, key_val)
        if xid >= len(self.key_arr) or xid < 0:
            return None
        if self.key_arr[xid] == key_val:
            with fsspec.open(self.bin_url, 'rb').open() as fp:
                fp.seek(self.pos1[xid])
                D = pickle.loads(fp.read(self.pos2[xid] - self.pos1[xid]))
                return dict(zip(self.columns, D))
        else:
            return None


def fetch_gaia_index():
    if si.gaiaIndex is not None:
        return
    base_url = ('https://portal.nersc.gov/project/desi/users/koposov/'
                'dr1/gaia_id_db/')
    path_base_path = os.path.dirname(os.path.abspath(__file__))
    parquet_fname = 'gaia-index-dr1-coadd_250319.parquet'
    parquet_url = f'{base_url}/' + parquet_fname
    bin_url = f'{base_url}/gaia-index-dr1-coadd_250319.bin'
    for i in range(2):
        try:
            with open(path_base_path + '/' + parquet_fname, 'rb') as fp:
                pqf = pq.ParquetFile(fp).read()
        except:  # noqa
            print('Downloading remote parquet file')
            try:
                with requests.get(parquet_url) as rp:
                    with open(path_base_path + '/' + parquet_fname,
                              'wb') as fp_out:
                        fp_out.write(rp.content)
                print('Successfully downloaded')
            except:  # noqa
                print('''Failed to download the gaia index file
You may want to update desi_retriever''')
                traceback.print_exc()
            continue

    with fsspec.open(bin_url, 'rb').open() as fp:
        header = 1000
        keys = pickle.loads(fp.read(header))

    D = {}
    for k in ['EDR3_SOURCE_ID', 'pos1', 'pos2']:
        D[k] = np.array(pqf[k])
    si.gaiaIndex = GaiaIndex(D['EDR3_SOURCE_ID'],
                             D['pos1'],
                             D['pos2'],
                             bin_url=bin_url,
                             columns=keys)


@lrudecorator(100)
def get_specs(gaia_edr3_source_id=None,
              tileid=None,
              night=None,
              fiber=None,
              hpx=None,
              targetid=None,
              expid=None,
              group_type='healpix',
              spec_type='coadd',
              dataset='iron',
              survey=None,
              program=None,
              spectrograph=None,
              mask=False,
              ivar=False,
              fibermap=False):
    """
    Get DESI spectra for a single object.
    Typically if you are getting a coadded object, you
    need to provide targetid, healpix, survey, program.

    Parameters
    ----------
    healpix: int
        The HEALPIX with the object ( Nside=64 nested)
    survey: string
        The name of the survey (i.e. main/sv1 etc)
    program: string
        The name of the program (i.e. bright/dark)
    targetid: int
        DESI target identifier (TARGETID)
    coadd: bool
         If true read coadded spectra
    fiber: int
    tileid: int
    night: int or string
         The night identifier (i.e. 20200220 or 'all' or 'deep' for coadds)
    expid: int
    mask: bool
         If true return the masks as well
    ivar: bool
         If true return the inverse variances
    fibermap: bool
         If true return the inverse fibermap

    Returns
    -------
    ret: list(dict)
        The list of dictionaries for each observation
        where each dictionary
        has keywords b_wavelength, r_wavelength, z_wavelength
        b_flux, b_mask, b_ivar

    """

    if spec_type not in ['coadd', 'cframe', 'spectra']:
        raise Exception('unknown')
    if group_type not in ['exposure', 'healpix', 'tiles/cumulative', 'tiles']:
        raise Exception('unknown')
    if group_type != 'healpix':
        if fiber is None:
            raise Exception(
                'Fiber must be specified as it is needed to identify the ' +
                'spectrograph')
        if spectrograph is None:
            spectrograph = fiber // 500

    if gaia_edr3_source_id is not None:
        fetch_gaia_index()
        res = si.gaiaIndex.search_id(gaia_edr3_source_id)
        if res is None:
            raise ValueError('object not found')
        survey = res['SURVEY']
        program = res['PROGRAM']
        hpx = res['hpx']
        targetid = res['TARGETID']

    if group_type == 'tiles/cumulative':
        night1 = f'thru{night}'
    else:
        night1 = night

    data_desi = 'https://data.desi.lbl.gov/public/dr1/'
    if group_type == 'tiles/cumulative':
        fname = f'{spec_type}-{spectrograph}-{tileid}-{night1}.fits'
        url = (f'{data_desi}/spectro/redux/{dataset}/tiles/cumulative/' +
               f'{tileid}/{night}/{fname}')
    elif group_type == 'tiles':
        fname = f'{spec_type}-{spectrograph}-{tileid}-{night1}.fits'
        url = (f'{data_desi}/spectro/redux/{dataset}/tiles/{tileid}/'
               f'{night}/{fname}')
    elif group_type == 'healpix':
        fname = f'{spec_type}-{survey}-{program}-{hpx}.fits'
        url = (f'{data_desi}/spectro/redux/{dataset}/healpix/{survey}/'
               f'{program}/{hpx//100}/{hpx}/{fname}')
        print(url)

    else:
        raise Exception('oops')
    return read_spectra(url,
                        targetid,
                        fiber=fiber,
                        expid=expid,
                        mask=mask,
                        ivar=ivar,
                        fibermap=fibermap)


@lrudecorator(100)
def get_rvspec_models(gaia_edr3_source_id=None,
                      tileid=None,
                      night=None,
                      fiber=None,
                      targetid=None,
                      expid=None,
                      hpx=None,
                      coadd=True,
                      survey=None,
                      program=None,
                      dataset='iron',
                      spec_type='coadd',
                      group_type='healpix',
                      run='240520'):
    """
    Get RVSpecfit models

    Parameters
    ----------
    tileid: int
    night: int
    fiber: int
    targetid: int
    expid: int (optional)
    coadd: bool
         If true read coadded spectra
    run: string
         The string identifying a software run
    dataset: the dataset fitted (i.e. andes/sv_daily)

    Returns
    -------
    ret: list(dict)
        The list of dictionaries where each dictionary
        has keywords b_wavelength, r_wavelength, z_wavelength
        b_model etc
    """
    if gaia_edr3_source_id is not None:
        fetch_gaia_index()
        res = si.gaiaIndex.search_id(gaia_edr3_source_id)
        if res is None:
            raise ValueError('object not found')
        survey = res['SURVEY']
        program = res['PROGRAM']
        hpx = res['hpx']
        targetid = res['TARGETID']
    if spec_type not in ['coadd', 'cframe', 'spectra']:
        raise Exception('unknown')
    if group_type != 'healpix':
        if fiber is None:
            raise ValueError(
                'Fiber must be specified as it is needed to identify the ' +
                'spectrograph')
        spectrograph = fiber // 500
    if group_type == 'cumulative':
        night1 = f'thru{night}'
    else:
        night1 = night

    data_desi = (f'https://data.desi.lbl.gov/public/dr1/vac/dr1/mws/'
                 f'{dataset}/v1.0/rv_output/{run}/')
    if group_type == 'tiles/cumulative':
        fname = f'rvmod_`{spec_type}-{spectrograph}-{tileid}-{night1}.fits'
        url = (f'{data_desi}/tiles/cumulative/' + f'{tileid}/{night}/{fname}')
    elif group_type == 'tiles':
        fname = f'rvmod_{spec_type}-{spectrograph}-{tileid}-{night1}.fits'
        url = (f'{data_desi}/tiles/{tileid}/'
               f'{night}/{fname}')
    elif group_type == 'healpix':
        fname = f'rvmod_{spec_type}-{survey}-{program}-{hpx}.fits'
        url = (f'{data_desi}/healpix/{survey}/'
               f'{program}/{hpx//100}/{hpx}/{fname}')
        print(url)

    return read_models(url, targetid, fiber=fiber)
