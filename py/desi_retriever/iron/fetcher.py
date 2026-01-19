import astropy.table as atpy
import copy
import httpio
import numpy as np
import astropy.io.fits as pyfits
import urllib3
import os
from pylru import lrudecorator, lrucache
import fsspec
import aiohttp

urllib3.disable_warnings()


class si:
    DESI_USER = None
    DESI_PASSWD = None
    cache = lrucache(100)
    gaiaIndex = None


from ..utils import (read_spectra, read_models, get_gaia_index,
                     get_desi_login_password, read_spectra_gz)

urllib3.disable_warnings()


class si:
    DESI_USER = None
    DESI_PASSWD = None
    cache = lrucache(100)
    gaiaIndex = None


class GaiaIndex:

    def __init__(self, gaiaid, targetid, survey, program, hpx, row):
        self.gaiaid = gaiaid
        self.targetid = targetid
        self.survey = survey
        self.program = program
        self.hpx = hpx
        self.row = row

    def search_id(self, gaiaid_new):
        xid = np.searchsorted(self.gaiaid, gaiaid_new)
        if xid >= len(self.gaiaid) or xid < 0:
            return None
        if self.gaiaid[xid] == gaiaid_new:
            return dict(targetid=self.targetid[xid],
                        survey=self.survey[xid],
                        program=self.program[xid],
                        hpx=self.hpx[xid],
                        row=self.row[xid])
        else:
            return None


def fetch_gaia_index():
    if si.gaiaIndex is not None:
        return
    user, pwd = get_desi_login_password()
    D = pyfits.getdata(
        fsspec.open(
            'https://data.desi.lbl.gov/desi/users/koposov/gaiaid_db/indexes/gaia-index-iron-coadd.fits',
            'rb',
            auth=aiohttp.BasicAuth(user, pwd)).open())
    si.gaiaIndex = GaiaIndex(D['EDR3_SOURCE_ID'], D['TARGETID'], D['survey'],
                             D['program'], D['hpx'], D['row'])


@lrudecorator(100)
def get_specs(gaia_edr3_source_id=None,
              tileid=None,
              night=None,
              fiber=None,
              hpx=None,
              targetid=None,
              expid=None,
              group_type='healpix',
              dataset='iron',
              spec_type='coadd',
              survey=None,
              subsurvey=None,
              program=None,
              spectrograph=None,
              mask=False,
              ivar=False,
              fibermap=False,
              nersc=False):
    """
    Get DESI spectra
    
    Parameters
    ----------
    tileid: int
    night: int or string
         The night identifier (i.e. 20200220 or 'all' or 'deep' for coadds)
    fiber: int
    targetid: int (optional)
    expid: int (optional)
    coadd: bool
         If true read coadded spectra
    mask: bool
         If true return the masks as well
    ivar: bool
         If true return the inverse variances

    Returns
    -------
    ret: list(dict)
        The list of dictionaries for each observation
        where each dictionary
        has keywords b_wavelength, r_wavelength, z_wavelength
        b_flux, b_mask, b_ivar

    """
    if nersc:
        user, pwd = None, None
    else:
        user, pwd = get_desi_login_password()

    if spec_type not in ['coadd', 'cframe', 'spectra']:
        raise Exception('unknown')
    if group_type not in ['exposure', 'healpix', 'tiles/cumulative', 'tiles']:
        raise Exception('unknown')
    if subsurvey is not None:
        print('Warning subsurvey is deprecated, use program keyword')
        program = subsurvey
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
        survey = res['survey']
        program = res['program']
        hpx = res['hpx']
        targetid = res['targetid']

    if group_type == 'tiles/cumulative':
        night1 = f'thru{night}'
    else:
        night1 = night

    if nersc:
        data_desi = '/global/cfs/cdirs/desi/'
    else:
        data_desi = 'https://data.desi.lbl.gov/desi/'
    if spec_type == 'spectra':
        fname_end = '.fits.gz'
    else:
        fname_end = '.fits'

    if group_type == 'tiles/cumulative':
        fname = f'{spec_type}-{spectrograph}-{tileid}-{night1}{fname_end}'
        url = (f'{data_desi}/spectro/redux/{dataset}/tiles/cumulative/' +
               f'{tileid}/{night}/{fname}')
    elif group_type == 'tiles':
        fname = f'{spec_type}-{spectrograph}-{tileid}-{night1}{fname_end}'
        url = (f'{data_desi}/spectro/redux/{dataset}/tiles/{tileid}/'
               f'{night}/{fname}')
    elif group_type == 'healpix':
        fname = f'{spec_type}-{survey}-{program}-{hpx}{fname_end}'
        url = (f'{data_desi}/spectro/redux/{dataset}/healpix/{survey}/'
               f'{program}/{hpx//100}/{hpx}/{fname}')
        print(url)
    else:
        raise Exception('oops')

    if spec_type == 'spectra' and not nersc:
        reader = read_spectra_gz
    else:
        reader = read_spectra

    return reader(url,
                  targetid,
                  fiber=fiber,
                  expid=expid,
                  mask=mask,
                  ivar=ivar,
                  fibermap=fibermap,
                  user=user,
                  pwd=pwd,
                  dataset=dataset)


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
                      subsurvey=None,
                      program=None,
                      spec_type='coadd',
                      group_type='healpix',
                      run='230211',
                      dataset='iron',
                      nersc=False):
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
    if nersc:
        user, pwd = None, None
    else:
        user, pwd = get_desi_login_password()
    if subsurvey is not None:
        print('Warning subsurvey is deprecated, use program keyword')
        program = subsurvey
    if gaia_edr3_source_id is not None:
        fetch_gaia_index()
        res = si.gaiaIndex.search_id(gaia_edr3_source_id)
        if res is None:
            raise ValueError('object not found')
        survey = res['survey']
        program = res['program']
        hpx = res['hpx']
        targetid = res['targetid']
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

    if nersc:
        data_desi = '/global/cfs/cdirs/desi/'
    else:
        data_desi = 'https://data.desi.lbl.gov/desi/'

    data_desi = (f'{data_desi}/science/mws/redux/'
                 f'{dataset}/rv_output/{run}/')
    if group_type == 'tiles/cumulative':
        fname = f'rvmod_{spec_type}-{spectrograph}-{tileid}-{night1}.fits'
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

    return read_models(url, targetid, fiber=fiber, user=user, pwd=pwd)
