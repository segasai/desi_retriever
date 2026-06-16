import os

import astropy.table as atpy
import urllib3
from pylru import lrudecorator

from ..utils import (get_desi_login_password, get_gaia_index, read_models,
                     read_spectra, read_spectra_gz)

urllib3.disable_warnings()

DEFAULT_NSIDE = 64
GAIA_PARQUET_FNAME = None
GAIA_BIN_FNAME = None


def _resolve_uniqpix(uniqpix=None, hpx=None, nside=DEFAULT_NSIDE):
    if uniqpix is not None:
        return int(uniqpix)
    if hpx is None:
        return None
    return int(4 * int(nside)**2 + int(hpx))


def _first_table_value(table, *names):
    for name in names:
        if name in table.colnames:
            return table[name][0]
    return None


def get_spectra_info_from_gaia_id(source_id, nersc=False):
    if GAIA_PARQUET_FNAME is None or GAIA_BIN_FNAME is None:
        raise NotImplementedError(
            'Matterhorn Gaia index filenames are not configured yet')

    gaia_index = get_gaia_index(GAIA_PARQUET_FNAME,
                                GAIA_BIN_FNAME,
                                cache_dir=os.path.dirname(
                                    os.path.abspath(__file__)),
                                nersc=nersc)

    res = gaia_index.search_id(source_id)
    if res is None:
        raise ValueError('object not found')
    print('Found %d spectra, returning 1st one' % (len(res['SURVEY'])))
    return atpy.Table(res)


@lrudecorator(100)
def get_specs(gaia_edr3_source_id=None,
              tileid=None,
              night=None,
              fiber=None,
              hpx=None,
              uniqpix=None,
              nside=DEFAULT_NSIDE,
              targetid=None,
              expid=None,
              group_type='pix',
              dataset='matterhorn',
              spec_type='coadd',
              survey=None,
              program=None,
              subsurvey=None,
              spectrograph=None,
              mask=False,
              ivar=False,
              fibermap=False,
              nersc=False):
    """
    Get DESI spectra for a single object.

    For Matterhorn coadds, the pixel path uses the MOC UNIQ encoding
    ``uniqpix = 4 * nside**2 + healpix``. Pass ``uniqpix`` directly, or pass
    ``hpx`` with ``nside``.

    Parameters
    ----------
    tileid: int
    night: int or string
         The night identifier (i.e. 20200220 or 'all' or 'deep' for coadds)
    fiber: int
    targetid: int (optional)
    expid: int (optional)
    nersc: bool
         If true we assume we are running on NERSC and we fetch local file
    group_type: str
        The group type (pix or tiles/cumulative)
    program: str
        The program name (i.e. backup/bright/dark)
    survey: string
        The name of the survey (i.e. main/sv1 etc)
    gaia_edr3_source_id: int
         The gaia edr3 source id
    hpx: int
        The HEALPix id, converted to uniqpix using nside
    uniqpix: int
        The MOC UNIQ pixel id used by Matterhorn
    nside: int
        The HEALPix NSIDE for converting hpx to uniqpix
    spec_type: str
        The type of the spectra (i.e. coadd/cframe/spectra)
    dataset: str
        The dataset (i.e. matterhorn/daily)
    coadd: bool
         If true read coadded spectra
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
        where each dictionary has keywords b_wavelength, r_wavelength,
        z_wavelength, b_flux, b_mask, b_ivar
    """
    if nersc:
        user, pwd = None, None
    else:
        user, pwd = get_desi_login_password()

    if spec_type not in ['coadd', 'cframe', 'spectra']:
        raise Exception('unknown')
    if group_type not in ['exposure', 'pix', 'tiles/cumulative', 'tiles']:
        raise Exception('unknown')
    if subsurvey is not None:
        print('Warning subsurvey is deprecated, use program keyword')
        program = subsurvey
    if group_type != 'pix':
        if fiber is None:
            raise Exception(
                'Fiber must be specified as it is needed to identify the ' +
                'spectrograph')
        if spectrograph is None:
            spectrograph = fiber // 500

    if gaia_edr3_source_id is not None:
        res = get_spectra_info_from_gaia_id(gaia_edr3_source_id, nersc=nersc)
        survey = res['SURVEY'][0]
        program = res['PROGRAM'][0]
        targetid = res['TARGETID'][0]
        uniqpix = _first_table_value(res, 'UNIQPIX', 'uniqpix')
        if uniqpix is None:
            hpx = _first_table_value(res, 'HPXPIXEL', 'HEALPIX', 'healpix',
                                     'hpx')
            nside_from_index = _first_table_value(res, 'HPXNSIDE', 'NSIDE',
                                                  'nside')
            if nside_from_index is not None:
                nside = nside_from_index

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
    elif group_type == 'pix':
        uniqpix = _resolve_uniqpix(uniqpix=uniqpix, hpx=hpx, nside=nside)
        if uniqpix is None:
            raise ValueError('uniqpix or hpx must be specified')
        fname = f'{spec_type}-{survey}-{program}-{uniqpix}{fname_end}'
        url = (f'{data_desi}/spectro/redux/{dataset}/{survey}/'
               f'{program}/{uniqpix//100}/{uniqpix}/{fname}')
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
                      uniqpix=None,
                      nside=DEFAULT_NSIDE,
                      coadd=True,
                      survey=None,
                      program=None,
                      spec_type='coadd',
                      group_type='pix',
                      run='260613',
                      model_type='rvmod',
                      dataset='matterhorn',
                      nersc=False):
    """
    Get RVSpecfit models.

    For Matterhorn coadds, the pixel path uses ``uniqpix`` rather than
    ``healpix``. Pass ``uniqpix`` directly, or pass ``hpx`` with ``nside``.

    Parameters
    ----------
    gaia_edr3_source_id: int
         The gaia edr3 source id
    tileid: int
    night: int
    fiber: int
    targetid: int
    hpx: int
        The HEALPix id, converted to uniqpix using nside
    uniqpix: int
        The MOC UNIQ pixel id used by Matterhorn
    nside: int
        The HEALPix NSIDE for converting hpx to uniqpix
    group_type: str
        The group type (pix or tiles/cumulative)
    program: str
        The program name (i.e. backup/bright/dark)
    expid: int (optional)
    coadd: bool
         If true read coadded spectra
    model_type: str
         The model type (i.e. rvmod or rvjmod)
    nersc: bool
         If true we assume we are running on NERSC and we fetch local file
    run: string
         The string identifying a software run
    dataset: the dataset fitted (i.e. andes/sv_daily)

    Returns
    -------
    ret: list(dict)
        The list of dictionaries where each dictionary has keywords
        b_wavelength, r_wavelength, z_wavelength, b_model etc
    """
    if nersc:
        user, pwd = None, None
    else:
        user, pwd = get_desi_login_password()
    if gaia_edr3_source_id is not None:
        res = get_spectra_info_from_gaia_id(gaia_edr3_source_id, nersc=nersc)
        survey = res['SURVEY'][0]
        program = res['PROGRAM'][0]
        targetid = res['TARGETID'][0]
        uniqpix = _first_table_value(res, 'UNIQPIX', 'uniqpix')
        if uniqpix is None:
            hpx = _first_table_value(res, 'HPXPIXEL', 'HEALPIX', 'healpix',
                                     'hpx')
            nside_from_index = _first_table_value(res, 'HPXNSIDE', 'NSIDE',
                                                  'nside')
            if nside_from_index is not None:
                nside = nside_from_index

    if spec_type not in ['coadd', 'cframe', 'spectra']:
        raise Exception('unknown')
    if group_type not in ['exposure', 'pix', 'tiles/cumulative', 'tiles']:
        raise Exception('unknown')
    if group_type != 'pix':
        if fiber is None:
            raise ValueError(
                'Fiber must be specified as it is needed to identify the ' +
                'spectrograph')
        spectrograph = fiber // 500
    if group_type == 'tiles/cumulative':
        night1 = f'thru{night}'
    else:
        night1 = night

    if nersc:
        data_desi = '/global/cfs/cdirs/desi/'
    else:
        data_desi = 'https://data.desi.lbl.gov/desi/'

    data_desi = data_desi + f'science/mws/redux/{dataset}/rv_output/{run}/'
    if group_type == 'tiles/cumulative':
        fname = (f'{model_type}_{spec_type}-{spectrograph}-{tileid}-{night1}'
                 '.fits')
        url = (f'{data_desi}/tiles/cumulative/' + f'{tileid}/{night}/{fname}')
    elif group_type == 'tiles':
        fname = (f'{model_type}_{spec_type}-{spectrograph}-{tileid}-{night1}'
                 '.fits')
        url = (f'{data_desi}/tiles/{tileid}/'
               f'{night}/{fname}')
    elif group_type == 'pix':
        uniqpix = _resolve_uniqpix(uniqpix=uniqpix, hpx=hpx, nside=nside)
        if uniqpix is None:
            raise ValueError('uniqpix or hpx must be specified')
        fname = f'{model_type}_{spec_type}-{survey}-{program}-{uniqpix}.fits'
        url = (f'{data_desi}/{survey}/'
               f'{program}/{uniqpix//100}/{uniqpix}/{fname}')
        print(url)

    return read_models(url, targetid, fiber=fiber, user=user, pwd=pwd)
