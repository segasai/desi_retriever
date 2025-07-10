import urllib3
import os
from pylru import lrudecorator
from ..utils import (read_spectra, read_models, get_gaia_index,
                     get_desi_login_password)

urllib3.disable_warnings()

GAIA_PARQUET_FNAME = 'gaia-index-loa-coadd_241127.parquet'
GAIA_BIN_FNAME = 'gaia-index-loa-coadd_241127.bin'


@lrudecorator(100)
def get_specs(gaia_edr3_source_id=None,
              tileid=None,
              night=None,
              fiber=None,
              hpx=None,
              targetid=None,
              expid=None,
              group_type='healpix',
              dataset='loa',
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
    Typically if you are getting a coadded object, you
    need to provide targetid, healpix, survey, program.

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
        The group type (healpix or tiles/cumulative)
    program: str
        The program name (i.e. backup/bright/dark)
    survey: string
        The name of the survey (i.e. main/sv1 etc)
    gaia_edr3_source_id: int
         The gaia edr3 source id
    hpx: int
        The healpix id
    spec_type: str
        The type of the spectra (i.e. coadd/cframe/spectra)
    dataset: str
        The dataset (i.e. loa/daily)
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
        where each dictionary
        has keywords b_wavelength, r_wavelength, z_wavelength
        b_flux, b_mask, b_ivar

    """
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
        gaia_index = get_gaia_index(GAIA_PARQUET_FNAME,
                                    GAIA_BIN_FNAME,
                                    cache_dir=os.path.dirname(
                                        os.path.abspath(__file__)),
                                    nersc=nersc)
        res = gaia_index.search_id(gaia_edr3_source_id)
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
    return read_spectra(url,
                        targetid,
                        fiber=fiber,
                        expid=expid,
                        mask=mask,
                        ivar=ivar,
                        fibermap=fibermap,
                        user=user,
                        pwd=pwd)


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
                      run='241119',
                      model_type='rvmod',
                      dataset='loa',
                      nersc=False):
    """
    Get RVSpecfit models

    Parameters
    ----------
    gaia_edr3_source_id: int
         The gaia edr3 source id

    tileid: int
    night: int
    fiber: int
    targetid: int
    hpx: int
        The healpix id
    group_type: str
        The group type (healpix or tiles/cumulative)
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
        The list of dictionaries where each dictionary
        has keywords b_wavelength, r_wavelength, z_wavelength
        b_model etc
    """
    if subsurvey is not None:
        print('Warning subsurvey keyword is deprecated, use program')
    user, pwd = get_desi_login_password()
    if gaia_edr3_source_id is not None:
        gaia_index = get_gaia_index(GAIA_PARQUET_FNAME,
                                    GAIA_BIN_FNAME,
                                    cache_dir=os.path.dirname(
                                        os.path.abspath(__file__)),
                                    nersc=nersc)

        res = gaia_index.search_id(gaia_edr3_source_id)
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
    elif group_type == 'healpix':
        fname = f'{model_type}_{spec_type}-{survey}-{program}-{hpx}.fits'
        url = (f'{data_desi}/healpix/{survey}/'
               f'{program}/{hpx//100}/{hpx}/{fname}')
        print(url)

    return read_models(url, targetid, fiber=fiber, user=user, pwd=pwd)
