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


def get_desi_login_password():
    if si.DESI_USER is None:
        config = os.environ['HOME'] + '/.desi_http_user'
        if not os.path.exists(config):
            raise Exception('''You need to specify the DESI_USER/DESI_PASSWD.
put them in $HOME/.desi_http_user like that
username:password
''')
        user, pwd = open(config).read().rstrip().split(':')
        si.DESI_USER, si.DESI_PASSWD = user, pwd
    return si.DESI_USER, si.DESI_PASSWD


class GaiaIndex:

    def __init__(self, gaiaid, targetid, survey, subsurvey, hpx, row):
        self.gaiaid = gaiaid
        self.targetid = targetid
        self.survey = survey
        self.subsurvey = subsurvey
        self.hpx = hpx
        self.row = row

    def search_id(self, gaiaid_new):
        xid = np.searchsorted(self.gaiaid, gaiaid_new)
        if xid >= len(self.gaiaid) or xid < 0:
            return None
        if self.gaiaid[xid] == gaiaid_new:
            return dict(targetid=self.targetid[xid],
                        survey=self.survey[xid],
                        subsurvey=self.subsurvey[xid],
                        hpx=self.hpx[xid],
                        row=self.row[xid])
        else:
            return None


def fetch_gaia_index():
    if si.gaiaIndex is not None:
        return
    D = pyfits.getdata(
        fsspec.open(
            'https://data.desi.lbl.gov/desi/users/koposov/gaiaid_db/indexes/gaia-fuji-coadd-hpx-index.fits',
            'rb',
            auth=aiohttp.BasicAuth(si.DESI_USER, si.DESI_PASSWD)).open())
    si.gaiaIndex = GaiaIndex(D['EDR3_SOURCE_ID'], D['TARGETID'], D['survey'],
                             D['subsurvey'], D['hpx'], D['row'])


def read_spectra(url, user, pwd, targetid, expid, fiber, mask, ivar):
    kw = dict(auth=(user, pwd), verify=False)
    block_size = 2880 * 10  # caching block
    with httpio.open(url, block_size=block_size, **kw) as fp:
        if url in si.cache:
            fp._cache = si.cache[url]
        hdus = pyfits.open(fp)

        ftab = atpy.Table(hdus['FIBERMAP'].data)

        if expid is not None:
            xind = ftab['EXPID'] == expid
        else:
            xind = np.ones(len(ftab), dtype=bool)
        if targetid is not None:
            xids = np.nonzero((ftab['TARGETID'] == targetid) & xind)[0]
        else:
            xids = np.nonzero((ftab['FIBER'] == fiber) & xind)[0]
        if len(xids) == 0:
            print('no spectra')
            return []

        waves = {}
        for arm in 'BRZ':
            waves[arm] = hdus[arm + '_WAVELENGTH'].data

        fluxes = {}
        for arm in 'BRZ':
            fluxes[arm] = hdus[arm + '_FLUX'].section

        masks = {}
        if mask:
            for arm in 'BRZ':
                masks[arm] = hdus[arm + '_MASK'].section

        ivars = {}
        if ivar:
            for arm in 'BRZ':
                ivars[arm] = hdus[arm + '_IVAR'].section

        rets = []
        for xid in xids:
            ret = {}
            for arm in 'BRZ':
                ret[arm.lower() + '_wavelength'] = waves[arm]
                ret[arm.lower() + '_flux'] = fluxes[arm][xid, :]
            if mask:
                for arm in 'BRZ':
                    ret[arm.lower() + '_mask'] = masks[arm][xid, :]
            if ivar:
                for arm in 'BRZ':
                    ret[arm.lower() + '_ivar'] = ivars[arm][xid, :]

            rets.append(ret)
        si.cache[url] = copy.copy(fp._cache)
        return rets


def read_models(url, user, pwd, targetid, fiber, expid=None):

    block_size = 2880 * 10  # caching block
    kw = dict(auth=(user, pwd), verify=False)

    with httpio.open(url, block_size=block_size, **kw) as fp:
        if url in si.cache:
            fp._cache = si.cache[url]
        hdus = pyfits.open(fp)
        ftab = atpy.Table(hdus['FIBERMAP'].data)

        if expid is not None:
            xind = ftab['EXPID'] == expid
        else:
            xind = np.ones(len(ftab), dtype=bool)
        if targetid is not None:
            xids = np.nonzero((ftab['TARGETID'] == targetid) & xind)[0]
        else:
            xids = np.nonzero((ftab['FIBER'] == fiber) & xind)[0]

        if len(xids) == 0:
            print('no spectra')
            return []

        waves = {}
        for arm in 'BRZ':
            waves[arm] = hdus[arm + '_WAVELENGTH'].data

        models = {}
        for arm in 'BRZ':
            models[arm] = hdus[arm + '_MODEL'].section

        rets = []

        for xid in xids:
            ret = {}
            for arm in 'BRZ':
                ret[arm.lower() + '_wavelength'] = waves[arm]
                ret[arm.lower() + '_model'] = models[arm][xid, :]
            rets.append(ret)
        si.cache[url] = copy.copy(fp._cache)
        return rets


@lrudecorator(100)
def get_specs(gaia_edr3_source_id=None,
              tileid=None,
              night=None,
              fiber=None,
              hpx=None,
              targetid=None,
              expid=None,
              group_type='healpix',
              dataset='fuji',
              spec_type='coadd',
              survey=None,
              subsurvey=None,
              spectrograph=None,
              mask=False,
              ivar=False):
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
    user, pwd = get_desi_login_password()

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
        survey = res['survey']
        subsurvey = res['subsurvey']
        hpx = res['hpx']
        targetid = res['targetid']

    if group_type == 'tiles/cumulative':
        night1 = f'thru{night}'
    else:
        night1 = night

    data_desi = 'https://data.desi.lbl.gov/desi/'
    if group_type == 'tiles/cumulative':
        fname = f'{spec_type}-{spectrograph}-{tileid}-{night1}.fits'
        url = (f'{data_desi}/spectro/redux/{dataset}/tiles/cumulative/' +
               f'{tileid}/{night}/{fname}')
    elif group_type == 'tiles':
        fname = f'{spec_type}-{spectrograph}-{tileid}-{night1}.fits'
        url = (f'{data_desi}/spectro/redux/{dataset}/tiles/{tileid}/'
               f'{night}/{fname}')
    elif group_type == 'healpix':
        fname = f'{spec_type}-{survey}-{subsurvey}-{hpx}.fits'
        url = (f'{data_desi}/spectro/redux/{dataset}/healpix/{survey}/'
               f'{subsurvey}/{hpx//100}/{hpx}/{fname}')
        print(url)

    else:
        raise Exception('oops')
    return read_spectra(url, user, pwd, targetid, expid, fiber, mask, ivar)


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
                      spec_type='coadd',
                      group_type='healpix',
                      run='220309',
                      dataset='fuji'):
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
    user, pwd = get_desi_login_password()
    if gaia_edr3_source_id is not None:
        fetch_gaia_index()
        res = si.gaiaIndex.search_id(gaia_edr3_source_id)
        if res is None:
            raise ValueError('object not found')
        survey = res['survey']
        subsurvey = res['subsurvey']
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

    data_desi = (f'https://data.desi.lbl.gov/desi/science/mws/redux/'
                 f'{dataset}/rv_output/{run}/')
    if group_type == 'tiles/cumulative':
        fname = f'rvmod_`{spec_type}-{spectrograph}-{tileid}-{night1}.fits'
        url = (f'{data_desi}/tiles/cumulative/' + f'{tileid}/{night}/{fname}')
    elif group_type == 'tiles':
        fname = f'rvmod_{spec_type}-{spectrograph}-{tileid}-{night1}.fits'
        url = (f'{data_desi}/tiles/{tileid}/'
               f'{night}/{fname}')
    elif group_type == 'healpix':
        fname = f'rvmod_{spec_type}-{survey}-{subsurvey}-{hpx}.fits'
        url = (f'{data_desi}/healpix/{survey}/'
               f'{subsurvey}/{hpx//100}/{hpx}/{fname}')
        print(url)

    return read_models(url, user, pwd, targetid, fiber)
