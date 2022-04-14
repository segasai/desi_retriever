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
            'https://data.desi.lbl.gov/desi/users/koposov/gaiaid_db/gaia-everest-coadd-hpx-index.fits',
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


@lrudecorator(100)
def get_specs(gaia_edr3_source_id=None,
              tileid=None,
              night=None,
              fiber=None,
              hpx=None,
              targetid=None,
              expid=None,
              coadd=True,
              group_type='healpix',
              dataset='everest',
              survey=None,
              subsurvey=None,
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

    if coadd:
        prefix = 'coadd'
    else:
        prefix = 'spectra'
    if coadd_type != 'healpix':
        if fiber is None:
            raise Exception(
                'Fiber must be specified as it is needed to identify the ' +
                'spectrograph')
        spectrograph = fiber // 500

    if coadd_type == 'cumulative':
        night1 = f'thru{night}'
    else:
        night1 = night
    if gaia_edr3_source_id is not None:
        fetch_gaia_index()
        res = si.gaiaIndex.search_id(gaia_edr3_source_id)
        if res is None:
            raise ValueError('object not found')
        survey = res['survey']
        subsurvey = res['subsurvey']
        hpx = res['hpx']
        targetid = res['targetid']
    if tileid is not None:
        url = f'https://data.desi.lbl.gov/desi/spectro/redux/{dataset}/tiles/{coadd_type}/{tileid}/{night}/{prefix}-{spectrograph}-{tileid}-{night1}.fits'
    elif hpx is not None:
        url = f'https://data.desi.lbl.gov/desi/spectro/redux/{dataset}/healpix//{survey}/{subsurvey}/{hpx//100}/{hpx}/{prefix}-{survey}-{subsurvey}-{hpx}.fits'
    return read_spectra(url, user, pwd, targetid, expid, fiber, mask, ivar)


@lrudecorator(100)
def get_rvspec_models(gaia_edr3_source_id=None,
                      tileid=None,
                      night=None,
                      fiber=None,
                      targetid=None,
                      expid=None,
                      coadd=True,
                      coadd_type='healpix',
                      run='210803',
                      dataset='everest'):
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
        if survey[:2] == 'sv':
            survey_short = 'sv'
    if coadd:
        prefix = 'rvmod_coadd'
    else:
        prefix = 'rvmod_spectra'
    if coadd_type != 'healpix':
        if fiber is None:
            raise ValueError(
                'Fiber must be specified as it is needed to identify the ' +
                'spectrograph')
        spectrograph = fiber // 500
    if coadd_type == 'cumulative':
        night1 = f'thru{night}'
    else:
        night1 = night
    if tileid is not None:
        url = f'https://data.desi.lbl.gov/desi/science/mws/redux/{dataset}/tiles/{coadd_type}/{tileid}/{night}/{prefix}-{spectrograph}-{tileid}-{night1}.fits'
    elif hpx is not None:
        url = f'https://data.desi.lbl.gov/desi/science/mws/redux/{dataset}/rv_output/{run}/healpix_{survey_short}_{subsurvey}/{hpx//100}/{hpx}/{prefix}-{survey}-{subsurvey}-{hpx}.fits'

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
