import astropy.table as atpy
import copy
import httpio
import numpy as np
import astropy.io.fits as pyfits
import urllib3
import os
from pylru import lrucache
import pyarrow.parquet as pq
import fsspec
import pickle
import traceback
import requests
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


def read_spectra(url,
                 targetid,
                 fiber=None,
                 expid=None,
                 mask=False,
                 ivar=False,
                 fibermap=False,
                 user=None,
                 pwd=None):
    """
    Read the spectra from the given url
    Parameters
    ----------
    url: str
        The url to the file (or file path)
    targetid: int
        The targetid of the object
    fiber: int
        The fiber number (optional)
    expid: int
        The exposure id (optional)
    mask: bool
        Whether to read the mask (optional)
    ivar: bool
        Whether to read the inverse variance (optional)
    fibermap: bool
        Whether to read the fibermap (optional)
    user: str
        The username for authentication (optional)
    pwd: str
        The password for authentication (optional)
    Returns
    -------
    Returns a list of dictionaries with the following keys:
        - fibermap: the fibermap of the object
        - b_wavelength: the wavelength of the blue arm
        - b_flux: the flux of the blue arm
        - b_mask: the mask of the blue arm (optional)
        - b_ivar: the inverse variance of the blue arm (optional)
        - r_wavelength: the wavelength of the red arm
        - r_flux: the flux of the red arm
        - r_mask: the mask of the red arm (optional)
        - r_ivar: the inverse variance of the red arm (optional)
        - z_wavelength: the wavelength of the z arm
        - z_flux: the flux of the z arm
        - z_mask: the mask of the z arm (optional)
        - z_ivar: the inverse variance of the z arm (optional)
    If there is only one object, the list will have only one element.
        """
    kw = dict(verify=False)
    if user is not None:
        kw['auth'] = (user, pwd)
    kw['block_size'] = 2880 * 10  # caching block
    local_mode = False
    if url[0] == '/':
        local_mode = True
    with open(url, 'rb') if local_mode else httpio.open(url, **kw) as fp:
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
            print('Warning no spectra was found with provided info')
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
            if fibermap:
                ret['fibermap'] = ftab[xid]
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
        if not local_mode:
            si.cache[url] = copy.copy(fp._cache)
        return rets


def read_models(url, targetid, fiber=None, expid=None, user=None, pwd=None):
    """
    Read the models from the given url

    Parameters
    ----------
    url: str
        The url to the file
    targetid: int
        The targetid of the object
    fiber: int
        The fiber number (optional)
    expid: int
        The exposure id (optional)
    user: str
        The username for authentication (optional)
    pwd: str
        The password for authentication (optional)
    """
    kw = dict(verify=False)
    kw['block_size'] = 2880 * 10  # caching block
    if user is not None:
        kw['auth'] = (user, pwd)
    local_mode = False
    if url[0] == '/':
        local_mode = True
    with open(url, 'rb') if local_mode else httpio.open(url, **kw) as fp:
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
        if not local_mode:
            si.cache[url] = copy.copy(fp._cache)
        return rets


class GaiaIndex:

    def __init__(self,
                 key_arr,
                 pos1,
                 pos2,
                 bin_url=None,
                 auth=None,
                 columns=None):
        self.key_arr = key_arr
        self.pos1 = pos1
        self.pos2 = pos2
        self.bin_url = bin_url
        self.auth = auth
        self.columns = columns

    def search_id(self, key_val):
        xid = np.searchsorted(self.key_arr, key_val)
        xid_right = np.searchsorted(self.key_arr, key_val, 'right')
        if xid >= len(self.key_arr) or xid < 0:
            return None
        if self.key_arr[xid] == key_val:
            Ds = []
            with fsspec.open(self.bin_url, 'rb', auth=self.auth).open() as fp:
                for cur_xid in range(xid, xid_right):
                    fp.seek(self.pos1[cur_xid])
                    D = pickle.loads(
                        fp.read(self.pos2[cur_xid] - self.pos1[cur_xid]))
                    Ds.append(D)
            ret = {}
            for i, c in enumerate(self.columns):
                ret[c] = []
                for D in Ds:
                    ret[c].append(D[i])
            return ret
        else:
            return None


def get_gaia_index(parquet_fname, bin_fname, cache_dir=None, nersc=False):
    if si.gaiaIndex is not None:
        return si.gaiaIndex
    path0 = '/desi/users/koposov/gaiaid_db/indexes/'
    if nersc:
        login, passwd = None, None
        auth = None
        base_url = f'/global/cfs/cdirs/{path0}'
    else:
        login, passwd = get_desi_login_password()
        auth = aiohttp.BasicAuth(login, passwd)
        base_url = f'https://data.desi.lbl.gov/{path0}'
    if cache_dir is None:
        cache_dir = os.path.dirname(os.path.abspath(__file__))
    local_parquet_fname = cache_dir + '/' + parquet_fname
    parquet_url = f'{base_url}/' + parquet_fname
    bin_url = f'{base_url}/{bin_fname}'
    for i in range(2):
        try:
            with open(local_parquet_fname, 'rb') as fp:
                pqf = pq.ParquetFile(fp).read()
        except:  # noqa
            print('Downloading remote parquet file')
            try:
                with open(local_parquet_fname, 'wb') as fp_out:
                    if nersc:
                        fp_out.write(open(parquet_url, 'rb').read())
                    else:
                        with requests.get(parquet_url,
                                          auth=(login, passwd)) as rp:
                            fp_out.write(rp.content)
                            print('Successfully downloaded')
            except:  # noqa
                print('''Failed to download the gaia index file
You may want to update desi_retriever''')
                traceback.print_exc()
            continue

    with fsspec.open(bin_url, 'rb', auth=auth).open() as fp:
        header = 1000
        keys = pickle.loads(fp.read(header))

    D = {}
    for k in ['EDR3_SOURCE_ID', 'pos1', 'pos2']:
        D[k] = np.array(pqf[k])
    si.gaiaIndex = GaiaIndex(D['EDR3_SOURCE_ID'],
                             D['pos1'],
                             D['pos2'],
                             bin_url=bin_url,
                             columns=keys,
                             auth=auth)
    return si.gaiaIndex
