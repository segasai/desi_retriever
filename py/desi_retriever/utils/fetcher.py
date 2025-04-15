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

urllib3.disable_warnings()


class si:
    cache = lrucache(100)


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
