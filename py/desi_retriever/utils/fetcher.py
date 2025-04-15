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
        si.cache[url] = copy.copy(fp._cache)
        return rets


def read_models(url, targetid, fiber=None, expid=None, user=None, pwd=None):

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
        si.cache[url] = copy.copy(fp._cache)
        return rets
