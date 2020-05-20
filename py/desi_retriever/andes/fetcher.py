import astropy.table as atpy
import httpio
import numpy as np
import astropy.io.fits as pyfits
import urllib3
import os
urllib3.disable_warnings()


class si:
    DESI_USER = None
    DESI_PASSWD = None


def get_desi_login_password():
    if si.DESI_USER is None:
        config = os.environ['HOME'] + '/.desi_http_user'
        user, pwd = open(config).read().rstrip().split(':')
        si.DESI_USER, si.DESI_PASSWD = user, pwd
    return si.DESI_USER, si.DESI_PASSWD


def get_specs(tileid=None,
              night=None,
              fiber=None,
              targetid=None,
              expid=None,
              coadd=False):
    """
    Get DESI spectra 
    
    Parameters
    ----------
    tileid: int
    night: int
    fiber: int
    targetid: int
    expid: int (optional)
    coadd: bool 
         If true read coadded spectra

    Returns
    -------
    ret: list(dict)
        The list of dictionaries where each dictionary
        has keywords b_wavelength, r_wavelength, z_wavelength
        b_flux etc

    """
    if coadd:
        prefix = 'coadd'
    else:
        prefix = 'spectra'
    spectrograph = fiber // 500
    url = f'https://data.desi.lbl.gov/desi/spectro/redux/andes/tiles/{tileid}/{night}/{prefix}-{spectrograph}-{tileid}-{night}.fits'
    user, pwd = get_desi_login_password()
    kw = dict(auth=(user, pwd), verify=False)

    with httpio.open(url, **kw) as fp:
        hdus = pyfits.open(fp)
        ftab = atpy.Table(hdus['FIBERMAP'].data)

        if expid is not None:
            xind = ftab['EXPID'] == expid
        else:
            xind = np.ones(len(ftab), dtype=bool)
        xids = np.nonzero((ftab['TARGETID'] == targetid) & xind)[0]
        if len(xids) == 0:
            print('no spectra')
            return []
        bwave = hdus['B_WAVELENGTH'].data
        rwave = hdus['R_WAVELENGTH'].data
        zwave = hdus['Z_WAVELENGTH'].data

        rets = []
        bdata = hdus['B_FLUX'].data
        rdata = hdus['R_FLUX'].data
        zdata = hdus['Z_FLUX'].data
        rets = []
        for xid in xids:
            bdat = bdata[xid, :]
            rdat = rdata[xid, :]
            zdat = zdata[xid, :]

            ret = dict(b_wavelength=bwave,
                       r_wavelength=rwave,
                       z_wavelength=zwave,
                       b_flux=bdat,
                       r_flux=rdat,
                       z_flux=zdat)
            rets.append(ret)
        return rets


def get_rvspec_models(tileid=None,
                      night=None,
                      fiber=None,
                      targetid=None,
                      expid=None,
                      coadd=False,
                      run='200507'):
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

    Returns
    -------
    ret: list(dict)
        The list of dictionaries where each dictionary
        has keywords b_wavelength, r_wavelength, z_wavelength
        b_model etc
    """

    if coadd:
        prefix = 'rvmod_spectra'
    else:
        prefix = 'rvmod_coadd'
    spectrograph = fiber // 500
    url = f'https://data.desi.lbl.gov/desi/science/mws/redux/andes/{run}/rv_output/{tileid}/{night}/{prefix}-{spectrograph}-{tileid}-{night}.fits'
    user, pwd = get_desi_login_password()
    kw = dict(auth=(user, pwd), verify=False)
    with httpio.open(url, **kw) as fp:
        hdus = pyfits.open(fp)
        ftab = atpy.Table(hdus['FIBERMAP'].data)

        if expid is not None:
            xind = ftab['EXPID'] == expid
        else:
            xind = np.ones(len(ftab), dtype=bool)
        xids = np.nonzero((ftab['TARGETID'] == targetid) & xind)[0]

        if len(xids) == 0:
            print('no spectra')
            return []
        bwave = hdus['B_WAVELENGTH'].data
        rwave = hdus['R_WAVELENGTH'].data
        zwave = hdus['Z_WAVELENGTH'].data

        rets = []
        bdata = hdus['B_MODEL'].data
        rdata = hdus['R_MODEL'].data
        zdata = hdus['Z_MODEL'].data
        rets = []
        for xid in xids:
            bdat = bdata[xid, :]
            rdat = rdata[xid, :]
            zdat = zdata[xid, :]

            ret = dict(b_wavelength=bwave,
                       r_wavelength=rwave,
                       z_wavelength=zwave,
                       b_model=bdat,
                       r_model=rdat,
                       z_model=zdat)
            rets.append(ret)
        return rets
