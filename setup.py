from __future__ import print_function
import os
from setuptools import setup
import glob
import subprocess


def get_revision():
    """
    Get the git revision of the code

    Returns:
    --------
    revision : string
        The string with the git revision
    """
    try:
        tmpout = subprocess.Popen(
            'cd ' + os.path.dirname(__file__) +
            ' ; git log -n 1 --pretty=format:%H -- setup.py',
            shell=True,
            bufsize=80,
            stdout=subprocess.PIPE).stdout
        revision = tmpout.read().decode()[:6]
        return revision
    except:
        return ''


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


VERSION = open('version.txt').read().rstrip()
VERSION1 = get_revision()
if VERSION1 != '':
    VERSION = VERSION + '+dev' + VERSION1
with open('py/desi_retriever/_version.py', 'w') as fd:
    print('version="%s" ' % VERSION, file=fd)
setup(
    name="desi_retriever",
    version=VERSION,
    author="Sergey Koposov",
    author_email="skoposov@ed.ac.uk",
    description="Fetch DESI spectra of individual objects by tile/night/fiber",
    license="BSD",
    keywords="DESI spectrum",
    url="https://github.com/segasai/desi_retriever",
    packages=[
        'desi_retriever',
        'desi_retriever/andes',
        'desi_retriever/blanc',
        'desi_retriever/denali',
        'desi_retriever/daily',
        'desi_retriever/fuji',
        'desi_retriever/iron',
        'desi_retriever/jura',
        'desi_retriever/loa',
        'desi_retriever/dr1',
        'desi_retriever/utils',
    ],
    #scripts = [fname for fname in glob.glob(os.path.join('bin', '*'))],
    package_dir={'': 'py'},
    #    package_data={'sqlutilpy':['tests']},
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    install_requires=open('requirements.txt').readlines(),
    classifiers=[
        "Development Status :: 3 - Alpha", "Topic :: Utilities",
        "License :: OSI Approved :: BSD License"
    ],
)
