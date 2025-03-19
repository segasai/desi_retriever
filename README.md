# desi_retriever

This is a tool that allows you to fetch desi spectra without manually downloading them.

You may also consider using SPARCL 
https://astrosparcl.datalab.noirlab.edu/
for another interface.


## Installation

To install just do 

`pip install -U desi_retriever@git+https://github.com/segasai/desi_retriever`

### Authentication

desi_retriever works with public DR1 data and private DESI data.
If you want to use non-public DESI data, you must set up the DESI user / password

You need to put in the standard DESI web login/password in the home folder in the file:
$HOME/.desi_http_user
in the following format
login:password

The login password are the same ones as used to access data.desi.lbl.gov 

## Usage

Then to fetch a given spectrum from coadd you can do 

```
import desi_retriever.dr1
SP = desi_retriever.dr1.get_specs(survey='sv1',
               program='dark',
               hpx=17683,
               targetid=39627652591521181)[0]
```

You could also get the RVSpecFit models by using 

`Dmod = desi_retriever.dr1.get_rvspec_models(survey='sv1',
               program='dark',
               hpx=17683,
               targetid=39627652591521181)[0]`

More examples are in jupyter notebooks in example folder.


