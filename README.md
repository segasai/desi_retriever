# desi_retriever

This is a tool that allows you to fetch desi spectra without manually downloading them 

## Installation

To install just do 

`pip install https://github.com/segasai/desi_retriever/archive/master.zip`

### Set up the DESI user / password

You need to put in the standard DESI web login/password in the home folder in the file:
$HOME/.desi_http_user
in the following format
login:password

The login password are the same ones as used to access data.desi.lbl.gov 

## Usage

Then to fetch a given spectrum from coadd you can do 
 
`D = desi_retriever.andes.get_specs(
    tileid=63075, night=20200219, fiber=0, targetid=35190987189916575, coadd=True)`
or to get list of multiple exposures

`D = desi_retriever.andes.get_specs(tileid=63075, night=20200219, fiber=0, targetid=35190987189916575)`

You could also get the RVSpecFit models by using 

`Dmod = desi_retriever.andes.get_rvspec_models(
    tileid=63075, night=20200219, fiber=0, targetid=35190987189916575,coadd=True)`

More examples are in jupyter notebooks in example folder
