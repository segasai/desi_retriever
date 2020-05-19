# desi_retriever

To install just do 

`pip install https://github.com/segasai/desi_retriever/archive/master.zip`

You also need to set-up the DESI user/password. You need to put the login password in the file 
$HOME/.desi_http_user in the format login:password

Then to fetch a given spectrum from coadd you can do 
 
`D = desi_retriever.andes.get_specs(
    tileid=63075, night=20200219, fiber=0, targetid=35190987189916575, coadd=True)`
or to get list of multiple exposures

`D = desi_retriever.andes.get_specs(tileid=63075, night=20200219, fiber=0, targetid=35190987189916575)`

You could also get the RVSpecFit models by using 

`Dmod = desi_retriever.andes.get_rvspec_models(
    tileid=63075, night=20200219, fiber=0, targetid=35190987189916575,coadd=True)`

More examples are in jupyter notebooks in example folder
