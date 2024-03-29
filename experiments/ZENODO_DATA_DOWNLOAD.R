# ZENODO DATA DOWNLOAD


library(rprojroot)
root <- rprojroot::find_root(".git/index"); 

# FIRST GO TO src/init.R and specifiy where you want to store your data (variable YOUR_DATA).
# perhaps outside of the git repo as the files are large.

source(file.path(root, "src/init.R"))

dir.create(ZEN_DATA)

library(zen4R)
options(timeout = max(1000, getOption("timeout")))
download_zenodo(doi="10.5281/zenodo.4357318" , 
              path = file.path(ZEN_DATA))

untar( file.path(ZEN_DATA,"data/reef-data.tar.gz"))