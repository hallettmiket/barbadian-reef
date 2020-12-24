# ZENODO DATA DOWNLOAD

YOUR_DATA <- "data"
# Modify to put these large files somewhere outside of your git repo?

library(rprojroot)
root <- rprojroot::find_root(".git/index"); 
dir.create(YOUR_DATA)
download.file("zenodolink" , 
              file.path(YOUR_DATA, "reef-data.tar.gz"))
untar( file.path(YOUR_DATA,"data/reef-data.tar.gz"))