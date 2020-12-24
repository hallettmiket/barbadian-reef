library(rprojroot)
library(xtable); library(ggplot2); library(vcd); library(MASS); library(FNN); library(rlang)
library(boot); library(ggrepel); library(ggthemes); library(ggpubr); library(pals)

source(file.path(root, "src/functions.R"))
source(file.path(root, "src/load.tree.R"))

REEF_DIR <- "/home/data/refined/reef/R/R"
mystat<- function( data, indices ) {return(mean(data[indices]))}


date <- "april.19"
vers <- "1.0"

original <- tree <- load.tree( path=REEF_DIR, vers="1.0", tidyup=TRUE, adjust = FALSE )
gs_adjustment <- load.adjustments( path=REEF_DIR, vers="1.0" )
 
