#
#
#
# percolate_tara function takes the processing done by include_tara() 
#   that assigns reads to the observations in the tara dataset
#   and percolates these up the tree.
#



percolate_tara <- function(parent, lb, ub) {
  
  c <- p2c( parent )
  ind <- tax_2_index( parent )
  
  #  if (length(level) < 100)
  #    cat("\nLevel ", level)
  
  
  if (length(c) > 0) {
    
    for (i in 1:length(c)) {
      tmp <- percolate_tara( c[i], lb, ub )
      res[ind, lb:ub] <<- res[ind, lb:ub] + tmp
    }
  }
  
  return( res[ind, lb:ub ] )
} # end of percolate

  


#
#
#
# include_tara function takes our NCBI taxonomy-based tree sructure and adds the 
#   tara oceans data into the data frame.
#
# Both ours and tara should be tibbles as input.
# 
# Output is a data.frame that looks like our data but adds the information for all the 
#   samples from tara.
#   It recursively percolates this informaiton up to the tree
#   assumes the root is in the first row of res
#
#   Assumes rows of res and ordered the same way as tree (original)
#
include_tara <- function( ours, tara ) {
  # Domain, Phylum, Class, Order, Family, Genus
  
  res <<- ours  # this will eventually be the output
  res$"OTU.rep" <<- NA
  lb <- ncol(res) + 1; ub <- lb + (ncol(tara) - 8)
  res[, lb:ub] <<- 0; colnames(res)[lb:ub] <<- colnames(tara)[8:ncol(tara)]
  
  confusing <- c()
  for (i in 1:nrow(tara)) {
    
    t <- tara[i,]
    found <- FALSE
    for (j in 6:1) {
      if ((!is.na(t[j])) & (t[j]!= "undef")) {
        idx_ours <- which(ours$name == as.character(t[j]))
 
        if (length(idx_ours) == 1) {
          res[idx_ours, "OTU.rep"] <<- as.character(t[7])
          res[idx_ours, lb:ub] <<- res[idx_ours, lb:ub] + as.numeric(t[8:ncol(t)])
          found <- TRUE; break
        } else if (length(idx_ours) > 1) {
          if (671 %in% idx_ours) {
            res[671, "OTU.rep"] <<- as.character(t[7])
            res[671, lb:ub] <<- res[671, lb:ub] + as.numeric(t[8:ncol(t)]) 
          } else if (10010 %in% idx_ours) {
            res[10010, "OTU.rep"] <<- as.character(t[7])
            res[10010, lb:ub] <<- res[10010, lb:ub] + as.numeric(t[8:ncol(t)]) 
          } else if (661 %in% idx_ours) {
            res[661, "OTU.rep"] <<- as.character(t[7])
            res[661, lb:ub] <<- res[661, lb:ub] + as.numeric(t[8:ncol(t)]) 
          } else if (9715 %in% idx_ours) {
            res[584, "OTU.rep"] <<- as.character(t[7])
            res[584, lb:ub] <<- res[584, lb:ub] + as.numeric(t[8:ncol(t)]) 
          } else if (840 %in% idx_ours) {
            res[840, "OTU.rep"] <<- as.character(t[7])
            res[840, lb:ub] <<- res[840, lb:ub] + as.numeric(t[8:ncol(t)]) 
          } else if (6042 %in% idx_ours) {
            res[671, "OTU.rep"] <<- as.character(t[7])
            res[671, lb:ub] <<- res[671, lb:ub] + as.numeric(t[8:ncol(t)]) 
          } else
          {
            cat("\n", idx_ours, tree[idx_ours, "name"])
            confusing <- c(confusing, t[7])
          }
        }
      }
    } # end j
    
     if (!found) { 
       print(t[7]); confusing <- c(confusing, t[7])
       # add anything not found to the root)
       res[1, lb:ub] <<- res[1, lb:ub] + as.numeric(t[8:ncol(t)])
    }

   } # end for i
  
  void <- percolate_tara( 1, lb, ub ) 
  return(list(res, confusing))
}


#
# Used by the heatmap function to select taxa for clustering
#
get_top_taxa <- function( m, lmt ) {
  
  #fav <- matrix(data= NA, nrow=lmt, ncol=ncol(m))
  return( sapply(1:ncol(m),  simplify=TRUE, FUN = function(col) {
    return( order(m[,col], decreasing=TRUE)[1:lmt] )
  }) )# end of apply
  
}

#
# Used to change nominal factors in the metadata for the tara oceans project
# into integers for heatmaps.
#
# remove some redundancies

# 10, environment, nominal
# 12, biome, nominal 
# 13, region, nominal 
# 14, MRGID, nominal 
prepare_metatara <- function(mt, barbados=TRUE) {
  rownames(mt) <- mt[,1]
 
  f <- c(); num <- c()
   for (i in 1:ncol(mt)) {
     if  (class(mt[,i])=="factor") { f <- c(f, i); mt[,i] <- as.character(mt[,i]) }
     if ((class(mt[,i])=='numeric') | (class(mt[, i])=='integer')) num <- c(num, i)
   }
   
  if (barbados) {
    mt <- prepare_barbados_sites( mt )
  }
  for (ff in f) { mt[,ff] <- factor(mt[,ff])}
  for (nn in num) { mt[,nn] <- as.numeric(mt[,nn])}
  
  
  # give Environment nicer shorter names
  mt[,10] <- factor(mt[,10], labels = c("DCM", "DCM + oxy min", "MES", "MES + oxy min", 
                                        "MIX", "SRF", "Barbados"))
   
  # now Region
  mt[, 13] <- factor(mt[,13], labels = c(
   "Indian", "Med Sea",    
   "North Atlantic", "North Pacific",
   "Red Sea",               "South Atlantic",
   "Southern",        "South Pacific", "Barbados" ))

  # now MRGID
  ttt <- matrix(unlist(strsplit(as.character(levels(mt[,14])), ")")), nrow=2 )[1,]
  ttt <- sapply(ttt, FUN = function(x) {  return(substr(x, 2, nchar(x))) })
  mt[, 14] <- factor(mt[,14], labels =  ttt )
    

  return(mt)
} # end of prepare_metatara

prepare_barbados_sites <- function(  mt ) {
  
  library(tidyverse)
  bellairs <- c(rep(NA, times=ncol(mt)))
  maycocks <- c(rep(NA, times=ncol(mt)))
  names(bellairs) <- colnames(mt)  
  names(maycocks) <- colnames(mt)


  bellairs[1] <- "br_bel"
  maycocks[1] <- "br_may"
  
  bellairs["Biome"] <- "Barbados"
  maycocks["Biome"] <- "Barbados"
  bellairs["Environment"] <- "Barbados"
  maycocks["Environment"] <- "Barbados"
  bellairs["Region"] <- "Barbados"
  maycocks["Region"] <- "Barbados"
  
  bellairs[8] <- 13.1130 #lat
  bellairs[9] <- 50.3829  #long
  maycocks[8] <- 13.1733
  maycocks[9] <- 50.3947
  
  bellairs[ "Size Fraction Up"] <- maycocks[ "Size Fraction Up" ] <- 1.6
  
  bellairs[ "Mean Date"] <- "30/01/18"; maycocks[ "Mean Date"] <- "31/01/18"
  bellairs[ "Mean Depth"] <- 5; maycocks[ "Mean Depth"] <- 5
  
  environ <- readRDS("/home/data/refined/reef/R/environ.RData" )
  options(pillar.sigfig =5)
  tmp <- environ %>% group_by(site) %>% summarise(across("Time":"Salinity(ppt)", ~mean(.x, na.rm = TRUE)))
  
  bellairs[ "Mean Temp"] <-  29.031
  maycocks[ "Mean Temp"] <-  29.134
  bellairs[ "Mean Salinity"] <- 35.033
  maycocks[ "Mean Salinity"] <- 34.783
  
  environ <- environ %>% mutate( adj_DO = `DO(mg/l)` * 1000 / 32 / 1.027 )  %>% 
                group_by(site) %>% summarise(m = mean(adj_DO))
  bellairs[ "Oxygen"] <- 172.93
  maycocks[ "Oxygen"] <- 183.79
  
  
  sea_chem <- readRDS("/home/data/refined/reef/R/sea_chemistry.RData" )
  
  #For nitrates. nitrites and Phosphate : mg/l to umol/l : following the conversion outlined in 
  # https://www.caryinstitute.org/sites/default/files/public/downloads/curriculum-project/4A1_Nitrogen_reading.pdf
  
  sea_chem <- sea_chem %>% mutate( adj_NO3_1 = (((`Nitrate#1(NO3-N)mg/l` / 1000 ) / 62.0049) / 0.000001 ) )
  sea_chem <- sea_chem %>% mutate( adj_NO3_2 = (((`Nitrate#2(NO3-N)mg/l` / 1000 ) / 62.0049) / 0.000001 ) )
  sea_chem <- sea_chem %>% mutate( adj_NO3_3 = (((`Nitrate#3(NO3-N)mg/l` / 1000 ) / 62.0049) / 0.000001 ) )
  inter <- sea_chem %>% group_by(site) %>% summarise(across("adj_NO3_1":"adj_NO3_3", ~mean(.x, na.rm = TRUE))) %>% rowwise() %>% mutate(m = mean(c(adj_NO3_1, adj_NO3_2, adj_NO3_3 )))
  bellairs[ "Mean Nitrates [umol/L]" ] <- as.numeric( inter %>% filter( site =="B") %>% select("m") )
  maycocks[ "Mean Nitrates [umol/L]" ] <- as.numeric( inter %>% filter( site =="M") %>% select("m") )

  sea_chem <- sea_chem %>% mutate( adj_NO2 = (((`Nitrites(NO2-N)mg/l` / 1000 ) / 46.0055) / 0.000001 ) )
  inter <- sea_chem %>% group_by(site) %>% summarise( m = mean(adj_NO2) )
  bellairs[ "NO2 [umol/L]" ] <- as.numeric( inter %>% filter( site =="B") %>% select("m"))
  maycocks[ "NO2 [umol/L]" ] <- as.numeric( inter %>% filter( site =="M") %>% select("m"))
    
  
  sea_chem <- sea_chem %>% mutate( adj_PO4 = (((`PO4_Mean` / 1000 ) / 94.971) / 0.000001 ) )
  inter <- sea_chem %>% group_by(site) %>% summarise( m = mean(adj_PO4) )
  bellairs[ "PO4 [umol/L]" ] <- as.numeric( inter %>% filter( site =="B") %>% select("m"))
  maycocks[ "PO4 [umol/L]" ] <- as.numeric( inter %>% filter( site =="M") %>% select("m"))
  

  mt <- rbind(maycocks,mt)
  mt <- rbind(bellairs,mt)
  rownames(mt)[1:2] <- c("br_bel", "br_may")
  return(mt)
}

