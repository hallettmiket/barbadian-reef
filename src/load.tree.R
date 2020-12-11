

load.tree <- function(path = REEF_DIR, vers="1.0", struct_type = ".RData", adjust = FALSE, tidyup = TRUE  ){
  
  load(file = file.path(path, paste0(paste0("tree_", vers), struct_type)))
  if (exists("tree"))  return(tree) else stop("Code assumes a data.frame of name tree.")
  return(tidayup(tree))
} # load.tree

load.adjustments <- function( path = "~", vers="1.0" ){
  load(file = file.path(path, paste0(paste0("genome_size_adjustment_", vers), ".RData")))
  if (exists("genome_size_adjustments"))  return(tree) else stop("Code assumes a list of linear models.")
  return(genome_size_adjustments)
}



tidyup <- function(  ) {
  
  # The following have already been incorporated into version 1.0
  
  # tree[t2i(1993653), "name"] <<- "1993653 Cand. Nitrosomarinus"
  # tree[t2i(1680826), "name"] <<- "1680826 Cand. Poseidoniales"
  # tree[t2i(570266), "name"] <<- "570266 Methanocella"
  # tree[t2i(1593364), "name"] <<- "1593364 Nitrosopelagicus"
  # tree[t2i(1980514), "name"] <<- "1980514 Halodesulfurarchaeum"
  # tree[t2i(1783275), "name"] <<- "1783275"
  # tree[t2i(1783276), "name"] <<- "1783276"
  # tree[t2i(1803821), "name"] <<- "1803821"
  # tree[t2i(1803822), "name"] <<- "1803822"
  # tree[t2i(2026739), "name"] <<- "2026739"
  # tree[t2i(2608109), "name"] <<- "2608240 Haptophyta"
  # tree[t2i(2611352), "name"] <<- "2611352 "
  # tree[t2i(554915), "name"] <<- "554915 Eumycetozoa"
  # tree[t2i(2611341), "name"] <<- "2611341 Parabasalia/Fornicata"
  # tree[t2i(554296), "name"] <<- "2608240 Apusomonadidae"
  # tree[t2i(2683617), "name"] <<- "2683617 Picozoa/Palpitomonas"
  # tree[t2i(2608240), "name"] <<- "2608240 Collodictyonidae"
  # tree[t2i(589304), "name"] <<- "589304 Dinophyceae envir samples"
  # tree[t2i(669373), "name"] <<- "669373 Dinophyceae"
  # tree[t2i(549779), "name"] <<- "549779 Mimiviridae"
  # tree[t2i(2559587), "name"] <<- "2559587 "
  # tree[t2i(548681), "name"] <<- "548681 Herpesviridae"
  # tree[t2i(687329), "name"] <<- "687329 Anelloviridae"
  # tree[t2i(1982239), "name"] <<- "1982239 Nipunavirus"
  # tree[t2i(2560199), "name"] <<- "2560199"
  # tree[t2i(1783272), "name"] <<- "1783272"
  # tree[t2i(1783270), "name"] <<- "1783270"
  # tree[t2i(1783257), "name"] <<- "1783257"
  # tree[t2i(1553900), "name"] <<- "1553900"
  # tree[t2i(580370), "name"] <<- "580370"
  # tree[t2i(2008785), "name"] <<- "2008785"
  # tree[t2i(1798711), "name"] <<- "1798711 Cyanobacteria"
  
  # end of version 1.0
  
} # end of tidyup