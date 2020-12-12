###############################################################
#
#
# Simpson et al. (2020) Methods 8 Taxonomic Analysis
#
#       Input: Bracken output {bellairs.bracken, maycocks.bracken} assigning reads to taxa.
#
#       Output: An R data frame called tree that combines information from NCBI's Taxonomy,
#               the counts of reads assigned to taxa from Bracken, and
#               various statistics that are needed for later analyses.
#               RData file stored in /reef/R/tree.latest.RData
#
#
################################################################

library(dplyr)
root <- rprojroot::find_root(".git/index"); setwd(root); 
source(file.path(root, "src/functions.R"))
SEQDATA <- "/home/data/raw/ReefData/WaterSamples"
KDB <- RES <-  "/home/data/refined/reef/kr_br"


# Load the Bracken data.
bellairs <- read.csv(paste0(RES, "temp_trimgl_krbr/bellairs_trimgl.bracken"), header=TRUE, sep = "\t",
                     col.names = c("name", "tax_id", "tax_lvl", "assigned_reads", "added_reads", "est_reads", "fraction"))
bellairs <- bellairs[ order( -bellairs$fraction, -bellairs$est_reads), ]
bellairs$site <- "b"

maycock  <- read.csv(paste0(RES, "temp_trimgl_krbr/maycocks_trimgl.bracken"), header=TRUE, sep = "\t",
                    col.names = c("name", "tax_id", "tax_lvl", "assigned_reads", "added_reads", "est_reads", "fraction"))
maycock <- maycock[ order( -maycock$fraction, -maycock$est_reads), ]
maycock$site <- "m"

all <- union( bellairs$tax_id, maycock$tax_id )
options(warn=-1)
m <- bind_rows(lapply( 1:length(all), FUN = function( t ) {
  b <- which( bellairs$tax_id == all[t]); m <- which( maycock$tax_id == all[t])
  
  if (length(b)!=0) { tmp <- bellairs[b, 1:3] } else
  { tmp <- maycock[m, 1:3] }
  
  if (length(b)==0) { tmp[4:7] <- 0 } else { tmp[4:7] <- bellairs[b, 4:7] }
  if (length(m)==0) { tmp[8:11] <- 0 } else { tmp[8:11] <- maycock[m, 4:7] }
  
  colnames(tmp) <-  c("name", "tax_id", "lvl", "b_assigned", "b_added", "b_est", "b_fraction", 
                      "m_assigned", "m_added", "m_est", "m_fraction")
  return(tmp)
} ))

m$diff <- m$b_fraction - m$m_fraction

######################################
# NCBI Taxonomy download
#   We assume that the taxonomy has been download during the Kraken processing
#   e.g. via  kraken2-build --download-taxonomy -use-ftp  --db $RES
#
######################################


nodes <- read.csv(paste0(RES, "/taxonomy/nodes.dmp"), header=FALSE, sep = "\t")
nodes <- nodes[, c(1,3,5,7,9)]
colnames(nodes) <- c( "tax_id", "parent", "rank", "embl_code", "division_id" )

nms <- read.csv(paste0(RES, "/taxonomy/names.dmp"), header=FALSE, sep = "\t")
nms <- nms[, c(1,3,5,7)]
colnames(nms) <- c( "tax_id", "name_txt", "unique name", "name_class" )

divisions <- read.csv(paste0(RES, "/taxonomy/division.dmp"), header=FALSE, sep = "\t")
divisions <- divisions[, c(1,3,5,7)]
colnames(divisions) <- c( "id", "code", "name", "comments" )

######################################


tree <- nodes;  tree$br_bel <- 0;       tree$br_may <- 0;       tree$bell_orig_est_reads <- 0; tree$bell_orig_fraction <- 0
                tree$br_bel_frac <-  0; tree$br_may_frac <-  0; tree$may_orig_est_reads <- 0; tree$may_orig_fraction <- 0

######################################
# Assigning Braken counts to leaves at Bellairs
######################################

b_notfound <-c()
for (i in 1:nrow(bellairs)) {
  t <- bellairs[i, "tax_id"]
  tt <- which(tree$tax_id == t)
  
  if (length(tt)==0) {
    cat("\n Not found Bellairs, ", i, t)
    b_notfound <- append(b_notfound, i)
  } else {
    tree[ tt, "br_bel"] <- tree[ tt, "bell_orig_est_reads"] <- bellairs[i, "est_reads"]
    tree[ tt, "br_bel_frac"] <- tree[ tt, "bell_orig_fraction"] <- bellairs[i, "fraction"]
  }
} # end of bellairs

######################################
# Assigning Braken counts to leaves at Maycocks
######################################

m_notfound <- c()
for (i in 1:nrow(maycock)) {
  t <- maycock[i, "tax_id"]
  tt <- which(tree$tax_id == t)
  
  if (length(tt)==0) {
    cat("\n Not found maycock, ", i, t)
    m_notfound <- append(m_notfound, i)
  } else {
    tree[ tt, "br_may"] <- tree[ tt, "may_orig_est_reads"] <- maycock[i, "est_reads"]
    tree[ tt, "br_may_frac"] <- tree[ tt, "may_orig_fraction"] <- maycock[i, "fraction"]
  }
}

######################################
# Set parent of root to NA to avoid
# infinite recursion.
######################################

tree[1, "parent"] <- NA

# percolate the reads at the leaves up through the tree to the root.
void <- percolate(1)

save(tree, file = "~/temp.percolate.tree.Rdata")
# prune away all nodes (taxa) without any reads assigned to them.

to_remove <- intersect( which(tree$br_bel==0), which(tree$br_may==0) )
tree <- tree[ -to_remove, ]

#########################
#
# add names directly to each row of the tree. 
# use the scientific name if it exists.
#
##########################

for (i in 1:nrow(tree)) {
  tmp <- which(nms$tax_id == tree[i, "tax_id"])
  if (length(tmp) == 0) next
  if (length(tmp) == 1) { tree[i, "name"] <- as.character(nms[tmp, "name_txt"]); next }
  tmp2 <- which( nms[tmp, "name_class"] == "scientific name")
  if (length(tmp2) == 0)  { tree[i, "name"] <- as.character(nms[tmp[1], "name_txt"]); next }
  if (length(tmp2) == 1) { tree[i, "name"] <- as.character(nms[tmp[tmp2], "name_txt"]); next }
  if (length(tmp2) > 1) { tree[i, "name"] <- as.character(nms[tmp[tmp2[1]], "name_txt"]); next }
  print(tree$name[i, "name"])
}
tree <- tree[ , c( 14, 1:13 )]


##################################################
# Now add frequencies.  src/function.R
##################################################

void <- local_frequencies(1)
void <- global_frequencies( 1 )

##################################################
# Now statistical tests. src/function.R
##################################################

void <- multinomial_tree_test(1)

void <- polarity_test( 1 )

##################################################
# Include the taxon to root path through the tree of life.  
#       src/functinon.R
##################################################

for (i in 1:nrow(tree)) {
  tree[i, "path"] <- unlist( paste( path2root(tree[i, "tax_id"])$name, collapse="."))
}

tree$DeltaFreq <- tree$Glob.Freq.Bel - tree$Glob.Freq.May

tree <- tree[, c(1:18, 23, 19:22)]

#save(tree, file = "/home/data/refined/reef/R/raw.tree.april.15.RData")
#write.csv(tree, file = "/home/data/refined/reef/R/raw.tree.april.15.csv")



