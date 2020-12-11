

######################################################
# Genome size versus read count comparisons
#
# tax_id2genome_size (t2g)
#
# Input:  tax_id: An NCBI tax_id
#         ncbi: a list of data.frames from the NCBI GENOME_REPORTS
# 
#
# Output: The genome size of the tax_id
#
#######################################################

t2g <- tax_id2genome_size <- function( tid, ncbi ){
  for (i in 1:length(ncbi)) {
    tmp <- which(ncbi[[i]]$TaxID == tid) 
    if (length(tmp) > 0) {
      return(mean(ncbi[[i]][tmp, "Size..Mb."], na.rm= TRUE ))
    }
  }
  return(NA)
} # end of t2g



######################################################
# Genome size versus read count comparisons
#
# percolate_genome_size
#
# Input:  
#         parent is the root of the tree to recurse down
#         ncbi: a list of data.frames from the NCBI GENOME_REPORTS
#         level is for diagnostics
# 
#
# Output: The genome size of the tax_id
#
#######################################################


pgs <- percolate_genome_size <- function( parent, ncbi, level= c(0) ) {
  c <- parentTax_2_childrenTax( parent )
  ind <- tax_2_index( parent )
  
 # if (length(level) < 100)
#    cat("\nLevel ", level)
  
  if (length(c) > 0) {   # not a leaf
    for (i in 1:length(c)) {
      void <- percolate_genome_size( c[i], ncbi, append(level, c(i, length(c)) ))
    }
    tmp <- mean( tree[t2i(c), "genome_size"], na.rm=TRUE)

    tree[ind, "genome_size"] <<- ifelse(is.nan(tmp), NA,  tmp)
    
  } else {     # leaf
    tree[ind, "genome_size"] <<- t2g( parent, ncbi )
  }
  
  return( NA )
} # end of percolate


