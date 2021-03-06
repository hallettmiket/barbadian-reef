## This contains functions and only functions for estimating differences in abundances
## between sites.
## The file is meant to be sourced.
# MH

library(dplyr)
library(reshape)


#######################################################
# functions to move around in the Tree of Life
# 
# index of parent: this function just adds a $parent field to the data frame so that each node knows its parents
#
#######################################################

indexOfParent <- function( tax ) {
  return( which(tree$tax_id == tree[which(tree$tax_id == tax ), "parent"] ) )
}

# for a parent (tax_id not index) it returns its children (tax_ids not indicies)
parentTax_2_childrenTax <- function( parent ) {
  return(  tree[ which(tree$parent == parent), "tax_id" ] )
}
p2c <- parentTax_2_childrenTax

# converts a tax_id into an index (row) of the data.frame
tax_2_index <- function( t ) {
  tmp <- c()
  for (i in 1:length(t)) {
 #   if (i %% 1000 == 0) print(i)
    tmp[i] <- which(tree$tax_id == t[i])
  }
  return(tmp) }
t2i <- tax_2_index

#######################################################

# After the bracken reads are assigned to taxa in the tree of life
#   this function "percolates" these numbers upward towards the root.
# For example, if a node p has children c1, c2, c3 with bracken counts b1, b2, b3,
#   the bracken count at p will be b1 + b2 + b3
#

percolate <- function( parent, level= c(0) ) {
  c <- parentTax_2_childrenTax( parent )
  ind <- tax_2_index( parent )
  
#  if (length(level) < 100)
#    cat("\nLevel ", level)
  
  if (length(c) > 0) {
    
    for (i in 1:length(c)) {
      tmp <- percolate( c[i], append(level, c(i, length(c)) ))
      tree[ind, "br_bel"] <<- tree[ind, "br_bel"] + tmp[1]
      tree[ind, "br_may"] <<- tree[ind, "br_may"] + tmp[2]
    }
  }
  
  return( tree[ind, c("br_bel", "br_may")] )
} # end of percolate



#######################################################
# Give back the subtree of the Tree of Life (tree) rooted at parent
#   we assume here the full tree of life is located at variable tree
# This function creates a new data.frame that is a subset of tree
#   containing only rows that are descendents of the parameter parent
#
#######################################################

induce_tree <- function( parent ) {
  
  c <- p2c( parent )
  ind <- t2i( parent )
  
  tmp <- tree[ind, ]
  if (length(c) > 0) {
    for (i in 1:length(c)) {
      tmp <- rbind( tmp, induce_tree(c[i]))
    }
  }
  return(tmp)
}

#######################################################
# Frequencies
#
#   local frequencies. for parent p with children c1, c2, c3 with bracken counts b1, b2, b3
#   the frequency for c1 is b1 / (b1+b2+b3)   (computed seprately for both the bellairs and
#   maycocks counts).
#
#######################################################

local_frequencies <- function( parent, level= c(0) ) {
  
  c <- p2c( parent );  c_i <- t2i(c); ind <- t2i( parent )
  
#  if (length(level) < 100)  cat("\nLevel ", level)
  
  if (length(c) > 0) {
    for (i in 1:length(c)) {
      tree[ c_i[i], "Local.Freq.Bel"] <<- tree[ c_i[i], "br_bel"] / sum( tree[ c_i, "br_bel"] )
      tree[ c_i[i], "Local.Freq.May"] <<- tree[ c_i[i], "br_may"] / sum( tree[ c_i, "br_may"] )
    } # i
    
    for (i in 1:length(c)) {
      local_frequencies( c[i], append(level, c(i, length(c)) ))
    } # i 
  }
  return( NULL )
} # end of local_frequencies


#   local frequencies. for parent p with children c1, c2, c3 with bracken counts b1, b2, b3
#   the frequency for c1 is b1 / (b1+b2+b3)


global_frequencies <- function( parent, B_tot = tree[1, "br_bel"], M_tot = tree[1, "br_may"] ) {
  
  c <- p2c( parent )
  ind <- t2i( parent )
  tree[ind, "Glob.Freq.Bel"] <<- tree[ind, "br_bel"] / B_tot
  tree[ind, "Glob.Freq.May"] <<- tree[ind, "br_may"] / M_tot
  
  if (length(c) > 0) {
    for (i in 1:length(c)) {
      global_frequencies(c[i], B_tot, M_tot)
    }
  }
  return(NULL)
}





########################################################
# Statistics
#
# multinomial_tree_test recursively goes through the phylogenetic tree and 
#   assigns the estimated p-value in tree$Multinom
#
# polarity_test  recursively goes through the phylogenetic tree and 
#   assigns the estimated p-value in Polarity (original) and Polarity.Adj (adjusted)
#
#####

multinomial_tree_test <- function(parent, level= c(0), silent = TRUE) {
  
  c <- p2c( parent )
  c_i <- t2i(c)
  ind <- t2i( parent )
  
  
#  if ((length(level) < 100) && !silent)
 #   cat("\nLevel ", level)
  
  tree[ ind, "Multinom"] <-  NA;   
  
  if (length(c) > 0) {
    M <- as.table(rbind(tree[c_i, "br_bel"], tree[c_i, "br_may"]))
    dimnames(M) <- list(site = c("B", "M"), subcat = tree[c_i, "name"])
    if (length(c)==1) { tree[ ind, "Multinom"] <<- 1 } else {
  #    print(parent)
   #   print(M)
      
      tree[ ind, "Multinom"] <<-  chisq.test(M)$p.value }
    
    for (i in 1:length(c)) {
      multinomial_tree_test( c[i], append(level, c(i, length(c)) ))
    } # i 
  }
} # end of multinomial_tree_test

polarity_test <- function(parent, level= c(0), silent = TRUE) {
  
  c <- p2c( parent )
  c_i <- t2i(c)
  ind <- t2i( parent )
#  if ((length(level) < 100) & !silent) cat("\nLevel ", level)


  if (length(c) > 0) {

    polarity.orig <- NA
    polarity.new <- NA
    if (length(c) > 2) {
    
      M <- as.table(rbind(tree[c_i, "br_bel"], tree[c_i, "br_may"]))
      if (sum(M[1,]) == 0) M[1,] <- 0 else M[1,] <- M[1,] / sum(M[1,])
      if (sum(M[2,]) == 0) M[2,] <- 0 else M[2,] <- M[2,] / sum(M[2,])
      X <- sum(M[1,] < M[2,])
      polarity.orig <- binom.test(x=X, n=ncol(M), p=0.5, alternative="two.sided" )$p.value 
      M <- M[ ,-which(M[1,]<M[2,])]

      
      if ((class(M) != "numeric") && (dim(M)[2]!=0)){
        # renormalize 
        M[1,] <- M[1,] / sum(M[1,]); M[2,] <- M[2,] / sum(M[2,])
    
        # and re-test
        if ((all(!is.na(M[1,]))) && (all(!is.na(M[2,])))) {
           X <- sum(M[1,] < M[2,])
          polarity.new <- binom.test(x=X, n=ncol(M), p=0.5, alternative="two.sided" )$p.value 
        }  else { polarity.new <- -1 }
      }

    } # end of internal lenght check c > 2
    tree[ind, "Polarity"]<<- polarity.orig 
    tree[ind, "Polarity.Adj"]<<- polarity.new 

    for (i in 1:length(c)) {
      polarity_test( c[i], append(level, c(i, length(c)) ), silent)
    } # i 
  }
} # end of polarity_test


#######################################################
# Find the path from a taxa to the root
#######################################################

path2root <- function( target ) {
  if (is.na(tree[t2i(target), "parent"])) { 
    return(tree[t2i(target), ])
  }
  else {
    return( rbind( path2root(tree[ t2i(target), "parent"]), tree[t2i(target),]))
  }
}

#######################################################
# Pretty print
#######################################################


mt <- make_table <- function( parent, precision = 3, relative_taxa = NA ) {
  c <- p2c( parent )
  c_i <- t2i(c)
  ind <- t2i( parent )
  partial <- tree[ c_i, 1:4]
  if (length(c) > 0) {
    for (i in 1:length(c)) {
      partial[ i, "Local.Freq.Bel"] <- sprintf(paste0(paste0("%.", precision), "f"), tree[ c_i[i], "Local.Freq.Bel"])
      partial[ i, "Local.Freq.May"] <- sprintf(paste0(paste0("%.", precision), "f"), tree[ c_i[i], "Local.Freq.May"])
      partial[ i, "log(BvsM)"] <- sprintf(paste0(paste0("%.", precision), "f"), 
                                          log( tree[ c_i[i], "Local.Freq.Bel"] / tree[ c_i[i], "Local.Freq.May"]))
      
      partial[ i, "Glob.Freq.Bel"] <- sprintf(paste0(paste0("%.", precision), "f"), tree[ c_i[i], "Glob.Freq.Bel"])
      partial[ i, "Glob.Freq.May"] <- sprintf(paste0(paste0("%.", precision), "f"), tree[ c_i[i], "Glob.Freq.May"])     
      partial[ i, "DeltaFreq"] <- sprintf(paste0(paste0("%.", precision), "f"), tree[ c_i[i], "DeltaFreq"])     
      
      if (!is.na(relative_taxa)) {
        partial[ i, "Rel.Freq.Bel"] <- sprintf(paste0(paste0("%.", precision), "f"), 
                                               tree[ c_i[i], "br_bel"] / tree[ t2i(relative_taxa), "br_bel"]
                                               )
        partial[ i, "Rel.Freq.May"] <-  sprintf(paste0(paste0("%.", precision), "f"), 
                                                     tree[ c_i[i], "br_may"] / tree[ t2i(relative_taxa), "br_may"]
                                              )
        partial[ i, "Rel.DeltaFreq"] <- sprintf(paste0(paste0("%.", precision), "f"), 
                                                as.numeric(partial[ i, "Rel.Freq.Bel"])  - as.numeric( partial[ i, "Rel.Freq.May"] ) 
                                                )     
        partial[ i, "Two-Portions"] <- sprintf(paste0(paste0("%.", precision), "f"), 
                                               prop.test(
                                                c(tree[ c_i[i], "br_bel"] , tree[ c_i[i], "br_may"]),
                                                c(tree[ t2i(relative_taxa), "br_bel"],tree[ t2i(relative_taxa), "br_may"])
                                               )$p.value
        )
      }
      
      partial[ i, "tmpFreq.Bel"] <-  tree[ c_i[i], "Local.Freq.Bel"]
      partial[ i, "tmpFreq.May"] <-  tree[ c_i[i], "Local.Freq.May"]
      
      
    }
    partial <- partial[ order( -partial$tmpFreq.Bel, partial$tmpFreq.May), 
                        -c((ncol(partial)-1): ncol(partial)) ]
  }
  rownames(partial) <- NULL
  colnames(partial)[1:4] <- c("Name",  "Tax. Id.", "Parent", "Rank")
  return( partial )
} # end of make_table


#######################################################
# simple newick
#######################################################

simple_newick <- function( parent ) {
  c <- p2c( parent );  ind <- t2i( parent )
  tmp <- ""
  if (length(c) == 0) {
    tmp2 <- unlist( paste(path2root( parent )$name, collapse=".") )
    return(tmp2) 
  } else {
    for (i in 1:length(c)) {
      tmp <- c( tmp, simple_newick(c[i]))
    }
  }
  return(tmp[2:length(tmp)])
} # end of simple newick

#######################################################
# reduce_tree
#######################################################

reduce_tree <- function( parent ) {
  tmp <- induce_tree(parent)
  return( tree[-t2i( tmp$tax_id ), ] )
}

#######################################################
# remove_update_tree: removes the subtypee at parent and 
#    recalcualtes all read counts at internal nodes with frequencies
#######################################################

remove_update_tree <- function( parent ) {
  tree <<- reduce_tree( parent )

  tree$br_bel <<- tree$bell_orig_est_reads;  tree$br_may <<- tree$may_orig_est_reads
  tree$br_bel_frac <<- tree$bell_orig_fraction; tree$br_may_frac <<- tree$may_orig_fraction
  
  void <- percolate( 1 )

  void <- local_frequencies(1)
  void <- global_frequencies( 1 )
  
  void <- multinomial_tree_test(1)
  void <- polarity_test( 1 )
  
  tree$DeltaFreq <<- tree$Glob.Freq.Bel - tree$Glob.Freq.May
  return(NULL)
}


#######################################################
# adjust_read_counts: takes a linear model (fitted with species
#   using log(reads) ~ log(genome_size) and adjusts the 
#   (un-logged) read counts br_bel and br_may for all species
#   in my_tree
#
#   mytree is a tree data.frame; f is a lm object
#
#######################################################

adjust_read_counts <- function( mytree, f ) {
  intercept <- f$coefficients[1]; slope <- f$coefficients[2]
  species_ind <- ((!is.na(mytree$genome_size)) & (mytree$rank=="species"))
  mytree$br_bel <- mytree$br_bel+1
  mytree$br_may <- mytree$br_may+1
  
  reads_adj_bel <- ifelse( species_ind, 
    log(mytree$br_bel) - (log(mytree$genome_size)*f$coefficients[2]+f$coefficients[1]),
    log(mytree$br_bel) - f$coefficients[1])
  
    
  reads_adj_may <- ifelse( species_ind, 
                           log(mytree$br_may) - (log(mytree$genome_size)*f$coefficients[2]+f$coefficients[1]),
                           log(mytree$br_may) - f$coefficients[1])
  mytree$"reads_adj_bel" <- reads_adj_bel
  mytree$"reads_adj_may" <- reads_adj_may
  return(mytree)
}



#######################################################
# Add a column to the tree d.f. isChild iff it is a child
#######################################################


recursive_isLeaf <- function( parent, level = c(0) ) {
  
  c <- p2c( parent )
  if (length(c) == 0) { tree[t2i(parent), "isLeaf"] <<- TRUE; return(NULL) }
  tree[t2i(parent), "isLeaf"] <<- FALSE
  
 # if (length(level) < 100) cat("\nLevel ", level); 
  
  c_i <- t2i(c)
  if (length(c) > 0) {
    for (i in 1:length(c)) {
      recursive_isLeaf( c[i],  append(level, c(i, length(c))) )
    } # i
  }
  return( NULL )
}




