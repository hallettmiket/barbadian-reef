
c2f <- count2frequency <- function( counts ) return( counts / sum(counts) )
f2c <- frequency2count <- function( freq, n ) return( as.vector(rmultinom( 1 , n, prob=freq )))

#######################################################
# augments is a vector the same or less than the length of counts; 
#   each non-zero represents the fraction of total (sum(counts)) that are
#   added to that component under a multinomial distribution (with adj freq vector)
#   eg counts == c( 1,1,1,1 ); augments <- c( 0.5, 0.5 ) so the new
#   frequency is c( 1 + 2, 1+ 2, 1, 1) /8 
#######################################################

prob_augment_counts <- function( counts, augments ) {
  new_augments <- vector( mode="numeric", length = length(counts))
  new_augments[1:length(augments)] <-  augments 
  new_counts <- counts + (sum(counts) * new_augments)
  return( f2c(c2f(  new_counts ), sum(new_counts) ) )}

#######################################################
# similar to prob_augments_counts but it does not use the mulitnomial but simply adds 
# the fraction to counts deterministically
#######################################################

augment_counts <- function( counts, augments ) {
  new_augments <- vector( mode="numeric", length = length(counts))
  new_augments[1:length(augments)] <-  augments 
  return( counts + (sum(counts) * new_augments) )}

#######################################################
# Returns the estimated p-value for component of the observed vector
# Each such X_i ~ Bin(n, p_i) where n is sum(observed)
#######################################################
test_each_component <- function( observed, freq, threshold = 0, verbose = FALSE) {
  
  big <- which( (observed/sum(observed)) > threshold)
  if (verbose) print(big)
  ans <- rep(1, times=length(observed))
  lapply( big, FUN = function(i) {
                   ans[i] <<- binom.test(x=observed[i], n=sum(observed), p=freq[i], alternative="two.sided" )$p.value })
  return(ans) 
}
                 

#######################################################
# For two observed count vectors (equal lenth), compute the
# freq vector assuming null hypothesis (from same distribution)
#######################################################
null_freq <- function( obs1, obs2 ) {
  if (length(obs1) != length(obs2)) stop("Count vectors must be of equal length")
  return(( obs1 + obs2 ) / (sum(obs1) + sum(obs2))) }


#######################################################
# decentralize: for two observed vectors, as long as the null hypothesis is rejected
#   identify components that are significantly different via the binomial test
#   and remove them.
#######################################################

# version 1  - it sucks because many low frequency guys have estimated pvalues of 0. so you end up
# removing a bunch of low frequency guys.

decentralize_v1 <- function( obs1, obs2, p = 0.05, verbose = FALSE ) {
  if (length(obs1) != length(obs2)) stop("Count vectors must be of equal length")
  
  removed_idx <- c(); futile <- FALSE;
  while ( (length(obs1) > 1) && (chisq.test(rbind(obs1, obs2))$p.value < p ) && (!futile)) {
    if (verbose) print(removed_idx)
    freq <- null_freq(obs1, obs2)
    out1 <- test_each_component(obs1, freq); out2 <- test_each_component(obs2, freq)
    min1 <- min(out1); idx1 <- which.min(out1); min2 <- min(out2); idx2 <- which.min(out2)
    if ((min1 > p) && (min2 > p)) { futile <- TRUE; next }
    
    ifelse( min1 < min2,  target <- idx1, target <- idx2 )
    removed_idx <- c(removed_idx, target)
    obs1 <- obs1[-target]; obs2 <- obs2[-target]
  } # end of while
  return(removed_idx)
}


# version 2  - this version only considers components that are significant at p-value p
# but then choose the entry that has the maximum frequency (obs counts to total counts at site)

decentralize_v2 <- function( nodes, p = 0.05, threshold = 0.00, verbose = FALSE ) {
  obs1 <- nodes$br_bel
  obs2 <- nodes$br_may

  removed_nodes <- NA
  removed_idx <- c(); futile <- FALSE;
  while ( (length(obs1) > 1) && (chisq.test(rbind(obs1, obs2))$p.value < p ) && (!futile)) {

    freq <- null_freq(obs1, obs2)
    out1 <- test_each_component(obs1, freq, threshold, verbose = TRUE)
    out2 <- test_each_component(obs2, freq, threshold, verbose = TRUE)
    sig1 <- which(out1 < p); sig2 <- which(out2 < p)
    
    if ((length(sig1) < 1) && (length(sig2) < 1)) { futile <- TRUE; next }
    
    max1 <- max(obs1[sig1]); max2 <- max(obs2[sig2])
    idx1 <- which.max(obs1[sig1]); idx2 <- which.max(obs2[sig2])
    
    ifelse( (max1/sum(obs1)) > (max2/sum(obs2)),  target <- sig1[idx1], target <- sig2[idx2] )
    removed_nodes <- rbind(removed_nodes, nodes[target,])
    nodes <- nodes[-target,]
    
    obs1 <- nodes$br_bel
    obs2 <- nodes$br_may
    
    if (verbose) { 
      cat("\n\n-----\nRemoved: "); print( removed_nodes[nrow(removed_nodes), c(1:3, 7:8)] ) 
      cat("\nLength of observed vector"); print(length(obs1))
      cat("\nTotal Count (Site 1):"); print(sum(obs1))
      cat("\nTotal Count (Site 2):"); print(sum(obs2))
      cat("\nGlobal multinomial:"); print(chisq.test(rbind(obs1, obs2))$p.value )
    }
       
  } # end of while
  print("Exit because:")
  print(length(obs1))
  print(futile)
  print(chisq.test(rbind(obs1, obs2))$p.value)
  return(removed_nodes)
}


# version 3  - KL divergence

decentralize_v3 <- function( nodes, p = 0.05, threshold = 0.00, verbose = FALSE ) {
  obs1 <- nodes$br_bel
  obs2 <- nodes$br_may
  all_obs <- list (list(obs1, obs2))
  
  removed_nodes <- NA
  removed_idx <- c(); futile <- FALSE;
  while ( (length(obs1) > 1) && (chisq.test(rbind(obs1, obs2))$p.value < p ) && (!futile)) {
    
    freq <- null_freq(obs1, obs2)
    out1 <- test_each_component(obs1, freq, threshold, verbose = TRUE)
    out2 <- test_each_component(obs2, freq, threshold, verbose = TRUE)
    sig1 <- which(out1 < p); sig2 <- which(out2 < p)
    
    if ((length(sig1) < 1) && (length(sig2) < 1)) { futile <- TRUE; next }
    
    max1 <- max(obs1[sig1]); max2 <- max(obs2[sig2])
    idx1 <- which.max(obs1[sig1]); idx2 <- which.max(obs2[sig2])
    
    ifelse( (max1/sum(obs1)) > (max2/sum(obs2)),  target <- sig1[idx1], target <- sig2[idx2] )
    removed_nodes <- rbind(removed_nodes, nodes[target,])
    nodes <- nodes[-target,]
    
    obs1 <- nodes$br_bel
    obs2 <- nodes$br_may
    
    if (verbose) { 
      cat("\n\n-----\nRemoved: "); print( removed_nodes[nrow(removed_nodes), c(1:3, 7:8)] ) 
      cat("\nLength of observed vector"); print(length(obs1))
      cat("\nTotal Count (Site 1):"); print(sum(obs1))
      cat("\nTotal Count (Site 2):"); print(sum(obs2))
      cat("\nGlobal multinomial:"); print(chisq.test(rbind(obs1, obs2))$p.value )
    }
    all_obs <- append(all_obs, list(list(obs1, obs2)))
    
  } # end of while
  print("Exit because:")
  print(length(obs1))
  print(futile)
  print(chisq.test(rbind(obs1, obs2))$p.value)
  plot_progress(all_obs)
  return(removed_nodes)
}


plot_progress <- function( obs, lmt = 300 ) {
  
  if (lmt > length(obs[[1]][[1]])) stop("Limit too high.")
  
  tmp1 <- lapply( obs, function(x) return(c2f(x[[1]]) ))
  tmp2 <- lapply( obs, function(x) return(c2f(x[[2]]) ))
  tmp3 <- list()
  
  for (i in 1:length(tmp1)) {
    o <- order(-abs(tmp1[[i]]))
    tmp1[[i]] <- tmp1[[i]][o][1:lmt]
    tmp2[[i]] <- tmp2[[i]][o][1:lmt]
    tmp3[[i]] <- tmp1[[i]] - tmp2[[i]]
  }
  

  diff <- do.call(rbind, tmp3 )


  df <- data.frame(  ); 
  for (i in 1:ncol(diff)) {
    for (j in 1:nrow(diff)) {
      df <- rbind( df, c( i, diff[j,i], j ))
    } }
  colnames(df) <- c("idx", "frac", "exper")
  
  ggplot(df, aes(idx, frac) ) +
    geom_line() + 
    facet_wrap(~exper) 

  
} #end of plot_progress



##############################################################################################################

#######################################################
# Routines for exploring the Tree of Life
#######################################################



#######################################################
# frequency distribution
#######################################################


freq_distribution <- function( parent, site ) {
  c <- p2c( parent )
  c_i <- t2i(c)

  if (length(c) <= 1) { return(NULL) }
  
  x <- (0:(length(c)-1)) / (length(c) - 1)
  tmp <- tree[ c_i, site ]
  y <- tmp[ order( -tmp ) ]
  
  return(list(x,y))

} # end of make_table


#######################################################
# get all frequency distributions
#######################################################


get_all_freq_distributions <- function( parent, site, level = c(0) ) {
  
  c <- p2c( parent )
  c_i <- t2i(c)
  ind <- t2i( parent )
  
  if (length(level) < 100)
    cat("\nLevel ", level); 
  
  if (length(c) > 2)
    master<<- append(master, list(freq_distribution( parent, site )))
  
  if (length(c) > 0) {
    for (i in 1:length(c)) {
      get_all_freq_distributions( c[i], site, append(level, c(i, length(c))) )
          } # i
  }
  
  return( NULL )
  
}

#######################################################
# Add a column to the tree d.f. isChild iff it is a child
#######################################################


recursive_isLeaf <- function( parent, level = c(0) ) {
  
  c <- p2c( parent )
  if (length(c) == 0) { tree[t2i(parent), "isLeaf"] <<- TRUE; return(NULL) }
  tree[t2i(parent), "isLeaf"] <<- FALSE
  
  if (length(level) < 100) cat("\nLevel ", level); 
  
  c_i <- t2i(c)
  if (length(c) > 0) {
    for (i in 1:length(c)) {
      recursive_isLeaf( c[i],  append(level, c(i, length(c))) )
    } # i
  }
  return( NULL )
}



