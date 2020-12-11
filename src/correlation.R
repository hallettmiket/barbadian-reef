Spearman_test <- function(parent, level= c(0)) {
  
  c <- p2c( parent )
  c_i <- t2i(c)
  ind <- t2i( parent )
  
  
  if (length(level) < 100)
    cat("\nLevel ", level)
  
  tree[ ind, "Spearman.Br"] <-  NA;   tree[ ind, "Spearman.Kr"] <-  NA; tree[ ind, "Spearman_CC.Br"] <-  NA;   tree[ ind, "Spearman_CC.Kr"] <-  NA; 
  
  if (length(c) > 0) {
    M <- as.table(rbind(tree[c_i, "br_b"], tree[c_i, "br_c"]))
    dimnames(M) <- list(site = c("B", "C"), subcat = tree[c_i, "name"])
    if (length(c)==1) { tree[ ind, "Spearman.Br"] <<- 1 } else {
      Sprman <- cor.test(M["B"],M["C"], method = "spearman")
      tree[ ind, "Spearman_CC.Br"] <<-  Sprman$estimate
      tree[ ind, "Spearman.Br"] <<-  Sprman$p.value }
      
      
    
    M <- as.table(rbind(tree[c_i, "kr_b_tot"], tree[c_i, "kr_c_tot"]))
    dimnames(M) <- list(site = c("B", "C"), subcat = tree[c_i, "name"])
    if (length(c)==1) {tree[ ind, "Spearman.Kr"] <<- 1 } else {
      SprmanKr <- cor.test(M["B"],M["C"], method = "spearman")
      tree[ ind, "Spearman_CC.Kr"] <<-  SprmanKr$estimate
      tree[ ind, "Spearman.Kr"] <<-  SprmanKr$p.value }
      
     
    
    for (i in 1:length(c)) {
      Spearman_test( c[i], append(level, c(i, length(c)) ))
    } # i 
  }
} # end of Spearman_test


kendall_test <- function(parent, level= c(0)) {
  
  c <- p2c( parent )
  c_i <- t2i(c)
  ind <- t2i( parent )
  
  
  if (length(level) < 100)
    cat("\nLevel ", level)
  
  tree[ ind, "Kendall.Br"] <-  NA;   tree[ ind, "Kendall.Kr"] <-  NA; tree[ ind, "Kendall_CC.Br"] <-  NA;   tree[ ind, "Kendall_CC.Kr"] <-  NA; 
  
  if (length(c) > 0) {
    M <- as.table(rbind(tree[c_i, "br_b"], tree[c_i, "br_c"]))
    dimnames(M) <- list(site = c("B", "C"), subcat = tree[c_i, "name"])
    if (length(c)==1) { tree[ ind, "Spearman.Br"] <<- 1 } else {
      kendall <- cor.test(tree[c_i, "br_b"],tree[c_i,"br_b"], method = "kendall")
      tree[ ind, "Kendall_CC.Br"] <<-  kendall$estimate
      tree[ ind, "Kendall.Br"] <<-  kendall$p.value }
    
    
    
    M <- as.table(rbind(tree[c_i, "kr_b_tot"], tree[c_i, "kr_c_tot"]))
    dimnames(M) <- list(site = c("B", "C"), subcat = tree[c_i, "name"])
    if (length(c)==1) {tree[ ind, "Kendall.Kr"] <<- 1 } else {
      kendallKr <- cor.test(tree[c_i, "kr_b_tot"],tree[c_i,"kr_c_tot"], method = "kendall")
      tree[ ind, "Kendall_CC.Kr"] <<-  kendallKr$estimate
      tree[ ind, "Kendall.Kr"] <<-  kendallKr$p.value }
    
    
    
    for (i in 1:length(c)) {
      kendall_test( c[i], append(level, c(i, length(c)) ))
    } # i 
  }
} # end of Kendall Test
