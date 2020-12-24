
freq_count_plot <- function( mytree, relative_to = 1, mytitle="Eukaryota", verbose = FALSE, 
                             lf_quant = 0.15, rf_quant = 0.89, lc_quant = 0.42, rc_quant = 0.99,
                             top_left = -20, top_even = -10, top_right = 15,
                             stretch = 1.0,
                             size = 3) {
  
  root_bel <-max(1,tree[t2i(relative_to), "br_bel"]) ; root_may <-max(1, tree[t2i(relative_to), "br_may"])

  mytree$total_counts <- log(mytree$br_bel+mytree$br_may)
  mytree$delta <- log( ((mytree$br_bel+1) / root_bel ) /  ((mytree$br_may+1) / root_may ) )

  results_freq <- boot( data = mytree$delta, statistic = mystat, R = 10000)
  ci <- boot.ci(results_freq)
  leftpt_freq <- ci$norm[2]; rightpt_freq <- ci$norm[3]
  
  tmp2 <- quantile( mytree$delta, c(lf_quant, rf_quant))
  tmp <- quantile( mytree$total_counts, c(lc_quant, rc_quant))
  leftpt_counts <- tmp[1]; rightpt_counts <- tmp[2]
  
  if (verbose) {
    cat("\n Counts at ancestor: ", root_bel, root_may)
    cat("\n Mean: ", mean(mytree$delta) )
    cat("\n Confidence Interval for mean: ", ci$norm)
    cat("\n Quantile on Frequencies:", tmp2 )
    cat("\n Quantile on counts:", tmp)
  }
  
  
  stretch_factor <-  stretch * (max(mytree$delta) -  min(mytree$delta))
  
  
  num_taxa <- length(unique(mytree$taxa))
  
  clrs <- glasbey(); 
  # make the colors a bit easier to read
  tmp <- clrs
  clrs[1] <- tmp[6]
  clrs[2] <- tmp[11]
  clrs[5] <- tmp[4]
  clrs[6] <- tmp[7]
  clrs[12] <- tmp[1]
  clrs[13] <- tmp[1]
  clrs[4] <- tmp[14]
  clrs[14] <- tmp[2]
  clrs[7] <- tmp[5]
  
  p <- ggplot(mytree, aes(y=total_counts, x=delta, color=taxa)) +
    geom_point(shape=20)+
    scale_color_manual(values=clrs) +
    theme( axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 8),
           axis.title.y = element_text(size = 14))+
    theme_tufte()+
    geom_rug(outside = TRUE, color="slategray2")+ 
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = c(0.1, 0.1), limits = c( min(mytree$total_counts)-3, max(mytree$total_counts)+3)) +
    scale_x_continuous(limits=c( min(mytree$delta)-stretch_factor, max(mytree$delta)+stretch_factor)) + 
    labs(title=mytitle, 
         x="Log change in frequency B versus M", y="Log of total number of counts") +
    geom_vline(xintercept=leftpt_freq, color="thistle3") +
    geom_vline(xintercept=rightpt_freq, color = "thistle3") +
    geom_text(aes(x=leftpt_freq, label="\n 95% CI:", y=2), colour="red", angle=90, text=element_text(size=3)) +
    
    geom_hline(yintercept=leftpt_counts, color="lightskyblue") + 
    geom_text(aes(y=leftpt_counts, label=paste0("\n Qnt:",lc_quant), x=4), colour="seagreen4", angle=0, text=element_text(size=3)) +
    
    geom_hline(yintercept=rightpt_counts, color = "lightskyblue") +
    geom_text(aes(y=rightpt_counts, label=paste0("Qnt:",rc_quant), x=4), colour="seagreen4", angle=0, text=element_text(size=3)) +
    
    geom_vline(xintercept=tmp2[1], color="lightskyblue") + 
    geom_text(aes(x=tmp2[1], label=paste0("Qnt:", lf_quant ), y=0), colour="seagreen4", angle=90, text=element_text(size=3)) +
    geom_vline(xintercept=tmp2[2], color = "lightskyblue")+
    geom_text(aes(x=tmp2[2], label=paste0("\n Qnt:", rf_quant ), y=0), colour="seagreen4", angle=90, text=element_text(size=3)) +
    
    geom_text_repel(
      max.iter=100000,
      aes(label=subset(mytree, ((delta < tmp2[1]) &
                                        (total_counts > leftpt_counts) ) )$name),
      size = size,
      data = subset(mytree, ((delta < tmp2[1]) & (total_counts > leftpt_counts) ) ),
      nudge_x       = top_left,   segment.size  = 0.1,    direction     = "y",    hjust         = 0.5, force=1
    ) +
    geom_text_repel(
      max.iter=100000,
      aes(label=subset(mytree, ((delta > tmp2[2])  & (total_counts > leftpt_counts) ) )$name),
      size = size,
      data = subset(mytree, ((delta > tmp2[2])  & (total_counts > leftpt_counts) ) ),
      nudge_x= top_right,    segment.size  = 0.1,    direction     = "y",
      force =5   #, hjust         = 0
    ) +
    geom_text_repel(
      max.iter=100000,
      force=5,
      aes(label=subset(mytree, ((delta > tmp2[1]) & (delta < tmp2[2]) & 
                                        (total_counts > rightpt_counts) ) )$name),
      size = size,
      data = subset(mytree, ((delta > tmp2[1]) & (delta < tmp2[2]) & 
                                     (total_counts > rightpt_counts) ) ), 
      segment.size =0.1, nudge_y = 4, nudge_x = top_even
     # nudge_x = top_even,
    #   angle        = 45, segment.size  = 0.0,  direction     = "x"   , hjust =0,
    #  force         = 5, # force_pull = 5,
    #  nudge_y = 12 #, vjust = 0.5
    ) 
  
   return(p)
}

 
make_figure <- function(tax_target, name, figurefile) {
  
  name <- gsub("[/ ]", "_", name)
  tree <- original # resets the tree just in case
  recursive_isLeaf( 1 )
  
  cat("\n", tax_target); mt(tax_target)  # target
  
  children <- list()
  for (i in 1:length(p2c(tax_target))) {
    tree <- original
    tmp <- induce_tree( p2c(tax_target)[i])  
    tmp$taxa <- tmp[1, "name"]
    children[[i]] <- tmp
  }
  tree <- original
  target_species <- do.call("rbind", children)
  target_species <- target_species[ which(target_species$rank == "genus"), ]
  
  
  relative_to <- tax_target
  p <- freq_count_plot( mytree=target_species, relative_to = relative_to, mytitle = "", verbose = TRUE,
                        lf_quant = 0.15, rf_quant = 0.89, lc_quant = 0.45, rc_quant = 0.99,
                        top_left = -6, top_even = 1, top_right = 6,
                        stretch = 0.8 ,
                        size = 3.0
  )
  p
#  ggsave( filename = paste0(paste0(paste0(name,".genus.rel_", relative_to)), ".png"), path = figurefile, device = "png", dpi = 300)
  
 # children <- list()
#  for (i in 1:length(p2c(tax_target))) {
#    tree <- original
#    tmp <- induce_tree( p2c(tax_target)[i])  
#    tmp$taxa <- tmp[1, "name"]
#    children[[i]] <- tmp
#  }
#  tree <- original
#  target_species <- do.call("rbind", children)
#  target_species <- target_species[ which(target_species$rank == "species"), ]
  
 # relative_to <- tax_target
#  p <- freq_count_plot( mytree=target_species, relative_to = relative_to, mytitle = paste0( name, " (species level, relative to root of target)"), verbose = TRUE,
 #                       lf_quant = 0.15, rf_quant = 0.89, lc_quant = 0.42, rc_quant = 0.99,
  #                      top_left = -20, top_even = -10, top_right = 15,
  #                      stretch = 1.0,
  #                      size = 2)
  #p
  #ggsave( filename = paste0(paste0(paste0(name,".species.rel_", relative_to)), ".png"), path = figurefile, device = "png", dpi = 300)
  
}


