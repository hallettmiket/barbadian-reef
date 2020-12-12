

select_freq_count_plot <- function( inmyspecies, relative_to = 1, mytitle="cleaning", verbose = FALSE, 
                             lf_quant = 0.22, rf_quant = 0.87, lc_quant = 0.5, rc_quant = 0.99,
                             top_left = -10, top_even = -10, top_right = 7,
                             stretch = 1,
                             size = 2) {
  
  require(pals)
  root_bel <-max(1,tree[t2i(relative_to), "br_bel"]) ; root_may <-max(1, tree[t2i(relative_to), "br_may"])
  
  myspecies <- tree[ inmyspecies, ]
  myspecies$total_counts <- log(myspecies$br_bel+myspecies$br_may)
  myspecies$delta <- log( ((myspecies$br_bel+1) / root_bel ) /  ((myspecies$br_may+1) / root_may ) )
  
  no_zero <- intersect(which(myspecies$br_bel > 0), which(myspecies$br_may > 0)) 
  delta_no_zero <- log( (myspecies[no_zero, "br_bel"] / root_bel ) /  (myspecies[no_zero, "br_may"] / root_may ) )
  
  results_freq <- boot( data = delta_no_zero, statistic = mystat, R = 10000)
  ci <- boot.ci(results_freq)
  leftpt_freq <- ci$norm[2]; rightpt_freq <- ci$norm[3]
  
  
  tmp2 <- quantile( myspecies$delta, c(lf_quant, rf_quant))
  tmp <- quantile( myspecies$total_counts, c(lc_quant, rc_quant))
  leftpt_counts <- tmp[1]; rightpt_counts <- tmp[2]
  
  if (verbose) {
    cat("\n Counts at ancestor: ", root_bel, root_may)
    cat("\n Mean: ", mean(myspecies$delta) )
    cat("\n Confidence Interval for mean: ", ci$norm)
    cat("\n Quantile on Frequencies:", tmp2 )
    cat("\n Quantile on counts:", tmp)
  }
  
  stretch_factor <-  stretch * (max(myspecies$delta) -  min(myspecies$delta))
  
  num_taxa <- length(unique(myspecies$taxa))
  
  clrs <- glasbey(); 
  # make the colors a bit easier to read
  tmp <- clrs[2]; clrs[2] <- clrs[6]; clrs[6] <- tmp
  tmp <- clrs[1]; clrs[1] <- clrs[12]; clrs[12] <- tmp
  tmp <- clrs[3]; clrs[3] <- clrs[12]; clrs[12] <- tmp
  
  
  # To use for line and point colors, add
  
  p <- ggplot(myspecies, aes(y=total_counts, x=delta)) +
    geom_point(shape=20 )+
    scale_color_manual(values=clrs) +
    theme( axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 8),
           axis.title.y = element_text(size = 14))+
    theme_tufte()+
    geom_rug(outside = TRUE, color="slategray2")+ 
    coord_cartesian(clip = "off") +
    scale_y_continuous(expand = c(0.1, 0.1), limits = c( min(myspecies$total_counts), max(myspecies$total_counts)+1)) +
    scale_x_continuous(limits=c( min(myspecies$delta)-stretch_factor, max(myspecies$delta)+stretch_factor)) + 
    labs(title=mytitle, 
         x="Log ratio of frequency B versus M", y="Log of total number of counts") +
    geom_vline(xintercept=leftpt_freq, color="thistle3") +
    geom_vline(xintercept=rightpt_freq, color = "thistle3") +
    geom_text(aes(x=leftpt_freq, label="\n 95% CI:", y=3), colour="red", angle=90, text=element_text(size=size)) +
    
    geom_hline(yintercept=leftpt_counts, color="lightskyblue") + 
    geom_text(aes(y=leftpt_counts, label=paste0("Qnt:",lc_quant), x=-7), colour="seagreen4", angle=0, text=element_text(size=size)) +
    
    geom_hline(yintercept=rightpt_counts, color = "lightskyblue") +
    geom_text(aes(y=rightpt_counts, label=paste0("Qnt:",rc_quant), x=-7), colour="seagreen4", angle=0, text=element_text(size=size)) +
    
    geom_vline(xintercept=tmp2[1], color="lightskyblue") + 
    geom_text(aes(x=tmp2[1], label=paste0("Qnt:", lf_quant ), y=2.4), colour="seagreen4", angle=90, text=element_text(size=size)) +
    geom_vline(xintercept=tmp2[2], color = "lightskyblue")+
    geom_text(aes(x=tmp2[2], label=paste0("Qnt:", rf_quant ), y=2.4), colour="seagreen4", angle=90, text=element_text(size=size)) +
    
    geom_text_repel(
      max.iter=100000,
      aes(label=subset(myspecies, ((delta < tmp2[1]) &
                                  (total_counts > leftpt_counts) ) )$name),
      size = size,
      data = subset(myspecies, ((delta < tmp2[1]) & (total_counts > leftpt_counts) ) ),
      nudge_x       = top_left,   segment.size  = 0.1,    direction     = "y",    hjust         = 0.5, force=1
      
    ) +
    geom_text_repel(
      max.iter=100000,
      aes(label=subset(myspecies, ((delta > tmp2[2])  & (total_counts > leftpt_counts) ) )$name),
      size = size,
      data = subset(myspecies, ((delta > tmp2[2])  & (total_counts > leftpt_counts) ) ),
      nudge_x= top_right,    segment.size  = 0.1,    direction     = "y" #, hjust         = 0
    ) +
    geom_text_repel(
      max.iter=100000,
      aes(label=subset(myspecies, ((delta > tmp2[1]) & (delta < tmp2[2]) & 
                                  (total_counts > rightpt_counts) ) )$name),
      size = size,
      data = subset(myspecies, ((delta > tmp2[1]) & (delta < tmp2[2]) & 
                               (total_counts > rightpt_counts) ) ), 
 
      nudge_x = top_even, segment.size  = 0.1,
      #nudge_y= 4, 
      angle        = 0  ,  direction     = "y" 
      
    ) 
  
  return(p)
}



