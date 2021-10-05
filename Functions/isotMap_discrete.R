source("Functions/STACYmap_PMIL_discrete.R")

isotMap_discrete <- function(isot_lyr, leg_name, point_lyr = NULL){
  Plot_lyr <- isot_lyr
  Plot_lyr[is.na(Plot_lyr)] = 100000
  Plot_lyr[Plot_lyr>0] <- Plot_lyr[Plot_lyr>0]+1 
  Plot_lyr[Plot_lyr<0] <- Plot_lyr[Plot_lyr<0]-1
  Plot_lyr[Plot_lyr>0] <- log10(Plot_lyr[Plot_lyr>0])
  Plot_lyr[Plot_lyr<0] <- - log10(abs(Plot_lyr[Plot_lyr<0]))
  Plot_lyr[abs(Plot_lyr)>=5] <- NA
  
  # layer <- point_lyr$layer
  # layer[is.na(layer)] = 1000
  # layer[layer>0] <- layer[layer>0]+1 
  # layer[layer<0] <- layer[layer<0]-1
  # layer[layer>0] <- log10(layer[layer>0])
  # layer[layer<0] <- - log10(abs(layer[layer<0]))
  # layer[abs(layer)>5] <- NA
  
  if(!is.null(point_lyr)){point_lyr$layer = - log10(abs(point_lyr$layer -1))}
  
  colorbar = rev(c(RColorBrewer::brewer.pal(11, 'BrBG')[c(2,4,5,6)],RColorBrewer::brewer.pal(11, 'Spectral')[7:11]))
  #colorbar = rev(RColorBrewer::brewer.pal(11, 'BrBG'))[c(1:8,10)]

  source("Functions/STACYmap_PMIL_discrete.R")
  
  if(!is.null(point_lyr)){
    plot <- STACYmap_discrete(gridlyr =  Plot_lyr, ptlyr = point_lyr,
                              legend_names = list(grid = leg_name), graticules = T, colorscheme = colorbar, 
                              limits_color = list("max" = log10(11), "min" = -log10(71)),
                              labels_discrete = c(10,5, 2,1,0,-1, -2,-5,-10,-20,"", -70), 
                              discrete_breaks = c(log10(11), log10(6),log10(3),log10(1.54),log10(1),-log10(1.54), -log10(3), -log10(6), -log10(11), -log10(21), -log10(36), -log10(71))) +
      theme(panel.border = element_blank(), legend.background = element_blank(), 
            axis.text = element_blank(), legend.text = element_text(size = 10),plot.title = element_text(hjust = 0.5))
  }else{
    plot <- STACYmap_discrete(gridlyr =  Plot_lyr,
                              legend_names = list(grid = leg_name), graticules = T, colorscheme = colorbar, 
                              limits_color = list("max" = log10(11), "min" = -log10(71)),
                              labels_discrete = c(10,5, 2,1,0,-1, -2,-5,-10,-20,"", -70), 
                              discrete_breaks = c(log10(11), log10(6),log10(3),log10(1.54),log10(1),-log10(1.54), -log10(3), -log10(6), -log10(11), -log10(21), -log10(36), -log10(71))) +
      theme(panel.border = element_blank(), legend.background = element_blank(), 
            axis.text = element_blank(), legend.text = element_text(size = 10),plot.title = element_text(hjust = 0.5))
  }
  
  
  
  #plot
  
  return(plot)
}





