source("Functions/STACYmap_PMIL_logscale.R")

isotMap <- function(isot_lyr, leg_name, point_lyr = NULL){
  Plot_lyr <- isot_lyr
  Plot_lyr[is.na(Plot_lyr)] = 1000
  Plot_lyr[Plot_lyr>0] <- Plot_lyr[Plot_lyr>0]+1 
  Plot_lyr[Plot_lyr<0] <- Plot_lyr[Plot_lyr<0]-1
  Plot_lyr[Plot_lyr>0] <- log10(Plot_lyr[Plot_lyr>0])
  Plot_lyr[Plot_lyr<0] <- - log10(abs(Plot_lyr[Plot_lyr<0]))
  Plot_lyr[abs(Plot_lyr)>5] <- NA
  
  # layer <- point_lyr$layer
  # layer[is.na(layer)] = 1000
  # layer[layer>0] <- layer[layer>0]+1 
  # layer[layer<0] <- layer[layer<0]-1
  # layer[layer>0] <- log10(layer[layer>0])
  # layer[layer<0] <- - log10(abs(layer[layer<0]))
  # layer[abs(layer)>5] <- NA
  
  point_lyr$layer = - log10(abs(point_lyr$layer -1))
  
  plot <- STACYmap_isot(gridlyr = Plot_lyr, ptlyr = point_lyr,
                        legend_names = list(grid =leg_name),
                        graticules = TRUE,
                        colorscheme = rev(RColorBrewer::brewer.pal(11, 'BrBG'))[c(1:8,10)],
                        breaks_isot = c(-log10(71),-log10(11),-log10(6), -log10(2), 0, log10(2), log10(6)),
                        labels_isot = c(-70, -10,-5, -1, 0, 1, 5),
                        min_grid = -log10(71), max_grid = log10(12))
  return(plot)
}

