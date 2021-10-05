library(plyr)
library(dplyr)
library(rgdal)
library(latex2exp)

source("Functions/STACYmap_PMIL.R")
source("Functions/projection_ptlyr.R")


coastline_map <- rgdal::readOGR(dsn ="Functions/naturalearth_10m_physical/ne_10m_coastline.shp", verbose = FALSE)
karst_map <- rgdal::readOGR(dsn = "Functions/naturalearth_10m_physical/karst_wgs.shp", verbose = FALSE)

karst_map_plot <- function(USED_SITES_d13C,
                           USED_SITES_rest, 
                           NOT_SITES, 
                      projection = as.character('+proj=robin +datum=WGS84'),
                      pt_size){
  
  karst_map2 <- karst_map %>% spTransform(., CRSobj = CRS(projection)) %>% fortify(.)
  coastline_map2 <- coastline_map %>% spTransform(., CRSobj = CRS(projection)) %>% fortify(.)
  
  USED_SITES_d13C_p <- projection_ptlyr(USED_SITES_d13C, projection)
  USED_SITES_rest_p <- projection_ptlyr(USED_SITES_rest, projection)
  NOT_SITES_p <- projection_ptlyr(NOT_SITES, projection)
  
  leg_name = c(expression(paste(delta^{plain(18)}, plain(O), " and ", delta^{plain(13)}, plain(C))), 
               expression(paste("only ",delta^{plain(18)}, plain(O))))
  
  plot <- STACYmap(coastline = TRUE, filledbg = TRUE) +
    geom_polygon(data = karst_map2, aes(x=long, y = lat, group = group), fill = '#A57F7C', color = NA) +
    new_scale_color() +
    #geom_point(data = NOT_SITES_p, aes(x = long, y = lat, shape = "4", color = '4'),# shape = 20,
    #           size = (pt_size), alpha = 0.7, show.legend = c(shape =TRUE)) +
    geom_point(data = USED_SITES_rest_p, aes(x = long, y = lat, shape = "3", color = '3'),# shape = 17,
               size = (pt_size), alpha = 0.9, show.legend = c(shape =TRUE)) +
    geom_point(data = USED_SITES_d13C_p, aes(x = long, y = lat, shape = "2", color = '2'),# shape = 17,
               size = (pt_size), alpha = 0.9, show.legend = c(shape =TRUE)) +
    # scale_color_manual(name = "SISAL v2 sites", labels = c("d18O and d13C", "only d18O", "other"), 
    #                     values = c('#DD2525', '#F98D11', '#000000')) +
    # scale_shape_manual(name = "SISAL v2 sites", labels = c("d18O and d13C", "only d18O", "other"), 
    #                    values = c(17, 17, 20)) +
    scale_color_manual(name = "SISAL v2 sites", labels = leg_name,
                       values = c(rgb(178, 24, 43, max = 255), rgb(33, 102, 172, max = 255))) +
    scale_shape_manual(name = "SISAL v2 sites", labels = leg_name,
                       values = c(17, 17)) +
    #geom_polygon(data = coastline_map2, aes(x=long, y=lat, group = group), color = 'black',  size = 0.2, fill = NA, alpha = 0.8) +
    theme(legend.position = c(0.02, 0.02),
          legend.justification = c(0, 0),
          legend.box = 'vertical',
          legend.box.background = element_blank(),
          legend.background = element_rect(colour = 'black'),
          panel.border = element_blank())
  
  rm(coastline_map2, karst_map2, NOT_SITES_p)
  rm(USED_SITES_rest_p, USED_SITES_d13C_p)
  
  return(plot)
}

