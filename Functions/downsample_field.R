source("Functions/SubsampleTimeseriesBlock_highresNA.R")
downsample_field <- function(fld,dim_lon,dim_lat){
  dim_fld_lon = dim(fld)[1]
  dim_fld_lat = dim(fld)[2]
  
  #first get longitude straight
  fld_lon_new = array(dim = c(dim_lon, dim_fld_lat))
  for(ii in 1:dim_fld_lat){
    fld_lon_new[,ii] = as.numeric(SubsampleTimeseriesBlock_highresNA(ts(data = fld[,ii], start = 0.5*360/dim_fld_lon, 
                                                                        end = (360-360/dim_fld_lon/2), 
                                                                        deltat = 360/dim_fld_lon),
                                                                     seq(0,360,length.out = dim_lon)))
    
  }
  #next, get latitide straight
  fld_lon_lat_new = array(dim = c(dim_lon, dim_lat))
  for(ii in 1:dim_lon){
    fld_lon_lat_new[ii,] = as.numeric(SubsampleTimeseriesBlock_highresNA(ts(data = fld_lon_new[ii,], start = 0.5*180/dim_fld_lat, 
                                                                            end = (180-180/dim_fld_lat/2), 
                                                                            deltat = 180/dim_fld_lat),
                                                                         seq(0,180,length.out = dim_lat)))
  }
  
  
  
  
  return(fld_lon_lat_new)
}