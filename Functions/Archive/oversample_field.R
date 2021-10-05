oversample_field <- function(fld,dim_lon,dim_lat){
  dim_fld_lon = dim(fld)[1]
  dim_fld_lat = dim(fld)[2]
  
  #first get longitude straight
  fld_lon_new = array(dim = c(dim_lon, dim_fld_lat))
  for(ii in 1:dim_fld_lat){
    fld_lon_new[,ii] = as.numeric(PaleoSpec::MakeEquidistant(t.x = seq(0,360,length.out = dim_fld_lon),
                                                             t.y = fld[,ii],
                                                             time.target = seq(0,360,length.out = dim_lon)))
  }
  #next, get latitide straight
  fld_lon_lat_new = array(dim = c(dim_lon, dim_lat))
  for(ii in 1:dim_lon){
    fld_lon_lat_new[ii,] = PaleoSpec::MakeEquidistant(t.x = seq(0,180,length.out = dim_fld_lat),
                                                      t.y = fld_lon_new[ii,],
                                                      time.target = seq(0,180,length.out = dim_lat))
  }
  
  return(fld_lon_lat_new)
}