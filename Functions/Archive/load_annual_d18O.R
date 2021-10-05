source("Functions/clear_data_matrix.R")

load_annual_d18O_field <- function(Model){
  if(Model == "HadCM3"){
    ncf <- ncdf4::nc_open("Data/Sim_data/HadCM3_d18O_850-1850.nc")
    ISOT_MONTHLY <- ncdf4::ncvar_get(ncf, "dO18")
    ncdf4::nc_close(ncf)
    
    ncf <- ncdf4::nc_open("Data/Sim_data/HadCM3_prec_850-1850.nc")
    PREC <- clear_data_matrix_neighbour(ncdf4::ncvar_get(ncf))*86400*360
    ncdf4::nc_close(ncf)
    ISOT_MONTHLY[PREC <= 0.001] = NA
    
    ISOT_MONTHLY = clear_data_matrix_neighbour(ISOT_MONTHLY)
  }
  if(Model == "ECHAM5"){
    ncf <- ncdf4::nc_open("Data/Sim_data/ECHAM5_d18O_850-1849.nc")
    ISOT_MONTHLY <- clear_data_matrix_neighbour(ncdf4::ncvar_get(ncf))
    ncdf4::nc_close(ncf)
  }
  if(Model == "CESM"){
    ncf <- ncdf4::nc_open("Data/Sim_data/CESM_d18O_850-1850.nc")
    ISOT_MONTHLY <- clear_data_matrix_neighbour(ncdf4::ncvar_get(ncf)[,96:1,1:12000])
    ncdf4::nc_close(ncf)
  }
  if(Model == "GISS"){
    ncf <- ncdf4::nc_open("Data/Sim_data/GISS_d18O_850-1849.nc")
    ISOT_MONTHLY <- clear_data_matrix_neighbour(abind::abind(ncdf4::ncvar_get(ncf)[73:144,90:1,],ncdf4::ncvar_get(ncf)[1:72,90:1,], along = 1))
    ncdf4::nc_close(ncf)
  }
  if(Model == "isoGSM"){
    ncf <- ncdf4::nc_open("Data/Sim_data/CCSM_prec_850-1849.nc")
    PREC <- clear_data_matrix_neighbour(ncdf4::ncvar_get(ncf)[,94:1,1:12000]*86400*360)
    ncdf4::nc_close(ncf)
    
    ncf <- ncdf4::nc_open("~/Dokumente/01_Promotion/06_Daten/11_CCSM/monthly/flxmon-prate1sfc.nc")
    ISOT2 <- clear_data_matrix_neighbour(ncdf4::ncvar_get(ncf)[,94:1,1:12000]*86400*360)
    ncdf4::nc_close(ncf)
    
    ISOT = (ISOT2/PREC -1)*1000
    ISOT[ISOT>50] = NA
    ISOT[ISOT< -100] = NA
    ISOT[PREC <= 0.001] = NA
    ISOT_MONTHLY = clear_data_matrix_neighbour(ISOT)
    rm(PREC,ISOT2, ISOT)
  }
  
  ISOT = array(dim = c(dim(ISOT_MONTHLY)[c(1,2)], dim(ISOT_MONTHLY)[3]/12))
  
  for(lon in 1:dim(ISOT_MONTHLY)[1]){
    if(lon%%10==0){print(lon)}
    for(lat in 1:dim(ISOT_MONTHLY)[2]){
      ISOT[lon,lat,] = zoo::rollapply(ISOT_MONTHLY[lon,lat,],3,mean, by = 12, na.rm = T)
    }
  }
  
  return(ISOT)
}
