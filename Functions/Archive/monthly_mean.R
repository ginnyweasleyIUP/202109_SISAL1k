monthlymean <- function(timeseries){
  years = floor(length(timeseries)/12)
  new_ts <- numeric(years)
  for(ii in 1:years){
    new_ts[ii] = mean(timeseries[(12*ii-11):(12*ii)], na.rm = T)
  }
  return(new_ts)
}

prec_weighting <- function(timeseries, prec_TS){
  years = floor(length(timeseries)/12)
  new_ts <- numeric(years)
  for(ii in 1:years){
    new_ts[ii] = sum(prec_TS[(12*ii-11):(12*ii)]*timeseries[(12*ii-11):(12*ii)], na.rm = T)/sum(prec_TS[(12*ii-11):(12*ii)], na.rm = T)
  }
  return(new_ts)
}

inf_weighting <- function(timeseries, prec_TS, evap_TS){
  years = floor(length(timeseries)/12)
  new_ts <- numeric(years)
  prec_inf = prec_TS-evap_TS
  prec_inf[prec_inf<0] = 0
  for(ii in 1:years){
    if(sum(prec_inf[(12*ii-11):(12*ii)], na.rm = T) == 0){
      new_ts[ii] = NA
    }else{
      new_ts[ii] = sum(prec_inf[(12*ii-11):(12*ii)]*timeseries[(12*ii-11):(12*ii)], na.rm = T)/sum(prec_inf[(12*ii-11):(12*ii)], na.rm = T)
    }
    
  }
  return(new_ts)
}