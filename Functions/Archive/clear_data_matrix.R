clear_data_matrix <- function(data, type){
  
  #temp, prec, isot
  max_value=array(c(70, 2.5e-4, 10, 40))
  min_value=array(c(-100, 0, -70, -40))
  
  n = 0
  
  for(lon in 1:dim(data)[1]){
    for(lat in 1:dim(data)[2]){
      for(ii in 1:dim(data)[3]){
        if(is.na(data[lon,lat,ii])){
          next
        } else if(data[lon,lat,ii] > max_value[type]){
          data[lon,lat,ii] <- NA
          n = n +1
        } else if(data[lon,lat,ii] < min_value[type]){
          data[lon,lat,ii] <- NA
          n = n + 1
        }
      }
    }
  }
  print(n)
  return(data)
}


clear_data_matrix_neighbour <- function(data, type = NULL){
  
  #this is a 3D neighbour approach...
  data[is.infinite(abs(data))] = NA
  if(type == "PREC"){
    data[data<0] = NA
    data[data>20000] = NA
  }
  if(type == "ISOT"){
    data[data>50] = NA
    data[data< -100] = NA
  }
  if(type == "TEMP"){
    data[data>  100] = NA
    data[data< -100] = NA
  }
  
  n = 0
  
  for(lon in 1:dim(data)[1]){
    print(lon)
    for(lat in 1:dim(data)[2]){
      for(time in 2:(dim(data)[3]-1)){
        if(is.na(data[lon,lat,time])){
          next
        }
        lon_before = lon-1
        lon_after = lon+1
        lat_before = lat-1
        lat_after = lat+1
        if(lon == 1){lon_before = dim(data[1])}
        if(lon == dim(data)[1]){lon_after = 1}
        if(lat == 1){lat_before = dim(data[2])}
        if(lat == dim(data)[2]){lat_after = 1}
          neighbours = c(data[lon_before,c(lat_before,lat,lat_after),(time-1):(time+1)],
                         data[lon_after,c(lat_before,lat,lat_after),(time-1):(time+1)],
                         data[lon,c(lat_after,lat_before), (time-1):(time+1)],
                         data[lon,lat,c((time-1),(time+1))])
          neighbours = na.omit(neighbours[!is.infinite(neighbours)])
          if(length(neighbours) < 2){
            data[lon,lat,time] = NA
            next
          }
          mean_neigh = mean(neighbours, na.rm = T)
          std_neigh = sd(neighbours, na.rm = T)

        if(data[lon,lat,time]>mean_neigh+5*std_neigh | data[lon,lat,time]<mean_neigh-5*std_neigh){
          #print(paste(lon,lat,time))
          data[lon,lat,time] = NA
          n = n+1
        }
      }
    }
  }
  print(n)
  return(data)
}
