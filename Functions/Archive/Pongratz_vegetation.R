#################################################
##2) SISAL TIME SERIES ##########################
#################################################

MODEL <- list()

ncf <- ncdf4::nc_open("/home/ginnyweasley/Dropbox/SISAL1k/Data/Pongratz_800-1992_vegtype000010.nc")
MODEL$VEG1 <- ncdf4::ncvar_get(ncf)
MODEL$lon <- ncdf4::ncvar_get(ncf, 'longitude')
MODEL$lat <- ncdf4::ncvar_get(ncf, 'latitude')
ncdf4::nc_close(ncf)

ncf <- ncdf4::nc_open("/home/ginnyweasley/Dropbox/SISAL1k/Data/Pongratz_800-1992_vegtype000009.nc")
MODEL$VEG2 <- ncdf4::ncvar_get(ncf)
ncdf4::nc_close(ncf)



# needs to be imported first, as then only the relevant cave sites will be extracted and calculated further
CAVES <- list()
print("...SISAL extracting")
source("Functions/SISAL_extracting.R")

data <- load_sisal_data(min_period = 600, min_dating = 3, min_d18O = 36, used_dates = F)

CAVES$entity_info <- data[[1]]

CAVES$site_info <- CAVES$site_info %>% filter(site_id %in% CAVES$entity_info$site_id)
for (ii in CAVES$entity_info$entity_id){
  name = paste0("ENTITY", ii)
  site <- CAVES$entity_info %>% filter(entity_id == ii) %>% distinct(site_id)
  CAVES$record_data[[name]] <- data[[2]] %>% filter(entity_id == ii) %>% distinct(entity_id, mineralogy, interp_age, d18O_measurement, d13C_measurement) %>%
    mutate(site_id = (site$site_id))
}

remove(data, site, ii, name, load_sisal_data)


#################################################
## 3) Extract data from Caves such that they are in a grid box that is the average of all surrounding
#################################################

print("...DATA extract")
source("Functions/extract_gridboxes_v2.R")

#for(Model in c("HadCM3", "ECHAM5", "CCSM", "CESM", "GISS")){
for (ii in 1:(length(CAVES$entity_info$entity_id))){
  lon_cave = CAVES$entity_info$longitude[ii]
  
  if(lon_cave<0){lon_cave = 360+lon_cave}
  
  lat_cave = CAVES$entity_info$latitude[ii]
  
  ratios <- extract_gridboxes_2(lon_cave, lat_cave, d.lon = 360/length(MODEL$lon), d.lat = 180/length(MODEL$lat))
  name <- paste0("ENTITY",CAVES$entity_info$entity_id[ii])
  
  for(var in c("VEG1", "VEG2")){
    CAVES$sim_data[[name]][[var]] <- rowSums(cbind(ratios$Q11*MODEL[[var]][ratios$Q11_lon, ratios$Q11_lat,],
                                                     ratios$Q12*MODEL[[var]][ratios$Q12_lon, ratios$Q12_lat,],
                                                     ratios$Q21*MODEL[[var]][ratios$Q21_lon, ratios$Q21_lat,],
                                                     ratios$Q22*MODEL[[var]][ratios$Q22_lon, ratios$Q22_lat,]), na.rm = T)
    
  }
}

remove(ratios, ii, lat_cave, lon_cave, name, extract_gridboxes_2, var)

# Yearly Data
print("...write csv file")
MONTHLY <- list()
for(entity in sort(CAVES$entity_info$entity_id)){
  data_new = array(dim = c(1292,5))
  data_new[,1] = CAVES$entity_info$site_id[CAVES$entity_info$entity_id == entity]
  data_new[,2] = entity
  data_new[,3] = seq(701,1992,1)
  data_new[,4] = CAVES$sim_data[[paste0("ENTITY", entity)]]$VEG1
  data_new[,5] = CAVES$sim_data[[paste0("ENTITY", entity)]]$VEG2

  colnames(data_new) = c("site_id", "entity_id", "year_BP", "VEG1", "VEG2")
  
  if(entity == 14){data = data_new}
  else{data = rbind(data, data_new)}
  
}
print("...almost done...file writing")
DATA_EXPORT_MONTHLY <- data
write.csv(DATA_EXPORT_MONTHLY, file = paste0("Data/SISAL1k_Pongratz.csv"), row.names = F)
