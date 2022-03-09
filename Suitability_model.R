library (raster)
library (plotly)
library(spatial.tools)
library(spatialEco)
library(whitebox)
library(quantreg)
library(rgrass7)
setwd("F:\\research_project\\Beta3")

if (!dir.exists('outputs')) dir.create('outputs') 
if (!dir.exists('intermediates')) dir.create('intermediates') 

#read rasters
dsm <- round(raster('inputs/dsm.tif'))
treecover <- raster('inputs/treecover.tif')
potapov_treeheight <- raster('inputs/potapov_treeheight.tif')
simard_treeheight <- raster('inputs/simard_treeheight.tif')
soils <- raster('inputs/soils.tif')
gl30 <- raster('inputs/gl_land_cover.tif')

dsm[is.na(dsm)] <- 0

treecover[is.na(treecover)] <- 0

potapov_treeheight[potapov_treeheight > 60 ] <- 0
potapov_treeheight[is.na(potapov_treeheight)] <- 0

simard_treeheight[simard_treeheight > 60 ] <- 0
simard_treeheight[is.na(simard_treeheight)] <- 0

gl30[gl30 ==255] <- 0


positions <- c(5,9,13)
## Ok. Stick with me. This section splits raster indexes 
index_y_start <- c(1,1,1,1,1,1)
index_y_end <- c(round(nrow(dsm)/2),round(nrow(treecover)/2),round(nrow(potapov_treeheight)/2),round(nrow(simard_treeheight)/2),round(nrow(soils )/2),round(nrow(gl30 )/2))
index_x_start <- c(1,1,1,1,1,1)
index_x_end <- c(round(ncol(dsm)/2),round(ncol(treecover)/2),round(ncol(potapov_treeheight)/2),round(ncol(simard_treeheight)/2),round(ncol(soils )/2),round(ncol(gl30 )/2))
index_name <- c('dsm', 'treecover', 'potapov_treeheight', 'simard_treeheight', 'soils', 'gl30')

############### Practically magic. Dont think about it too much. Just run it. it works ############

for (a in 1: length(index_y_start)){

assign(paste0("index_",index_name[a]), data.frame(  index_y_start[a], index_y_end[a], #TOP LEFT
                          index_x_start[a], index_x_end[a],
                   
                          index_y_start[a], index_y_end[a], ##TOP RIGHT
                          index_x_end[a],  ncol(get(index_name[a])),
                   
                          index_y_end[a], nrow(get(index_name[a])),        ##BOTTOM LEFT
                          index_x_start[a], index_x_end[a],
                   
                          index_y_end[a],  nrow(get(index_name[a])),       ##BOTTOM RIGHT
                          index_x_end[a],  ncol(get(index_name[a]))  ) )   

}


############# Believe it or not, we can crop the rasters in a forloop. Saves... code. ####################
for ( position in positions ){

for( b in 1: length(index_name)){
assign( paste0(index_name[b],'_cropped'), crop(get(index_name[b]), extent(get(index_name[b]),
                                                    get(paste0("index_",index_name[b]))[1,position], 
                                                    get(paste0("index_",index_name[b]))[1,position+1],
                                                    get(paste0("index_",index_name[b]))[1,position+2],
                                                    get(paste0("index_",index_name[b]))[1,position+3])))
  print(index_name[b])
}


#RASTER RESAMPLING 
gl30_cropped  <- resample(gl30_cropped , dsm_cropped, method =  'ngb')
#soils_cropped  <- resample(soils_cropped , dsm_cropped, method =  'ngb')
treecover_cropped <- resample(treecover_cropped, dsm_cropped, method =  'ngb')
potapov_treeheight_cropped <- resample(potapov_treeheight_cropped, dsm_cropped, method = 'ngb')
simard_treeheight_cropped <- resample(simard_treeheight_cropped, dsm_cropped, method = 'ngb')


#TREE OFFSET METHOD

pot_binary_map <- reclassify(potapov_treeheight_cropped, c(0,5, 0, 5, Inf, 1))

binaryMap <- reclassify(treecover_cropped, c(0, 80 , 0,
                                     80, 100,   2))

comb <- binaryMap + pot_binary_map

simard_mask <- comb
simard_mask[simard_mask != 1] <- 0

simard_treeheight_masked <- simard_mask * simard_treeheight_cropped

treeheight <- potapov_treeheight_cropped
treeheight[simard_mask == 1 ] <- simard_treeheight_masked[simard_mask == 1]

comb[comb >=2 ] <- 1
treeheight <- comb * treeheight

treeheight_gaussian <- focal(treeheight, w = gaussian.kernel(2,5))

dsm_class <- reclassify(dsm_cropped, c(-Inf, 30, 0, 30 , Inf , 1))

dsm_class[dsm_class == 1] <- 0.9
dsm_class[dsm_class == 0] <- 0.6
dsm_class_gauss <-  focal(dsm_class, w = gaussian.kernel(2,5))

dem <- dsm_cropped - (dsm_class_gauss * treeheight_gaussian)
dem[dem < 0 ] <- 0
writeRaster(dem, 'intermediates/dem.tif', overwrite=TRUE)

rm(dsm_cropped)
rm(comb)
rm(binaryMap)
rm(dsm_class)
rm(dsm_class_gauss)
rm(treecover_cropped)
rm(treeheight)
rm(treeheight_gaussian)
rm(potapov_treeheight_cropped)
rm(simard_mask)
rm(simard_treeheight_cropped)
rm(simard_treeheight_masked)
rm(pot_binary_map)
rm(dem)
gc()

#SMOOTHING METHOD
wbt_lee_sigma_filter('F:\\research_project\\Beta3\\intermediates\\dem.tif', "F:\\research_project\\Beta3\\intermediates\\dems.tif", 
                     filterx = 5, filtery = 5, sigma = 3 , m = 9)
dems <- raster('intermediates/dems.tif')



slope_degrees <- terrain(dems, opt="slope", unit='degrees', neighbors=8)
slope_degrees[is.na(slope_degrees)] <- 0
slope_percent <- calc(slope_degrees, fun = function(x){round(tan(x*pi/180)*100,2)})
#slope_percent[slope_percent >=100] <- 100

writeRaster(slope_percent, paste0("outputs/slopes_", position,".tif"), overwrite = TRUE)

slopes_percent_class <- reclassify(slope_percent, c(0,15,0,10, Inf, 1))
slopes_percent_class <- aggregate(slopes_percent_class, fact = c(30,30), fun = median)
slopes_percent_class <- focal(slopes_percent_class, w=matrix(1,3,3) , fun = modal )
slopes_percent_class <- disaggregate(slopes_percent_class, fact = c(30,30))


soils_cropped[is.na(soils_cropped)] <- 0

merge_class <- soils_cropped
merge_class[slopes_percent_class == 1] <- 7
writeRaster(merge_class, 'intermediates/merge_class.tif', overwrite = TRUE)


depressions <- gl30_cropped
depressions[gl30_cropped == 10] <- 1
depressions[gl30_cropped == 20] <- 0
depressions[gl30_cropped == 30] <- 1
depressions[gl30_cropped == 40] <- 1
depressions[gl30_cropped == 50] <- 1
depressions[gl30_cropped == 60] <- 0
depressions[gl30_cropped == 70] <- 1
depressions[gl30_cropped == 80] <- 1
depressions[gl30_cropped == 90] <- 1
depressions[gl30_cropped == 100] <- 1
depressions[depressions >  100] <- 0

hydro <- gl30_cropped
hydro[hydro != 60] <- 0
hydro[hydro == 60 ] <- 1

hydro_gauss <- focal(hydro, w= gaussian.kernel(sigma=0.5, n =3) )
writeRaster(hydro_gauss, 'intermediates/hydro_gauss.tif', overwrite= TRUE)
dem_h <- dems - 5*(hydro_gauss)
writeRaster(dem_h, paste0('outputs/dem_h_', position,'.tif'), overwrite=TRUE)
writeRaster(depressions, 'intermediates/depressions.tif', overwrite=TRUE)
rm(dems)
rm(gl30_cropped)
rm(hydro_gauss)
rm(merge_class)
rm(slope_degrees)
rm(slope_percent)
rm(slopes_percent_class)
rm(soils_cropped)
gc()


###### INTERMEDIATE INPUT FILES.
loc <- initGRASS("C:/Program Files/GRASS GIS 7.8", 
                 home=paste0(getwd(),'/Grass'), location = "area_1", 
                 gisDbase = paste0(getwd(),'/Grass'),
                 mapset = "PERMANENT", override = TRUE)

execGRASS("g.proj", flags = c("c"), epsg = 4326)
execGRASS("r.import", flags = c("overwrite"), parameters=list(input=paste0(getwd(),paste0('/outputs/dem_h_', position,'.tif')), output="dem_h"))
execGRASS("r.import", flags = c("overwrite"), parameters=list(input=paste0(getwd(),'/intermediates/depressions.tif'), output="depressions"))
execGRASS("g.region", parameters=list( raster="dem_h"))
execGRASS('r.watershed', flags=c( "overwrite","s"), parameters = list(elevation="dem_h", depression= 'depressions', accumulation='flow_accum'  ))
execGRASS("r.out.gdal", flags=c("overwrite"), parameters= list(input = 'flow_accum', output=paste0(getwd(),"/intermediates/flow_accum.tif"), format="GTiff"))

flow_accum <- raster("intermediates/flow_accum.tif")
drainage_values <- read.csv('drainage/drain_stats.csv')
drainage_weights <- raster('intermediates/merge_class.tif')

for (drainage_class in drainage_values$classmajority ){
drainage_weights[drainage_weights == drainage_class ] <- 1000/drainage_values$median[drainage_values$classmajority == drainage_class ]
}

hydro_flow <- hydro * 1000
hydro_flow[hydro_flow == 0] <- 1

weighted_flow_accum <- flow_accum * drainage_weights
weighted_flow_accum <- weighted_flow_accum * hydro_flow
weighted_flow_accum <- (weighted_flow_accum * ((depressions-1)*(-1)) )
weighted_flow_accum[is.null(weighted_flow_accum)] <- 0
weighted_flow_accum[weighted_flow_accum > 5000] <- 5000

writeRaster(weighted_flow_accum, 'intermediates/weighted_flow_accumulation.tif',overwrite= TRUE)

execGRASS("r.import", flags = c("overwrite"), parameters=list(input=paste0(getwd(),'/outputs/dem_h_', position,'.tif'), output="dem_h"))
execGRASS("r.import", flags = c("overwrite"), parameters=list(input=paste0(getwd(),'/intermediates/weighted_flow_accumulation.tif'), output="weighted_flow_accum"))
execGRASS("g.region", parameters=list( raster="weighted_flow_accum"))
execGRASS('r.stream.extract', flags=c( "overwrite", "verbose"), parameters = list(elevation="dem_h", accumulation='weighted_flow_accum', stream_length= 5,threshold=1000,
                                                                       stream_raster="streams", direction="drain_dir"  ))

execGRASS("r.out.gdal", flags=c("overwrite"), parameters= list(input = 'streams',  output=paste0(getwd(),"/outputs/streams_", position,".tif"), format="GTiff") )
execGRASS("r.out.gdal", flags=c("overwrite"), parameters= list(input = 'drain_dir', output=paste0(getwd(),"/intermediates/drain_dir.tif"), format="GTiff") )

execGRASS("g.region", parameters=list( raster="streams"))
execGRASS('r.stream.distance', flags=c( "overwrite"), parameters = list(stream_rast = "streams", elevation="dem_h", direction="drain_dir" , difference = "hand" ))

execGRASS("r.out.gdal", flags=c("overwrite"), parameters= list(input = 'hand', output=paste0(getwd(),"/outputs/hand_", position,".tif"), format="GTiff") )

rm(weighted_flow_accum)
rm(hydro_flow)
rm(hydro)
rm(drainage_values)
rm(drainage_weights)
rm(depressions)
gc()
#SUITABILITY MAP METHOD

slopes <-raster(paste0('outputs/slopes_', position,'.tif'))
hand <-  raster(paste0('outputs/hand_', position,'.tif'))
srtm <- raster(paste0('outputs/dem_h_', position,'.tif'))
FCM <- raster('inputs/potapov_treeheight.tif')


#slopes
slopes_suitable <- aggregate(slopes, fact = c(3,3), fun = mean)
slopes_suitable <- reclassify(slopes_suitable, c(0, 10, 100,
                                                 10, 15, 75,
                                                 15, 20, 50,
                                                 20, 25, 25,
                                                 25, 30, 10,
                                                 30, Inf, 0)  )
writeRaster(slopes_suitable, paste0("outputs/slopes_suitabilityIndex_", position,".tif"), overwrite  = TRUE)

#hand
hand <- resample(hand, srtm, method='ngb')
hand_suitable <- aggregate(hand, fact = c(3,3), fun = mean)
hand_suitable <- reclassify(hand_suitable, c(-Inf, 0, 0,
                                                0,  0.9,   0,
                                              0.9,  2,    10,
                                                2,  5,    25,
                                                5, 10,    50,
                                               10, 20,    75,
                                               20, Inf , 100)  )
hand_suitable[is.na(hand_suitable)] <- 0
writeRaster(hand_suitable, paste0("outputs/hand_suitabilityIndex_", position,".tif"), overwrite  = TRUE)

#srtm
srtm_suitable <- aggregate(srtm, fact = c(3,3), fun = mean) 
srtm_suitable <- reclassify(srtm_suitable, c(-Inf, 0, 0,
                                             0,  2,   0,
                                             2,  5,   10,
                                             5,  10,  25,
                                             10, 25,  75,
                                             25, 400,  100,
                                             400, Inf, 0)  )

writeRaster(srtm_suitable,paste0("outputs/srtm_suitabilityIndex_", position,".tif"), overwrite  = TRUE)


#fcm
FCM <- resample( FCM, slopes, method="bilinear")
fcm_suitable <-  aggregate(FCM, fact = c(3,3), fun = modal) 
fcm_suitable <- reclassify(fcm_suitable,   c(0,  5,   0,
                                             5,  10,  10,
                                             10, 20,  25,
                                             20, 25,  50,
                                             25, 60, 100,
                                             60, Inf, 0)  )

writeRaster(fcm_suitable, paste0("outputs/fcm_suitabilityIndex_", position,".tif"), overwrite  = TRUE)


rabmethod <- function(x,...){
  
  si <- x[1] * sqrt(x[2]/100 * x[3]/100 * x[4]/100) 
  
  return(si)
}

sqrmethod <- function(x,...){
  
  rmin <- min(x)
  
  y <- x[x != which.min(x)]
  
  si <- rmin * sqrt(y[1]/100 * y[2]/100 * y[3]/100) 
  
  return(si)
}

storiemethod <- function(x,...){
  
  si <- x[1] * x[2]/100 * x[3]/100 * x[4]/100
  
  return(si)
}

#DJANEL CROISANT GEROOKTE zalm

rasterstack <- stack(slopes_suitable, srtm_suitable, hand_suitable, fcm_suitable)

rabiaSuitabilityIndex <- aggregate(rasterstack, fact= c(1,1,4), fun = rabmethod)
sqrSuitabilityIndex <- aggregate(rasterstack, fact= c(1,1,4), fun = sqrmethod)
storieSuitabilityIndex <- aggregate(rasterstack, fact= c(1,1,4), fun = storiemethod)

writeRaster(rabiaSuitabilityIndex, paste0("outputs/rabia_suitabilityIndex_", position,".tif"), overwrite  = TRUE)

writeRaster(sqrSuitabilityIndex, paste0("outputs/sqr_suitabilityIndex_", position,".tif"), overwrite  = TRUE)

writeRaster(storieSuitabilityIndex, paste0("outputs/storie_suitabilityIndex_", position,".tif"), overwrite  = TRUE)

rm(FCM)
rm(fcm_suitable)
rm(hand_suitable)
rm(hand)
rm(slopes_suitable)
rm(slopes)
rm(srtm_suitable)
rm(srtm)
rm(rasterstack)
rm(rabiaSuitabilityIndex)
rm(sqrSuitabilityIndex)
rm(storieSuitabilityIndex)
gc()

}

