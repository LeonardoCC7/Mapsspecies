######################################################### 
#### Evolutionary Biogeography and Spatial Modelling ####
#### Hands-on session on Conservation Biogeography   ####
#########################################################

# Questions: 
# 1 - Which bioclimatic variables seem to have more influence in explaining the 
#     predicted occurrence of each species;
# 2 - Discuss model performance results for each of the study species;
# 3 - Debate the loss/ gain in area of habitat suitability for both species in 50 and 70 years, 
#     under CO2 emissions of 2.6 (optimistic scenario) and 8.5 (pessimistic scenario);
# 4 - Which species is currently more protected by the Protected Areas (PAs) 
#     network? And which will be more protected in 50 and 70 years?
# 5 - Considering all results above and the literature available, discuss possible 
#     actions to enhance current and future protection of these endangered species

# Set Working directory
# Go to "Session" > "Set Working Directory" > "Choose Directory..."

# Install packages needed to run the hands-on session
install.packages(c("sdm","dplyr","tidyr","rgbif","maps","mapview","wdpar","ggmap"))
installAll(sdm)

# load packages
require(sdm)
require(dplyr)
require(tidyr)
library(rgbif)
library(maps)
library(mapview)
library(wdpar)
library(ggmap)

#2 Russia (Carnivores)----
sp1 = dismo::gbif("Panthera","tigris",download= T,geo = T,removeZeros = T)
sp2 = dismo::gbif("Ursus","thibetanus",download= T,geo = T,removeZeros = T)

table(sp1$basisOfRecord)
sp1 = sp1 %>%
  filter(basisOfRecord %in% c("HUMAN_OBSERVATION","LIVING_SPECIMEN","OCCURRENCE","PRESERVED_SPECIMEN"))
nrow(sp1)
table(sp2$basisOfRecord)
sp2 = sp2 %>%
  filter(basisOfRecord %in% c("HUMAN_OBSERVATION","LIVING_SPECIMEN","OCCURRENCE","PRESERVED_SPECIMEN"))
nrow(sp2)

sp1 = sp1 %>%
  filter(country %in% "Russia")
sp2 = sp2 %>%
  filter(country %in% "Russia")

sp1 = sp1 %>%
  select(lon,lat)
sp2 = sp2 %>%
  select(lon,lat)

sp1 = sp1 %>%
  drop_na()
sp2 = sp2 %>%
  drop_na()

# Remove duplicates
sp1 = unique(sp1)
sp2 = unique(sp2)

# Check sample size
dim(sp1)
dim(sp2)

# Download WorldClim bio climatic environmental variables
bio_raw = raster::getData("worldclim",var="bio",res=5)
# NOTE - Import downloaded RasterStack from your computer, to avoid an error latter while running the vifstep function
#rm(bio_raw)
#bio_raw = raster::stack(list.files(path = "/Users/vitorhpaiva/wc5",pattern=".bil$",full.names=T))
plot(bio_raw[[1]])
points(sp1, col = "blue")
points(sp2, col = "red""
geo_area = raster::drawExtent() # Draw extent directly on the map, then run the next line
bio_crop = raster::crop(bio_raw,geo_area)
# Confirm that all locs are inside your study area
plot(bio_crop[[1]])
points(sp1, col = "blue")
points(sp2, col = "red")

sp1$pres_sp1 = 1
sp2$pres_sp2 = 1

coordinates(sp1) = c("lon","lat")
proj4string(sp1) = raster::projection(raster::raster())
coordinates(sp2) = c("lon","lat")
proj4string(sp2) = raster::projection(raster::raster())

bio_val_sp1 = raster::extract(bio_crop,sp1)
bio_val_sp2 = raster::extract(bio_crop,sp2)

bio_vif_sp1 = usdm::vifstep(bio_val_sp1)
bio_vif_sp2 = usdm::vifstep(bio_val_sp2)

bio_sp1 = usdm::exclude(bio_crop,bio_vif_sp1)
bio_sp2 = usdm::exclude(bio_crop,bio_vif_sp1)

#### SDM ####

# Creating Background
d_sp1 = sdmData(pres_sp1~.,
                sp1,
                predictors = bio_sp1,
                bg = list(method = "gRandom", n = 1000))
d_sp2 = sdmData(pres_sp2~.,
                sp2,
                predictors = bio_sp2,
                bg = list(method = "gRandom", n = 1000))

# Model fitting
getmethodNames()
m_sp1 = sdm(pres_sp1~.,
            d_sp1,
            c('glm','brt','rf','fda'),
            replication = "sub",
            test.p = 30,
            n = 3)
m_sp1
m_sp2 = sdm(pres_sp2~.,
            d_sp2,
            c('glm','brt','rf','fda'),
            replication = "sub",
            test.p = 30,
            n = 3)
m_sp2

# Check model statistics and performance
gui(m_sp1)
gui(m_sp2)

# Model Prediction
#p_sp1 = predict(bio_sp1,
#                m_sp1)
#p_sp2 = predict(bio_sp2,
#                m_sp2)

color = colorRampPalette(c("#3E49BB",
                           "#3498DB",
                           "yellow",
                           "orange",
                           "red",
                           "darkred"))

color2 = colorRampPalette(c("red",
                            "orange",
                            "yellow",
                            "grey",
                            "green",
                            "blue"))

# Ensemble predictions ----
en_Today_sp1 = ensemble(m_sp1, 
                  bio_sp1, 
                  filename='en_Today_sp1.img',
                  setting=list(method='weighted',
                               stat='tss',
                               opt=2),
                  overwrite=T)
en_Today_sp2 = ensemble(m_sp2, 
                  bio_sp2, 
                  filename='en_Today_sp2.img',
                  setting=list(method='weighted',
                               stat='tss',
                               opt=2),
                  overwrite=T)

# Inspect ensemble models prediction on dynamic maps
mapview(en_Today_sp1,
        col.regions = color(200),
        na.color = "transparent") + sp1
mapview(en_Today_sp2,
        col.regions = color(200),
        na.color = "transparent") + sp2

mapview(en_Today_sp1,
        col.regions = color(200),
        na.color = "transparent") +
  mapview(en_Today_sp2,
        col.regions = color(200),
        na.color = "transparent") + sp2 + sp1


#_____________________________________________
# Future projections - RCP 2.6 - 50 years ----
#_____________________________________________

biof = raster::getData('CMIP5', var='bio', res=5, rcp=26, model='CN', year=50)

names(biof)
plot(biof[[2]])

biof_crop = raster::crop(biof,geo_area)
plot(biof_crop[[2]])

names(biof_crop) = names(bio_raw)

enf_50_sp1 = ensemble(m_sp1, 
                   biof_crop,
                   filename='enf_50_sp1.img',
                   setting=list(method='weighted',
                                stat='tss',
                                opt=2),
                   overwrite=T)
plot(enf_50_sp1)
enf_50_sp2 = ensemble(m_sp2,
                   biof_crop,
                   filename='enf_50_sp2.img',
                   setting=list(method='weighted',
                                stat='tss',
                                opt=2),
                   overwrite=T)
plot(enf_50_sp2)

cl = colorRampPalette(c('#3E49BB','#3498DB','yellow','orange','red','darkred'))

par(mfrow=c(1,2)) # more plots in the same frame (1 row, 2 columns)
plot(en_Today_sp1, col=cl(200), main="Current")
plot(enf_50_sp1, col=cl(200), main="Projected (RCP 2.6, 50 years)")

par(mfrow=c(1,2)) # more plots in the same frame (1 row, 2 columns)
plot(en_Today_sp2, col=cl(200), main="Current")
plot(enf_50_sp2, col=cl(200), main="Projected (RCP 2.6, 50 years)")

proj4string(sp1) = projection(en_Today_sp1)
proj4string(sp2) = projection(en_Today_sp2)

p_en_Today_sp1 = en_Today_sp1
p_en_Today_sp2 = en_Today_sp2
p_en_Today_sp1[] = ifelse(en_Today_sp1[] >= 0.5, 1, 0)
p_en_Today_sp2[] = ifelse(en_Today_sp2[] >= 0.5, 1, 0)

p_enf_50_sp1 = enf_50_sp1
p_enf_50_sp2 = enf_50_sp2
p_enf_50_sp1[] = ifelse(enf_50_sp1[] >= 0.5, 1, 0)
p_enf_50_sp2[] = ifelse(enf_50_sp2[] >= 0.5, 1, 0)

# Continuous 
#mapview(en_Today_sp1,col.regions=cl(200)) +
  #mapview(enf_50_sp1,col.regions=cl(200)) +
  #sp1
#mapview(en_Today_sp2,col.regions=cl(200)) +
#mapview(enf_50_sp2,col.regions=cl(200)) +
 # sp2

mapview(en_Today_sp1,col.regions=cl(200)) +
  mapview(enf_50_sp1,col.regions=cl(200)) +
  mapview(en_Today_sp2,col.regions=cl(200)) +
  mapview(enf_50_sp2,col.regions=cl(200)) +
  sp2 + sp1


ch_sp1 = enf_50_sp1 - en_Today_sp1
cl2 = colorRampPalette(c('red','orange','yellow','gray','green','blue'))
plot(ch_sp1,col=cl2(200),main="Habitat suitability variation sp1")
ch_sp2 = enf_50_sp2 - en_Today_sp2
cl2 = colorRampPalette(c('red','orange','yellow','gray','green','blue'))
plot(ch_sp2,col=cl2(200),main="Habitat suitability variation sp2")

ev_sp1 = getEvaluation(m_sp1, stat = c("AUC","TSS","threshold"), opt=2)
mean(ev_sp1$threshold) # Use the value obtained here two lines down
pa_sp1 = raster(en_Today_sp1)
pa_sp1[] = ifelse(en_Today_sp1[] >= 0.5887408, 1, 0)
cl = colorRampPalette((c("red","grey80","blue")))
plot(pa_sp1,col=cl(3),main="Predicted Current presences sp1")
paf_sp1 = raster(enf_50_sp1)
paf_sp1[] = ifelse(enf_50_sp1[] >= 0.5887408, 1, 0)
cl = colorRampPalette((c("red","grey80","blue")))
plot(paf_sp1,col=cl(3),main="Predicted Future presences sp1")

ev_sp2 = getEvaluation(m_sp2, stat = c("AUC","TSS","threshold"), opt=2)
mean(ev_sp2$threshold) # Use the value obtained here two lines down
pa_sp2 = raster(en_Today_sp2)
pa_sp2[] = ifelse(en_Today_sp2[] >= 0.7343112, 1, 0)
cl = colorRampPalette((c("red","grey80","blue")))
plot(pa_sp2,col=cl(3),main="Predicted Current presences sp2")
paf_sp2 = raster(enf_70_sp2)
paf_sp2[] = ifelse(enf_70_sp2[] >= 0.7343112, 1, 0)
cl = colorRampPalette((c("red","grey80","blue")))
plot(paf_sp2,col=cl(3),main="Predicted Future presences sp2")

par(mfrow=c(1,2)) # more plots in the same frame (1 row, 2 columns)
chp_sp1 = paf_sp1 - pa_sp1
plot(chp_sp1,col=c('red','gray','blue'), main="Variation Predicted presences sp1")
chp_sp2 = paf_sp2 - pa_sp2
plot(chp_sp2,col=c('red','gray','blue'), main="Variation Predicted presences sp2")

#_________________________________________________________________
# Inspect the coverage of PAs in relation to our study species----

# Read data obtained from Protected Planet (https://www.protectedplanet.net/en) 
wdpa_latest_version()
webdriver::install_phantomjs()

# download protected area data for multiple of countries
raw_pa_data =
  "Russia" %>%
  lapply(wdpa_fetch, wait = TRUE,
         download_dir = rappdirs::user_data_dir("wdpar")) %>%
  bind_rows()

pa_data = wdpa_clean(raw_pa_data)

#mapview(en_Today_sp1,col.regions=cl(200)) +
 # mapview(enf_50_sp1,col.regions=cl(200)) +
  #mapview(raw_pa_data,col.regions=cl(200)) +
  #sp1

#mapview(en_Today_sp2,col.regions=cl(200)) +
 # mapview(enf_50_sp2,col.regions=cl(200)) +
  #mapview(pa_data,col.regions=cl(200)) +
  #sp2

mapview(en_Today_sp1,col.regions=color(200)) +
  mapview(enf_50_sp1,col.regions=color(200)) + sp1 +
  mapview(en_Today_sp2,col.regions=color(200)) +
  mapview(enf_50_sp2,col.regions=color(200)) + sp2 +
  mapview(raw_pa_data,col.regions="green")
#_____________________________________________
# Future projections - RCP 8.5 - 70 years ----
#_____________________________________________

biof = raster::getData('CMIP5', var='bio', res=5, rcp=85, model='CN', year=70) #### IDEA FOR LECTURES - change this to 50 years and 70 years + CHANGE RCP scenario!! ####

names(biof)
plot(biof[[2]])

biof_crop = raster::crop(biof,geo_area)
plot(biof_crop[[2]])

names(biof_crop) = names(bio_raw)

enf_70_sp1 = ensemble(m_sp1, 
                   biof_crop,
                   filename='enf_70_sp1.img',
                   setting=list(method='weighted',
                                stat='tss',
                                opt=2),
                   overwrite=T)
plot(enf_70_sp1)
enf_70_sp2 = ensemble(m_sp2,
                   biof_crop,
                   filename='enf_70_sp2.img',
                   setting=list(method='weighted',
                                stat='tss',
                                opt=2),
                   overwrite=T)
plot(enf_70_sp2)

cl = colorRampPalette(c('#3E49BB','#3498DB','yellow','orange','red','darkred'))

par(mfrow=c(1,2)) # more plots in the same frame (1 row, 2 columns)
plot(en_Today_sp1, col=cl(200), main="Current")
plot(enf_70_sp1, col=cl(200), main="Projected (RCP 8.5, 70 years)")

par(mfrow=c(1,2)) # more plots in the same frame (1 row, 2 columns)
plot(en_Today_sp2, col=cl(200), main="Current")
plot(enf_70_sp2, col=cl(200), main="Projected (RCP 8.5, 70 years)")

proj4string(sp1) = projection(en_Today_sp1)
proj4string(sp2) = projection(en_Today_sp2)

p_en_Today_sp1 = en_Today_sp1
p_en_Today_sp2 = en_Today_sp2
p_en_Today_sp1[] = ifelse(en_Today_sp1[] >= 0.5, 1, 0)
p_en_Today_sp2[] = ifelse(en_Today_sp2[] >= 0.5, 1, 0)

p_enf_70_sp1 = enf_70_sp1
p_enf_70_sp2 = enf_70_sp2
p_enf_70_sp1[] = ifelse(enf_70_sp1[] >= 0.5, 1, 0)
p_enf_70_sp2[] = ifelse(enf_70_sp2[] >= 0.5, 1, 0)

mapview(en_Today_sp1,col.regions=cl(200)) +
  mapview(enf_70_sp1,col.regions=cl(200)) +
  sp1
mapview(en_Today_sp2,col.regions=cl(200)) +
  mapview(enf_70_sp2,col.regions=cl(200)) +
  sp2

mapview(en_Today_sp1,col.regions=cl(200)) +
  mapview(enf_70_sp1,col.regions=cl(200)) +
  sp1
mapview(en_Today_sp2,col.regions=cl(200)) +
  mapview(enf_70_sp2,col.regions=cl(200)) +
  sp2 + sp1

ch_sp1 = enf_70_sp1 - en_Today_sp1
cl2 = colorRampPalette(c('red','orange','yellow','gray','green','blue'))
plot(ch_sp1,col=cl2(200),main="Habitat suitability variation sp1")
ch_sp2 = enf_70_sp2 - en_Today_sp2
cl2 = colorRampPalette(c('red','orange','yellow','gray','green','blue'))
plot(ch_sp2,col=cl2(200),main="Habitat suitability variation sp2")

ev_sp1 = getEvaluation(m_sp1, stat = c("AUC","TSS","threshold"), opt=2)
mean(ev_sp1$threshold) # Use the value obtained here two lines down
pa_sp1 = raster(en_Today_sp1)
pa_sp1[] = ifelse(en_Today_sp1[] >= 0.5887408, 1, 0)
cl = colorRampPalette((c("red","grey80","blue")))
plot(pa_sp1,col=cl(3),main="Predicted Current presences sp1")
paf_sp1 = raster(enf_70_sp1)
paf_sp1[] = ifelse(enf_70_sp1[] >= 0.5887408, 1, 0)
cl = colorRampPalette((c("red","grey80","blue")))
plot(paf_sp1,col=cl(3),main="Predicted Future presences sp1")

ev_sp2 = getEvaluation(m_sp2, stat = c("AUC","TSS","threshold"), opt=2)
mean(ev_sp2$threshold) # Use the value obtained here two lines down
pa_sp2 = raster(en_Today_sp2)
pa_sp2[] = ifelse(en_Today_sp2[] >= 0.7343112, 1, 0)
cl = colorRampPalette((c("red","grey80","blue")))
plot(pa_sp2,col=cl(3),main="Predicted Current presences sp2")
paf_sp2 = raster(enf_70_sp2)
paf_sp2[] = ifelse(enf_70_sp2[] >= 0.7343112, 1, 0)
cl = colorRampPalette((c("red","grey80","blue")))
plot(paf_sp2,col=cl(3),main="Predicted Future presences sp2")

par(mfrow=c(1,2)) # more plots in the same frame (1 row, 2 columns)
chp_sp1 = paf_sp1 - pa_sp1
plot(chp_sp1,col=c('red','gray','blue'), main="Variation Predicted presences sp1")
chp_sp2 = paf_sp2 - pa_sp2
plot(chp_sp2,col=c('red','gray','blue'), main="Variation Predicted presences sp2")

#_________________________________________________________________
# Inspect the coverage of PAs in relation to our study species----

#mapview(en_Today_sp1,col.regions=cl(200)) +
 # mapview(enf_70_sp1,col.regions=cl(200)) +
  #mapview(pa_data,col.regions=cl(200)) +
  #sp1

#mapview(en_Today_sp2,col.regions=cl(200)) +
 # mapview(enf_70_sp2,col.regions=cl(200)) +
  #mapview(pa_data,col.regions=cl(200)) +
  #sp2

 mapview(en_Today_sp1,col.regions=color(200)) +
   mapview(enf_50_sp1,col.regions=color(200)) +
   mapview(enf_70_sp1,col.regions=color(200)) +sp1 +
   mapview(en_Today_sp2,col.regions=color(200)) +
   mapview(enf_50_sp2,col.regions=color(200)) + 
   mapview(enf_70_sp2,col.regions=color(200))+ sp2 +
   mapview(raw_pa_data,col.regions="green")

# End of script