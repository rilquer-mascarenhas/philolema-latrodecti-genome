setwd('/home/rilquer/lgeo_plat/')

require(flexsdm)
require(terra)
require(tidyverse)
require(sf)

## Reading environmental data ####
#files <- list.files('data/env_data/',pattern = '_v1.4.0.nc',full.names = T)
#env <- lapply(1:length(files),function(i){
#  var <- rast(files[i])
#  return(var[[800]])
#})
#env <- do.call(what = c,env)
#names(env) <- gsub('_800','',names(env))

## Alternatively, using worldclim with 2.5 min resolution
require(pastclim)
set_data_path('data/env_data/worldclim/')
bio <- c("bio01","bio04","bio05","bio06","bio07","bio08",
         "bio09","bio10","bio11","bio12","bio13","bio14",
         "bio15","bio16","bio17","bio18","bio19")
env <- region_slice(
  time_ce = 1985,
  bio_variables = bio,
  dataset = "WorldClim_2.1_2.5m",
  ext = ext(c(-65,-31,-40,10))
)

## Getting samples data ####
#samples <- read_csv('data/samples_info/thom_data/thom_samples_chosen.csv') %>%
#  select(sra_accession,Species,Longitude,Latitude) %>%
#  rename(code = 'sra_accession',binomial = 'Species',
#         longitude = 'Longitude',latitude = 'Latitude') %>%
#  add_row(read_csv('data/samples_info/mydata/samples.csv') %>%
#  select(code,binomial,longitude,latitude))

# `occ_download` for all Latrodectus
#https://docs.ropensci.org/rgbif/articles/getting_occurrence_data.html

### GBIF download ####
require(rgbif)
d_gbif=occ_download(
  type="and",
  pred("taxonKey", 2157920),
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus","PRESENT"), 
  pred_not(pred_in("basisOfRecord",c("FOSSIL_SPECIMEN"))),
  pred_or(  
    pred_lt("coordinateUncertaintyInMeters",10000),
    pred_isnull("coordinateUncertaintyInMeters")
  ),
  format = "SIMPLE_CSV"
)
occ_download_wait('0008802-260221153910048')
lat_occ <- occ_download_get('0008802-260221153910048') %>%
  occ_download_import()

lat_occ <- lat_occ %>% filter(taxonRank == 'SPECIES') %>% filter(species != '')
write_csv(lat_occ,'data/latrodecus_2157920.csv')

require(ggplot2)
require(ggspatial)
require(rnaturalearth)
## Map check
world <- ne_countries(scale = "medium", returnclass = "sf")
theme_set(theme_bw()) #Setting theme
spp <- unique(lat_occ$species)
for (x in spp) {
  coords <- lat_occ %>% filter(species == x) %>% select(decimalLongitude,decimalLatitude) %>% drop_na()
  ggplot(data = world) +
    geom_sf(fill= "ghostwhite", size = 0.1)+
    annotation_scale(location = "br", width_hint = 0.5) +
    annotation_north_arrow(location = "br", which_north = "true",
                           pad_x = unit(0.3, "in"), pad_y = unit(0.5, "in"),
                           style = north_arrow_fancy_orienteering) +
    geom_point(data = coords, aes(x = decimalLongitude,y = decimalLatitude),
               size = 1, color = 'black',alpha=0.6)+
    scale_x_continuous(name = "Longitude")+
    scale_y_continuous(name = "Latitude")+
    labs(title = x)+
    theme(legend.position = 'none',
          panel.grid.major = element_line(color = gray(.9),linetype = "dashed", linewidth = 0),
          panel.background = element_rect(fill = "aliceblue"),
          plot.title = element_text(face = "italic"))
  ggsave(paste0('output/spp_maps/',tolower(gsub(' ','_',x)),'.png'),
         width = 8, height = 6,dpi = 600)
}
    
### Map of L. geometricus invasion####
require(ggplot2)
require(ggspatial)
require(rnaturalearth)
world <- ne_countries(scale = "medium", returnclass = "sf")
theme_set(theme_bw()) #Setting theme

lgeo_occ <- read_csv('data/latrodectus_2157920.csv') %>% 
  filter(species == 'Latrodectus geometricus') %>% 
  filter(basisOfRecord != 'HUMAN_OBSERVATION') %>% 
  select(decimalLongitude,decimalLatitude,year) %>% drop_na()

ggplot(data = world) +
  geom_sf(fill= "ghostwhite", size = 0.1)+
  geom_point(data = lgeo_occ %>% filter(year <= 1920), aes(x = decimalLongitude,y = decimalLatitude),
             size = 2,alpha=0.8,color='brown')+
  theme(panel.grid.major = element_line(color = gray(.9),linetype = "dashed", linewidth = 0),
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
ggsave('output/lgeo_maps/lgeo_upto_1920.png',width = 5, height = 2.5,dpi = 600)

ggplot(data = world) +
  geom_sf(fill= "ghostwhite", size = 0.1)+
  geom_point(data = lgeo_occ %>% filter(year <= 1950), aes(x = decimalLongitude,y = decimalLatitude),
             size = 2,alpha=0.8,color='brown')+
  theme(panel.grid.major = element_line(color = gray(.9),linetype = "dashed", linewidth = 0),
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
ggsave('output/lgeo_maps/lgeo_upto_1950.png',width = 5, height = 2.5,dpi = 600)

ggplot(data = world) +
  geom_sf(fill= "ghostwhite", size = 0.1)+
  geom_point(data = lgeo_occ %>% filter(year <= 1990), aes(x = decimalLongitude,y = decimalLatitude),
             size = 2,alpha=0.8,color='brown')+
  theme(panel.grid.major = element_line(color = gray(.9),linetype = "dashed", linewidth = 0),
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
ggsave('output/lgeo_maps/lgeo_upto_1990.png',width = 5, height = 2.5,dpi = 600)

ggplot(data = world) +
  geom_sf(fill= "ghostwhite", size = 0.1)+
  geom_point(data = lgeo_occ %>% filter(year <= 2026), aes(x = decimalLongitude,y = decimalLatitude),
             size = 2,alpha=0.8,color='brown')+
  theme(panel.grid.major = element_line(color = gray(.9),linetype = "dashed", linewidth = 0),
        panel.background = element_rect(fill = "aliceblue"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
ggsave('output/lgeo_maps/lgeo_upto_2026.png',width = 5, height = 2.5,dpi = 600)


## Environmental filtering ####
## This step removes NA points, which deals with points that are outside the
## geographic area we set (erroneous points + C. pareola and P. mystaceus that are too far).
## This is done before setting calibration area to avoid creating erroneous calibration polygons.
spp <- sort(unique(lat_occ$species))

occs_filt <- lapply(occs,function(occs) {
  occs$id <- 1:nrow(occs) # adding unique id to each row
  return(occs %>%
           occfilt_env(data = .,x = "decimalLongitude",y = "",
                       id = "id",nbins = 8,env_layer = env) %>%
           left_join(occs, by = c("id", "x", "y")))
})
names(occs_filt) <- names(occs)

## Setting calibration area ####
ca <- lapply(occs_filt,function(x) {
  return(calib_area(data = x,x = "x",y = "y",method = c("buffer", width = 200000),crs = crs(env)))
})

## Map check
for (i in 1:length(occs_filt)) {
  # Plotting raster
  ggplot(data=data.frame(crds(env[[1]]),bio01=as.data.frame(env[[1]])))+
    geom_tile(aes(x=x,y=y,fill=bio01))+
    coord_fixed()+scale_fill_gradient2(low = '#4575b4',mid= '#ffffbf', high = '#d73027',
                                       midpoint = min(as.data.frame(env[[1]]))+((max(as.data.frame(env[[1]]))-min(as.data.frame(env[[1]])))/2))+
    geom_sf(data=st_as_sf(ca[[i]]),fill='lightgrey',alpha=0.3)+
    #Adding points
    geom_point(data=occs_filt[[i]],aes(x=x,y=y),alpha=0.5)+
    labs(title = spp[i])+
    theme(legend.position = 'none',
          panel.grid.major = element_line(color = gray(.9),linetype = "dashed", size = 0),
          panel.background = element_rect(fill = "aliceblue"),
          plot.title = element_text(face = "italic"))
  ggsave(paste0('output/enms/occ_points_maps/',tolower(gsub(' ','_',spp[i])),'.png'),
         width = 7, height = 9,dpi = 600)
}

## Spatial partition ####
## Using part_senv to decide best clustering based on env similarity, spatial autocorrelation and
## number of presence points in each partition.
## That method is not ideal for C. cearae, S. cearensis and X. atlanticus, because they have
## few points.
## For C. cearae, we will use part_sband to divide points into two latitudinal bands: a northern
## band with 4 points (Ceará and coast) and a southern band with remaining points (diamantina).
## For S. cearensis and X. atlanticus, the best approach seems to be leave-one-out
set.seed(10)
part <- lapply(1:length(occs_filt),function(i) {
  occs_filt[[i]]$pr_ab <- rep(1,nrow(occs_filt[[i]]))
  if (spp[i] %in% c('Conopophaga cearae','Sclerurus cearensis','Xiphorhynchus atlanticus')) {
    if (spp[i] == 'Conopophaga cearae') {
      # Latitudinal band for C. cearae
      return(occs_filt[[i]] %>% part_sband(data = .,pr_ab = "pr_ab",
                                           x = 'x', y = 'y',
                                           env_layer = env,type='lat',
                                           n_part = 2, min_occ = 3))
    } else {
      # Leave-one-out for S. cearensis and X. atlanticus
      return(occs_filt[[i]] %>% part_random(data = .,pr_ab = "pr_ab",
                                            method = c(method = 'loocv')))
    }
  } else {
    # Optimal partitions to Kmeans clustering for everyone else
    return(occs_filt[[i]] %>% part_senv(data = .,pr_ab = "pr_ab",
                                        x = 'x', y = 'y',
                                        env_layer = env,
                                        max_n_groups = (nrow(occs_filt[[i]])-2)))
  }
})

occs_filt <- lapply(1:length(part),function(i){
  if (is.null(nrow(part[[i]]))) {
    return(part[[i]]$part)
  } else {
    return(part[[i]])
  }
})

## Background and pseudo-absences ####
## Background points
set.seed(10)
bg <- lapply(1:length(occs_filt), function(i) {
  npart <- sort(unique(occs_filt[[i]]$.part))
  bg <- lapply(npart, function(x) {
    pts <- sample_background(
      data = occs_filt[[i]],
      x = "x",
      y = "y",
      n = sum(occs_filt[[i]]$.part == x) * 10,
      method = "random",
      rlayer = env[[1]],
      calibarea = ca[[i]])
    pts$.part <- rep(x,nrow(pts))
    return(pts)
  }) %>% bind_rows()
  return(bg)
})

## Pseudo-absence
set.seed(10)
psa <- lapply(1:length(occs_filt), function(i) {
  psa <- sample_pseudoabs(
    data = occs_filt[[i]],
    x = "x",
    y = "y",
    n = nrow(occs_filt[[i]]),
    method = "random",
    rlayer = env[[1]],
    calibarea = ca[[i]]
  )
  # Randomizing partition vector to psa points
  psa$.part <- sample(occs_filt[[i]]$.part)
  return(psa)
})

# Bind presences and pseudo-absences
occ_pa <- lapply(1:length(occs_filt),function(i){
  return(bind_rows(occs_filt[[i]], psa[[i]]))
})

# Extracting environmental data for presence-absence and background
occ_pa <- lapply(occ_pa, function(x){
  return(x %>% sdm_extract(
    data = .,
    x = "x",
    y = "y",
    env_layer = env,
    filter_na = TRUE))
})

bg <- lapply(bg,function(x){
  return(x %>%
           sdm_extract(
             data = .,
             x = "x",
             y = "y",
             env_layer = env,
             filter_na = TRUE
           ))
})

saveRDS(occ_pa,'RData/enms/flexsdm/occ_pa.rds')
saveRDS(bg,'RData/enms/flexsdm/bg.rds')

## Model tuning ####
# Maxent
model <- vector('list',length(occ_pa))
for (i in 1:length(model)) {
  model[[i]] <- tune_max(
    data = occ_pa[[i]],
    response = "pr_ab",
    predictors = names(env),
    background = bg[[i]],
    partition = ".part",
    grid = expand.grid(
      regmult = seq(0.5, 5, 0.5),
      classes = c("l", "lq", "lqp","lqhp","lqhpt")
    ),
    thr = c("max_sens_spec"),
    metric = "TSS",
    clamp = TRUE,
    pred_type = "cloglog"
  )
}
saveRDS(model,paste0('RData/enms/flexsdm/maxmodel_results_',Sys.Date(),'.rds'))

## Model output ####

# Reading model
model <- readRDS('RData/enms/flexsdm/maxmodel_results_2025-03-31.rds')

### Table of model results ####

lapply(model,function(x){return(x$performance)}) %>% bind_rows() %>%
  add_column(spp,.before = 'regmult') %>%
  select(spp,regmult,classes,n_presences,n_absences,OR_mean,OR_sd,TSS_mean,TSS_sd) %>%
  write_csv('RData/enms/flexsdm/models.csv')

### Response curves ####
# Here, we plot response curves by making predictions on the calibrated area
# 1) First, we mask climate data to the calibration area
# 2) Then, we retrieve climatic values for all points in the calibration area
# 3) To get response curve for a variable, we set all variables to their mean
# value except for the focal variable.
# 4) Then we make predictions on that modified dataset and plot the predictions
# along with the values of the variable
require(corrr)
bio_cor <- vector('list',length(model))
## Code for all response curves per species
resp_curves <- vector('list',length(model))
for (i in 1:length(model)) {
  # Masking environment
  env_spp <- mask(env,ca[[i]]) %>% crop(ext(ca[[i]]))
  
  # Getting clim values from all raster
  climdata <- as.points(env_spp) %>% as.data.frame()# Transform to spatial points
  
  # Creating correlation matrix for further investigation
  bio_cor[[i]] <- correlate(climdata) %>%
    pivot_longer(starts_with('bio'),names_to = 'var',values_to = 'rho') %>% 
    drop_na() %>% mutate(species = spp[i])
  
  # Making predictions for each variable
  pred <- lapply(1:length(bio),function(j) {
    org <- climdata[,bio[j]] # Saving original values of focal variable
    climdata <- climdata %>% mutate_all(mean) # Setting all variables to mean
    climdata[,bio[j]] <- org # Adding back the original values of focal variable
    pred <- predict(model[[i]]$model,climdata) %>% as.data.frame() %>%
      rename(pred = 'V1') # Renaming column to pred
    # Min-max normalization
    if (max(pred$pred)!=min(pred$pred)) { # First checking if response is not flat
      pred <- pred %>%
        mutate(pred = (pred-min(pred))/max((pred)-min(pred)))
    } else { # If it's flat, we just set pred to 0.5
      pred <- pred %>%
        mutate(pred = 0.5)
    }
    pred <- pred %>%
      mutate(values = climdata[,bio[j]], # Adding original values
             var = bio[j], # Adding name of variable
             spp = spp[i]) %>% # Adding name of species
      as_tibble()
    return(pred)
  })
  resp_curves[[i]] <- pred %>% bind_rows()
}

names(resp_curves) <- spp
# Saving correlation info
bio_cor <- bio_cor %>% bind_rows()
saveRDS(bio_cor,'RData/enms/bio_cor.rds')
saveRDS(resp_curves,'RData/enms/resp_curves.rds')



# Per species
lapply(1:length(resp_curves),function(i){
  resp_curves[[i]] %>% 
    ggplot(aes(x = values,y=pred))+geom_line()+
    facet_wrap(~var,ncol=4,scales = 'free')+
    labs(x="Bioclimatic values",y="Predicted suitability",
         title = spp[i])+
    theme_bw()+
    theme(panel.grid = element_line(size=0.2),
          axis.text = element_text(size=11),
          axis.title = element_text(size=15),
          strip.text = element_text(size=11),
          plot.title = element_text(face = "italic"))
  ggsave(paste0('output/enms/res_curves_per_spp/',tolower(gsub(' ','_',spp[i])),'.png'),
         width = 10,height = 9,dpi=600)
})

# Labels
lbl <- c('Bio 01',
         'Bio 04',
         'Bio 05',
         'Bio 06',
         'Bio 07',
         'Bio 08',
         'Bio 09',
         'Bio 10',
         'Bio 11',
         'Bio 12',
         'Bio 13',
         'Bio 14',
         'Bio 15',
         'Bio 16',
         'Bio 17',
         'Bio 18',
         'Bio 19')
names(lbl) <- c('bio01','bio04','bio05','bio06','bio07','bio08','bio09',paste0('bio',10:19))

# Northern species
resp_curves %>% bind_rows() %>%
  filter(spp %in% c("Chiroxiphia pareola","Conopophaga cearae","Hemitriccus mirandae",
                    "Platyrinchus mystaceus","Sclerurus cearensis","Xiphorhynchus atlanticus")) %>% 
  # Reducing species name
  mutate(spp = paste0(str_sub(str_split_i(spp,' ',1),1,1),'. ',str_split_i(spp,' ',2))) %>% 
  ggplot(aes(x = values,y=pred))+geom_line()+
  facet_grid(spp~var,scales = 'free',switch='both',
             labeller = labeller(var = lbl))+
  labs(x="Bioclimatic values",y="Species")+
  theme_bw()+
  theme(panel.grid = element_line(size=0.2),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=15),
        strip.text = element_text(size=12,face='italic'),
        plot.title = element_text(face = "italic"))
ggsave('output/enms/northernspp_res_curves.png',
       width = 19,height = 10,dpi=600)

# Southern species
resp_curves %>% bind_rows() %>%
  filter(!spp %in% c("Chiroxiphia pareola","Conopophaga cearae","Hemitriccus mirandae",
                     "Platyrinchus mystaceus","Sclerurus cearensis","Xiphorhynchus atlanticus")) %>% 
  # Reducing species name
  mutate(spp = paste0(str_sub(str_split_i(spp,' ',1),1,1),'. ',str_split_i(spp,' ',2))) %>% 
  ggplot(aes(x = values,y=pred))+geom_line()+
  facet_grid(spp~var,scales = 'free',switch='both',
             labeller = labeller(var = lbl))+
  labs(x="Bioclimatic values",y="Species")+
  theme_bw()+
  theme(panel.grid = element_line(size=0.2),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=15),
        strip.text = element_text(size=12,face='italic'),
        plot.title = element_text(face = "italic"))
ggsave('output/enms/southernnspp_res_curves.png',
       width = 24,height = 18,dpi=600)

# All curves
resp_curves %>% bind_rows() %>%
  # Reducing species name to acronyms
  mutate(spp = paste0(tolower(str_sub(str_split_i(spp,' ',1),1,3)),'_',str_sub(str_split_i(spp,' ',2),1,3))) %>% 
  ggplot(aes(x = values,y=pred))+geom_line()+
  facet_grid(spp~var,scales = 'free')+
  labs(x="Bioclimatic values",y="Species")+
  theme_bw()+
  theme(panel.grid = element_line(size=0.2),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size=15),
        strip.text = element_text(size=9),
        plot.title = element_text(face = "italic"))
ggsave('output/enms/all_res_curves.png',
       width = 15,height = 12,dpi=600)

# Extracting slope and max of curves for plotting

# We are extracting slope using the differences
# between consecutive points in the x and y axes.
# We use `diff()` to calculate these differences for x and y
# and then use dy / dx to calculate how much y differs at
# every change of x. In cases where x repeats consecutively,
# diff would be 0 and the division yields NA, so we remove
# NA values to disconsider differences between equal x values.
slope <- lapply(1:length(resp_curves),function(i){
  resp_curves[[i]] %>% group_by(var) %>% arrange(values) %>%
    mutate(dy = c(diff(pred),0),
           dx = c(diff(values),0)) %>%
    mutate(slope = dy / dx) %>%
    drop_na() %>%
    summarize(slope = mean(slope)) %>%
    mutate(spp = spp[i],
           slope = case_when(slope == 0 ~ NA, .default = slope))
}) %>% bind_rows()

# For max values, we retrieved the value of x axis where y
# reaches its peak. This will give a hint at the optimum
# bioclimatic value for the species.

slope_max <- lapply(1:length(resp_curves),function(i){
  d <- resp_curves[[i]] %>% group_by(var) %>%
    summarize(max_suit = max(pred)) %>%
    # We retrieve the slope of bioclim for that species
    # to correct for this later
    left_join(slope %>% filter(spp == names(resp_curves)[i]),by = 'var')
  max_bio = apply(d,1,function(x){
    var <- resp_curves[[i]] %>% filter(var == x[[1]])
    return(mean(var$values[which(var$pred==as.numeric(x[[2]]))]))
  })
  res <- d %>% mutate(max_bio) %>%
    # Changing to 0 cases where slope = 0
    mutate(max_bio = case_when(is.na(slope) ~ NA, .default = max_bio))
  return(res)
}) %>% bind_rows()

# Plotting
# Slope
require(RColorBrewer)
c <- brewer.pal(11, 'RdYlGn')
slope_max %>%
  mutate(spp = paste0(tolower(str_sub(str_split_i(spp,' ',1),1,3)),'_',str_sub(str_split_i(spp,' ',2),1,3))) %>% 
  ggplot(aes(x = var,y=fct_rev(spp),fill=slope))+geom_tile()+
  scale_fill_gradient2(name = 'Slope',
                       low = "#A50026", mid= "#FFFFBF", high = "#006837",
                       na.value="azure")+
  labs(x="Bioclimatic values",y="Species")+
  theme_bw()+
  theme(panel.grid = element_line(size=0.2),
        axis.text.x = element_text(size=15,vjust = .5,angle = 90),
        axis.text.y = element_text(size=15),
        axis.ticks = element_blank(),
        axis.title = element_text(size=15))
ggsave('output/enms/res_curves/curves_slope.png',
       width = 9,height = 9,dpi=600)

# Max Bio
require(RColorBrewer)
c <- brewer.pal(11, 'RdYlGn')
slope_max %>%
  mutate(spp = paste0(tolower(str_sub(str_split_i(spp,' ',1),1,3)),'_',str_sub(str_split_i(spp,' ',2),1,3))) %>% 
  ggplot(aes(x = var,y=fct_rev(spp),fill=max_bio))+geom_tile()+
  scale_fill_gradient2(name = 'Value of max suit.',
                       low = "#A50026", mid= "#FFFFBF", high = "#006837",
                       na.value="azure")+
  labs(x="Bioclimatic values",y="Species")+
  theme_bw()+
  theme(panel.grid = element_line(size=0.2),
        axis.text.x = element_text(size=15,vjust = .5,angle = 90),
        axis.text.y = element_text(size=15),
        axis.ticks = element_blank(),
        axis.title = element_text(size=15))
ggsave('output/enms/res_curves/max_bio.png',
       width = 9,height = 9,dpi=600)

# Table with bioclim names
bio_names <- data.frame(var = bio,
                        var_names = c('Annual Mean Temperature','Temperature Seasonality','Max Temperature of Warmest Month',
                                      'Min Temperature of Coldest Month','Temperature Annual Range (BIO5-BIO6)',
                                      'Mean Temperature of Wettest Quarter','Mean Temperature of Driest Quarter',
                                      'Mean Temperature of Warmest Quarter','Mean Temperature of Coldest Quarter',
                                      'Annual Precipitation','Precipitation of Wettest Month','Precipitation of Driest Month',
                                      'Precipitation Seasonality (Coefficient of Variation)','Precipitation of Wettest Quarter',
                                      'Precipitation of Driest Quarter','Precipitation of Warmest Quarter','Precipitation of Coldest Quarter'))

# Checking curves data
curves_info <- slope_max %>% group_by(var) %>%
  summarize(n = 19 - sum(is.na(slope)),
            sd_slope = sd(na.omit(slope)),
            sd_max = sd(na.omit(max_bio)),
            max_slope = max(na.omit(slope)),
            min_slope = min(na.omit(slope)),
            max_maxbio = max(na.omit(max_bio)),
            min_maxbio = min(na.omit(max_bio))) %>% 
  arrange(desc(n)) %>% left_join(bio_names)
saveRDS(curves_info,'RData/enms/curves_info.rds')

## Save data on response curves
## Long-pivotting to create final dataset with one column named `stat` to store both slope and maxbio
## in one column
curves_stats <- slope_max %>% pivot_longer(c('slope','max_bio'),names_to = 'stat',values_to = 'values')
saveRDS(curves_stats,'RData/predictors/resp_curves_stats.rds')