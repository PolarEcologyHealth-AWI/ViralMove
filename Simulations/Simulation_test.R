##########################################################
#### Titel: 
#### Author: Marlena Lohse, Simeon Lisovski
#### Date: 05.02.2024
##########################################################

#### 1. Data/Result folder
{
  if(any(Sys.info()=="PGeo07m013.dmawi.de")) {
    data_folder <- "/Volumes/projects/bioing/data/PathogenTransport/ViralMove_data"
  }
  if(any(Sys.info()=="linux4")) {
    data_folder <- "/bioing/data/PathogenTransport/ViralMove_data"
  }
}

#### 2. Packages and functions
{
  library(Rcpp)
  library(gitcreds)
  library(RcppArmadillo)
  library(RcppProgress)
  library(sf); sf::sf_use_s2(FALSE)
  library(geosphere) ## To be removed
  library(tidyr)
  library(dplyr)
  library(tibble)
  library(lubridate)
  
  source("OptimSDP/R/OptimSDP.R", echo=FALSE)
  sourceCpp('OptimSDP/src/OptimSDP_infection.cpp', showOutput = FALSE)
}

#### 3. Load site parameters and empirical 
{
  load(glue::glue("{data_folder}/Data/Map/eaafMap.rda"))
  load(glue::glue("{data_folder}/Data/breedTab_revision.rda"))
  load(glue::glue("{data_folder}/Data/mudflatTab.rda"))
  load(glue::glue("{data_folder}/Data/tempTab_revision.rda"))
  load(glue::glue("{data_folder}/Results/empTrackListNW.rda"))
  load(glue::glue("{data_folder}/Data/phenology_tracks.RData"))
  }

#### 4. Species
{
  sps <- c("Godwit", "RedKnot", "CurlewSandpiper", "RedNeckedStint")
  empTrackList <- empTrackListNW
  
  ## species
  spParms <- setNames(lapply(c(250, 105, 55, 25), sizeParams), sps)
  breedTab <- breedTab %>% filter(species%in%sps) %>% st_transform(4326) %>%
    dplyr::select(-dep, -arr, -arr_breed) %>%
    left_join(phen %>% dplyr::select(-Species), by = join_by(id==ID))
  
  spCols   <- c("darkblue", "brown3", "darkgoldenrod2", "yellow2")
}


#### 5. Simulation (past, present)
  
sp <- 'Godwit'
subBreedTab <- breedTab %>% filter(species == sp) %>% slice(5)

ind <- 1
      
      #######################
      ### site parameters ###
      #######################
      StartEnd_cell <- tibble(lon = c(subBreedTab[ind,] %>% pull("lon_start"), st_coordinates(subBreedTab)[ind,1]),
                                              lat = c(subBreedTab[ind,] %>% pull("lat_start"), st_coordinates(subBreedTab)[ind,2])) %>%
          st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326) %>% rowwise() %>%
          mutate(index = which.min(st_distance(geometry, mudflatTab %>% st_transform(4326)))) 
      
      
      
      size     <- (mudflatTab %>% st_centroid() %>% st_transform(4326) %>%
                     mutate(areaHist = ifelse(is.na(histArea_inner), currArea_outer, histArea_outer),
                            areaCurr = currArea_outer,
                            areaMang = mangArea_outer,
                            lake     = ifelse(is.na(lake_area), FALSE, TRUE),
                            pmax     = pmax(areaCurr, areaHist)) %>% 
                     dplyr::select(pmax, areaHist, areaCurr, areaMang, lake) %>%
                     rownames_to_column(var = "index") %>% mutate(index = as.integer(index))) %>%
        filter((areaHist>0 | areaCurr>0 | areaMang>0 | lake)) %>%
        filter(st_coordinates(.)[,1] > 104 | st_coordinates(.)[,1] < -150) %>% 
        arrange(st_distance(geometry, StartEnd_cell[1,])) %>%
        bind_rows(StartEnd_cell[2,] %>% mutate(areaHist = NA, areaCurr = NA, areaMang = NA, lake = FALSE) %>% dplyr::select(names(.))) %>%
        relocate(geometry, .after = last_col()) %>%
        suppressWarnings() 
      
      distM  <- st_distance(size, by_element = F)/1000
      bearM  <- abs(distm(st_coordinates(size), fun = bearing))
      
      ### Revision
      reward <- c(subBreedTab[ind,] %>% pull(start_hist), subBreedTab[ind,] %>% pull(start_curr),  subBreedTab[ind,] %>% pull(start_future))
      
      #######################
      ### sdp Objects #######
      #######################
      sdpObjects <- makeSDPobjects(
        list(
          direction = 1,
          species   = subBreedTab$species[ind],
          spParms   = spParms[[which(names(spParms)==sp)]],
          minT      = as.numeric(format(as.POSIXct("2012-01-01"), "%j")),
          maxT      = subBreedTab[ind,] %>% pull(start_hist) + 15,
          MaxX      = 100,
          B0        = 3,
          w         = 0.028, ## 0.028 for nm
          xc        = 10, ## 10 for nm
          WindAssist = 0,
          WindProb   = 1,
          ZStdNorm   = c(-2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,  1.5,  2.0,  2.5),
          PStdNorm   = c(0.0092, 0.0279, 0.0655, 0.1210, 0.1747, 0.2034, 0.1747, 0.1210, 0.0655, 0.0279, 0.0092),
          reward     = reward,
          decError   = 250000,
          pred    = c(1e-3, 1e-4, 1e-6, 2, 2),
          sites   = size,
          inScale = c(0.25, 1),
          penalty = c(0,0.0012, 0),
          latF    = 1,
          dist  = distM,
          bear  = bearM,
          tTab  = tempTab[size$index,,] 
        ))
      
        model   <- bwdIteration(sdpObjects[[1]])

        image(list(x = 1:dim(model@Results$FitnessMatrix)[1],
                   y = 1:dim(model@Results$FitnessMatrix)[2],
                   z = model@Results$FitnessMatrix[,,101,1]))
        
        image(list(x = 1:dim(model@Results$DecisionMatrix)[1],
                   y = 1:dim(model@Results$DecisionMatrix)[2],
                   z = model@Results$DecisionMatrix[,,101,1,1]))
        
        
        simu    <- tryCatch(fwdSimulation(model, 100, start_t = 1, start_site = 1, start_x = c(30,50)), error = function(e) NULL)
        condProfile(simu, model)
        matplot(model@Results$FitnessMatrix[,400,], type = 'o', col = 'grey80', pch = 16)
        simNetwork(simu, model, crds_ind = mudflatTab %>% st_centroid() %>% st_coordinates() %>% suppressWarnings(), plot = T)
