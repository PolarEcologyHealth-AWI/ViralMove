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
  
  source("Analysis/OptimSDP/R/OptimSDP.R", echo=FALSE)
  sourceCpp('Analysis/OptimSDP/src/OptimSDP.cpp', showOutput = FALSE)
}

#### 3. Load site parameters and empirical tracks
{
  load(glue::glue("{data_folder}/Data/Map/eaafMap.rda"))
  load(glue::glue("{data_folder}/Data/breedTab_revision.rda"))
  load(glue::glue("{data_folder}/Data/mudflatTab.rda"))
  load(glue::glue("{data_folder}/Data/tempTab_revision.rda"))
  load(glue::glue("{data_folder}/Results/empTrackList.rda"))  
  load(glue::glue("{data_folder}/Data/phenology_tracks.RData"))
}

#### 4. Species
{
  sps <- c("Godwit", "RedKnot", "CurlewSandpiper", "RedNeckedStint")
  empTrackList <- empTrackList[sapply(empTrackList, function(x) x[[3]]$species[1]) %in% sps]
  
  ## species
  spParms <- setNames(lapply(c(250, 105, 55, 25), sizeParams), sps)
  breedTab <- breedTab %>% filter(species%in%sps) %>% st_transform(4326) %>%
    dplyr::select(-dep, -arr, -arr_breed) %>%
    left_join(phen %>% dplyr::select(-Species), by = join_by(id==ID))
  
  spCols   <- c("darkblue", "brown3", "darkgoldenrod2", "yellow2")
}

penSeq <- seq(0, 0.002, length = 6)
intSeq <- matrix(c(0.25, 0.5, 0.75, 1, 1.25, 1.5), nrow = 3)
latSeq <- seq(1, 2, length = 4)

sim_list <- expand.grid(pen = 1:length(penSeq),
                   int = 1:nrow(intSeq),
                   lat = 1:length(latSeq))

for(run in 1:nrow(sim_list)) {

  #### 5. Simulation (past, present)
  dir <- 2
  
  if(dir == 2){breedTab <- breedTab %>% filter (!is.na(Arrival))
    } else {}
  
  allSpSim <- lapply(names(spParms), function(sp) {
    
    subBreedTab <- breedTab %>% filter(species == sp)
    
    indSim <- parallel::mclapply(1:nrow(subBreedTab), function(ind) {
      
      #######################
      ### site parameters ###
      #######################
      StartEnd_cell <- if (dir == 1) { tibble(lon = c(subBreedTab[ind,] %>% pull("lon_start"), st_coordinates(subBreedTab)[ind,1]),
                                              lat = c(subBreedTab[ind,] %>% pull("lat_start"), st_coordinates(subBreedTab)[ind,2])) %>%
          st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326) %>% rowwise() %>%
          mutate(index = which.min(st_distance(geometry, mudflatTab %>% st_transform(4326)))) 
      } else { tibble(lon = c(st_coordinates(subBreedTab)[ind,1], subBreedTab[ind,] %>% pull("lon_start")),
                      lat = c(st_coordinates(subBreedTab)[ind,2], subBreedTab[ind,] %>% pull("lat_start"))) %>%
          st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326) %>% rowwise() %>%
          mutate(index = which.min(st_distance(geometry, mudflatTab %>% st_transform(4326)))) 
      }
      
      
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
          species = subBreedTab$species[ind],
          spParms = spParms[[which(names(spParms)==sp)]],
          minT    =  if (dir ==1) {as.numeric(format(as.POSIXct("2012-01-01"), "%j"))
          } else {yday(subBreedTab[ind,] %>% pull(Departure_breed)) - 15},
          maxT    = if (dir ==1) {subBreedTab[ind,] %>% pull(start_hist) + 15
          }else {yday(subBreedTab[ind,] %>% pull(Arrival)) + 15},
          MaxX    = 100,
          B0      = 3,
          w       = 0.028,
          xc      = 10,
          WindAssist = 0,
          WindProb   = 1,
          ZStdNorm   = c(-2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,  1.5,  2.0,  2.5),
          PStdNorm   = c(0.0092, 0.0279, 0.0655, 0.1210, 0.1747, 0.2034, 0.1747, 0.1210, 0.0655, 0.0279, 0.0092),
          reward     = reward,
          decError   = 250000,
          pred    = c(1e-3, 1e-4, 1e-6, 2, 2),
          sites   = size,
          inScale = intSeq[sim_list[run,2],],
          penalty = c(0,penSeq[sim_list[run,1]], 0),
          latF    = latSeq[sim_list[run,3]],
          dist  = distM,
          bear  = bearM,
          tTab  = tempTab[size$index,,] 
        ))
      
      parallel::mclapply(1:2, function(x) {
        model   <- bwdIteration(sdpObjects[[x]])
        simu    <- tryCatch(fwdSimulation(model, 100, start_t = 1, start_site = 1, start_x = c(30,50)), error = function(e) NULL)
        simNetwork(simu, model, crds_ind = mudflatTab %>% st_centroid() %>% st_coordinates() %>% suppressWarnings(), plot = F)
      }, mc.cores = 2)
      
    }, mc.cores = 8)
    
    list(
      Reduce("+", lapply(indSim, function(ind) ind[[1]][[1]])[sapply(indSim, function(ind) class(ind[[1]][[1]])[1])=="matrix"]),
      Reduce("+", lapply(indSim, function(ind) ind[[2]][[1]])[sapply(indSim, function(ind) class(ind[[2]][[1]])[1])=="matrix"]),
      do.call("rbind", lapply(indSim, function(ind) tryCatch(ind[[1]][[2]], error = function(e) NULL))),
      do.call("rbind", lapply(indSim, function(ind) tryCatch(ind[[2]][[2]], error = function(e) NULL))),
      do.call("rbind", lapply(indSim, function(ind) tryCatch(ind[[1]][[3]], error = function(e) NULL))),
      do.call("rbind", lapply(indSim, function(ind) tryCatch(ind[[2]][[3]], error = function(e) NULL))),    do.call("rbind", lapply(1:length(indSim), function(ind) tryCatch(as_tibble(indSim[[ind]][[1]][[4]]) %>% setNames(c("id", "time", "site", "x")) %>% mutate(id = glue::glue("{ind}_{id}")), error = function(e) NULL))),
      do.call("rbind", lapply(1:length(indSim), function(ind) tryCatch(as_tibble(indSim[[ind]][[2]][[4]]) %>% setNames(c("id", "time", "site", "x")) %>% mutate(id = glue::glue("{ind}_{id}")), error = function(e) NULL))),
      do.call("rbind", lapply(1:length(indSim), function(ind) rbind(tryCatch(tibble(t = 1, f = indSim[[ind]][[1]][[5]]), error = function(e) NULL), 
                                                                    tryCatch(tibble(t = 2, f = indSim[[ind]][[2]][[5]]), error = function(e) NULL)))))
    
  })
  save(allSpSim, file = glue::glue("{data_folder}/Results/Southward/allSpSim_{penSeq[sim_list[run,1]]}_int_{intSeq[sim_list[run,2],2]}_lat_{latSeq[sim_list[run,3]]}_SW.rda"))
  

  #### 6. Diagnostic plot(s)
  
  # source("Analysis/Figure_Script.R", echo= FALSE)
  # load(glue::glue("{data_folder}/Results/Southward/allSpSim_{penSeq[sim_list[run,1]]}_int_{intSeq[sim_list[run,2],2]}_lat_{latSeq[sim_list[run,3]]}_SW.rda"))
  # plotMigrationData (allSpSim, empTrackList, spCols, breedTab, eaafMap, mudflatTab, spParms, 
  #                    glue::glue("fig_pen_{penSeq[sim_list[run,1]]}_int_{intSeq[sim_list[run,2],2]}_lat_{latSeq[sim_list[run,3]]}.pdf")) 

}
