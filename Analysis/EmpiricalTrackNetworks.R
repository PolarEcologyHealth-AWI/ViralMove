### Empirical tracks

library(sf); sf::sf_use_s2(FALSE)
library(geosphere)
library(tidyverse)
source("Analysis/OptimSPD/R/OptimSPD.R", echo=FALSE)

load("Data/Map/eaafMap.rda")
load("breedTab.rda")
load("mudflatTab.rda")

## species
spParms <- setNames(lapply(c(250, 144, 110, 90, 50, 20), sizeParams), 
                    c("Godwit", "GreatKnot", "RedKnot", "Sanderling", "CurlewSandpiper", "RedNeckedStint"))
breedTab <- breedTab %>% filter(species%in%names(spParms)) %>% st_transform(4326)

###############################
### Arrival at breeding site ##
###############################

# arr <- do.call("c", lapply(1:nrow(breedTab), function(sp) {
#         sp <- breedTab[sp,]
#         tab <- read.csv(glue::glue("Data/Tracks/{sp$species}/{sp$id}_Grouped_movementSummary.csv")) %>%
#           as_tibble()
#         if(sp$species%in%c("Godwit", "RedKnot")) {
#           tab$Arrival <- as.POSIXct(tab$Arrival.50., tz = "GMT")
#           tab$Departure <- as.POSIXct(tab$Departure.50., tz = "GMT")
#         }
#         if(sp$species=="GreatKnot") {
#             tab$Arrival <- as.POSIXct(tab$Arrival, tz = "GMT")
#             tab$Departure <- as.POSIXct(tab$Departure, tz = "GMT")
#             tab <- tab %>% rename(Type = type)
#         }
#         if(sp$species=="Sanderling") {
#           tab$Arrival <- as.POSIXct(tab$Arrival, tz = "GMT")
#           tab$Departure <- as.POSIXct(tab$Departure, tz = "GMT")
#           tab <- tab %>%  rename(Type = type)
#         }
#         if(sp$species=="CurlewSandpiper") {
#           tab$Arrival <- as.POSIXct(tab$StartTime, tz = "GMT")
#           tab$Departure <- as.POSIXct(tab$EndTime, tz = "GMT")
#         }
#         if(sp$species=="RedNeckedStint") {
#           tab$Arrival <- as.POSIXct(tab$StartTime, tz = "GMT")
#           tab$Departure <- as.POSIXct(tab$EndTime, tz = "GMT")
#         }
#         if(any(tab$Type==2)) {
#           tab[min(which(tab$Type==2)),] %>% pull(Arrival)
#         } else (tab[max(which(tab$Type==1)),] %>% pull(Departure))+ 24*60*60
# }))
# dist <- do.call("c", lapply(1:nrow(breedTab), function(sp) {
#   spT <- breedTab[sp,]
#   tab <- read.csv(glue::glue("Data/Tracks/{spT$species}/{spT$id}_Grouped_movementSummary.csv")) %>%
#     as_tibble()
#   if(spT$species%in%c("GreatKnot", "Sanderling")) {
#     tab$lon = tab$Lon
#     tab$lat = tab$Lat
#     tab <- tab %>% rename(Type = type)
#   }
#   if(spT$species%in%c("Godwit", "CurlewSandpiper", "RedNeckedStint", "RedKnot")) {
#     tab$lon = tab$Lon.50.
#     tab$lat = tab$Lat.50.
#   }
#   if(any(tab$Type==2)) {
#     as.numeric(sum((tab[1:min(which(tab$Type==2)),] %>%
#              st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326) %>% st_distance())[cbind(1:(min(which(tab$Type==2)-1)), 2:min(which(tab$Type==2)))])/1000)
#   } else as.numeric(sum((tab[1:max(which(tab$Type==1)),] %>% select(lon, lat) %>% rbind(spT$geometry %>% st_coordinates() %>% as_tibble() %>% setNames(c("lon", "lat"))) %>%
#                            st_as_sf(coords = c("lon", "lat")) %>% st_set_crs(4326) %>% st_distance())[cbind(1:(max(which(tab$Type==1)-1)), 2:max(which(tab$Type==1)))])/1000)
# }))

# breedTab <- breedTab %>% bind_cols(tibble(arr_breed = as.POSIXct(arr, tz = "GMT"))) %>% relocate(arr_breed, .after = arr) %>%
#   mutate(migDur = as.numeric(arr_breed-dep)) %>% relocate(migDur, .after = arr_breed)
# save(breedTab, file = "breedTab.rda")
load("breedTab.rda")

breedTab %>% st_drop_geometry() %>% group_by(species) %>% summarise_at(c("migDur", "dist"), list(median, min, max))

########################
### Species Networks ###
########################

spNetMat <- lapply(breedTab %>% group_split(species), function(sp) {
  # sp <- (breedTab %>% group_split(species))[[6]]
  Reduce("+", parallel::mclapply(1:nrow(sp), function(i) {
    tab <- read.csv(glue::glue("Data/Tracks/{sp$species[i]}/{sp$id[i]}_Grouped_movementSummary.csv")) %>%
      as_tibble()
    if(sp$species[i]%in%c("Godwit", "RedKnot")) {
      tab$Arrival.50. <- as.POSIXct(tab$Arrival.50.)
      tab$Departure.50. <- as.POSIXct(tab$Departure.50.)
      tab$Days <- abs(as.numeric(difftime(tab$Departure.50., tab$Arrival.50., units = "days")))
    }
    if(sp$species[i]=="GreatKnot") { 
      if(any(names(tab)%in%"days")) {
        tab <- tab %>% mutate(Days = ifelse(is.na(days), 0, days)) %>%
          rename(Type = type, Lon.50. = Lon, Lat.50. = Lat)
      } else {
        tab$Arrival <- as.POSIXct(tab$Arrival)
        tab$Departure <- as.POSIXct(tab$Departure)
        tab$Days <- abs(as.numeric(difftime(tab$Arrival, tab$Departure, units = "days")))
        tab <- tab %>% mutate(Days = ifelse(is.na(Days), 0, Days)) %>%  rename(Type = type, Lon.50. = Lon, Lat.50. = Lat)
      }
    }
    if(sp$species[i]=="Sanderling") {
      tab$Arrival <- as.POSIXct(tab$Arrival)
      tab$Departure <- as.POSIXct(tab$Departure)
      tab$Days <- abs(as.numeric(difftime(tab$Arrival, tab$Departure, units = "days")))
      tab <- tab %>% mutate(Days = ifelse(is.na(Days), 0, Days)) %>%  rename(Type = type, Lon.50. = Lon, Lat.50. = Lat)
    }
    
    sum <- bind_rows(tibble(days = 0, X = sp[i,] %>% pull("lon_start"), Y = sp[i,] %>% pull("lat_start")),
                     tab %>% filter(Type==1 & Days>0) %>% select(Days, Lon.50., Lat.50.) %>% setNames(c("days", "X", "Y")),
                     cbind(days = 0, st_coordinates(sp[i,])) %>% as_tibble()) %>% st_as_sf(coords = c("X", "Y")) %>% st_set_crs(4326) %>% st_transform(st_crs(eaafMap$map)) %>% rowwise() %>%
      mutate(site = which.min(st_distance(geometry, mudflatTab %>% st_centroid() %>% suppressWarnings()))) %>% st_drop_geometry() %>% select(site, days)
    
    trans <- matrix(0, nrow = nrow(mudflatTab), ncol = nrow(mudflatTab))
    trans[cbind(sum[-nrow(sum),1], sum[-1,1]) %>% as.matrix()] <- 1
    diag(trans)[sum$site] <- sum$days
    
    trans  
  }, mc.cores = parallel::detectCores()-1))
  
})
pr       <- mudflatTab %>% st_centroid() %>% st_coordinates() %>% suppressWarnings() 

spCols   <- c("chartreuse4", "darkblue", "cadetblue", "brown3", "yellow2", "darkgoldenrod2")
opar <- par(mfrow = c(2,3), mar = c(.5,.5,.5,.5))
for(i in c(2,3,4,1,6,5)) {
  plot(eaafMap$grid, col = "grey80", lty = 3)
  plot(eaafMap$map, col = "grey90", border = "grey60", add = T)
  plot(eaafMap$bbox, add = T)
  
  sites <- diag(spNetMat[[i]]) 
  trans <- spNetMat[[i]]; diag(trans) <- 0
  
  transT <- cbind(rep(1:length(sites), length(sites)), rep(1:length(sites), each = length(sites)), c(trans))
  transT <- transT[transT[,3]>0, ]
  segments(pr[transT[,1],1], pr[transT[,1],2], pr[transT[,2],1], pr[transT[,2],2], lwd = approx(c(0, 18), c(1, 7), transT[,3])$y, col = adjustcolor("grey20", alpha.f = 0.7))
  
  points(pr, pch = 21, cex = approx(range(sites), c(0, 8), sites)$y, bg = adjustcolor(spCols[i], alpha.f = 0.7))
}      
par(opar)


#################################
### Species Relative Site Use ###
#################################

spStrategy <- lapply(breedTab %>% group_split(species), function(sp) {
  do.call("rbind", parallel::mclapply(1:nrow(sp), function(i) {
    tab <- read.csv(glue::glue("Data/Tracks/{sp$species[i]}/{sp$id[i]}_Grouped_movementSummary.csv")) %>%
      as_tibble() 
    if(sp$species[i]%in%c("Godwit", "RedKnot")) {
      tab$Arrival.50. <- as.POSIXct(tab$Arrival.50.)
      tab$Departure.50. <- as.POSIXct(tab$Departure.50.)
      tab$Days <- abs(as.numeric(difftime(tab$Departure.50., tab$Arrival.50., units = "days")))
    }
    if(sp$species[i]=="GreatKnot") { 
      if(any(names(tab)%in%"days")) {
        tab <- tab %>% mutate(Days = ifelse(is.na(days), 0, days)) %>%
          rename(Type = type, Lon.50. = Lon, Lat.50. = Lat)
      } else {
        tab$Arrival <- as.POSIXct(tab$Arrival)
        tab$Departure <- as.POSIXct(tab$Departure)
        tab$Days <- abs(as.numeric(difftime(tab$Arrival, tab$Departure, units = "days")))
        tab <- tab %>% mutate(Days = ifelse(is.na(Days), 0, Days)) %>%  rename(Type = type, Lon.50. = Lon, Lat.50. = Lat)
      }
    }
    if(sp$species[i]=="Sanderling") {
      tab$Arrival <- as.POSIXct(tab$Arrival)
      tab$Departure <- as.POSIXct(tab$Departure)
      tab$Days <- abs(as.numeric(difftime(tab$Arrival, tab$Departure, units = "days")))
      tab <- tab %>% mutate(Days = ifelse(is.na(Days), 0, Days)) %>%  rename(Type = type, Lon.50. = Lon, Lat.50. = Lat)
    }
    tab %>% filter(Type == 1 & Days > 0) %>% select(Days, Lon.50., Lat.50.) %>% 
      mutate(species = sp$species[i], id = sp$id[i]) %>% arrange(desc(Days)) %>%
      mutate(siteInd = suppressWarnings(do.call("c",lapply(split(., 1:nrow(.)), function(x) {
        st_as_sf(x, coords = c("Lon.50.", "Lat.50."), crs = 4326) %>% st_transform(st_crs(mudflatTab)) %>%
          st_distance(mudflatTab %>% st_centroid()) %>% which.min()
      })))) %>% select(-Lon.50., -Lat.50.) %>% relocate(siteInd, .after = species) %>%
      mutate(perc  = round((Days/sum(Days))*100,1)) %>% rownames_to_column(var = "site") %>% mutate(site = as.numeric(site)) %>% select(-Days)
  }))
})

opar <- par(mfrow = c(2,3), mar = c(.5,.5,.5,.5))
for(i in c(2,3,4,1,6,5)) {
  plot(spStrategy[[i]][,c("site", "perc")], ylim = c(0, 100), xlim = c(1, 8), type = "n")
  spStrategy[[i]] %>% group_split(id) %>% lapply(function(x) lines(x$site, x$perc, type = "o", pch = 16, cex = 0.6, 
                                                                   col = adjustcolor(spCols[i], alpha.f = 0.2))) %>% invisible()
  points(spStrategy[[i]] %>% group_by(site) %>% summarise(med_perc = median(perc, na.rm = T)),
         type = "o", pch = 21, col = spCols[i], lwd = 1.5)
  mx <- max(spStrategy[[i]]$site)
}      
par(opar)


### Results summary
## Yellow Sea Region
centrYS <- st_point(c(122.160166, 36.302840)) %>% st_sfc(crs = 4326) %>% st_buffer(12)
plot(eaafMap$map)
plot(centrYS %>% st_transform(st_crs(eaafMap$map)), add = T, col = "orange")

centrSA <- st_point(c(138.476402, -11.909262)) %>% st_sfc(crs = 4326) %>% st_buffer(12)
plot(eaafMap$map)
plot(centrSA %>% st_transform(st_crs(eaafMap$map)), add = T, col = "orange")


mudflatTab %>% st_area() %>% median()*1e-6

## Great Knot
i = 3
inYS <- which(st_intersects(mudflatTab, centrYS %>% st_transform(st_crs(mudflatTab)), sparse = FALSE))
  
spStrategy[[i]] %>%  mutate(inYS = siteInd %in% inYS) %>% group_split(id) %>%
  sapply(function(x) sum(x$perc[x$inYS])) %>% median()


## Red knot
i = 4
inSA <- which(st_intersects(mudflatTab, centrSA %>% st_transform(st_crs(mudflatTab)), sparse = FALSE))

spStrategy[[i]] %>%  mutate(inSA = siteInd %in% inSA) %>% group_split(id) %>%
  sapply(function(x) sum(x$perc[x$inSA])) %>% median()
###



empTrackList <- lapply(c(2,3,4,1,6,5), function(x) {
  sites <- diag(spNetMat[[x]]) 
  trans <- spNetMat[[x]]; diag(trans) <- 0
  
  transT   <- cbind(rep(1:length(sites), length(sites)), rep(1:length(sites), each = length(sites)), c(trans))
  transT   <- transT[transT[,3]>0, ]
  list(lapply(1:nrow(transT), function(t) st_linestring(pr[transT[t,1:2],])) %>% st_sfc() %>% st_set_crs(st_crs(eaafMap$map)) %>%
    st_sf() %>% mutate(trans = transT[,3]) %>% select(trans, geometry),
    st_as_sf(pr %>% as.data.frame(), coords = c("X", "Y")) %>% mutate(ts = sites) %>% st_set_crs(st_crs(eaafMap$map)) %>%
    select(ts, geometry) %>% filter(ts>0),
    spStrategy[[x]])
})

save(empTrackList, file = "Results/empTrackList.rda")
