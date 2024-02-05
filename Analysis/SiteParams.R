library(sf); sf::sf_use_s2(FALSE)
library(raster)
library(tidyverse)
library(stars)

## Basic flyway map
proj <- "+proj=laea +lat_0=15 +lon_0=162 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"
# xlim <- c(60, 210)
# ylim <- c(-55, 82) 454,55
# 
# poly <- st_as_sf(as(extent(c(xlim, ylim)), "SpatialPolygons")) %>% st_set_crs(4326) %>% st_transform(CRS(proj)) %>% st_bbox() %>%
#   st_as_sfc() %>% st_set_crs(proj)
# 
# map <- read_sf("Data/Map/ne_50m_land/ne_50m_land.shp") %>%
#   st_difference(st_as_sf(as(extent(c(-10,10,-90,90)), "SpatialPolygons")) %>% st_set_crs(4326)) %>%
#   st_shift_longitude() %>% st_geometry() %>% st_transform(crs = proj) %>% st_buffer(0) %>% st_intersection(poly)
# 
# grid <- read_sf("Data/Map/ne_10m_graticules_20/ne_10m_graticules_20.shp", ) %>%
#   st_shift_longitude() %>% st_geometry() %>% st_transform(crs = CRS(proj)) %>% st_intersection(poly)
# 
# eaafMap <- list(map = map, grid = grid, bbox = poly)
# save(eaafMap, file = "Data/Map/eaafMap.rda")
load("Data/Map/eaafMap.rda")

opar <- par(mar = c(0.1,0.1,0.1,0.1))
plot(eaafMap$grid, col = "grey80", lty = 3)
plot(eaafMap$map, col = "grey70", border = "grey70", add = T)
# plot(st_cast(eaafMap$map, "LINESTRING") %>% st_buffer(20000) %>% st_intersection(eaafMap$bbox %>% st_buffer(-30000)), col = "orange", border = NA, add = T)
plot(eaafMap$bbox, add = T, border = "grey60")
par(opar)

hex <- st_make_grid(eaafMap$bbox, cellsize = 400000, square = FALSE) %>% st_set_crs(proj)

plot(hex %>% st_intersection(eaafMap$bbox), border = "grey80", add = T)
plot(eaafMap$bbox, add = T, border = "grey60")
plot(hex, add = T)
# text(hex%>%st_centroid()%>%st_coordinates(), as.character(1:length(hex)), cex = 0.6)
# save(hex, file = "hex.rda")
# load("hex.rda")

# ggplot() +
#   geom_sf(data = eaafMap$map) +
#   geom_sf(data = hex, fill = "transparent") +
#   geom_sf(data = eaafMap$bbox, fill = "transparent") +
#   # geom_text(data = st_coordinates(hex %>% st_centroid())[,1:2] %>% as_tibble() %>% rownames_to_column(var = "id"),
#   #           mapping = aes(x = X, y = Y, label = id), size = 1.9) +
#   geom_sf(data = hex[303], fill = "darkred", alpha = 0.7) +
#   geom_sf(data = hex[355], fill = "orange", alpha = 0.7) +
#   geom_sf(data = hex[574], fill = "darkblue", alpha = 0.7) +
#   theme_light()




###############################################
## Intertidal mudflats                       ##
### Historical mudflat extent (Yellow Sea)   ##
###############################################

hist <- read_sf("Data/Murray_etal_2014/East_Asia_Tidal_Flat_1952_1964.shp") %>%     
  st_transform(st_crs(eaafMap$map)) %>% st_union()

areaHist <- hex %>% st_as_sf() %>% rownames_to_column(var = "id") %>% mutate(id = as.numeric(id)) %>% mutate() %>% st_intersection(hist %>% st_buffer(0)) %>%
  mutate(area = as.numeric(st_area(.))/1000) %>% st_drop_geometry() %>% as_tibble() %>% full_join(
    st_centroid(hex) %>% st_buffer(400000) %>% st_as_sf() %>% rownames_to_column(var = "id") %>% mutate(id = as.numeric(id)) %>% st_intersection(hist %>% st_buffer(0)) %>%
      mutate(area_bfr = as.numeric(st_area(.))/1000) %>% st_drop_geometry() %>% as_tibble()
  ) %>% full_join(hex %>% st_as_sf() %>% rownames_to_column(var = "id") %>% mutate(id = as.numeric(id)) %>% rename(geometry = x)) %>% arrange(id) %>% st_as_sf()

# plot(eaafMap$map, col = NA, border = "grey60", add = T)


###############################################
### Current mudflat extent                   ##
###############################################

fls  <- list.files("Data/Murray_etal_2018/Extracted polygons", pattern = "*extracted_union.shp", full.names = T)

# hex_example <- hex[303]
# bf <- st_centroid(hex_example) %>% st_buffer(400000) %>% st_as_sf()
# 
#   hist_p <- hist %>% st_intersection(bf)
#   
#   pols <- lapply(fls, function(f) {
#     mud <- read_sf(f) %>% st_union() %>% st_transform(st_crs(hex)) %>% st_intersection(bf)
#     if(length(mud)>0) mud else NULL
#   })
#     
#   
#   ggplot() +
#     geom_sf(data = bf) +
#     geom_sf(data = pols %>% Reduce("rbind",.) %>% st_as_sfc() %>% st_set_crs(st_crs(hex)), fill = "darkred", color = "transparent") +
#     geom_sf(data = hex_example, fill = "transparent") +
#     theme_light()

# mudTabs <- do.call("cbind", parallel::mclapply(fls, function(f) {
#   mud <- read_sf(f) %>% st_union() %>% st_transform(st_crs(hex))
#   areaHist %>% select(id) %>% st_intersection(mud) %>%
#     mutate(innter_area = as.numeric(st_area(.))/1000) %>% st_drop_geometry() %>% as_tibble() %>% full_join(
#       st_centroid(areaHist) %>% select(id) %>% st_buffer(400000) %>% st_intersection(mud) %>%
#         mutate(full_area = as.numeric(st_area(.))/1000) %>% st_drop_geometry() %>% as_tibble()
#     ) %>% full_join(tibble(id = 1:length(hex))) %>% arrange(id) %>% select(innter_area, full_area)
# }, mc.cores = 3))
# 
# areaCurr <- mudTabs %>% setNames(paste0(rep(c("a", "b"), ncol(.)/2), seq_along(.))) %>%
#   mutate(currArea_inner = rowSums(select(., starts_with("a")), na.rm = T),
#          currArea_outer = rowSums(select(., starts_with("b")), na.rm = T))


###############################################
### All mudflat extent                       ##
###############################################

# mudflatTab <- areaHist %>% rename(histArea_inner = area, histArea_outer = area_bfr) %>%
#   bind_cols(areaCurr %>% select(c(currArea_inner, currArea_outer))) %>% select(histArea_inner, histArea_outer, currArea_inner, currArea_outer)

# ggplot(eaafMap$map) + geom_sf() +
# geom_sf(mudflatTab, mapping = aes(geometry = geometry, fill = currArea_outer), col = NA)


###############################################
### Mangroves                                ##
###############################################

# mang <- read_sf("/Volumes/slisovski/RemoteSensedData/MangroveDistribution/01_Data/14_001_WCMC010_MangroveUSGS2011_v1_3.shp") %>% st_transform(proj) %>%
#   st_buffer(0) %>% st_intersection(st_bbox(mudflatTab) %>% st_as_sfc())
# 
# mangArea <- mudflatTab %>% rownames_to_column(var = "id") %>% mutate(id = as.numeric(id)) %>% select(id) %>%
#   st_intersection(mang %>% select(AREA_KM2)) %>% st_drop_geometry() %>% group_by(id) %>% summarise(mangArea_inner = sum(AREA_KM2)) %>%
#   full_join(mudflatTab %>% st_centroid() %>% rownames_to_column(var = "id") %>% mutate(id = as.numeric(id)) %>% select(id) %>% st_buffer(400000) %>%
#               st_intersection(mang %>% select(AREA_KM2)) %>% st_drop_geometry() %>% group_by(id) %>% summarise(mangArea_outer = sum(AREA_KM2))) %>%
#   full_join(tibble(id = 1:nrow(mudflatTab))) %>% arrange(id)
# 
# ### join with mudflat table
# mudflatTab <- mudflatTab %>% mutate(mangArea_inner = ifelse(is.na(mangArea$mangArea_inner), 0, mangArea$mangArea_inner),
#                                     mangArea_outer = ifelse(is.na(mangArea$mangArea_outer), 0, mangArea$mangArea_outer)) %>%
#   relocate(geometry, .after = last_col())
# 
# ggplot(eaafMap$map) + geom_sf() +
#   geom_sf(mudflatTab, mapping = aes(geometry = geometry, fill = mangArea$mangArea_outer), col = NA)


###############################################
### Lakes                                    ##
###############################################

# mudflatTab <- mudflatTab %>% mutate(lake_area = mudflatTab %>% rownames_to_column(var = "id") %>% mutate(id = as.numeric(id)) %>% select(id) %>% st_intersection(
#       read_sf("Data/Map/ne_10m_lakes/ne_10m_lakes.shp" )%>% st_transform(proj) %>%
#       st_buffer(0) %>% st_intersection(st_bbox(mudflatTab) %>% st_as_sfc()) %>% 
#       filter(st_coordinates(st_centroid(.) %>% st_transform(4326))[,2]<65 & 
#                (st_coordinates(st_centroid(.) %>% st_transform(4326))[,1]>100 |
#                   st_coordinates(st_centroid(.) %>% st_transform(4326))[,1]<0)) %>%
#       mutate(area = as.numeric(st_area(.))) %>% filter(area >= quantile(.$area, probs = 0.9)) %>% st_combine()) %>%
#   mutate(lake_area = as.numeric(st_area(.))/1000) %>% st_drop_geometry() %>% as_tibble() %>% 
#   full_join(tibble(id = 1:nrow(mudflatTab)))  %>% arrange(id) %>% pull(lake_area)) %>% relocate(geometry, .after = last_col())
# 
# save(mudflatTab, file = "mudflatTab.rda")
load("mudflatTab.rda")
# 
# ggplot(eaafMap$map) + geom_sf() +
#   geom_sf(mudflatTab, mapping = aes(geometry = geometry, fill = lake_area), col = NA)

###############################################
### Temperature (daily)                      ##
###############################################

# path <- "/Volumes/slisovski/RemoteSensedData/NCEP/Surface_temp/"
# 
# flsTemp <- tibble(path = list.files(path)) %>% mutate(year = as.numeric(unlist(lapply(strsplit(path, ".nc"), function(x) substring(x, 12, 15))))) %>%
#   filter(year%in%c(1950:2070) | year%in%c(2000:2020))
# 
#    ## template
#    tmp_r <- brick(paste0(path, flsTemp$path[1]))
#    tmp_i <- rasterize(mudflatTab$geometry %>% st_transform(4326) %>% st_shift_longitude() %>% as("Spatial"), tmp_r[[1]])
# 
# tempArray <- abind::abind(parallel::mclapply(flsTemp$path, function(x) {
# 
#   (brick(paste0(path, x)) %>% st_as_stars() %>% st_set_crs(4326) %>% st_as_sf() %>% st_centroid() %>%
#     st_transform(st_crs(hex)) %>%
#     mutate(cell = apply(st_intersects(., hex, sparse = FALSE), 1, function(t) ifelse(any(t), which(t), NA))) %>%
#     filter(!is.na(cell)) %>% st_drop_geometry() %>% group_split(cell) %>%
#     lapply(., function(z) apply(z[,-ncol(z)], 2, function(y) quantile(as.numeric(y), probs = 0.4, na.rm = T))) %>% Reduce("rbind", .))[,1:364]

  # dim(tt)
  #
  # tmp_r <- brick(paste0(path, flsTemp$path[1]))
  #
  # proj4string(tmp_r) <- "+proj=longlat +ellps=WGS84 +pm=-360 +datum=WGS84 +no_defs"
  # tt <- data.frame(lon = coordinates(tmp_r)[,1], lat = coordinates(tmp_r)[,2], temp = tmp_r[[1]][])
  # tt$lon[tt$lon>180] <-  - (360 - tt$lon[tt$lon>180])
  # rast <- raster::rasterize(tt[,1:2], raster(xmn = -179.9, xmx = 179.9, ymn = -90, ymx = 90, res = 2.5), field = tt$temp)
  # tempGr <- st_make_grid(eaafMap$bbox, cellsize = 150000) %>% st_transform("+proj=longlat")
  # extr   <- raster::extract(rast, (tempGr %>% st_centroid() %>% st_coordinates())[,1:2])
  # tempPol <- tempGr %>% st_as_sf() %>% mutate(temp = extr) %>% st_transform(proj)
  # # plot(eaafMap$grid, col = "grey80", lty = 3)
  # # plot(eaafMap$map, col = "grey70", border = "grey70", add = T)
  # # plot(tempPol, add = T, border = NA, lwd = 0.1)
  # # plot(eaafMap$bbox, add = T, border = "grey60")
  #
  # tmp_t <- tibble(index = tmp_i[]) %>% bind_cols(tmp_r[]) %>% filter(!is.na(index)) %>% pivot_longer(cols = starts_with("X")) %>%
  #   group_by(index, name) %>% summarise(temp = quantile(value, probs = 0.4, na.rm = T)) %>%
  #   pivot_wider(names_from = name, values_from = temp) %>% dplyr::select(-index) %>% as.matrix()
  # tmp_t[,1:364]

# }, mc.cores = 7), along = 3)

# save(tempArray, file = "tempArray_new.rda")
# load("tempArray_new.rda")

# ## daily values (two periods)
# tempTab <- abind::abind(parallel::mclapply(1:dim(tempArray)[1], function(x) {
# 
#   tmp01 <- tibble(y = rep(flsTemp$year, each = dim(tempArray)[2]), doi = rep(1:364, dim(tempArray)[3]), temp = c(tempArray[x,,]) - 273.15)
#   tmp02 <- tmp01 %>% bind_rows(tmp01 %>% mutate(doi = doi-364)) %>% bind_rows(tmp01 %>% mutate(doi = doi+364))
# 
#   # plot(tmp02$doi, tmp02$temp, pch = 16, cex = 0.25)
#   mod0 <- predict(loess(temp~doi, data = tmp02 %>% filter(y > 1960 & y < 1970), span = 0.2), newdata = data.frame(doi = 1:365))
#   mod1 <- predict(loess(temp~doi, data = tmp02 %>% filter(y > 2010 & y < 2020), span = 0.2), newdata = data.frame(doi = 1:365))
#   
#   # tab <-  tibble(doy = tmp02$doi, year = tmp02$y, loess = FALSE, temperature = tmp02$temp) %>%
#   #           mutate(per = if_else(year < 1970, "1960-1970", "2010-2020")) %>% dplyr::select(-year) %>%
#   #             bind_rows(tibble(doy = 1:365, loess =TRUE, temperature = mod0, per = "1960-1970") %>%
#   #                 bind_rows(tibble(doy = 1:365, loess =TRUE, temperature = mod1, per = "2010-2020")))
# 
#   # ggplot(tab) +
#   #   geom_point(tab %>%filter(!loess, doy %in% c(1:365)), mapping = aes(x = doy, y = temperature), size = 1) +
#   #   geom_line(tab  %>% filter(loess, doy %in% c(1:365)), mapping = aes(x = doy, y = temperature), linewidth = 1.4, col = "orange") +
#   #   # geom_point(mapping = aes(x = doy, y = temperature), size = 1) +
#   #   facet_wrap(~per) +
#   #   theme_light()
#   
#   
#   
#   # lines(1:365, mod0, lwd = 2, col = "orange")
#   # lines(1:365, mod1, lwd = 2, col = "magenta")
# 
#   array(c(as.numeric(mod0), as.numeric(mod1)), dim = c(1, 365, 2))
# 
# }, mc.cores = 8), along = 1)
# 
# tempTab <- abind::abind(tempTab, apply(tempTab, 1:2, function(x)
#   predict(lm(temp~years, data = tibble(years = c(1960, 2010), temp = c(x))), newdata = data.frame(years = 2060))),
#   along = 3)

# save(tempTab, file = "tempTab.rda")
load("tempTab.rda")


### CIMP4 AWI SPP2 Scenario
cimp  <- (read_stars("Data/CIMP6_SPP2/tas_day_AWI-CM-1-1-MR_ssp245_r1i1p1f1_gn_20600101-20651231_v20190529.nc") %>%
           st_set_crs(4326)) %>% suppressMessages() %>% suppressWarnings()
dates   <- st_get_dimension_values(cimp, 3)
breedID <- unique(apply(st_intersects(breedTab %>% st_transform(st_crs(mudflatTab)), mudflatTab, sparse = F), 1, which))

tempCimp <- matrix(nrow = nrow(mudflatTab), ncol = 365)

for(x in 1:nrow(mudflatTab)) {
  
  cat(sprintf('\n pixel %d if %d', x, nrow(mudflatTab))) 
  
  if((!is.na(sum(mudflatTab[x,] %>% st_drop_geometry(), na.rm = T)) & 
     sum(mudflatTab[x,] %>% st_drop_geometry(), na.rm = T) > 0 &
     all(is.na(tempCimp[x,]))) || x %in% breedID) {
  
  poly <- mudflatTab[x,] %>% st_geometry() %>% st_transform(st_crs(cimp)) %>% st_shift_longitude()
  rst  <- st_extract(cimp, poly) %>% suppressMessages() %>% suppressWarnings()
  arr  <- apply(rst[[1]], 1:2, function(y) c(y))
  
  # tab <- tibble(date = rep(dates, dim(arr)[2]*dim(arr)[3]), temp = c(unlist(arr)) - 273.15) %>% filter(!is.na(temp)) %>%
  #   mutate(doy = as.numeric(format(date, "%j"))) %>% dplyr::select(doy, temp) %>%
  #   group_by(doy) %>% summarise(temp = median(temp))
  
  tab <- tibble(date = dates, temp = c(unlist(arr)) - 273.15) %>% filter(!is.na(temp)) %>%
    mutate(doy = as.numeric(format(date, "%j"))) %>% dplyr::select(doy, temp) %>%
    group_by(doy) %>% summarise(temp = median(temp))
  tabMod    <- rbind(tab %>% mutate(doy = doy - 365), tab, tab %>% mutate(doy = doy + 365))
  loessTemp <- predict(loess(temp~doy, data = tabMod, span = 0.2), newdata = data.frame(doy = 1:365))
  tempCimp[x,] <- as.numeric(loessTemp)
  invisible(gc())
  }
}
# save(tempCimp, file = "tempCimp.rda")
# 
# indNA <- which(apply(tempCimp, 1, function(x) any(!is.na(x))))
# indPlot <- indNA[150]
# plot(tempTab[indPlot,,3], ylim = range(c(tempTab[indPlot,,3], tempCimp[indPlot,])))
# points(tempCimp[indPlot,], pch = 16, col = "orange")
# 
# tempTab[,,3] <- tempCimp  
# save(tempTab, file = "tempTab_revision.rda")

load("tempTab_revision.rda")

###############################################
### Snow melt timing                         ##
###############################################
# library(bbmle)
# gaussMLE <- function(day, size, prob, thresh) {
# 
#   tab <- data.frame(day = day, size = size, p = prob)
# 
#   gauss.curve <- function(parms, intv = 1) {
#     t <- seq(1, 366, intv)
#     parms <- as.list(parms)
#     fit1 <- 1 - exp(-((parms$a1 - t[1:(which(t==floor(parms$a1)))])/parms$a4)^parms$a5)
#     fit2 <- 1 - exp(-((t[which(t==floor(parms$a1)):length(t)]-parms$a1)/parms$a2)^parms$a3)
#     c(fit1, fit2[-1])
#   }
# 
#   gauss.loglik <- function(a1, a2, a3, a4, a5) {
#     fit <- gauss.curve(parms = list(a1=a1, a2=a2, a3=a3, a4=a4, a5=a5), 1)
#     fit <- ifelse(fit>0.999, 1-(1e-5), ifelse(fit<0.001, 1e-5, fit))
#     # cat(paste(c(a1, a2, a3, a4, a5), sep = "  "), "\r")
#     -sum(dbinom(x = round(tab[,3]*tab[,2],0), size = rep(100, length(fit)), prob = fit[day], log=TRUE), na.rm=T)
#   }
# 
#   mle <- suppressWarnings(mle2(gauss.loglik, method="L-BFGS-B",
#                                start=list(a1 = 225, a2 = 40,  a3 = 9,  a4 = 40, a5 = 9),
#                                lower=list(a1 = 50,  a2 = 5,   a3 = 0.5,a4=5,  a5 = 0.5),
#                                upper=list(a1 = 225, a2 =  Inf,  a3 =  Inf, a4 =  Inf, a5 =  Inf),
#   ))
# 
#   t <- seq(1, 366, 1)
#   fit <- gauss.curve(coef(mle), intv = 1)
# 
#   start <- t[min(which(fit<thresh))]
# 
#   list(result = data.frame(week = t, fit = fit), start = start, end = end)
# }
## end maximum likelihood function

geoTab <- readxl::read_excel("Output/geolocTab.xlsx") %>% mutate(migDur = as.numeric(difftime(arr, dep, units = "days")))

ggplot(geoTab %>% group_by(species) %>% summarise(migDur = median(migDur))) +
  geom_point(mapping = aes(x = as.factor(species), y = migDur))

# View(geoTab)
# 
# library(ncdf4)
# 
# projSnow <- "+proj=stere +lat_0=90 +lon_0=10"
# ncDat <- nc_open("~/Google Drive/My Drive/GeoDat/nhsce_v01r01_19661004_20220103.nc")
# t <- ncvar_get(ncDat, "time")
# z <- ncvar_get(ncDat, "snow_cover_extent")
# ylat <- ncvar_get(ncDat, "longitude")
# xlon <- ncvar_get(ncDat, "latitude")
# nc_close(ncDat)

# ## date
# start <- as.POSIXct("1966-10-03", "GMT")
# date  <- start + ((t-6)*24*60*60)
# doy   <- as.numeric(format(date, "%j"))
#
# 
# dat <- matrix(c(z), ncol = dim(z)[3], byrow = F)
# 
# pl <- (st_as_sf(data.frame(lat = c(xlon), lon = c(ylat), index = 1:(dim(z)[2]*dim(z)[1])), coords = c("lon", "lat")) %>%
#                   st_set_crs(4326) %>% st_transform(projSnow)) %>% mutate(snow = dat[,122]) %>% dplyr::select(snow) %>%
#   st_transform("+proj=longlat") %>% st_coordinates()
# 
# 
# rast <- raster::rasterize(pl[,1:2], raster(xmn = -179.9, xmx = 179.9, ymn = -90, ymx = 90, res = 5), field = dat[,122])
# rast[is.na(rast[])] <- 0
# tempGr <- st_make_grid(eaafMap$bbox, cellsize = 150000) %>% st_transform("+proj=longlat")
# extr   <- raster::extract(rast, (tempGr %>% st_centroid() %>% st_coordinates())[,1:2])
# tempPol <- tempGr %>% st_as_sf() %>% mutate(temp = extr) %>% st_transform(proj)
# plot(eaafMap$grid, col = "grey80", lty = 3)
# plot(eaafMap$map, col = "grey70", border = "grey70", add = T)
# plot(tempPol, add = T, border = NA, lwd = 0.1, col = ifelse(tempPol$temp==1, "darkblue", "grey80"))
# plot(eaafMap$bbox, add = T, border = "grey60")
# 
# 
# extr   <- raster::extract(rast, (tempGr %>% st_centroid() %>% st_coordinates())[,1:2])
# tempPol <- tempGr %>% st_as_sf() %>% mutate(temp = extr) %>% st_transform(proj)
# plot(eaafMap$grid, col = "grey80", lty = 3)
# plot(eaafMap$map, col = "grey70", border = "grey70", add = T)
# plot(tempPol, add = T, border = NA, lwd = 0.1)
# plot(eaafMap$bbox, add = T, border = "grey60")

# rastInd <- (st_as_sf(data.frame(lat = c(xlon), lon = c(ylat), index = 1:(dim(z)[2]*dim(z)[1])), coords = c("lon", "lat")) %>%
#               st_set_crs(4326) %>% st_transform(projSnow)) %>%
#   filter(apply(dat, 1, function(x) sum(x, na.rm = T)>0))

# crds_breed <- st_as_sf(geoTab, coords = c("lon_end", "lat_end"), crs = 4326) %>%
#   st_transform(projSnow) %>% rowwise() %>% mutate(NHWCindex = rastInd$index[which.min(st_distance(geometry, rastInd))])
# 
# sm <-  do.call("rbind", lapply(crds_breed$NHWCindex, function(x) {
# 
# #   ## example
# #   # expl <- crds_breed[7,]
# #   #   hex_id <- as.numeric(unlist(expl %>% st_transform(st_crs(hex)) %>% st_intersects(hex)))
# #   # x <- expl$NHWCindex
# #
#   tmp <- tibble(year = as.numeric(format(date, "%Y")), doy = doy, y = dat[x,])
#   # with(tmp %>% filter(year%in%c(1975:1985)), plot(doy, y))
#   mod1 <- suppressWarnings(with(tmp %>% filter(year%in%c(1975:1985)), gaussMLE(day = doy, size = rep(100, length(doy)), prob = y, thresh = 0.75)))
#   # lines(mod1$result)
#   mod2 <- suppressWarnings(with(tmp %>% filter(year%in%c(2000:2020)), gaussMLE(day = doy, size = rep(100, length(doy)), prob = y, thresh = 0.75)))
#   # lines(mod2$result, col  ="orange")
# 
#   # ggplot(tmp) +
#   #   geom_point(aes(x = doy, y = y), size = 1, alpha = 0.6) +
#   #   geom_line(data = tibble(doy = mod1$result$week, pred = mod1$result$fit, gr = "1950s") %>%
#   #               bind_rows(tibble(doy = mod2$result$week, pred = mod2$result$fit, gr = "2010s")),
#   #             mapping = aes(x = doy, y = pred, col = gr), linewidth = 1.5) +
#   #   scale_color_manual(values = c("orange", "darkblue"), name = "") +
#   #   geom_point(data = tibble(doy = c(mod1$start, mod2$start), pr = rep(0.75, 2), gr = c("1950s", "2010s")),
#   #              mapping = aes(x = doy, y = pr, col = gr), size = 4)+
#   #   xlab("Day of the year") + ylab("Snow (1 = snow covered, 0 = snow free)") +
#   #   theme_light()
# 
# 
#   tibble(start_hist = mod1$start, start_curr = mod2$start)
#   }))
# 
# breedTab <- crds_breed %>% bind_cols(sm)

# plot(breedTab %>% select(start_hist, start_curr) %>% st_drop_geometry() %>%
#           pivot_longer(cols = c("start_hist", "start_curr")) %>% mutate(name = as.factor(name)))
# 
# 
# ### Start 2060s

# snowTempTab <- tibble(ID = apply(st_intersects(breedTab %>% st_transform(st_crs(mudflatTab)), mudflatTab, sparse = F), 1, which),
#                       start_hist = breedTab$start_hist, start_curr = breedTab$start_curr) %>%
#                mutate(temp_hist    = apply(., 1, function(x) tempTab[x[1], x[2], 1]), temp_curr = apply(., 1, function(x) tempTab[x[1], x[3], 2])) %>%
#                mutate(start_future = apply(., 1, function(x) min(which(tempTab[as.numeric(x[1]), , 3] >= median(as.numeric(x[4:5])))))) %>%
#                mutate(temp_future  = apply(., 1, function(x) tempTab[x[1], x[6], 3])) %>% dplyr::select("ID", "start_hist", "start_curr", "start_future",
#                                                                                                         "temp_hist", "temp_curr", "temp_future")
# 
# ggplot(snowTempTab[,1:4] %>% setNames(c("ID", "1960s", "2010s", "2060s")) %>% pivot_longer(-ID), aes(x = value, fill = as.factor(name))) +
#   geom_histogram(show.legend = FALSE) +
#   scale_fill_manual(values = c("darkred", "darkblue", "orange")) +
#   ylab("Frequency") + xlab("Snowmelt [doy]") +
#   facet_wrap(~name, nrow = 3) +
#   theme_light()
# 
# 
# breedTab <- breedTab %>% mutate(start_future = snowTempTab$start_future)
# 
# save(breedTab, file = "breedTab_revision.rda")
load("breedTab_revision.rda")



