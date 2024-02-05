setClass(
  "SDPMig",
  slots = c(Init    = "list",
            Species = "list",
            Sites   = "list",
            Results = "list")
)

sizeParams     <- function(lbm) {
  
  out <- list(lbm   = lbm,
              mbm   = approx(c(15, 300), c(23.05, 643.4), xout = lbm)$y)
  
  ### Fuel deposition rate
  out$X1gX  <- (out$mbm-out$lbm)/100
  out$X1xkJ <- out$X1gX*32.77
  out$FDR   <- 2.37*(out$lbm/1000)^-0.27
  out$FDRg  <- (out$FDR*(out$lbm))/100
  out$FDRx  <- out$FDRg/out$X1gX
  out$daysToRefuel <- (out$mbm-out$lbm)/out$FDRg
  
  ### Metabolic rate
  out$DER    <- 912*(lbm/1000)^0.704
  out$BMR    <- 5.06*(lbm/1000)^0.729
  out$cond   <- 0.0067*(lbm)^0.4384
  out$Kesm   <- 41-(out$BMR/(1.2*out$cond))
  out$W      <- 57.3*(lbm/1000)^0.81
  
  ### Flight parameters
  x1  <- c(250, 144, 110, 20)
  dst <- c(14000, 7500, 6500, 3000) 
  fit1 <- lm(log(dst)~log(x1))
  pr   <- exp(predict(fit1))
  out$FlightCap <- exp(predict(fit1, newdata = data.frame(x1 = lbm))) + 1000
  
  # Energy exp ~ temp
  t  <- seq(-35, 45, length = 100)
  te <- 41 - (out$BMR/(1.2*out$cond))
  b = out$BMR - (-out$cond)*te
  
  tm <- data.frame(tm = seq(-35, te, length = 100), W = -out$cond*seq(-25, te, length = 100) + b)
  tm <- rbind(tm, data.frame(tm = seq(te, 45, length = 100), W = rep(out$BMR, 100)))
  
  out$EEFnc  <- suppressWarnings(approxfun(x = tm[,1], y = 2*tm[,2]*86.4, rule = 3))
  
  out$speed  <- 1440
  out$c      <- out$FlightCap / (1-(1+99/100)^-0.5)
  
  return(out)
}
Mat3Array      <- function(x) {
  
  rows = length(x)
  cols = length(x[[1]])
  d3 = length(x[[1]][[1]])
  
  out <- array(dim = c(rows, cols, d3))
  
  for(r in 1:rows) {
    for(c in 1:cols) {
      for(d.1 in 1:d3) {
        out[r, c, d.1] <- x[[r]][[c]][[d.1]]
      }
    }
  }
  out 
}
CalcTR         <- function(x, t) {
  TR <- approx(parms$xFTReward, parms$yFTReward, t, rule = 3)$y
  s  <- parms$w * (x - parms$xc)
  
  if(x == 0) { SR <- 0 } else {
    SR <- (((exp(s) - exp(-s)) / (exp(s) + exp(-s))) + 1)/2
  }
  TR * SR + parms$B0
}

makeSDPobjects <- function(paramTracks, parms) {
  
  mclapply(paramTracks, function(parTab) {
    
    # plot(st_as_sf(parTab$sites, coords = c("X", "Y"), crs = st_crs(map$map)) %>% mutate(temp = c(apply(parTab$tTab, 1, median), NA)) %>% select(temp))
    # plot(st_as_sf(parTab$sites, coords = c("X", "Y"), crs = st_crs(map$map)) %>% mutate(temp = parTab$dist[,1]) %>% select(temp))
    # plot(st_as_sf(parTab$sites, coords = c("X", "Y"), crs = st_crs(map$map)) %>% mutate(temp = parTab$bear[,1]) %>% select(temp))
    
    sites <- parTab$sites %>% as_tibble() %>% mutate(id = 1:nrow(parTab$sites)) %>%
      mutate(area0 = ifelse(areaHist==0, areaCurr, areaHist) + areaMang * parms["Mang"],
             qLake = quantile(area0, probs = parms["qLake"], na.rm = T),
             areaHist  = ifelse(lake, qLake, area0)) %>%
      mutate(area0 = areaCurr + areaMang * parms["Mang"],
             qLake = quantile(area0, probs = parms["qLake"], na.rm = T),
             areaCurr  = ifelse(lake, qLake, area0)) %>% select(-areaMang, -lake, -area0, -qLake) %>%
      mutate(pmax = pmax(areaHist, areaCurr))
    sites$pmax[1] <- max(sites$pmax) 
      
    sites <- filter(sites, pmax>quantile(pmax, probs = parms["Cut"], na.rm = T) | is.na(pmax)) %>%
      ## gain
      mutate(Qmat_hist = approx(range(pmax, na.rm = T), c(0.5, parms["Qmax"]), areaHist, rule = 2)$y,
             Qmat_curr = approx(range(pmax, na.rm = T), c(0.5, parms["Qmax"]), areaCurr, rule = 2)$y,
             gain_hist = (Qmat_hist * parTab$spParms$FDRx + parTab$spParms$EEFnc(parTab$spParms$Kesm)/parTab$spParms$X1xkJ) * 
               approx(seq(0,90, length = 100), as.numeric(parms["latF"])^c(seq(0,1, length = 100)), abs(Y))$y,
             gain_curr = (Qmat_curr * parTab$spParms$FDRx + parTab$spParms$EEFnc(parTab$spParms$Kesm)/parTab$spParms$X1xkJ) * 
               approx(seq(0,90, length = 100), as.numeric(parms["latF"])^c(seq(0,1, length = 100)), abs(Y))$y)
    
    ## expenditure
    expend <- t(apply(parTab$tTab[sites$id[-nrow(sites)],] + as.numeric(parms["Tfac"]), 1, function(x) parTab$spParms$EEFnc(x)/parTab$spParms$X1xkJ))  
    # plot(st_as_sf(sites, coords = c("X", "Y"), crs = st_crs(map$map)) %>% mutate(temp = c(apply(expend, 1, median), NA)) %>% select(temp))
    # plot(st_as_sf(sites, coords = c("X", "Y"), crs = st_crs(map$map)) %>% mutate(temp = sites$gain_curr) %>% select(temp))
    # plot(st_as_sf(sites, coords = c("X", "Y"), crs = st_crs(map$map)) %>% mutate(temp = parTab$dist[sites$id,sites$id][1,]) %>% select(temp))
    # plot(st_as_sf(sites, coords = c("X", "Y"), crs = st_crs(map$map)) %>% mutate(temp = parTab$bear[sites$id,sites$id][1,]) %>% select(temp))
    
    new(
      "SDPMig",
      Init  = list(
        MinT   = parTab$minT,
        MaxT   = parTab$maxT,
        NSites = nrow(sites)-1,
        MaxX   = parTab$MaxX
      ),
      Species = list(
        B0    = parTab$B0,
        w     = parTab$w,
        xc    = parTab$xc,
        c     = parTab$spParm$c,
        speed = parTab$spParm$speed,
        WindAssist = parTab$WindAssist,
        WindProb   = parTab$WindProb,
        ZStdNorm = parTab$ZStdNorm,
        PStdNorm = parTab$PStdNorm,
        xFTReward  = c(0, c(unlist(c(parTab$reward[1]-1, parTab$reward, parTab$reward[2]+1)), parTab$maxT) - parTab$minT),
        yFTReward  = c(0,0,2,2,0,0),
        decError   = parTab$decError
      ),
      Sites = list(
        crds  = as.matrix(sites[,c("X", "Y")]),
        dist  = parTab$dist[sites$id,sites$id],
        bear  = parTab$bear[sites$id,sites$id],
        b0    = sites$pred1,
        b1    = sites$pred2,
        b2    = sites$pred3,
        pred_a1 =  sites$pred_a1[1],
        pred_a2 =  sites$pred_a2[1],
        expend  =  expend,
        gain    =  if(parms["Time"]==0) sites$gain_hist else sites$gain_curr
      ),
      Results = list(
        FitnessMatrix     = NA,
        DecisionMatrix    = NA,
        ProbMatrix        = NA
      )
    )
    
  }, mc.cores = 4)
}

##########################
### Backward iteration ###
##########################

bwdIteration <- function(obj) {
  
  Init(obj@Init$MinT, 
       obj@Init$MaxT, 
       obj@Init$NSites, 
       obj@Init$MaxX,
       obj@Species$w,
       obj@Species$xc,
       obj@Species$B0,
       obj@Sites$b0,
       obj@Sites$b1,
       obj@Sites$b2,
       obj@Sites$pred_a1,
       obj@Sites$pred_a2,
       obj@Species$c,
       obj@Species$speed,
       obj@Species$WindAssist,
       obj@Species$WindProb,
       obj@Species$ZStdNorm,
       obj@Species$PStdNorm,
       as.numeric(obj@Species$xFTReward),
       obj@Species$yFTReward,
       obj@Species$decError,
       as.matrix(obj@Sites$dist),
       as.matrix(obj@Sites$bear),
       obj@Sites$gain,
       as.matrix(obj@Sites$expend))
  
  out <- BackwardIteration()
  obj@Results$FitnessMatrix <- Mat3Array(out[[1]])
  
  # # matplot(obj@Results$FitnessMatrix[,1109,], type= "l", col = "grey80", lwd = 1, lty = 1)
  # # plot(raster(obj@Results$FitnessMatrix[,,101]))

  DM <- array(dim = c(dim(obj@Results$FitnessMatrix)[c(2,1,3)], 2))
  DM[,,,1] <- Mat3Array(out[[2]])
  DM[,,,2] <- Mat3Array(out[[3]])
  obj@Results$DecisionMatrix <- DM

  PM <- array(dim = c(dim(obj@Results$FitnessMatrix)[c(2,1,3)], 2))
  PM[,,,1] <- Mat3Array(out[[4]])
  PM[,,,2] <- Mat3Array(out[[5]])
  obj@Results$ProbMatrix <- PM

  obj
}


############################
#### Forward Simulation#####
############################

fwdSimulation <- function(model, NrInd, start_t, start_site, start_x) {
  
  InitSim(model@Init$MinT, 
          model@Init$MaxT, 
          model@Init$NSites, 
          model@Init$MaxX,
          model@Species$w,
          model@Species$xc,
          model@Species$B0,
          model@Sites$b0,
          model@Sites$b1,
          model@Sites$b2,
          model@Sites$pred_a1,
          model@Sites$pred_a2,
          model@Species$c,
          model@Species$speed,
          model@Species$WindAssist,
          model@Species$WindProb,
          model@Species$ZStdNorm,
          model@Species$PStdNorm,
          model@Species$xFTReward,
          model@Species$yFTReward,
          model@Species$decError,
          model@Sites$dist,
          model@Sites$bear,
          model@Sites$gain,
          model@Sites$expend)

  
  x <- round(runif(NrInd, start_x[1], start_x[2]),0)
  
  if(length(start_site)>1 & length(start_site)<start_x[1]) {
    stop("start_site must have same length as numbers of individuals or a single site.")
  }
  if(length(start_site)==1) start_site <- rep(start_site, NrInd)
  
  SimOut = array(dim = c(length(x), 6, dim(model@Results$FitnessMatrix)[1]))
  
  ### First entry
  for(i in 1:dim(SimOut)[1]) {
    SimOut[i, ,start_t] <- c(start_t, start_site[i], x[i], 0, 0, 0)
  }

  
  ## SimOut: 1 = time, 2 = site, 3 = x, 4 = decision, 5 = flying, 6 = dead {
  for(time in 1:(dim(SimOut)[3]-1)) {
    
    for(ind in 1:dim(SimOut)[1]) {
      
      ## Not dead, not arrived, not flying
      if(!SimOut[ind, 6, time] & 
         sum(SimOut[ind, 2, ] >= nrow(model@Sites$crds), na.rm = T)<1 & !SimOut[ind, 5, time]) {
        
        ## Decision
        if(runif(1) <  model@Results$ProbMatrix[SimOut[ind, 2, time], time, SimOut[ind, 3, time], 1]) {
          decision  <- model@Results$DecisionMatrix[SimOut[ind, 2, time], time, SimOut[ind, 3, time], 1]
        } else {
          decision  <- model@Results$DecisionMatrix[SimOut[ind, 2, time], time, SimOut[ind, 3, time], 2]
        }
        
        ## Action
        if(decision>=0) { ## Flying
          
          fl_help = simFlying(decision, time-1, SimOut[ind, 2, time]-1, SimOut[ind, 3, time])
          
          nextt = fl_help[1] + 1
          if(nextt<=time) nextt <- time+1
          if(nextt>dim(SimOut)[3]) time <- dim(SimOut)[3]
          
          nextx = fl_help[2] + 1
          if(nextx <  0) {
            nextx = 0
            dead  = 1 } else  dead = 0
          if(nextx > model@Init$MaxX) nextx = model@Init$MaxX
          
          SimOut[ind,,nextt] = c(nextt, decision+1, nextx, NA, 0, dead)
          if(nextt>(time+1)) SimOut[ind,5:6,(time+1):(nextt-1)] = cbind(1,0)
          if(SimOut[ind, 6, nextt])  SimOut[ind, 6, nextt:dim(SimOut)[3]] = 1 ## if dead make dead till the end
          
          if(SimOut[ind, 2, nextt]==nrow(model@Sites$crds)) {
            SimOut[ind,2:6, nextt:dim(SimOut)[3]] <- SimOut[ind, 2:6, nextt]
            SimOut[ind,1,   nextt:dim(SimOut)[3]] <- seq(nextt, dim(SimOut)[3])
          }
          
        } else { ## Feeding
          
          fo_help = simForaging(abs(decision+1.0), time-1, SimOut[ind, 2, time]-1, SimOut[ind, 3, time])
          
          newx = fo_help[1]+1
          dead = fo_help[2]
          if(newx<=0) dead = 1
          
          if(newx > model@Init$MaxX) newx = model@Init$MaxX
          
          SimOut[ind,,time+1] = c(time+1, SimOut[ind, 2, time], newx, abs(decision+1.0), 0, dead)
          
          if(SimOut[ind, 6, time+1])  SimOut[ind, 6, (time+1):dim(SimOut)[3]] = 1 ## if dead make dead till the end
          
        }
        
      }
      
    } ## Ind loop
  } ## time loop
  
  SimOut       <- SimOut[,-5,] 
  SimOut[,2,] <- SimOut[,2,]-1
  
  SimOut
}

condProfile <- function(simu, model) {
  trkS  <- simu[,2,]+1
  trkX  <- simu[,3,]
  trkD  <- simu[,5,]
  
  cls <- rainbow(max(trkS, na.rm = T))[sample(1:max(trkS, na.rm = T))]
  
  opar <- par(mfrow = c(1,1))
  matplot(model@Init$MinT:model@Init$MaxT, t(trkX), type = "n", 
         xaxt = "n", xlab  = "", ylab = "Body condition", las = 1)
  for(i in 1:dim(simu)[1]) {
    if(all(trkD[i,]==0)){
    tmp <- cbind(model@Init$MinT:model@Init$MaxT, trkS[i,], trkX[i,])[!is.na(trkS[i,]),]
    points(tmp[,1], tmp[,3], type = "l", col = "grey90")
    points(tmp[,1], tmp[,3], pch = 16, cex = 0.5, col = cls[tmp[,2]])
    }
  }
  axis(1, at = seq(model@Init$MinT,model@Init$MaxT, by = 15), labels = format(as.POSIXct("2012-01-01")+seq(model@Init$MinT,model@Init$MaxT, by = 15)*24*60*60, "%b-%d"))
  par(opar)
}

pltNetwork <- function(simu, model, map = eaafMap) {

  sitesCrds <- model@Sites$crds
  dead      <- apply(simu[,5,], 1, sum, na.rm = T)
  trkS      <- simu[dead<1,2,]+1
 
  diagT <-   apply(trkS, 1, function(x) {
    tt <- as.data.frame(table(x[!is.na(x)]))
    m  <- merge(data.frame(Var1 = 1:nrow(sitesCrds)), tt, all.x = T)
    m[nrow(m),2] <- NA
    ifelse(is.na(m[,2]), 0, m[,2])
  })
  tmp.out   <- matrix(0, ncol = nrow(sitesCrds), nrow = nrow(sitesCrds)) 
  diag(tmp.out) <- apply(diagT, 1, sum)
  
  
  transT <- as.data.frame(table(apply(do.call("rbind", apply(trkS, 1, function(x) {
    tmp01 <- cbind(x[!is.na(x)][-sum(!is.na(x))], x[!is.na(x)][-1])
    tmp01[tmp01[,1]!=tmp01[,2] & !is.na(tmp01[,1]) & !is.na(tmp01[,2]),]
  })), 1, function(y) paste(y, collapse = "_"))))
  trans <- cbind(t(apply(transT, 1, function(z) as.numeric(strsplit(z[1], "_")[[1]]))), transT[,2])
  
  
  lwd.f <- approxfun(c(0, max(trans[,3])), c(0.1, 5), rule = 3)
  
  opar <- par(mfrow = c(1,1))
  plot(map$map, col = "grey90", border = "grey50")
  pr <- project(model@Sites$crds, proj = st_crs(map$map)$input)
  segments(pr[trans[,1],1], pr[trans[,1],2], pr[trans[,2],1], pr[trans[,2],2], lwd = lwd.f(trans[,3]), col = "grey20")
  
  site01 <- cbind(sitesCrds, d = diag(tmp.out))
  cex.f <- approxfun(x = c(1, max(site01[-1,3])), y = c(0.1, 15), rule = 3)
  
    cex <- cex.f(site01[,3])
    cex[model@Sites$start] <- 2
    
    pch = rep(21, nrow(site01))
    pch[model@Sites$start] <- 23
    
  points(project(as.matrix(site01[,1:2]), proj = st_crs(map$map)$input), cex = cex, pch = pch, col  = "grey10",
         bg = c("white", rep(adjustcolor("orange", alpha.f = 0.55), nrow(site01)-1)), lwd = 1.4)
  points(pr[1,1], pr[1,2], pch = 21, bg = "white")
  par(opar)
  
  cat(sum(dead<1)/length(dead))
  
}

plotBatchNetwork <- function(simu, model, map = eaafMap) {
  require(gridExtra)
  
  ###
  st_crds <- st_as_sf(as.data.frame(model[[1]]@Sites$crds), coords = c("X", "Y"), crs = 4326) %>% st_transform(st_crs(map$map)) %>%
    mutate(area = model[[1]]@Sites$gain)

  mp <- ggplot(map$map) +
    geom_sf() + theme_bw() +
    geom_sf(st_crds, mapping = aes(geometry = geometry, col = area, size = area), shape = 16, show.legend = F)
  
  maps <- mclapply(1:length(simu), function(ind) {
    
    subSimu <- simu[[ind]]
    
    if(!is.null(subSimu)) {
      
      crds_sf <- st_as_sf(as.data.frame(model[[ind]]@Sites$crds), coords = c("X", "Y"), crs = 4326) %>% st_transform(st_crs(map$map))
      
      dead      <- apply(subSimu[,5,], 1, sum, na.rm = T)
      trkS      <- subSimu[dead<1,2,]+1
      
      fitness   <- apply(subSimu[dead<1, c(2,3),], 1, function(x) {
        if(any(x[1,!is.na(x[1,])]==(nrow(crds_sf)-1))) {
          time = min(which(x[1,]==(nrow(crds_sf)-1)))
          cond = x[2,time]
          model[[ind]]@Results$FitnessMatrix[time, nrow(crds_sf), cond]
        } else NA
      })
      
      
      diagT <-  apply(apply(trkS, 1, function(x) {
        tt <- as.data.frame(table(x[!is.na(x)]))
        m  <- merge(data.frame(Var1 = 1:nrow(crds_sf)), tt, all.x = T)
        m[nrow(m),2] <- NA
        ifelse(is.na(m[,2]), 0, m[,2])
      }), 1, sum); diagT[1] <- 0
      
      crds_sf <- crds_sf %>% mutate(ts = diagT)
      
      # mp +
      #   geom_sf(crds_sf, mapping = aes(geometry = geometry, size = ts))
      
      transT <- as.data.frame(table(apply(do.call("rbind", apply(trkS, 1, function(x) {
        tmp01 <- cbind(x[!is.na(x)][-sum(!is.na(x))], x[!is.na(x)][-1])
        tmp01[tmp01[,1]!=tmp01[,2] & !is.na(tmp01[,1]) & !is.na(tmp01[,2]),]
      })), 1, function(y) paste(y, collapse = "_"))))
      
      trans <- cbind(t(apply(transT, 1, function(z) as.numeric(strsplit(z[1], "_")[[1]]))), transT[,2])
      
      trans_sf <- st_as_sfc(lapply(1:nrow(trans), function(x) st_linestring(crds_sf[trans[x,1:2],] %>% st_coordinates())), crs = st_crs(crds_sf))
      
      mapOut <- mp + 
        geom_sf(trans_sf, mapping = aes(geometry = trans_sf, size = trans[,3]), col = adjustcolor("grey35", alpha.f = 0.6), show.legend = FALSE) +
        geom_sf(crds_sf %>% filter(ts>0), mapping = aes(geometry = geometry, size = ts), fill = adjustcolor("orange", alpha.f = 0.7), shape = 21, show.legend = FALSE) +
        geom_point(data = data.frame(st_as_sfc(list(st_point(as.numeric(geoTab[ind, c("lon_end", "lat_end")]))), crs = 4326) %>%
                    st_transform(st_crs(trans_sf)) %>% st_coordinates), mapping = aes(x = X, y = Y), shape = 23, size = 5) +
        geom_point(data = data.frame(st_as_sfc(list(st_point(as.numeric(geoTab[ind, c("lon_start", "lat_start")]))), crs = 4326) %>%
                    st_transform(st_crs(trans_sf)) %>% st_coordinates), mapping = aes(x = X, y = Y), shape = 23, size = 5) +
        scale_x_continuous(limits = st_bbox(eaafMap$map)[c(1,3)], expand = c(0, 0)) +
        scale_y_continuous(limits = st_bbox(eaafMap$map)[c(2,4)], expand = c(0, 0)) +
        labs(title = geoTab$species[ind],
             subtitle = glue::glue("survival: {sum(dead<1)/length(dead)}, success = {round(median(fitness),3)}."))
      
      mapOut
      
    } else {
      
      mapOut <- map + scale_x_continuous(limits = st_bbox(eaafMap$map)[c(1,3)], expand = c(0, 0)) +
        scale_y_continuous(limits = st_bbox(eaafMap$map)[c(2,4)], expand = c(0, 0)) +
        labs(title = geoTab$species[ind],
             subtitle = glue::glue("survival: {NA}, success = {NA}."))
      
    }
  }, mc.cores = 4)
  
  do.call("grid.arrange", c(maps, ncol=3))
  
}

rangeSDP <- function(x){
  t <- (x-min(x))/(max(x)-min(x))
  approxfun(c(0, max(t)), c(0,2))(t)
}

depth <- function(this) ifelse(is.list(this), 1L + max(sapply(this, depth)), 0L)


transTrack <- function(sf_track, start, end) {
  
  sf_stat_tracks <- sf_track %>% mutate(days = as.numeric(date2 - date1)) %>% filter(days>1.5)
  
  out <- merge(data.frame(tm = start:end), 
               data.frame(tm = c(as.numeric(format(sf_stat_tracks$date1, "%j")),
                                 as.numeric(format(sf_stat_tracks$date2, "%j"))),
                          lon = rep(st_coordinates(sf_stat_tracks)[,1], 2), 
                          lat = rep(st_coordinates(sf_stat_tracks)[,2], 2),
                          id  = rep(1:2, each = nrow(sf_stat_tracks))), by = "tm", all.x = T)
  
  int <- is.na(out$id) & c(diff(zoo::na.approx(out$id, rule = 3)),0)<0
  
  st_as_sf(as.data.frame(apply(out[!int,2:3], 2, zoo::na.approx, rule = 3)), coords = c("lon", "lat"), crs = st_crs(sf_track))
  
}
