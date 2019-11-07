##' Construct a sampler to draw n bivariate Normal deviates.
##'
##' Constructs a function that draws n deviates from n bivariate
##' Normal distributions with covariances determined by \code{S} which
##' must be either 2x2 covariance matrix or an nx2x2 array of
##' covariance matrices.
##' @title Bivariate Normal Samples
##' @param S a covariance matrix or an array of covariance matrices
##' @param s a scale factor applied to S
##' @param n number of deviates to draw.
##' @return A function that draws bivariate Normal deviates with mean
##' given by its first argument.
##' @export
bvnorm <- function(S,s=1,n=1) {
  if(length(dim(S))==2)
    S <- array(rep(S,each=n),c(n,2,2))
  else
    n <- dim(S)[1]
  
  A <- matrix(0,n,2)
  B <- double(n)
  for(k in 1:n) {
    L <- chol(s^2*S[k,,])
    A[k,] <- L[1,]
    B[k] <- L[2,2]
  }
  function(x) {
    x <- x+A*rnorm(n)
    x[,2] <- x[,2]+B*rnorm(n)
    x
  }
}



invChanges <- function(twl, calib, fixed = NULL, location = NULL, cpu = detectCores()-1) {
  
  require(doParallel)
  
  if(!is.null(fixed)) {
    if(is.null(location)) stop("location cannot be NULL if fixed range has been specified.")
    
    if(!is.null(fixed[1])) ind <- which(twl$Twilight<=as.POSIXct(fixed[1], tz = "GMT"))
    if(!is.null(fixed[1])) ind <- c(ind, which(twl$Twilight>=as.POSIXct(fixed[2], tz = "GMT")))
    
    twl$Twilight[ind] <- twilight(twl$Twilight[ind], lon = location[1], lat = location[2], twl$Rise[ind], zenith = calib[1], iters = 3)
  }
  
  
  
  tmp <- apply(cbind(rep(1, 1000)), 1, function(x) {
    tmp1 <- rgamma(500, calib[3], calib[4])*ifelse(twl[x,2], 1, -1)  + (rgamma(1, calib[3], calib[4])*ifelse(twl[x,2], -1, 1))
    tmp1/60
  })
  f1 <- fitdistr(c(tmp), "normal")$estimate
  f2 <- fitdistr(apply(cbind(rnorm(1000, -5, f1[2]), rnorm(1000, 5, f1[2])), 1, mean), "normal")$estimate

  parms <- matrix(c(f1, f2), ncol = 2, byrow = T)
  
  nsim = 1e+05
  s <- seq(-5, 5, by =10/120)
  
  cl <- parallel::makeCluster(cpu)
  registerDoParallel(cl)
  # clusterExport(cl, c("parms"))
  
  ovl <- apply(cbind(s), 1, function(x) {
    
    post  <-  hist(rnorm(100000, 0, parms[1,2]), breaks = seq(-10, 10, length.out = 1000), plot = F)$density
    prior <-  hist(rnorm(100000, x, parms[1,2]), breaks = seq(-10, 10, length.out = 1000), plot = F)$density
    
    mindens <- apply(cbind(post, prior), 1, min)
    o1 <- sum(mindens) * (10 - -10)/1000
    
    post  <-  hist(rnorm(100000, 0, parms[2,2]), breaks = seq(-10, 10, length.out = 1000), plot = F)$density
    prior <-  hist(rnorm(100000, x, parms[2,2]), breaks = seq(-10, 10, length.out = 1000), plot = F)$density
    
    mindens <- apply(cbind(post, prior), 1, min)
    o2 <- sum(mindens) * (10 - -10)/1000
    
    c(o1, o2)
  })
  
  ovlF <- list(Tw = approxfun(s, ovl[1,], rule = 3), No = approxfun(s, ovl[2,], rule = 3))  
  
  spl1 <- split(twl[,c("Twilight", "Rise")], twl$Rise)
  
  ovl1 <- lapply(spl1, function(x) {
    ind <- list(cbind(1:nrow(x), 0:(nrow(x)-1)), 
                cbind(2:(nrow(x)+1), 0:(nrow(x)-1)),
                cbind(1:nrow(x), -1:(nrow(x)-2)),
                cbind(3:(nrow(x)+2), -1:(nrow(x)-2)))
    do.call("cbind", lapply(ind, function(y) {
      foreach(h=1:nrow(y), .export = c("x", "ovlF"), .combine = "c") %dopar% {
        z <- y[h,]
        if(all(z>=1 & z<nrow(x))) {
          ovlF[[1]](as.numeric(abs(difftime(x[z[2],1], x[z[1],1], units = "hours"))-(24*abs(diff(z)))))
        } else 1
        }
    }))
  })
  ovl1[[1]][apply(ovl1[[1]],1,function(x) ifelse(all(x<0.01), TRUE, FALSE)),] <- 1
  ovl1[[2]][apply(ovl1[[2]],1,function(x) ifelse(all(x<0.01), TRUE, FALSE)),] <- 1
  
  gl <- export2GeoLight(twl)
  nm <- gl[,1] + difftime(gl[,2], gl[,1], units = "secs")/2
  
  spl2 <- split(nm, gl[,3])
  
  ovl2 <- lapply(spl2, function(x) {
    ind <- list(cbind(1:length(x), 0:(length(x)-1)), 
                cbind(2:(length(x)+1), 0:(length(x)-1)),
                cbind(1:length(x), -1:(length(x)-2)),
                cbind(3:(length(x)+2), -1:(length(x)-2)))
    do.call("cbind", lapply(ind, function(y) {
      foreach(h=1:nrow(y), .export = c("x", "ovlF"), .combine = "c") %dopar% {
        z <- y[h,]
        if(all(z>=1 & z<length(x))) {
          ovlF[[2]](as.numeric(abs(difftime(x[z[2]], x[z[1]], units = "hours"))-(24*abs(diff(z)))))
        } else 1
      }
    }))
  })
  ovl2[[1]][apply(ovl2[[1]],1,function(x) ifelse(all(x<0.01), TRUE, FALSE)),] <- 1
  ovl2[[2]][apply(ovl2[[2]],1,function(x) ifelse(all(x<0.01), TRUE, FALSE)),] <- 1
  
  
  stopCluster(cl)
  
  res <- list(Rise = cbind(data.frame(Twilight = twl$Twilight[twl$Rise],  Rise = twl$Rise[twl$Rise],  ovl1[[which(names(spl1)=="TRUE")]])),
              Set  = cbind(data.frame(Twilight = twl$Twilight[!twl$Rise], Rise = twl$Rise[!twl$Rise], ovl1[[which(names(spl1)=="FALSE")]])),
              Noon = cbind(data.frame(Time = spl2[[which(names(spl2)=="1")]], ovl2[[which(names(spl2)=="1")]])),
              Moon = cbind(data.frame(Time = spl2[[which(names(spl2)=="2")]], ovl2[[which(names(spl2)=="2")]])))

  res
  
}



extractMovements <- function(res, offset, threshold = 0.6, days = 2, na.break = NULL, exclude = NULL, plot = TRUE) {
  
  twl <- do.call("rbind", res[c(1,2)])[order(do.call("rbind", res[c(1,2)])[,1]),c(1:2)]
  
  if(length(threshold)==1) threshold <- rep(threshold, 4)
  threshold[is.na(threshold)] <- -1
  
  opar <- par(mfrow = c(4,1), mar = c(2,3,1,4), oma = c(3,3,1,3))
  
  layout(matrix(c(1:10), ncol = 2, byrow = T), widths = c(0.8, 0.2))
  
  ind1 <- apply(res[[1]][,c(-c(1:2))], 1, function(x) all(x<=threshold[1]))
  
  plot(res[[1]]$Twilight, hourOffset(as.hour(res[[1]]$Twilight), offset%%24), pch = 21, 
       bg = adjustcolor("firebrick", alpha.f = 0.2), yaxt = "n", ylab = "")
  points(res[[1]]$Twilight[ind1], hourOffset(as.hour(res[[1]]$Twilight), offset%%24)[ind1], pch = 16, col = "firebrick")
  
  par(new = TRUE)
  plot(res[[1]]$Twilight, rep(1, nrow(res[[1]])), type = "n", ylim = c(0,1), las = 1, xaxt = "n", yaxt = "n", ylab = "")
  abline(h = threshold[1], lty = 2, col = adjustcolor("grey80", alpha.f = 0.8))
  segments(res[[1]]$Twilight, apply(res[[1]][,c(-c(1:2))], 1, min, na.rm = T),
           res[[1]]$Twilight, apply(res[[1]][,c(-c(1:2))], 1, max, na.rm = T), col = adjustcolor("grey50", alpha.f = 0.5))
  points(res[[1]]$Twilight,   apply(res[[1]][,c(-c(1:2))], 1, max, na.rm = T), pch = 16, cex = 0.4)
  axis(4, las = 1)
  
  h <- hist(do.call("c", res[[1]][,c(-c(1:2))]), breaks = seq(0, 1, length = 11), plot = FALSE)
  barplot(h$counts, horiz = T, col = ifelse(h$mids<=threshold[1], "firebrick", "grey90"))

  ind2 <- apply(res[[2]][,c(-c(1:2))], 1, function(x) all(x<=threshold[2]))
  plot(res[[2]]$Twilight, hourOffset(as.hour(res[[2]]$Twilight), offset%%24), pch = 21, 
       bg = adjustcolor("cornflowerblue", alpha.f = 0.2), yaxt = "n", ylab = "")
  points(res[[2]]$Twilight[ind2], hourOffset(as.hour(res[[2]]$Twilight), offset%%24)[ind2], pch = 16, col = "cornflowerblue")
  
  par(new = TRUE)
  plot(res[[2]]$Twilight, rep(1, nrow(res[[2]])), type = "n", ylim = c(0,1), las = 1, xaxt = "n", yaxt = "n", ylab = "")
  abline(h = threshold[2], lty = 2, col = adjustcolor("grey80", alpha.f = 0.8))
  segments(res[[2]]$Twilight, apply(res[[2]][,c(-c(1:2))], 1, min, na.rm = T),
           res[[2]]$Twilight, apply(res[[2]][,c(-c(1:2))], 1, max, na.rm = T), col = adjustcolor("grey50", alpha.f = 0.5))
  points(res[[2]]$Twilight,   apply(res[[2]][,c(-c(1:2))], 1, max, na.rm = T), pch = 16, cex = 0.4)
  axis(4, las = 1)
  
  h <- hist(do.call("c", res[[2]][,c(-c(1:2))]), breaks = seq(0, 1, length = 11), plot = FALSE)
  barplot(h$counts, horiz = T, col = ifelse(h$mids<=threshold[2], "firebrick", "grey90"))
  
  
  ind3 <- apply(res[[3]][,-1], 1, function(x) all(x<=threshold[3]))
  plot(res[[3]]$Time, hourOffset(as.hour(res[[3]]$Time)), pch = 21, 
       bg = adjustcolor("orange", alpha.f = 0.2), yaxt = "n", ylab = "")
  points(res[[3]]$Time[ind3], hourOffset(as.hour(res[[3]]$Time))[ind3], pch = 16, col = "orange")
  par(new = TRUE)
  plot(res[[3]]$Time, rep(1, nrow(res[[3]])), type = "n", ylim = c(0,1), las = 1, xaxt = "n", yaxt = "n", ylab = "")
  abline(h = threshold[3], lty = 2, col = adjustcolor("grey80", alpha.f = 0.8))
  segments(res[[3]]$Time, apply(res[[3]][,c(-c(1))], 1, min, na.rm = T),
           res[[3]]$Time, apply(res[[3]][,c(-c(1))], 1, max, na.rm = T), col = adjustcolor("grey50", alpha.f = 0.5))
  points(res[[3]]$Time,   apply(res[[3]][,c(-c(1))], 1, max, na.rm = T), pch = 16, cex = 0.4)
  axis(4, las = 1)
  
  h <- hist(do.call("c", res[[3]][,-1]), breaks = seq(0, 1, length = 11), plot = FALSE)
  barplot(h$counts, horiz = T, col = ifelse(h$mids<=threshold[3], "firebrick", "grey90"))
  
  
  ind4 <- apply(res[[4]][,-1], 1, function(x) all(x<=threshold[4]))
  plot(res[[4]]$Time, hourOffset(as.hour(res[[4]]$Time)), pch = 21, 
       bg = adjustcolor("darkgreen", alpha.f = 0.2), yaxt = "n", ylab = "")
  points(res[[4]]$Time[ind4], hourOffset(as.hour(res[[4]]$Time))[ind4], pch = 16, col = "darkgreen")
  par(new = TRUE)
  plot(res[[4]]$Time, rep(1, nrow(res[[4]])), type = "n", ylim = c(0,1), las = 1, xaxt = "n", yaxt = "n", ylab = "")
  abline(h = threshold[4], lty = 2, col = adjustcolor("grey80", alpha.f = 0.8))
  segments(res[[4]]$Time, apply(res[[4]][,c(-c(1))], 1, min, na.rm = T),
           res[[4]]$Time, apply(res[[4]][,c(-c(1))], 1, max, na.rm = T), col = adjustcolor("grey50", alpha.f = 0.5))
  points(res[[4]]$Time,   apply(res[[4]][,c(-c(1))], 1, max, na.rm = T), pch = 16, cex = 0.4)
  axis(4, las = 1)
  
  
  h <- hist(do.call("c", res[[4]][,-1]), breaks = seq(0, 1, length = 11), plot = FALSE)
  barplot(h$counts, horiz = T, col = ifelse(h$mids<=threshold[4], "firebrick", "grey90"))
  
  if(!is.null(na.break)) {
    ind5.0 <- which(diff(as.numeric(twl$Twilight)/60/60)>median(diff(as.numeric(twl$Twilight)/60/60))*na.break)+c(0, 1)
    ind5 <- rep(FALSE, nrow(twl)); ind5[ind5.0] <- TRUE
  } else ind5 <- rep(FALSE, length = nrow(res$Noon))

  tm <- as.POSIXct(c(as.numeric(min(twl$Twilight))-500, as.numeric(max(twl$Twilight))+500,                  ## Start/End
                     as.numeric(res[[1]]$Twilight[ind1])-500, as.numeric(res[[2]]$Twilight[ind2])-500,      ## Rise/Set 
                     as.numeric(res[[3]]$Time[ind3])-17*60*60, as.numeric(res[[4]]$Time[ind4])-17*60*60,    ## Moon/Noon
                     as.numeric(twl$Twilight[ind5])+ c(+500, -500)) ,                                      
                     origin = "1970-01-01", tz = "GMT")
  
  gr <- cut(as.numeric(twl[,1]), unique(as.numeric(tm[order(tm)])), labels = FALSE)
  
  tbl <- as.data.frame(table(gr))
  
    gr[gr%in%tbl[tbl[,2]<=(days*2),1]] <- NA
  
    
  if(is.numeric(exclude)) {
    
    tmp    <- split(data.frame(twl$Twilight, gr), f = gr)
    gr_tmp <- as.vector(unlist(lapply(tmp, function(x) {
      if(nrow(x)>1 & nrow(x)>exclude*2) {
        tmp01 <- x[,2]
        tmp01[c(1:exclude, (nrow(x)-exclude+1):nrow(x))] <- NA
        tmp01
      } else {rep(NA, nrow(x))}
    })))
    
    tmpOut <- rbind(data.frame(Twilight = twl$Twilight[!is.na(gr)], Gr = gr_tmp),
                    data.frame(Twilight = twl$Twilight[is.na(gr)], Gr = NA))
    
    
    gr <- tmpOut[order(tmpOut$Twilight),2]
    # gr <- na.approx(gr_tmp, rule = 3)
  }
    
    
  out <- c(1)
  for(i in 2:nrow(twl)) out <- c(out, ifelse(is.na(gr[i]) | (!is.na(gr[i]) & is.na(gr[i-1])) | gr[i]!=gr[i-1], out[i-1]+1, out[i-1]))
    
  cols <- rainbow(max(out, na.rm = T))[sample(1:max(out, na.rm = T))]
  
  plot(twl$Twilight,  hourOffset(as.hour(twl$Twilight)), 
       type  = "n", xlab = "", ylab = "", yaxt = "n")
  r <- range(hourOffset(as.hour(twl$Twilight)))
  for(i in unique(out[!is.na(out)])){
    if(length(twl$Twilight[!is.na(out) & out==i])>1) {
    rect(min(twl$Twilight[!is.na(out) & out==i]), r[1], 
             max(twl$Twilight[!is.na(out) & out==i]), r[2], lty = 2, col = "grey90", border = "grey40") 
    }
  }
  points(twl$Twilight,  
         hourOffset(as.hour(twl$Twilight)), 
         pch = 21, bg = ifelse(out%in%as.data.frame(table(out))[which(as.data.frame(table(out))[,2]==1),1], "grey99", cols[out]), yaxt = "n")
 
  par(opar)
  
  out
}
  


gammaSens        <- function(rise, set, calib, sr.proposal, ss.proposal, range, n.thin, n.iters) {
  
  sr <- solar(rise)
  ss <- solar(set)
  
  ## Initialize chain from best estimate of location if stationary
  p0 <- as.vector(thresholdEstimate(rise, set, calib[2]))
  sr.p <- p0
  ss.p <- p0
  
  ## Initialize cached solar times
  sr.solar <- twilightSolartime(sr,sr.p[1L],sr.p[2L],TRUE,calib[2])
  ss.solar <- twilightSolartime(ss,ss.p[1L],ss.p[2L],FALSE, calib[2])
  
  ## Ensure approximation consistency
  sr$solarTime <- sr.solar
  ss$solarTime <- ss.solar
  
  ## Initialize cached log posterior
  sr.logp <- dgamma(0, calib[3], calib[4], log=TRUE)
  ss.logp <- dgamma(0, calib[3], calib[4], log=TRUE)
  
  P.sr <- matrix(0,2L,n.iters)
  P.ss <- matrix(0,2L,n.iters)
  for(k1 in 1:n.iters) {
    for(k2 in 1:n.thin) {
      
      ## Propose new sunrise location
      sr.p1 <- sr.proposal(sr.p)
      sr.solar1 <- twilightSolartime(sr,sr.p1[1L],sr.p1[2L],TRUE,calib[2])
      if(is.finite(sr.solar1) && gcdist(sr.p1,ss.p) < range) {
        ## When proposal in range compute time error
        sr.delta <- 4*(sr$solarTime-sr.solar1)
        if(sr.delta>0) {
          ## Metropolis rule
          sr.logp1 <- dgamma(sr.delta,calib[3],calib[4],log=TRUE)
          if(sr.logp1-sr.logp > log(runif(1))) {
            ## Accept proposal
            sr.p <- sr.p1
            sr.solar <- sr.solar1
            sr.logp <- sr.logp1
          }
        }
      }
      
      ## Propose new sunset location
      ss.p1 <- ss.proposal(ss.p)
      ss.solar1 <- twilightSolartime(ss,ss.p1[1L],ss.p1[2L],FALSE,calib[2])
      if(is.finite(ss.solar1) && gcdist(sr.p,ss.p1) < range) {
        ## When proposal in range compute time error
        ss.delta <- 4*(ss.solar1-ss$solarTime)
        if(ss.delta>0) {
          ## Metropolis rule
          ss.logp1 <- dgamma(ss.delta,calib[3],calib[4],log=TRUE)
          if(ss.logp1-ss.logp > log(runif(1))) {
            ## Accept proposal
            ss.p <- ss.p1
            ss.solar <- ss.solar1
            ss.logp <- ss.logp1
          }
        }
      }
      
    }
    ## Record locations at sr/ss
    P.sr[,k1] <- sr.p
    P.ss[,k1] <- ss.p
  }
  
  list(p0=p0,rise=t(P.sr),set=t(P.ss))
}
groupSensitivity <- function(twl, gr, calib, range = 50, n.thin = 10, n.iters = 1000, cpu = detectCores()-1) {
  
  require(doSNOW)
  require(abind)
  
  ## Great circle distance (km)
  gcdist <- function(a,b) {
    rad <- pi/180
    6378.137*acos(pmin.int(
      cos(rad*a[2L])*cos(rad*b[2L])*
        cos(rad*(b[1L]-a[1L]))+sin(rad*a[2L])*
        sin(rad*b[2L]),
      1))
  }
  
  p0 <- thresholdPath(twl$Twilight, twl$Rise, zenith = calib[1])$x
  
  
  cl <- parallel::makeCluster(cpu)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = nrow(twl), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  out1 <- abind(foreach(i = 1:nrow(twl), .packages = c("SGAT", "geosphere"), .export = c("bvnorm", "gammaSens", "gcdist"), .options.snow = opts) %dopar% {
    
    if(i>1 && gr[i]==gr[i-1]) {
      
      ## Solar properties at ss/sr
      if(twl$Rise[i]) {
        rise <- twl$Twilight[i]
        set  <- twl$Twilight[i-1]
      } else {
        rise <- twl$Twilight[i-1]
        set  <- twl$Twilight[i]
      }  
      
      ## Initial proposals
      sr.proposal <- bvnorm(diag(c(0.3,0.3)^2))
      ss.proposal <- bvnorm(diag(c(0.3,0.3)^2))
      
      fit <- gammaSens(rise, set, calib, sr.proposal, ss.proposal, range, n.thin, n.iters)
          
        if(!is.nan(fit$p0[2])) {
          
          ## Update proposals
          sr.proposal <- tryCatch(bvnorm(S=cov(fit$rise[!is.nan(fit$rise[,2]),]),s=0.4),
                                  error = function(e) bvnorm(diag(c(0.3,0.3)^2)))
          ss.proposal <- tryCatch(bvnorm(S=cov( fit$set [!is.nan(fit$set[,2]),]),s=0.4),
                                  error = function(e) bvnorm(diag(c(0.3,0.3)^2)))
          
          
          fit <- gammaSens(rise, set, calib, sr.proposal, ss.proposal, range, n.thin, n.iters)
          
          
          midPoint(fit$rise[!is.nan(fit$rise[,1]),],
                   fit$set[!is.nan(fit$set[,1]),])
        } else {
          matrix(p0[i,], ncol = 2, nrow = n.iters, byrow = T)
        }
        } else {
          matrix(p0[i,], ncol = 2, nrow = n.iters, byrow = T)
        }
    
  }, along = 3)
  close(pb)
  stopCluster(cl) 
  
  out <- aperm(out1, c(3,2,1))
  
  list(twl = twl, gr = gr, x = out)
}
groupSummary     <- function(fit,  probs=0.95) {
  
  stat <- function(x) c(mean=mean(x, na.rm = T),sd=sd(x, na.rm = T),quantile(x,prob=c(0.5,(1-probs)/2,1-(1-probs)/2), na.rm = T))
  
  sm <- data.frame(Tmin = as.POSIXct(tapply(fit$twl$Twilight, fit$gr, min), origin = "1970-01-01", tz = "GMT"), 
                   Tmax = as.POSIXct(tapply(fit$twl$Twilight, fit$gr, max), origin = "1970-01-01", tz = "GMT"),
                   Group = unique(fit$gr))
  
  lon = do.call("rbind", tapply(fit$x[,1,], rep(fit$gr, dim(fit$x)[3]), stat))
    colnames(lon) <- paste("Lon",colnames(lon),sep=".")
  lat <- do.call("rbind", tapply(fit$x[,2,], rep(fit$gr, dim(fit$x)[3]), stat))
    colnames(lat) <- paste("Lat",colnames(lat),sep=".")

  as.data.frame(cbind(sm, lon, lat))
  
}
mergeGroups      <- function(fit, n.sample = 1000, threshold = 0.6, plot = TRUE) {
  
  gr <- fit$gr
  
  GrInd  <- unique(gr[gr%in%gr[duplicated(gr)]])
  indTab <- cbind(GrInd[-length(GrInd)], GrInd[-1], NA) 
  
  repeat{
    for(i in 1:nrow(indTab)) {
      
      s1 <- apply(fit$x[gr==indTab[i,1],,], 2, I)
        s1[,1] <- ifelse(s1[,1]>180, -180+(s1[,1]-180), s1[,1])
      s2 <- apply(fit$x[gr==indTab[i,2],,], 2, I)
        s2[,1] <- ifelse(s2[,1]>180, -180+(s2[,1]-180), s2[,1])
        
      d1 <- c(distm(s1[sample(1:nrow(s1), n.sample),]),
              distm(s2[sample(1:nrow(s2), n.sample),]))
      d2 <- distm(s1[sample(1:nrow(s1), n.sample),],
                  s2[sample(1:nrow(s2), n.sample),])
  
      indTab[i,3] <- suppressWarnings(birdring::overlap(d1, d2, from = 0, to = max(c(d1, d2), na.rm = T)))
      
      if(indTab[i,3]>=threshold) {
        
        gr[min(which(gr==indTab[i,1])):max(which(gr==indTab[i,2]))] <- indTab[i,1]
        gr[(which(diff(gr)>1)+1):length(gr)] <-  gr[(which(diff(gr)>1)+1):length(gr)]-(max(diff(gr))-1)
        
        GrInd  <- unique(gr[gr%in%gr[duplicated(gr)]])
        indTab <- cbind(GrInd[-length(GrInd)], GrInd[-1], NA) 
        
      } else {
        if(i==nrow(indTab)) break
      }
    }
    
    if(all(!is.na(indTab[,3]))) break
  }    
  
  
  if(plot) {
    opar <- par(mfrow = c(2,1))
  
    gr.old <- fit$gr
    gr.old[!gr.old%in%gr.old[duplicated(gr.old)]] <- NA
    
    cols <- rainbow(max(gr.old, na.rm = T))[sample(1:max(gr.old, na.rm = T))]
    
    plot(fit$twl$Twilight,  hourOffset(as.hour(fit$twl$Twilight)), 
         type  = "n", xlab = "", ylab = "", yaxt = "n")
    r <- range(hourOffset(as.hour(fit$twl$Twilight)))
    for(i in unique(gr.old[!is.na(gr.old)])){
      rect(min(fit$twl$Twilight[!is.na(gr.old) & gr.old==i]), r[1], 
           max(fit$twl$Twilight[!is.na(gr.old) & gr.old==i]), r[2], lty = 2, col = "grey90", border = "grey40") 
    }
    points(fit$twl$Twilight,  
           hourOffset(as.hour(fit$twl$Twilight)), 
           pch = 21, bg = cols[gr.old], yaxt = "n")
    
    
    
    gr.new <- gr
    gr.new[!gr.new%in%gr.new[duplicated(gr.new)]] <- NA
    
    
    plot(fit$twl$Twilight,  hourOffset(as.hour(fit$twl$Twilight)), 
         type  = "n", xlab = "", ylab = "", yaxt = "n")
    r <- range(hourOffset(as.hour(fit$twl$Twilight)))
    for(i in unique(gr.new[!is.na(gr.new)])){
      rect(min(fit$twl$Twilight[!is.na(gr.new) & gr.new==i]), r[1], 
           max(fit$twl$Twilight[!is.na(gr.new) & gr.new==i]), r[2], lty = 2, col = "grey90", border = "grey40") 
    }
    points(fit$twl$Twilight,  
           hourOffset(as.hour(fit$twl$Twilight)), 
           pch = 21, bg = cols[gr.new], yaxt = "n")
    
    par(opar)
  }    


list(gr = gr, ovlp = indTab)  
  
}
