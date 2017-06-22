require(xts)
require(zoo)

create_dataset <- function(season,ar_ord,years) {
  # import data
  dat <- read.table('soldat-clim.txt',header=FALSE,sep=',',colClasses='numeric')
  # convert MATLAB datenums to POSIX dates
  dat[,1] <- as.POSIXct(strftime(as.POSIXct((dat[,1]-719529)*86400,origin='1970-1-1',tz='UTC'),
                                 format='%Y-%m-%d %H:%M',tz='UTC',usetz=FALSE),tz='UTC')
  # convert data to xts object
  dat <- xts(dat[,-1],order.by=dat[,1])
  indexTZ(dat) <- 'EST'
  # calculate kt and kb by averaging GHI and GNI over hours
  yearind <- years-1900
  dat <- dat[.indexyear(dat) %in% yearind,]
  dat.agg <- period.apply(dat[,c(3,4,5)],endpoints(dat,'hours'),colMeans)
  dat.all <- merge(dat.agg,dat[,2],join='left')
  names(dat.all) <- c('DNI','GHI','G0','cosz')
  dat.all$kt <- dat.all$GHI/dat.all$G0
  dat.all$kb <- (dat.all$DNI*dat.all$cosz)/dat.all$G0
  # do some quality control for data that is erronous
  dat.all <- dat.all[dat.all[,'G0'] >= 0,]
  dat.all[dat.all[,'kt'] < 0 & dat.all[,'kt'] > -0.1,'kt'] <- 0
  dat.all[dat.all[,'kb'] < 0 & dat.all[,'kb'] > -0.1,'kb'] <- 0
  dat.all[dat.all[,'kt'] > 1,'kt'] <- NA
  dat.all[dat.all[,'kb'] > 1,'kb'] <- NA
  
  # set months for data based on season
  if (season == "summer") {
    months <- c(5,6,7)
    max.length <- 24
  } else if (season == "winter") {
    months <- c(11,0,1)
    max.length <- 18
  } else if (season == "spring") {
    months <- c(2,3,4)
    max.length <- 22
  } else if (season == "fall") {
    months <- c(8,9,10)
    max.length <- 20
  } else {
    months <- 0:11
    max.length <- 20
  }
  
  # construct observation matrix y
  y <- sapply(
    lapply(
      lapply(split(dat.all[.indexmon(dat.all) %in% months,c('kt','kb')],'days'),t)
      ,matrix,ncol=1),
    'length<-',max(lengths(split(dat.all[.indexmon(dat.all) %in% months,c('kt','kb')],'days'))))
  
  # gather cosine of zenith angle
  cosz <- lapply(split(dat.all[.indexmon(dat.all) %in% months,'cosz'],'days'),'as.numeric')
                 
  cosz <- sapply(cosz,'length<-',max(lengths(cosz)))
  
  # get rid of days with NAs mid-day (this needs to be dealt with later as missing data)
  cosz <- cosz[1:(max.length/2),which(apply(y[1:max.length,],2,function(x) all(!is.na(x))))]
  
  y <- y[1:max.length,which(apply(y[1:max.length,],2,function(x) all(!is.na(x))))]
  
  # move data into interior (0,1) and apply logit transform
  eps <- 1e-3
  y[which(y < 1e-3)] <- eps
  y[which(y > (1-1e-3))] <- 1-eps
  y <- log(y/(1-y))
  
  # prepare covariate matrix
  X <- array(0,dim=c(nrow(y)-(ar_ord*2),(ar_ord+2)*2,ncol(y)))
  for (i in 1:ncol(y)) {
    X1 <- na.omit(as.matrix(lag(zoo(y[seq(1,nrow(y),2),i]),-(1:ar_ord))))
    X2 <- na.omit(as.matrix(lag(zoo(y[seq(2,nrow(y),2),i]),-(1:ar_ord))))
    X[seq(1,2*nrow(X1),2),1:(ar_ord+2),i] <- as.matrix(cbind(1+mat.or.vec(nr=nrow(X1),nc=1),
                                                             cosz[-c(1:ar_ord*2),i],
                                                             as.matrix(X1),deparse.level=0))
    X[seq(2,2*nrow(X1),2),(ar_ord+3):((ar_ord+2)*2),i] <- as.matrix(cbind(1+mat.or.vec(nr=nrow(X1),nc=1),
                                                                          cosz[-c(1:ar_ord*2),i],
                                                                          as.matrix(X2),deparse.level=0))
  }
#  y.init <- y[c(1:(ar_ord*2)),]
  y <- y[-c(1:(ar_ord*2)),]

  return(list(y,X))
  
}