Sys.setenv("TF_NUM_INTEROP_THREADS"=12)
Sys.setenv("TF_NUM_INTRAOP_THREADS"=12)
Sys.setenv("ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"=12)
################################################################################
library( ANTsRNet )
library( ANTsR )
library( patchMatchR )
library( tensorflow )
library( keras )
library( reticulate )
library( ggplot2 )
library( tfdatasets )
################################################################################
predicted2segmentation <- function( xList, domainImageList, nClasses ) {
  refDomainImage = domainImageList[[1]]
  refDomainImage1chan = splitChannels( refDomainImage )[[1]]
  nx = length( xList )
  avgProbImage = refDomainImage1chan * 0
  for ( nmodels in 1:nx ) {
    x = as.array( tf$squeeze( xList[[nmodels]] ) )
    xdim = dim( x )
    nvoxels = prod( head(xdim,2) )
    probImage = makeImage( splitChannels(domainImageList[[nmodels]])[[1]]*0+1, x[,] ) %>%
      iMath("Normalize")
    avgProbImage = avgProbImage + resampleImageToTarget( probImage, refDomainImage1chan ) / nmodels
  }
  return( avgProbImage )
}

plotColor <- function(imgList, scale=TRUE, vectors=NULL, points=NULL, paths=NULL) {
  
  if (class(imgList) == "antsImage") {
    imgList = list(imgList, imgList, imgList)
  }
  
  direction = antsGetDirection( imgList[[1]] )
  
  # max in all images
  maxi = 1.0
  if ( scale )
  {
    maxi = max( unlist( lapply( imgList, function(x) { max(x) } ) ) )
  }
  
  rgbList = lapply( imgList, function(x) { apply(t(as.matrix(x)),2,rev) / maxi })
  
  col <- rgb(rgbList[[1]], rgbList[[2]], rgbList[[3]])
  
  d = dim(rgbList[[1]])
  
  x = rep(1:d[2],each=d[1])
  y = rep(1:d[1], d[2])
  pts = antsTransformIndexToPhysicalPoint( imgList[[1]], cbind(x,y) )
  
  dat = data.frame(x=pts[,1], y=pts[,2], col=col)
  x1 = min(pts[,1])
  x2 = max(pts[,1])
  y1 = min(pts[,2])
  y2 = max(pts[,2])
  
  g = ggplot(dat) + geom_raster(aes(x=x, y=y, fill=col), hjust=0, vjust=0, alpha=1) + theme(legend.position="none", aspect.ratio=1,text=element_blank(),axis.ticks=element_blank(), panel.grid=element_blank() ) + scale_fill_manual(values=as.character(levels(factor(col))) )
  
  g = g + coord_cartesian( xlim=c(x1,x2), ylim=c(y1,y2) )
  if ( direction[1,1] > 0 ) {
    g = g + scale_x_continuous( lim=c(x1,x2) )
  }
  else {
    g = g + scale_x_reverse( lim=c(x2,x1) )
  }
  if ( direction[2,2] > 0 ) {
    g = g + scale_y_continuous( lim=c(y1,y2) )
  }
  else {
    g = g + scale_y_reverse( lim=c(y2,y1) )
  }
  
  if ( !is.null(points) ) {
    pdat = data.frame( x=points[,1], y=points[,2], id=factor(1:dim(points)[1]) )
    g = g + geom_point( data=pdat, aes(x=x, y=y, colour=id ))
  }
  
  if ( !is.null(paths) ) {
    g = g + geom_path(data=paths, aes(x=x,y=y,group=id,colour=id))
  }
  
  if ( !is.null(vectors) ) {
    xvec = as.vector( t(as.matrix(vectors[[1]])) )
    yvec = as.vector( -t(as.matrix(vectors[[2]])) )
    vpts = antsTransformIndexToPhysicalPoint( imgList[[1]], cbind(x+0.5,y+0.5) )
    
    mag = sqrt(xvec*xvec + yvec*yvec)
    elim = which(mag < 0.01)
    if (length(elim) > 0 ) {
      xvec = xvec[-elim]
      yvec = yvec[-elim]
      vpts = vpts[-elim,]
    }
    vdat = data.frame(x=vpts[,1]-xvec, y=vpts[,2]-yvec, xend=vpts[,1]+xvec, yend=vpts[,2]+yvec)
    g = g + geom_segment(data=vdat, aes(x=x,y=y,xend=xend,yend=yend), colour="white", alpha=0.5)
  }
  
  suppressWarnings(print(g))
}

read.fcsv<-function( x, skip=3 ) {
  df = read.table( x, skip=skip, sep=',' )
  colnames( df ) = c("id","x","y","z","ow","ox","oy","oz","vis","sel","lock","label","desc","associatedNodeID")
  #  df$y = df$y * (-1)
  return( df )
}

toLuminance <- function( x ) {
  antsAverageImages( splitChannels( x ), verbose=FALSE )
}

polarX <- function(X) {
  x_svd <- svd(X)
  P <- x_svd$u %*% diag(x_svd$d) %*% t(x_svd$u)
  Z <- x_svd$u %*% t(x_svd$v)
  if (det(Z) < 0)
    Z = Z * (-1)
  return(list(P = P, Z = Z, Xtilde = P %*% Z))
}

randAff <- function( loctx,  txtype = "Rigid", sdAffine,
                     idparams, fixParams, seeder ) {
  idim = 2
  set.seed( seeder )
  noisemat = stats::rnorm(length(idparams), mean = 0, sd = sdAffine)
  if (txtype == "Translation")
    noisemat[1:(length(idparams) - idim )] = 0
  idparams = idparams + noisemat
  idmat = matrix(idparams[1:(length(idparams) - idim )],
                 ncol = idim )
  idmat = polarX(idmat)
  if (txtype == "Rigid")
    idmat = idmat$Z
  if (txtype == "Affine")
    idmat = idmat$Xtilde
  if (txtype == "ScaleShear")
    idmat = idmat$P
  if ( rnorm(1,0,1) < 0 ) { # controls frequency of flipping
    flipper = diag( 2 )
    flipper[1,1] = -1
    idmat = idmat %*% flipper
  }
  if ( rnorm(1,0,1) < 0 & FALSE ) { # controls frequency of flipping
    flipper = diag( 2 )
    flipper[2,2] = -1
    idmat = idmat %*% flipper
  }
  idparams[1:(length(idparams) - idim )] = as.numeric(idmat)
  setAntsrTransformParameters(loctx, idparams)
  setAntsrTransformFixedParameters( loctx, fixParams )
  return(loctx)
}

randomRotateImage <- function( image, sdAff=0.1, seeder ) {
  fixedParams = getCenterOfMass( image * 0 + 1 )
  loctx <- createAntsrTransform(precision = "float",
                                type = "AffineTransform", dimension = image@dimension  )
  setAntsrTransformFixedParameters(loctx, fixedParams)
  idparams = getAntsrTransformParameters( loctx )
  setAntsrTransformParameters( loctx, idparams )
  setAntsrTransformFixedParameters(loctx, fixedParams)
  loctx = randAff( loctx, sdAffine=sdAff, txtype = 'Rigid',
                   idparams = idparams, fixParams = fixedParams, seeder = seeder )
  imageR = applyAntsrTransformToImage( loctx, image, image,
                                       interpolation = "nearestNeighbor" )
  return( list( imageR, loctx ) )
}

generateDataBinSeg <- function( imgIn, segger,
                                batch_size = 16, mySdAff=0.15, visualize = FALSE,
                                inference=FALSE, ychan=1 ) {
  # choice of this parameter can have a strong effect on outcome
  # below, we use a guess at an automated scaling approach
  aboutthelow = 512
  if ( ! inference ) aboutthelow = sample( c( 256, 512, 128, 192 ), 1)
  mySubSam = ( tail( dim( imgIn ), 1 ) / aboutthelow )
  subSampling = round(mySubSam)
  imgSub = resampleImage( imgIn, dim( imgIn )/subSampling, useVoxels=T)
  temp = splitChannels( imgSub )[1:3]
  nChannels = length( temp )
  myLum = toLuminance( imgSub )
  for ( k in 1:nChannels ) {
    temp[[k]] = ANTsRNet::padImageByFactor(temp[[k]], 16 )
  }
  imgSub = mergeChannels( temp )
  bodySeg = resampleImageToTarget( segger, temp[[1]], interpType='linear' )
  X = array( dim = c( batch_size, dim( imgSub  ), nChannels ) )
  Y = array( dim = c( batch_size, dim( imgSub  ), ychan ) )
  for ( k in 1:batch_size ) {
    myseed = as.integer( sample(1:10000000,1) )
    splitter = splitChannels( imgSub )
    for ( j in 1:nChannels ) {
      if ( !inference ) {
        rr = randomRotateImage( splitter[[j]], sdAff=mySdAff, seeder = myseed )
        loctx = rr[[2]]
        if ( j == 1 ) {
          bodySegRot = applyAntsrTransformToImage( loctx, bodySeg, bodySeg, interpolation='linear' )
        }
        rr = iMath( rr[[1]], "Normalize" )
      } else {
        rr = iMath( splitter[[j]], "Normalize" )
        bodySegRot = bodySeg
      }
      X[k,,,j] = as.array( rr  )
    }
    if ( ychan > 1 ) Y[k,,,2] = as.array( 1.0 - bodySegRot )
    Y[k,,,1] = as.array( bodySegRot )
  }
  list( X=X-0.5, Y=Y, imgSub=imgSub )
}

fishInferenceSeg<- function( imgIn, mdl, visualize=FALSE ) {
  tL = toLuminance
  temp = splitChannels(imgIn)[[1]]*0+1
  gg = generateDataBinSeg( imgIn, temp, batch_size=1, mySdAff=0, inference=TRUE )
  predictions = tf$nn$sigmoid( mdl( gg$X ) )
  heatP = as.antsImage( as.array(predictions[1,,,1] ) )
  temp2 = splitChannels(gg$imgSub)[[1]]*0+1
  heatP = antsCopyImageInfo2( heatP, temp2 ) %>% iMath("Normalize")
  flll = thresholdImage(heatP ,0.5,Inf)  %>% morphology("close",2) %>%
    iMath("FillHoles") %>% labelClusters(50)
  if ( visualize & var(heatP) > 1e-3  )  {
    layout(matrix(1:2,nrow=1))
    plot( tL(gg$imgSub), iMath(heatP,"Normalize") * 255, alpha=0.5, window.overlay=c(24,255))
    plot( tL(gg$imgSub), flll, alpha=0.5 )
  }
  heatP = resampleImageToTarget( heatP, temp )
  flll = resampleImageToTarget( flll, temp, interpType='nearestNeighbor' )
  return( list( softThreshold=heatP, hardThreshold=flll ) )
}


generateData <- function( imgSub, segSub, ptsSub,
                          batch_size = 16, mySdAff=0.15, visualize = FALSE,
                          subSampling, inference=FALSE, species, allSpecies )
{
  nPoints = nrow(ptsSub)
  whichSpecies = which( species == allSpecies  )
  speciesBin = rep( 0, length( allSpecies ) )
  speciesBin[ whichSpecies ] = 1
  
  # calculate reorientation first
  reofi = reorientImage( segSub, c( 1, 0 ) )
  temp = splitChannels( imgSub )
  mytx = c( reofi$txfn )
  segSub = antsApplyTransforms( segSub, segSub, mytx, interpolator='nearestNeighbor' )
  segCrop = cropImage( segSub, segSub )
  temp = splitChannels( imgSub )[1:3]
  nChannels = length( temp )
  for ( k in 1:nChannels ) {
    temp[[k]] = antsApplyTransforms( temp[[k]], temp[[k]], mytx )
    temp[[k]] = cropImage( temp[[k]], segSub ) 
    lodim = dim( temp[[k]] )[2]/128
    if ( lodim > 0 ) {
      targetdim = round(  dim( temp[[k]] )/lodim  )
       temp[[k]] = resampleImage( temp[[k]], targetdim, useVoxels = TRUE )
      } 
    temp[[k]] = ANTsRNet::padImageByFactor( temp[[k]], 16 )
    }

  ptsSub = antsApplyTransformsToPoints( 2, ptsSub, rev(mytx), whichtoinvert=c(TRUE,TRUE) )
  imgSub = mergeChannels( temp )
  print(dim(imgSub))
  X = array( dim = c( batch_size, dim( imgSub  ), nChannels ) )
  Xm = array( dim = c( batch_size, dim( imgSub  ), nPoints ) )
  mycc = array( dim = c( batch_size, dim( imgSub  ), imgSub@dimension ) )
  mymasks = array( dim = c( batch_size, dim( imgSub  ), 1 ) )
  ptsrot = array( dim = c( batch_size, nPoints, imgSub@dimension ) )
  mySp = array( 0, dim = c( batch_size, dim( imgSub  ), length( allSpecies ) ) )
  for ( k in 1:batch_size ) {
    myseed = Sys.time()
    splitter = splitChannels( imgSub )
    for ( j in 1:nChannels ) {
      if ( !inference ) {
        rr = randomRotateImage( splitter[[j]], sdAff=mySdAff, seeder = myseed )
        loctx = rr[[2]]
        rr = iMath( rr[[1]], "Normalize" )
      } else {
        rr = iMath( splitter[[j]], "Normalize" )
      }
      if ( j == 1 ) { # get coord conv results and point results
        coordconver = coordinateImages( rr*0+1 )
        mycc[k,,,1] = as.array( coordconver[[1]] )
        mycc[k,,,2] = as.array( coordconver[[2]] )
        if ( !inference ) {
          loctxInv = invertAntsrTransform( loctx )
          pointsR = applyAntsrTransform( loctxInv, data.matrix(ptsSub), dataType = 'point' )
          ptsrot[k,,] = data.matrix( pointsR )
        } else {
          ptsrot[k,,] = data.matrix( ptsSub )
        }
      }
      if ( imgSub@dimension == 2 ) {
        X[k,,,j] = as.array( rr  )
        mySp[k,,,whichSpecies] = mySp[k,,,whichSpecies] + as.array( rr )/3.0
      }
      if ( imgSub@dimension == 3 ) {
        X[k,,,,j] = as.array( rr  )
        mySp[k,,,,whichSpecies] = mySp[k,,,,whichSpecies] + as.array( rr  )/3.0
      }
    }
    if ( ! inference )
      for ( pp in 1:nPoints ) {
        temp = thresholdImage( segSub, pp, pp )
        temp = smoothImage( temp, 10.0 ) %>% iMath("Normalize" )
        temp = randomRotateImage( temp, sdAff=mySdAff, seeder = myseed  )[[1]]
        if ( visualize & max(temp) > 0 & FALSE ) {
          plot( rr * 255, temp * 255, doCropping=FALSE )
        }
        if ( var( temp ) == 0 ) message(paste("WARN:",pp,"VAR0"))
        if ( imgSub@dimension == 2 ) Xm[k,,,pp]    = as.array( temp  )
        if ( imgSub@dimension == 3 ) Xm[k,,,,pp]   = as.array( temp  )
      }
  }
  # gg$ximg, gg$cc, gg$hout, bigger$yptRot
  list( ximg=X,  cc=mycc, hout=Xm, yptRot = ptsrot,
        masks=mymasks, imgSub=imgSub, bodySeg=NA,
        ptsSub=ptsSub, reorientation=reofi, species=mySp )
}


fishInferenceLM <- function( imgIn, segIn, mdl, truePoints, heatThresh = 0.65,
                             visualize=FALSE, species, allSpecies ) {
  tL = toLuminance
  mpi = NULL
  if ( missing( truePoints ) ) haveGT=FALSE else haveGT=TRUE
  if ( !haveGT ) truePoints = matrix( rnorm(24*2), nrow=24 )
  gg = generateData( imgIn, segIn, truePoints, batch_size=1, mySdAff=0,
                     inference=TRUE, species=species, allSpecies=allSpecies )
  truPoints = gg$ptsSub
  if ( length( mdl$inputs ) == 3 )
    predictions = mdl( list( gg$ximg, gg$species, gg$cc ) )
  if ( length( mdl$inputs ) == 2 )
    predictions = mdl( list( gg$ximg, gg$cc ) )
  guessPoints = as.matrix( as.array( predictions[[2]] )[1,,] )
  guessPoints2 = guessPoints * 0
  heatP = as.array( predictions[[1]] )
  coords = patchMatchR::coordinateImages( toLuminance( gg$imgSub )*0+1 )
  errs = errsRaw = rep( NA, tail(dim(heatP),1) )
  heatList = list()
  heatArr = array( dim=c(1,dim(gg$imgSub),24 ) )
  for ( k in 1:tail(dim(heatP),1) ) {
    heatList[[k]] = as.antsImage( heatP[1,,,k]  ) %>% antsCopyImageInfo2( gg$imgSub )
    heatList[[k]] = iMath( as.antsImage( heatP[1,,,k]  ), "Normalize" ) %>% antsCopyImageInfo2( gg$imgSub )
    temp = thresholdImage( heatList[[k]], heatThresh, Inf )
    heatArr[1,,,k] = as.array( heatList[[k]] * temp )
    #    plot( toLuminance(gg$imgSub), heatList[[k]]*255, doCropping=F, window.overlay=c(5,255), alpha=0.5 )
    selheat = heatList[[k]] > quantile( heatList[[k]], heatThresh )[1]
    selheat = heatList[[k]] >  heatThresh
    #    selheat = temp == 1
    xcoordvec = coords[[1]][ selheat ]
    xcoordwt = heatList[[k]][ selheat ]
    xcoord = sum(xcoordvec*xcoordwt/sum(xcoordwt))
    ycoordvec = coords[[2]][ selheat ]
    ycoordwt = heatList[[k]][ selheat ]
    ycoord = sum(ycoordvec*ycoordwt/sum(ycoordwt))
    guessPoints2[k,] = c( xcoord, ycoord )
    if ( haveGT )
      errsRaw[k] = mean(abs(as.numeric(truPoints[k,])-guessPoints[k,] ))
    errs[k] = mean(abs(as.numeric(truPoints[k,])-guessPoints2[k,] ))
  }
  if ( haveGT ) {
    print(paste("MeanErr:",mean(errsRaw),"MeanErrPost-Hoc:",mean(errs)))
  }
  if ( visualize & haveGT ) {
    mpiTr = makePointsImage( truPoints, toLuminance(gg$imgSub)*0+1, radius = 3 )
    mpi = makePointsImage( guessPoints2, toLuminance(gg$imgSub)*0+1, radius = 3 )
    mpiB = makePointsImage( guessPoints, toLuminance(gg$imgSub)*0+1, radius = 3 )
    layout( matrix(1:3,nrow=1))
    plot( toLuminance(gg$imgSub), mpiTr, doCropping=FALSE )
    plot( toLuminance(gg$imgSub), mpi, doCropping=FALSE  )
    plot( toLuminance(gg$imgSub), mpiB, doCropping=FALSE  )
  }
  if ( visualize & ! haveGT ) {
    mpi = makePointsImage( guessPoints2, toLuminance(gg$imgSub)*0+1, radius = 3 )
    mpiB = makePointsImage( guessPoints, toLuminance(gg$imgSub)*0+1, radius = 3 )
    layout( matrix(1:3,nrow=1))
    plot( toLuminance(gg$imgSub), doCropping=FALSE )
    plot( toLuminance(gg$imgSub), mpi, doCropping=FALSE  )
    plot( toLuminance(gg$imgSub), mpiB, doCropping=FALSE  )
  }
  return( list( pointImage = mpi, points = guessPoints2, img=gg$imgSub, heatArr=heatArr ) )
  #  return( )
  # gg$ximg, gg$cc, gg$hout, bigger$yptRot
  #  list( ximg=X,  cc=mycc, masks=mymasks,subSampling=subSampling )
}

groupVar='Family'
demog=read.csv( "9nov20/LM_training_all_GLIN_9nov20.csv" )
names(demog) = c("id"  ,    "Family",  "Genu"  ,  "species", "Nothing")
segmdl = load_model_hdf5( "models/binary_fish_seg.h5" )
lmMdl = load_model_hdf5( "models/species_guided_unet_pretrain.h5" )
# in theory, we could glue these together and do full end-to-end training
# with appropriate augmentation to improve results
unseenDemog = read.csv("unseen/unseen.csv")
unseens = sample( c(
  Sys.glob("unseen/images/*JPG"),  Sys.glob("unseen/images/*jpg") ), replace=F )
unseens = sample(  Sys.glob("unseen/images/*lat*jpg"), replace=F )
for ( fn in unseens ) {
  identifier = basename( fn ) %>% tools::file_path_sans_ext()
  selector = grep( identifier, unseenDemog$id )
  if ( length( selector ) == 1 ) {
    localSpecies = unseenDemog$Family[selector]
    print( paste( fn, localSpecies ) )
    uimg = antsImageRead( fn )
    segger = fishInferenceSeg( uimg, segmdl, FALSE )$hardThreshold
    biggestseg = iMath( segger, "GetLargestComponent" )
    usplit = splitChannels(uimg)
    for ( u in 1:3 )
      usplit[[u]] = cropImage( usplit[[u]], biggestseg )
    plotColor(usplit)
    finf = fishInferenceLM(
      imgIn = uimg,
      segIn = biggestseg,
      mdl = lmMdl,
      heatThresh = 0.80,
      species = localSpecies,
      allSpecies = unique( demog[,groupVar] ),
      visualize = TRUE )
    
    Sys.sleep( 5 )
  } else  print( paste("Cannot find unique", fn, "via identifier", identifier, "in demog" ) )
}
