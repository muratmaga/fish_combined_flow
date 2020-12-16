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
  probList = list()
  for ( nc in 1:nClasses ) {
    avgProbImage = refDomainImage1chan * 0
    for ( nmodels in 1:nx ) {
      x = as.array( tf$squeeze( xList[[nmodels]] ) )
      xdim = dim( x )
      nvoxels = prod( head(xdim,2) )
      probImage = makeImage( splitChannels(domainImageList[[nmodels]])[[1]]*0+1, x[,,nc] )
      avgProbImage = avgProbImage + resampleImageToTarget( probImage, refDomainImage1chan ) / nmodels
    }
    probList[[nc]] = avgProbImage
  }
  nvoxels = prod( dim( refDomainImage1chan ) )
  pmat = imageListToMatrix( probList, refDomainImage1chan * 0 + 1 )
  segvec = apply( pmat, MARGIN=2, FUN=which.max )
  seg = makeImage( refDomainImage1chan*0+1, segvec )
  return( seg )
}


mvkmFishIso <-function( x, initSeg, nc=3, off = 2 ) {
  xx = splitChannels( x )
  for ( k in 1:length( xx ) ) {
    xx[[k]] = n3BiasFieldCorrection( xx[[k]] * initSeg, 4 ) %>% denoiseImage()
  }
  mk = xx[[1]] * 0 + 1
  cl <- kmeans( cbind( xx[[1]][mk==1],  xx[[2]][mk==1],  xx[[3]][mk==1] ) ,  nc,
                nstart = 10, algorithm='Hartigan-Wong' )
  segimg = makeImage( mk, cl$cluster )
  segimg2 = makeImage( mk, cl$cluster )
  segimg2[1:4,1:4]=0
  plot( toLuminance( x ), segimg2, window.overlay=c(1,nc) )
  rawtbl=(table(cl$cluster))
  mxdimx = head( dim(segimg),1 )
  mxdimy = tail( dim(segimg),1 )
  # corners are background and so is the largest component
  bkglabs = c( as.numeric(names(rawtbl)[which.max(rawtbl)]),
               segimg[1:off,1:off],
               segimg[(mxdimx-off):(mxdimx),1:off],
               segimg[1:off,(mxdimy-off):(mxdimy)],
               segimg[(mxdimx-off):(mxdimx),(mxdimy-off):(mxdimy)] )
  mytbl = table( bkglabs )
  print( mytbl )
  for ( bk in sort(unique(bkglabs)) )
    segimg[ segimg == bk ] = 0
  segimg = thresholdImage( segimg, 1, Inf ) %>% iMath("GD",1) %>%
    labelClusters( minClusterSize = 25 )
  maskvals = sort( unique( segimg[ initSeg == 1 ] ) )
  temp = maskImage(segimg,segimg,maskvals[maskvals>0], binarize=TRUE) %>% iMath("FillHoles")
  finalseg = iMath( temp, "GetLargestComponent") %>% morphology("close",1)
  return( finalseg  )
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
  loctx = randAff( loctx, sdAffine=sdAff, txtype = 'Affine',
                   idparams = idparams, fixParams = fixedParams, seeder = seeder )
  imageR = applyAntsrTransformToImage( loctx, image, image,
                                       interpolation = "nearestNeighbor" )
  return( list( imageR, loctx ) )
}


generateData <- function( imgIn, ptsIn,
                          batch_size = 16, mySdAff=0.15, visualize = FALSE,
                          subSampling, inference=FALSE, species, allSpecies ) {
  whichSpecies = which( species == allSpecies  )
  speciesBin = rep( 0, length( allSpecies ) )
  speciesBin[ whichSpecies ] = 1
  # choice of this parameter can have a strong effect on outcome
  # below, we use a guess at an automated scaling approach
  
  mySubSam = ( tail( dim( imgIn ), 1 ) / 192 )
  if (mySubSam < 1 ) mySubSam = 1
  subSampling = round(mySubSam)
  print( dim( imgIn )/subSampling )
  imgSub = resampleImage( imgIn, dim( imgIn )/subSampling, useVoxels=T)
  
  antsSetSpacing( imgSub, c(1,1) )
  temp = splitChannels( imgSub )[1:3]
  nChannels = length( temp )
  myLum = toLuminance( imgSub )
  ptsSub = ptsIn/subSampling
  segger = makePointsImage( ptsSub, myLum*0+1, radius = 3 )
  if ( visualize ) {
    plot( myLum, segger )
    print( dim(segger) )
  }
  nPoints = nrow(ptsIn)
  for ( k in 1:nChannels ) {
    temp[[k]] = ANTsRNet::padImageByFactor(temp[[k]], 16 )
  }
  imgSub = mergeChannels( temp )
  reofi = reorientImage( antsAverageImages( splitChannels( imgSub ), verbose=FALSE ) , c( 1, 0 ) )
  temp = splitChannels( imgSub )
  mytx = c( reofi$txfn )
  for ( k in 1:nChannels ) {
    temp[[k]] = antsApplyTransforms( temp[[k]], temp[[k]], mytx )
  }
  imgSub = mergeChannels( temp )
  segSub = antsApplyTransforms( temp[[1]], segger, mytx, interpolator='nearestNeighbor' )
  ptsSub = antsApplyTransformsToPoints( 2, ptsSub, rev(mytx), whichtoinvert=c(TRUE,TRUE) )
  if ( visualize ) {
    plot( toLuminance( imgSub ), segSub )
    plot( toLuminance( imgSub ), makePointsImage( ptsSub, toLuminance( imgSub )*0+1, radius = 3 ) )
  }
  X = array( dim = c( batch_size, dim( imgSub  ), nChannels ) )
  kMeansK = nrow(ptsIn)
  Xm = array( dim = c( batch_size, dim( imgSub  ), kMeansK ) )
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
      for ( pp in 1:kMeansK ) {
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
        masks=mymasks,
        subSampling=subSampling, imgSub=imgSub, bodySeg=NA,
        ptsSub=ptsSub, reorientation=reofi, species=mySp )
}


fishInference <- function( imgIn, mdl, truePoints, heatThresh = 0.65,
                           visualize=FALSE, species, allSpecies ) {
  tL = toLuminance
  mpi = NULL
  if ( missing( truePoints ) ) haveGT=FALSE else haveGT=TRUE
  if ( !haveGT ) truePoints = matrix( rnorm(24*2), nrow=24 )
  gg = generateData( imgIn, truePoints, batch_size=1, mySdAff=0,
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
    selheat = heatList[[k]] > quantile( heatList[[k]], heatThresh )[1]
    selheat = heatList[[k]] >  heatThresh
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
}

groupVar='Family'
demog=read.csv( "9nov20/LM_training_all_GLIN_9nov20.csv" )
segmdl = load_model_hdf5( "models/unet_seg_dec_2020.h5" )
lmMdl = load_model_hdf5( "models/species_guided_unet_pretrain.h5" )
# with appropriate augmentation to improve results
unseenDemog = demog # read.csv("unseen/unseen.csv")
stop("Get all the raw fish data and run it through here to get training segmentations")
unseens = c( Sys.glob("9nov20/images/*JPG"),  Sys.glob("9nov20/images/*jpg") )
for ( fn in unseens ) {
  identifier = basename( fn ) %>% tools::file_path_sans_ext()
  selector = grep( identifier, unseenDemog$id )
  ofn = paste0( "9nov20/segmentation/",identifier,"_seg.nii.gz" )
  if ( length( selector ) == 1 & ! file.exists( ofn ) ) {
    localSpecies = unseenDemog$Family[selector]
    print( paste( fn, localSpecies ) )
    uimg = antsImageRead( fn )
    lodim = c( 256, 256 ) # should train this way
    lodim2 = c( 384, 384 ) # not this way
    lodim3 = c( 512, 512 ) # not this way
    subSampling = round( tail( dim( uimg ), 1 ) / 256 )
    if( subSampling < 1) subSampling = 1
    reduced_img = resampleImage( uimg, dim( uimg )/subSampling, useVoxels=T)
    print( reduced_img )
    #    reduced_img = resampleImage( uimg, lodim, useVoxels = TRUE, interpType = 'linear')
    reduced_img2 = resampleImage( uimg, lodim2, useVoxels = TRUE, interpType = 'linear')
    reduced_img3 = resampleImage( uimg, lodim3, useVoxels = TRUE, interpType = 'linear')
    splitter = splitChannels( reduced_img )
    splitter2 = splitChannels( reduced_img2 )
    splitter3 = splitChannels( reduced_img3 )
    for ( j in 1:3 ) { # channels
      splitter[[j]] = ANTsRNet::padImageByFactor( splitter[[j]], 16 )
      splitter2[[j]] = ANTsRNet::padImageByFactor( splitter2[[j]], 16 )
      splitter3[[j]] = ANTsRNet::padImageByFactor( splitter3[[j]], 16 )
    }
    X <- array( data = NA, dim = c( 1, dim( splitter[[j]] ), 3 ) )
    X2 <- array( data = NA, dim = c( 1, dim( splitter2[[j]] ), 3 ) )
    X3 <- array( data = NA, dim = c( 1, dim( splitter3[[j]] ), 3 ) )
    for ( j in 1:3 ) { # channels
      X[1,,, j] <- as.array( splitter[[j]] ) #populate with images
      X2[1,,, j] <- as.array( splitter2[[j]] ) #populate with images
      X3[1,,, j] <- as.array( splitter3[[j]] ) #populate with images
    }
    reduced_img = mergeChannels( splitter )
    reduced_img2 = mergeChannels( splitter2 )
    reduced_img3 = mergeChannels( splitter3 )
    seg = predicted2segmentation(
      ( list( segmdl(X), segmdl(X2), segmdl(X3) )),
      ( list( reduced_img, reduced_img2, reduced_img3 )), 4 )
    if ( TRUE ) {
      layout(matrix(1:3,nrow=1))
      plot( toLuminance( reduced_img ), seg, window.overlay=c(2,4), alpha=0.5)
      fishEst = thresholdImage(seg,2,2) %>% iMath("MD",2)
      if ( length( grep("lat",fn) ) == 1 ) {
        dd = dim( fishEst )
        fishEst[1:dd[1],round(0.66*dd[2]):dd[2]]=0
      }
      seg = mvkmFishIso( reduced_img, fishEst, 10, off = 5 )
      plot(toLuminance(reduced_img),seg,window.overlay=c(0.5,1.5))
    }
    usplit = splitChannels(uimg)
    segger = resampleImageToTarget( seg, usplit[[1]] )
    print( paste( ofn, "done" ))
    antsImageWrite( segger, ofn )
    if ( FALSE ) {
      for ( u in 1:3 )
        usplit[[u]] = cropImage( usplit[[u]], segger ) # %>% iMath("PadImage",16)
      finf = fishInference(
        imgIn = mergeChannels(usplit),
        mdl = lmMdl,
        heatThresh = 0.80,
        species = localSpecies,
        allSpecies = unique( demog[,groupVar] ),
        visualize = TRUE )
    }
  } else  print( paste("Cannot find unique", fn, "via identifier", identifier, "in demog" ) )
}

