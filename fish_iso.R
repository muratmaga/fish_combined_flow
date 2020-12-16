setwd('/fish')
Sys.setenv("CUDA_VISIBLE_DEVICES"=0)
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
library(Rvcg)
################################################################################
#from Morpho
write.fcsv <- function(x,filename=dataname,description=NULL) {
  dataname <- deparse(substitute(x))
  if (!grepl("*.fcsv$",filename))
    filename <- paste0(filename,".fcsv")
  cat("# Markups fiducial file version = 4.4\n# CoordinateSystem = 0\n# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\n",file=filename)
  ptdim <- ncol(x)
  ptn <- nrow(x)
  if (ptdim == 2)
    x <- cbind(x,0)
  else if (ptdim != 3)
    stop("only 2D or 3D point clouds")
  rowids <- paste0("vtkMRMLMarkupsFiducialNode_",1:ptn)
  associatedNodeID <- rep("vtkMRMLVectorVolumeNode12",ptn)
  buffer <- rep(0,ptn)
  buffer <- cbind(buffer,0,0,1,1,1,0)
  desc <- rep("",ptn)
  names <- rownames(x)
  if (is.null(names))
    names <- paste0("F_",1:ptn)
  outframe <- data.frame(rowids,x,buffer,names,desc,associatedNodeID)
  write.table(format(outframe, scientific = F, 
                     trim = T), file = filename, sep = ",", append = TRUE, 
              quote = FALSE, row.names = FALSE, col.names = FALSE, 
              na = "")
}

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

mvkm <-function( x, nc=3, bkgThresh=200, off = 4 ) {
  xx = splitChannels( x )
  mk = xx[[1]] * 0 + 1
  cl <- kmeans( cbind( xx[[1]][mk==1],  xx[[2]][mk==1],  xx[[3]][mk==1] ) ,  nc,
                nstart = 10, algorithm='Hartigan-Wong' )
  segimg = makeImage( mk, cl$cluster )
  mxdimx = head( dim(segimg),1 )
  mxdimy = tail( dim(segimg),1 )
  bkglabs = c(
    segimg[1:off,], segimg[,1:off], segimg[(mxdimx-off):(mxdimx),], segimg[,(mxdimy-off):(mxdimy)] )
  mytbl = table( bkglabs )
  #  myord = order( table( bkglabs ), decreasing = TRUE )
  bkglabs = as.numeric( names(mytbl)[ mytbl > bkgThresh ] )
  for ( bk in bkglabs )
    segimg[ segimg == bk ] = 0
  segimg = thresholdImage( segimg, 1, Inf ) %>%
    morphology("close",5) %>%
    iMath("GetLargestComponent")
  segimg
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

#demog = read.csv( "bgnn_5museum_trainingset.csv")
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
  #  bodySeg = mvkm( imgSub, 8, 200 )
  if ( visualize ) {
    plot( myLum, segger )
    print( dim(segger) )
  }
  nPoints = nrow(ptsIn)
  for ( k in 1:nChannels ) {
    #    temp[[k]] = cropImage( temp[[k]] * bodySeg, iMath( bodySeg, "MD", 4 ) )
    temp[[k]] = # iMath(temp[[k]] , "PadImage", 64 ) %>%
      ANTsRNet::padImageByFactor(temp[[k]], 16 )
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
  #  bodySeg = antsApplyTransforms( temp[[1]], bodySeg, mytx, interpolator='nearestNeighbor' )
  ptsSub = antsApplyTransformsToPoints( 2, ptsSub, rev(mytx), whichtoinvert=c(TRUE,TRUE) )
  if ( visualize ) {
    #    plot( toLuminance( imgSub ), bodySeg )
    plot( toLuminance( imgSub ), segSub )
    plot( toLuminance( imgSub ), makePointsImage( ptsSub, toLuminance( imgSub )*0+1, radius = 3 ) )
  }
  # } else {
  #  reofi = NULL
  #  segSub = resampleImageToTarget( segger, splitChannels(imgSub)[[1]], interpType='nearestNeighbor' )
  #  bodySeg = resampleImageToTarget( bodySeg, splitChannels(imgSub)[[1]], interpType='nearestNeighbor' )
  # }
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
        #        bodySegRot = applyAntsrTransform( loctx, bodySeg, bodySeg, interpolation='nearestNeighbor' )
        rr = iMath( rr[[1]], "Normalize" )
      } else {
        rr = iMath( splitter[[j]], "Normalize" )
        #        bodySegRot = bodySeg
      }
      if ( j == 1 ) { # get coord conv results and point results
        #        mymasks[k,,,1] = as.array( bodySegRot )
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
demog=read.csv("9nov20/LM_training_all_GLIN_9nov20.csv")
# THIS WAS YOUR PROBLEM
names(demog) = c("id"  ,    "Family",  "Genu"  ,  "species", "Nothing")
segmdl = load_model_hdf5( "models/unet_seg_dec_2020.h5")
lmMdl = load_model_hdf5("models/species_guided_unet_pretrain.h5")
# in theory, we could glue these together and do full end-to-end training
# with appropriate augmentation to improve results

unseenDemog = read.csv("unseen/unseen.csv")
unseens = sample( c(
  Sys.glob("unseen/images/*JPG"),  Sys.glob("unseen/images/*jpg") ), replace=F )
#unseens = sample(  Sys.glob("unseen/images/*lat*jpg"), replace=F )
if ( any( (unique(unseenDemog[,groupVar]) %in%  unique(demog[,groupVar])) == FALSE ) )
     stop("Unseen data not in training data.")

for ( fn in unseens ) {
  identifier = basename( fn ) %>% tools::file_path_sans_ext()
  selector = grep( identifier, unseenDemog$id )
  if ( length( selector ) == 1 ) {
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
    # segger = morphology(seg,"dilate",4) %>% iMath("GetLargestComponent")
    for ( u in 1:3 )
      #      usplit[[u]] = cropImage( usplit[[u]], segger )
      usplit[[u]] = cropImage( usplit[[u]], segger ) # %>% iMath("PadImage",16)
    plotColor(usplit)
    finf = fishInference(
      imgIn = mergeChannels(usplit),
      mdl = lmMdl,
      heatThresh = 0.80,
      species = localSpecies,
      allSpecies = unique( demog[,groupVar] ),
      visualize = TRUE )
    write.csv(finf$points, paste0('/fish/predicted_landmarks/', identifier, '.csv'),row.names = F)
    
    Sys.sleep( 5 )
  } else  print( paste("Cannot find unique", fn, "via identifier", identifier, "in demog" ) )
}
