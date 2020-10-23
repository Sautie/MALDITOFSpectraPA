library("MALDIquant")
library("MALDIquantForeign")
library("stringr")
library("dplyr")

#' @Description Detect_Peaks is a pipeline with two stages: 1) import and processing of Maldi_Tof spectra and 2) peak detection 
#'              details: 1) import, check of quality, transformation, smoothing, baseline removing, normalization, and alignment of Maldi_Tof spectra
#'              details: 2) peak detection, binning, filtering and merging with metadata ("Filered_Meta.csv")
#' @param       Dp: folder containing the ".mzXML" files, filein: metadata input csv file, fileout: output csv file, 
#' @param       keyCode: code for dataframe joining (default:"Identifiant_MALDI"), 
#' @param       speColumn: name of isolate column ("Taxonomie"),
#' @param       t: spectrum transformation, smooth: smoothing method, baseline: bs removing method, normalization: normalization algorithm
#' @param       Iter: iterations for baseline removing (100) SN_R:signal_to_noise ratio (2), 
#' @param       minFreq: the minimum peak frequency for spectrum selection (0.25), 
#' @param       align: boolean Boolean parameter to indicate whether or not spectrum alignment is performed 
#' @return      merged dataframes and csv file
#' @examples  df_Peaks<-Detect_Peaks("newdata", "Filered_Meta.csv", "Filered_Meta_Peaks.csv"), 
#'            df_Peaks_0<-Detect_Peaks("newdata", "Filered_Meta.csv", "Filered_Meta_Peaks_0.csv", minFreq=0, align=FALSE)
#'            df_PeaksTH<-Detect_Peaks("newdata", "Filered_Meta.csv", "Filered_Meta_Peaks.csv", baseline="TopHat")
Detect_Peaks <- function(Dp, filein, fileout, keyCode="Identifiant_MALDI",speColumn="Taxonomie",
                           t="sqrt", smooth="SavitzkyGolay", baseline="SNIP", normalization="TIC", 
                           Iter=100, SN_R=2, minFreq=0.25, align=TRUE){

if (!dir.exists(Dp))  { warning("this directory does not exist in the current working directory") }
  else {
    print("spectra importing...")
    spectra<- importMzXml(Dp, verbose=FALSE) 
    print("spectra importing finished...")
    
  if (any(sapply(spectra, isEmpty)))  { warning("at least one empty spectrum")  }
  if (!all(sapply(spectra, isRegular))) { warning("irregular spectra")  }
  spectra <- trim(spectra)
  
  if (t=="sqrt")
    spectra <- transformIntensity(spectra, method="sqrt")
  else if (t=="log")
    spectra <- transformIntensity(spectra, method="log")
  
  if (smooth=="SavitzkyGolay")
      spectra <- smoothIntensity(spectra, method="SavitzkyGolay", halfWindowSize=10)
  else if (smooth=="MovingAverage")
      spectra <- smoothIntensity(spectra, method="MovingAverage",halfWindowSize=2)
  else if (smooth=="WMovingAverage")
      spectra <- smoothIntensity(fiedler2009subset, method="MovingAverage", halfWindowSize=2, weighted=TRUE)
  
  if (baseline=="SNIP")
        spectraRB <- removeBaseline(spectra, method="SNIP", iterations=Iter)
  else if (baseline== "TopHat")
        spectraRB <- removeBaseline(spectra, method="TopHat")
  else if (baseline== "ConvexHull")
        spectraRB <- removeBaseline(spectra, method="ConvexHull")
  
  if (normalization=="TIC") spectraRB <- calibrateIntensity(spectraRB, method="TIC")
  else if (normalization=="PQN") spectraRB <- calibrateIntensity(spectraRB, method="PQN")
  
  # spectraRB <- alignSpectra(spectraRB)
  
  if (align) suppressWarnings(spectraRB <- alignSpectra(spectraRB, halfWindowSize=20, SNR=SN_R,  tolerance=20, warpingMethod="lowess"))
  
  print("start detecting peaks......")
  peaks <- detectPeaks(spectraRB, SNR=SN_R, halfWindowSize=20)
  peaks <- binPeaks(peaks)
  if (minFreq>0) peaks <- filterPeaks(peaks, minFrequency=minFreq)
  print("spectra processing completed...")
 
  for (i in 1:length(peaks)) {
    sx<-str_locate(peaks[[i]]@metaData$file,".mzXML")
    sd<-str_locate(peaks[[i]]@metaData$file, Dp)
    code<-str_sub(peaks[[i]]@metaData$file, sd[2]+2, sx[1]-1)
    # peaks[[i]]@metaData$file<-code
    if (i==1)
      codes<-code
    else 
      codes<-c(codes, code)
  }
  
  featureMatrix <- intensityMatrix(peaks, spectraRB)
  rownames(featureMatrix) <- codes
  featureMatrix <-cbind(featureMatrix, codes)
  colnames(featureMatrix)[ncol(featureMatrix)]<-keyCode
  df<-as.data.frame(featureMatrix)
  
  df2<-read.csv(filein, stringsAsFactors=FALSE)  #sep=";", 
  colnames(df2)[1]<-keyCode
  colnames(df2)[2]<-speColumn
  print(dim(df))
  print(dim(df2))
  df_p <- inner_join(df2, df, by = keyCode)
  gen<-unlist(strsplit(df_p[1, 2], split = " "))[1]
   for (row in 2:nrow(df_p)) {
     gen <- c(gen, unlist(strsplit(df_p[row, 2], split = " "))[1])
   }
   df_p$Genre<-gen
   df_p<-df_p[,c(1,2,ncol(df_p),4:ncol(df_p)-1)]
   nr <- paste(df_p$Taxonomie,df_p$Identifiant_MALDI, sep="_")
   rownames(df_p)<-nr
   write.csv(df_p , fileout)
  
  return(df_p)
  }
}
