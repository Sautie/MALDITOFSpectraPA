library("MALDIquant")
library("MALDIquantForeign")
library("stringr")
library("dplyr")


Generate_Peaks_Meta <- function(Dp="newdata", filein="Filered_Meta1.csv", fileout="Filered_Meta_Peaks.csv",
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
        spectraRB <- removeBaseline(spectra, method="SNIP", iterations=100)
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
  colnames(featureMatrix)[ncol(featureMatrix)]<-"Identifiant_MALDI"
  df<-as.data.frame(featureMatrix)
  
  df2<-read.csv(filein, sep=";", stringsAsFactors=FALSE)
  colnames(df2)[1]<-"Identifiant_MALDI"
  colnames(df2)[2]<-"Taxonomie"

  df_m <- inner_join(df2, df, by = 'Identifiant_MALDI')
  gen<-unlist(strsplit(df_m[1, 2], split = " "))[1]
   for (row in 2:nrow(df_m)) {
     gen <- c(gen, unlist(strsplit(df_m[row, 2], split = " "))[1])
   }
   df_m$Genre<-gen
   df_m<-df_m[,c(1,2,ncol(df_m),4:ncol(df_m)-1)]
   nr <- paste(df_m$Taxonomie,df_m$Identifiant_MALDI, sep="_")
   rownames(df_m)<-nr
   write.csv(df_m , fileout)
  
  return(df_m)
  }
}
