library("dplyr")
library("readxl")

#'AgFilter
#' @noRd internal function 
AgFilter <- function(data_Scores, fu="select", y=2){
    p<-length(data_Scores)
    if (fu=="means") {
       num<-rowMeans(data_Scores[,y:p]) 
    } else if (fu=="medians"){
       num<-rowMedians(data_Scores[,y:p])     
    } else if (fu=="select") {
       num<-data_Scores[,y]
    } else if (fu=="max") {
       num<-apply(data_Scores[,y:p],1,max)
    } else if (fu=="min") {
        num<-apply(data_Scores[,y:p],1,min)
    }
  return(num)
}
#' PMetaData
#' @description PMetaData extracts from csv files the names of the identified isolates and transforms them according to the values of the scores#'            
#' @param  
#'        Dirp: directory (and path) where the csv files are located, fileout: name of the output file where Maldi identifiers and species names 
#'        are found, fu: this parameter determines whether only one of the scores (y, default value=2) is considered, or instead, the mean, median, minimum or 
#'        maximum of the scores corresponding to each of the isolates, "means", "medians", "select", "max" and "min", id: determine whether unidentified 
#'        isolates are labeled with "NRI" (default value) or removed , sc: columns containing the scores  spNames:columns containing the names of the identified species
#' @return dataframe (and csv file) with three columns: Maldi identifiers, isolate identifications and the corresponding scores
#' @examples Rn<-PMetaData("Bees project 2019-2020"), R<-PMetaData("Bees project 2019-2020", fileout="exp0.csv", id="0") 
#'           Rm<-PMetaData("C:/Users/14389/Documents/Bees project 2019-2020", fu="medians")
#Dirp="C:/Users/14389/Documents/Bees project 2019-2020"   fileout="exp.csv"
PMetaData <- function(Dirp, fileout="exp.csv", fu="select", y=2, id="NRI", sc=c(1,8,13,18,23,28,33,38,43, 48,53), spNames=c(1,5,10,15,20,25,30,35,40,45,50)){
  if (!dir.exists(Dirp))  {
    warning("this directory does not exist in the current working directory")
  }
  else {
    all_data_frames <-
      lapply(
        list.files(
          Dirp,
          pattern = ".csv",
          recursive = TRUE,
          full.names = TRUE
        ),
        read.csv,
        header = FALSE,
        stringsAsFactors = FALSE,
        sep = ';'
      )
    single_data_frame <- Reduce(rbind, all_data_frames)
    sc <- c(1, 8, 13, 18, 23, 28, 33, 38, 43, 48, 53)
    spNames <- c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
    data_Scores <- single_data_frame[, sc]
    data_Names <- single_data_frame[, spNames]
    v0 <- rep("v", length(data_Scores))
    names(data_Scores) <- paste0(v0, sc)
    vn <- rep("v", length(data_Names))
    names(data_Names) <- paste0(vn, spNames)
    data_Scores <- subset(data_Scores, v1 != "BTS")
    data_Names <- subset(data_Names, v1 != "BTS")
    data_Scores <- subset(data_Scores, v1 != "CTL")
    data_Names <- subset(data_Names, v1 != "CTL")
    num <- AgFilter(data_Scores)
    ind <- length(num)
    i = 1
    isolate <- rep("NRI", ind)
    for (row in 1:nrow(data_Names)) {
      Pri <- grepl("not reliable identification", data_Names[row, ])
      Priv <- which(!Pri)
      if (num[i] >= 2) {
        if (length(Priv) > 1) {
          isolate [i] <- data_Names[i, Priv[2]]
        }
      } else if ((num[i] < 2) && (num[i] >= 1.7)) {
        if (length(Priv) > 1)  {
          rt <- strsplit(data_Names[i, Priv[2]], " ")
          isolate [i] <- paste(rt[[1]][1], "sp.")
        }
      }
      i <- i + 1
    }
    Maldi_code <- data_Scores[, 1]
    df <- data.frame(Maldi_code, isolate, num,  stringsAsFactors = FALSE)
    if (id!="NRI") df<-df[df$isolate!="NRI", ] 
    write.csv(df, fileout, row.names = FALSE)
  }
return(df)
}


#' @description Merges two dataframes, one is the output of the funct. PMetaData1 and the other comes from a xlsx file. Merging between
#'              both dataframes is done through the Maldi code. Therefore it is important to ensure error-free transcription of 
#'              the Maldi codes. 
#' @param   
#'          df: input dataframe (output of PMetaData1), Dp: folder containing the metadata files, filein: xlsx files for metadata 
#'          fileout: csv file for both joined dataframes, keyCode: "Maldi_code" (Default)
#' @return  both dataframes joined in a new dataframe and an output csv file
#' @examples dfm<-MergeMetaData(Rn, "Bees project 2019-2020","Bees metadata.xlsx", "Filered_Meta.csv")
         
MergeMetaData <- function(df, Dp, filein, fileout, keyCode="Maldi_code"){
  if (!dir.exists(Dp))  { 
    warning("this directory does not exist in the current working directory") } 
  else {
  filed<-file.path(Dp, filein)
  df2 <- read_excel(filed)
  colnames(df2)[2]<-keyCode
  df_m <- inner_join(df, df2, by = keyCode)
  write.csv(df_m , fileout, row.names = FALSE)
  }
  return(df_m)
}

#' DropCMetadata
#' @Description remove variables from the metadata csv file
#' @param   dm: input dataframe, Dp: folder containing the metadata files, fileout: csv file for both joined dataframes
#' @return  output dataframe and csv file
#' @examples dfc<-DropCMetadata(dfm, "Bees project 2019-2020","DCFilered_Meta.csv", c("Date d'analyse","Date de récolte")) 

DropCMetadata<- function(dm, Dp, fileout, Co) {
  if (!dir.exists(Dp))  { 
    warning("this directory does not exist in the current working directory") }
  else {
  gm<-dm[ , !(names(dm) %in% Co)]
  write.csv(gm, fileout, row.names = FALSE)
  }
  return(gm)
}
  
  


