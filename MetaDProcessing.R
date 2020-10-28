library("dplyr")
library("readxl")
library("matrixStats")

#' @noRd AgFilter is an internal function for column selection ("select" , default) or aggregation based on "means", "medians", "max" and "min".
#'       of scores. the function of aggregation is determined by the fu value, The value "select"
#'       is for selecting one of the score column determined by y. (y=2, default).
AgFilter <- function(data_Scores, fu="select", y=2){
    p<-length(data_Scores)
    if (fu=="means") {
       num<-rowMeans(data.matrix(data_Scores[,y:p]))
    } else if (fu=="medians"){
       num<-rowMedians(data.matrix(data_Scores[,y:p]))
    } else if (fu=="select") {
       num<-data_Scores[,y]
    } else if (fu=="max") {
       num<-apply(data_Scores[,y:p],1,max)
    } else if (fu=="min") {
        num<-apply(data_Scores[,y:p],1,min)
    }
  return(num)
}

#' @title PMetaData
#' @description extracts from csv files the taxonomic identifications of the isolates and transforms them according to the values of the scores
#' @param  Dirp: directory (and path) where the csv files are located,
#' @param fileout: name of the output file where Maldi identifiers and species names  are found,
#' @param fu: this parameter determines whether only one of the scores (y, default value=2) is considered, or instead, the mean, median, minimum or
#'            maximum of the scores corresponding to each of the isolates, "means", "medians", "select", "max" and "min",
#' @param id: determine whether unidentified isolates are labeled with "NRI" (default value) or removed ,
#' @param sc: columns containing the scores  spNames:columns containing the names of the identified species
#' @return     dataframe (and csv file) with three columns: Maldi identifiers, isolate identifications and the corresponding scores
#' @examples Rn<-PMetaData("Bees project 2019-2020"), Rn<-PMetaData("Bees project 2019-2020", fileout="expMeans.csv", fu="means")
#'           Rn<-PMetaData("Bees project 2019-2020", fileout="expMedians.csv", fu="medians")
#            Rn<-PMetaData("Bees project 2019-2020", fileout="expMin.csv", fu="min"), Rn<-PMetaData("Bees project 2019-2020", fileout="expMax.csv", fu="max")
#'@export
PMetaData <- function(Dirp, fileout="expMax.csv", fu="max", y=2, id="NRI", sc=c(1,8,13,18,23,28,33,38,43, 48,53), spNames=c(1,5,10,15,20,25,30,35,40,45,50)){
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
    num <- AgFilter(data_Scores, fu, y)
    ind <- length(num)
    i = 1
    isolate <- rep(id, ind)
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
    if (id!="NRI") df<-df[df$isolate!=id, ]
    write.csv(df, fileout, row.names = FALSE)
  }
return(df)
}

#' @title  MergeMetaData
#' @description Merges two dataframes, one is the output of the funct. PMetaData1 and the other comes from a xlsx file. Merging is done through the Maldi code. Therefore it is important to ensure error-free transcription of Maldi codes.
#' @param  df: input dataframe (output of PMetaData1), Dp: folder containing the metadata files,
#' @param  filein: xlsx files for metadata
#' @param  fileout: csv file for both joined dataframes,
#' @param  keyCode: key code used for Dataframe merging, "Maldi_code" (Default)
#' @return  both dataframes joined in a new dataframe and an output csv file
#' @examples dfm<-MergeMetaData(Rn, "Bees project 2019-2020","Bees metadata.xlsx", "Filtered_Meta.csv")
#' @export

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


#' @title DropCMetadata
#' @description removes variables from the metadata csv file
#' @param   dm: input dataframe,
#' @param   fileout: csv file for both joined dataframes
#' @return  output dataframe and csv file
#' @examples dfc<-DropCMetadata(dfm,"DCFilered_Meta.csv", c("Date d'analyse","Date de récolte"))
#' @export

DropCMetadata<- function(dm, fileout, Co) {
  gm<-dm[ , !(names(dm) %in% Co)]
  write.csv(gm, fileout, row.names = FALSE)
  return(gm)
}




