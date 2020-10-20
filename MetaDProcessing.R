# "C:/path/to/somewhere/else"
# pat<-"C:/Users/14389/Documents/Bees project 2019-2020"
# list.files(pat, pattern = ".csv", recursive = TRUE)
# cr<- read.csv(file = '191009-1309-1011018996.csv', stringsAsFactors = FALSE, sep = ';')
# cr<- read.csv(file = 'Bees project 2019-2020/190912-1312-1011019334.csv', header = FALSE, stringsAsFactors = FALSE, sep = ';')

# state <- function(x){
#  return(ifelse (x>=2, 2, ifelse ( ((x<2) && (x>=1.7)), 1, 0) ))
# }
library("dplyr")
library("readxl")

AgFilter <- function(fu="select", y=2){
    if (fu=="means") {
       num<-rowMeans(data_Scores[,2:11]) 
    } else if (fu=="medians"){
       num<-rowMedians(data_Scores[,2:11])     
    } else if (fu=="select") {
       num<-data_Scores[,y]
    } else if (fu=="max") {
       num<-apply(data_Scores[,2:11],1,max)
    } else if (fu=="min") {
        num<-apply(data_Scores[,2:11],1,min)
    }
  return(num)
}

Recoding<-function() {
  ind<-length(num)
  base<-rep("NRI", ind)
  for (row in 1:nrow(data_Names)) {
  Pri <- grepl("not reliable identification", data_Names[row,]) 
  Priv<-which(!Pri)
  if (num[row]>=2) {
   if (length(Priv)>1) { base[row]<- data_Names[row, Priv[2]]  }
   } else if  ((num[row]<2) && (num[row]>=1.7)) {
   if (length(Priv)>1)  { 
      rt<-strsplit(data_Names[row, Priv[2]], " ") 
      base[row]<- paste(rt[[1]][1],"spp")  }  
     }
    }
return(base)
}

PMetaData1 <- function(Dirp,fileout="exp.csv",  fu="select", y=2){
  wd<-getwd()
  pat<-file.path(wd, Dirp)
  if (!dir.exists(pat))  { 
    warning("this directory does not exist in the current working directory") }
  else {
    all_data_frames <- lapply(list.files(pat, pattern = ".csv", recursive = TRUE, full.names = TRUE), read.csv, header = FALSE, stringsAsFactors = FALSE, sep = ';' )
    single_data_frame <- Reduce(rbind, all_data_frames)
    data_Scores<-single_data_frame[,c(1,8,13,18,23,28,33,38,43, 48,53)]
    data_Names<-single_data_frame[,c(1,5,10,15,20,25,30,35,40,45,50)]
    names(data_Scores)<- c("v1","v8","v13", "v18","v23","v28","v33","v38","v43","v48","v53")
    names(data_Names)<- c("v1", "v5", "v10", "v15", "v20","v25","v30","v35","v40","v45","v50")
    data_Scores<-subset(data_Scores, v1!="BTS")
    data_Names<-subset(data_Names, v1!="BTS")
    data_Scores<-subset(data_Scores, v1!="CTL")
   data_Names<-subset(data_Names, v1!="CTL")
# y<-2
   num<-AgFilter(fu, y)
# if (fu=="means") {
#  num<-rowMeans(data_Scores[,2:11]) 
#   } else if (fu=="medians"){
#  num<-rowMedians(data_Scores[,2:11])     
#   } else if (fu=="select") {
#  num<-data_Scores[,y]
#    } else if (fu=="max") {
#  num<-apply(data_Scores[,2:11],1,max)
#   } else if (fu=="min") {
#  num<-apply(data_Scores[,2:11],1,min)
#   }
   base<-Recoding()
# ind<-length(num)
# i=1
# base<-rep("NRI", ind)
# for (row in 1:nrow(data_Names)) {
#   Pri <- grepl("not reliable identification", data_Names[row,]) 
#   Priv<-which(!Pri)
#   if (num[i]>=2) {
#     if (length(Priv)>1) { base [i]<- data_Names[i, Priv[2]]  }
#   } else if  ((num[i]<2) && (num[i]>=1.7)) {
#     if (length(Priv)>1)  { 
#       rt<-strsplit(data_Names[i, Priv[2]], " ") 
#       base [i]<- paste(rt[[1]][1],"spp")  }
#   }
#   i<-i+1
# }
# fileout<-"exp.csv"
code<-data_Scores[, 1] 
filed<-file.path(Dirp, fileout)
df<-data.frame(code, base, stringsAsFactors=FALSE)
write.csv(df, filed, row.names = FALSE)
  }
return(df)
}


MergeMetaData <- function(df, Dp="Bees project 2019-2020", filein="Bees metadata.xlsx", fileout="Filered_Meta.csv"){
  colnames(df)[1]<-"Identifiant MALDI"
  if (!dir.exists(Dp))  { 
    warning("this directory does not exist in the current working directory") }
  else {
  filed<-file.path(Dp, filein)
  df2 <- read_excel(filed)
  colnames(df2)[2]<-"Identifiant MALDI"
  df_m <- inner_join(df, df2, by = 'Identifiant MALDI')
  filed<-file.path(Dp, fileout)
  write.csv(df_m , filed, row.names = FALSE)
  }
  return(df_m)
}

DropCMetadata<- function(dm, Dp="Bees project 2019-2020", fileout="DCFilered_Meta.csv", Co=c("Date d'analyse","Date de récolte")) {
  if (!dir.exists(Dp))  { 
    warning("this directory does not exist in the current working directory") }
  else {
  filed<-file.path(Dp, fileout)
  gm<-dm[ , !(names(dm) %in% Co)]
  write.csv(gm, filed, row.names = FALSE)
  }
  return(gm)
}
  
  


