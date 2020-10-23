library("stringr") 
library("readxl")

source("MetaDProcessing.R")

Dirp<-"Bees project 2019-2020"
newdir<-"newdata"


#' @Description  RenameXml copy and renames the mzXML files that are not in the folders: "BTS", "BTS_Validation", "CTL", "Autocalibration"
#' @param      
#'            Dirp: folder containing the mzXML files, newdir: the folder to which the mzXML files are copied
#' @examples  CopyRenameXml("Bees project 2019-2020", "newdata")

CopyRenameXml <- function(Dirp, newdir) {
  if (!dir.exists(Dirp))  {
    warning("this directory does not exist in the current working directory")
  }
  else {
    while (dir.exists(newdir)) {
      rn <- floor(runif(1, min = 0, max = 10))
      newdir <- paste(newdir, rn, sep = "")
    }
    dir.create(newdir)
    vmz <-
      list.files(
        Dirp,
        pattern = ".mzXML",
        recursive = TRUE,
        full.names = TRUE
      )
    print("running....")
    for (v in vmz) {
      vvf <- unlist(strsplit(v, "/"))
      Dirp <- paste("\\<", Dirp, "\\>" , sep = "")
      pw0 <- which(grepl(Dirp, vvf))
      pw <- pw0 + 2
      if (!(vvf[pw] %in% c("BTS", "BTS_Validation", "CTL", "Autocalibration"))) {
        newfile <- paste(vvf[pw], ".mzXML", sep = "")
        cvvf <- c(newdir, newfile)
        newname <- str_c(cvvf, collapse = "/")
        file.copy(v, newname)
      }
    }
  }
}

#' @Description  XmlNCleaning  copies all mzXML files with at least n lines to a new folder or delete those with less than n lines,
#'               (the mzXML files have 26 lines) this way incomplete or truncated files are removed...the number of characters of line n=18 
#'               is also verified
#' @param  
#'             Dirp:folder containing the Maldi_Tof spectrum files, newdir:folder where the chosen files are copied, 
#'             op: "delete"(default) or "copy", n number of lines (n=18 default), number of characters of line 8 (440518 default)  
#'            
#' @examples  XmlNCleaning("newdata5")  XmlNCleaning("newdata5", "newdata", op="copy"),XXmlNCleaning("newdata5", "newdata", op="copy", 26)

XmlNCleaning<- function(Dirp, newdir, op="delete", n=18, pl=440518) {
  if (!dir.exists(Dirp))  { warning("this directory does not exist in the current working directory") }
  else {
    if (op!="delete"){
    while (dir.exists(newdir)){ 
      rn<-floor(runif(1, min=0, max=10))
      newdir<-paste(newdir, rn, sep = "")
    }
    dir.create(newdir)
    }
    vmz<-list.files(Dirp,recursive = TRUE, full.names = TRUE)
    print("running....")
    for (v in vmz) {
      d<-readLines(v, n)
      t<-length(d)
      if (n>=18) j<-nchar(d[18])
      if (op=="delete") { 
        if (n<18){ if (t<n) file.remove(v)}
        else {
          if ((t<n)||(j<pl))  file.remove(v)
        }
        }
      else { 
        if (n<18){
        if (t==n){
          vm<- unlist(strsplit(v, "/"))
          newfile<-paste(newdir,vm[length(vm)],sep = "/")
          file.copy(v,newfile)
          }
        }
          else
          { 
            if ((t==n)&&(j>=pl)) {
              vm<- unlist(strsplit(v, "/"))
              newfile<-paste(newdir,vm[length(vm)],sep = "/")
              file.copy(v,newfile)
              }
          }
          
        }
      }
    }
  }


#' @Description  SlowXml26Cleaning  copies all mzXML files with 26 lines to a new folder or delete those with less than 26 lines,
#'               slow and old version of XmlNCleaning 
#' @param  
#'             Dirp:folder containing the Maldi_Tof spectrum files, newdir:folder where the chosen files are copied, 
#'             op: "delete"(default) or "copy", metadata xlsx file    
#'            
#' @examples  SlowXml26Cleaning("newdata5", "newdata")  SlowXml26Cleaning("newdata5", "newdata", op="copy")

SlowXml26Cleaning<- function(Dirp, newdir, op="delete") {
  if (!dir.exists(Dirp))  { warning("this directory does not exist in the current working directory") }
  else {
    if (op!="delete"){
    while (dir.exists(newdir)){ 
      rn<-floor(runif(1, min=0, max=10))
      newdir<-paste(newdir, rn, sep = "")
    }
    dir.create(newdir)
    }
    vmz<-list.files(Dirp,recursive = TRUE, full.names = TRUE)
    print("running....")
    for (v in vmz) {
      opened<- file(v, open = "r")
      readop <- paste("wc -l ", v, " | awk '{ print $1 }'", sep="")
      suppressWarnings(n <- system(command=readop, intern=TRUE))
      t<-as.numeric(gsub(".*?([0-9]+).*", "\\1",  n[5])) 
      close(opened)
      if (op=="delete") { if (t<26) {
        file.remove(v)} }
      else { 
        if (t==26) {
          vm<- unlist(strsplit(v, "/"))
          newfile<-paste(newdir,vm[length(vm)],sep = "/")
          file.copy(v,newfile)
        }
      }
    }
  }
}



#' @Description XmlChoosing copies to a new folder those mzXML files whose Maldi identifiers are found in the metadata  xlsx file 
#'               or deletes the mzXML files that have no spectrum identifier in the metadata file
#' @param 
#'             Dirp:folder containing the Maldi_Tof spectrum files, newdir:folder to copy the chosen files, 
#'             metadata xlsx file
#'

XmlChoosing <- function(Dirp, newdir, filein) {
  if (!dir.exists(Dirp))  {
    warning("this directory does not exist in the current working directory")
  }
  else {
    
      while (dir.exists(newdir)) {
        rn <- floor(runif(1, min = 0, max = 10))
        newdir <- paste0(newdir, rn)
      }
      dir.create(newdir)
    filed <- file.path(Dirp, filein)
    df2 <- read_excel(filed)
    dfc <- df2[, "Identifiant MALDI"]
    vmz <-list.files(Dirp, pattern = ".mzXML", recursive = TRUE, full.names = TRUE)
    for (v in vmz) {
      vm <- unlist(strsplit(v, "/"))
      Dirp <- paste("\\<", Dirp, "\\>" , sep = "")
      pw0 <- which(grepl(Dirp, vm))
      pw <- pw0 + 2
      namew <- paste("\\<",vm[pw] , "\\>" , sep = "")
      pww <- length(which(grepl(namew, dfc)))
        if (pww > 0){
          newf <- paste(vm[pw], ".mzXML", sep = "")
          newfile <- paste(newdir, newf, sep = "/")
          file.copy(v, newfile)
        }
      
    }
    
  }
}
  


