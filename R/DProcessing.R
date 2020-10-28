library("stringr")
# library("readxl")


#' @title CopyRenameXml
#' @description CopyRenameXml is a function aimed at copying and renaming the mzXML files not being in the folders: "BTS", "BTS_Validation", "CTL", "Autocalibration"
#' @param   Dirp: folder containing the mzXML files,
#' @param   newdir: the folder to which the mzXML files are copied
#' @examples  CopyRenameXml("Bees project 2019-2020", "newdata")
#'@export

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

#' @title XmlNCleaning
#' @description copies all mzXML files with at least n lines to a new folder or delete those with less than n lines, (the mzXML files have 26 lines) this way incomplete or truncated files are removed...the number of characters of line n=18 is also verified
#' @param      Dirp: folder containing the Maldi_Tof spectrum files,
#' @param      newdir: folder where the chosen files are copied,
#' @param      op: "delete"(default value) or "copy",
#' @param      n: number of lines (n=18, default value),
#' @param      number: n of characters of line 8 (440518, default value)
#' @examples  XmlNCleaning("newdata")  XmlNCleaning("newdata5", "newdata", op="copy"),XXmlNCleaning("newdata5", "newdata", op="copy", n=26)
#'            its recommended to run this function with op="delete" after RenameXml,thus cleaning the output directory of the latter function
#'@export
#'
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


#' @title SlowXml26Cleaning
#' @description copies all mzXML files with 26 lines to a new folder or delete those with less than 26 lines, slow and old version of XmlNCleaning
#' @param  Dirp: folder containing the Maldi_Tof spectrum files,
#' @param  newdir: folder where the chosen files are copied,
#' @param  op: "delete"(default value) or "copy", metadata xlsx file
#' @examples  SlowXml26Cleaning("newdata5", "newdata"),  SlowXml26Cleaning("newdata5", "newdata", op="copy")
#'@export
#'
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



