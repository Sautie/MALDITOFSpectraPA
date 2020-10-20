library(stringr) 
library("readxl")

source("prueba.R")
Dirp<-"Bees project 2019-2020"
newdir<-"newdata"

RenameXml<- function(Dirp) {

  if (!dir.exists(Dirp))  { warning("this directory does not exist in the current working directory") }
  else {
 vmz<-list.files(Dirp,  pattern = ".mzXML", recursive = TRUE, full.names = TRUE)
 for (v in vmz) {
  vvf<-unlist(strsplit(v, "/"))
  Dirp<-paste("\\b", Dirp, "\\b" ,sep = "")
  pw<-which(grepl(Dirp,vvf) )+2
  vvf[length(vvf)]<-paste(vvf[pw],".mzXML",sep = "")
  newname<-str_c(vvf, collapse = "/")
  file.rename(v, newname)
}
  }
}


CopyRenameXml<- function(Dirp, newdir) {
   if (!dir.exists(Dirp))  { warning("this directory does not exist in the current working directory") }
  else {
  while (dir.exists(newdir)){ 
    rn<-floor(runif(1, min=0, max=10))
    newdir<-paste(newdir, rn, sep = "")
  }
dir.create(newdir)
vmz<-list.files(Dirp, pattern = ".mzXML", recursive = TRUE, full.names = TRUE)
for (v in vmz) {
  vw<-strsplit(v, "/")
  vvf<-unlist(vw)
  Dirp<-paste("\\<", Dirp, "\\>" ,sep = "")
     # bvf<-str_detect(vvf, Dirp)
  pw0<-which(grepl(Dirp,vvf) )
  pw<-pw0+2
  if (!(vvf[pw] %in% c("BTS","BTS_Validation","CTL", "Autocalibration")) ) {
  newfile<-paste(vvf[pw],".mzXML",sep = "")
  cvvf<-c(vvf[1:(pw0-1)],newdir,newfile)
  newname<-str_c(cvvf, collapse = "/")  
  file.copy(v, newname) }
}
  }
}

XmlCleaning<- function(Dirp, newdir, op) {
  if (!dir.exists(Dirp))  { warning("this directory does not exist in the current working directory") }
 else {
 while (dir.exists(newdir)){ 
  rn<-floor(runif(1, min=0, max=10))
  newdir<-paste(newdir, rn, sep = "")
  }
 dir.create(newdir)
 patnew<-file.path(wd, newdir)
 vmz<-list.files(Dirp,recursive = TRUE, full.names = TRUE)
 r<-0
 for (v in vmz) {
   con<- file(v, open = "r")
   t<-0
   while ((length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 1)&&(t<4)) {t<-t+1} 
   if (op=="delete") { if (t==1) {
     file.remove(v)} }
   else { 
     if (t>1) { 
       vw<-strsplit(v, "/")
       newfile<-paste(patnew,vm(length(vm)),sep = "/")
       copy.file(v,newfile)
     }
   }
 }
  }
}

XmlChoosing <- function(Dirp, newdir, op, filein="Bees metadata.xlsx") {
   if (!dir.exists(Dirp))  { warning("this directory does not exist in the current working directory") }
  else {
    while (dir.exists(newdir)){ 
      rn<-floor(runif(1, min=0, max=10))
      newdir<-paste(newdir, rn, sep = "")
    }
    dir.create(newdir)
    patnew<-file.path(wd, newdir)
    filed<-file.path(Dirp, filein)
    df2 <- read_excel(filed)
    dfc<- df2[, "Identifiant MALDI"]
    vmz<-list.files(Dirp, pattern = ".mzXML", recursive = TRUE, full.names = TRUE)
    for (v in vmz) {
      vw<-unlist(strsplit(v, "/"))
      name<-unlist(strsplit(vm[length(vm)], "."))
      namew<-paste("\\<", name[1], "\\>" ,sep = "")
      pw0<-length(which(grepl(namew,dfc) ))
      if (op=="delete") { if (pw0==0) {file.remove(v)} }
      else {
        if (pw0>0)
          newfile<-paste(patnew,vm(length(vm)),sep = "/")
          copy.file(v,newfile)
              }
    }
  
  }
}
  
  
  


