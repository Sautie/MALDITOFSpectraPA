library("pvclust")


# default values dft<-BHclus(df)
# dft<-BHclus(df, vcol="all")
BHclus <- function(df_m, meth="ward.D2", d="euclidean", nb=100, vcol="genre", pgen="Lactobacillus", 
                   pspe="Pediococcus pentosaceus", fl="Erica cinerea", ruch="VI",  fig=TRUE){
 
   if (vcol=="genre") {
     df_m<-df_m[df_m$Genre==pgen, ]
     namefile<-paste(pgen, "txt", sep = ".")
   }
  else if (vcol=="species") {
    df_m<-df_m[df_m$Taxonomie==pspe, ] 
    namefile<-paste(pspe, "txt", sep = ".")
  }
  else if (vcol=="flower") {
    df_m<-df_m[df_m$Nutrition==fl, ]
    namefile<-paste(fl, "txt", sep = ".")
  }
  else if (vcol=="Ruche"){
    df_m<-df_m[df_m$Ruche==ruch, ] 
    namefile<-paste(ruch, "txt", sep = ".")
  }

set.seed(873) 
pv <- pvclust(t(df_m),
              method.hclust=meth,
              method.dist=d, nboot=nb)
if (fig) {
  plot(pv, hang = -1, cex = 0.5)
  pvrect(pv)
}
clus<- pvpick(pv)
clus<-as.matrix(clus)
write.table(clus, namefile)
return (clus)
}




dd <- dist(scale(USArrests), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
plot(as.phylo(hc), edge.color = "black", tip.color = "black" ,cex = 0.6, label.offset = 0.5)




