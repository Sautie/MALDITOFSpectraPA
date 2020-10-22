
library("NbClust")
library("factoextra")
library("fpc")



OptClusters <- function(df_m, meth="ward.D2", dist="euclidean", varCat1="Genre", value="Lactobacillus", minc=2, maxc=20){
  # meth: "ward.D","ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid",  "kmeans"
  # d: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"
  # varCat1: "Identifiant_MALDI"	"Taxonomie"	"Genre"	"Date.d.analyse"	"Origine"	"Ruche"	"Nutrition"	"Date.de.récolte"	"Lieu.de.la.ruche"
  
  ni <- 1
  nf <- 9
  if (varCat1 %in% colnames(df_m)) {
    if (value == "All")
      chosenrows <- df_m$varCat1
    else {
      if (value1 %in% df_m$varCat1) {
        chosenrows <- (df_m$varCat1 == value)
        if (length(which(chosenrows)) > 3)  {
          df_S <- df_m[chosenrows, ]
          df_S <- df_S[, -ni:-nf]
          for (j in 1:ncol(df_S)) {
            df_S[, j] <- as.numeric(df_S[, j])
          }
        }
        else
          stop("too small number of individuals...try another value")
      }
      else
        stop("this value is not found ...try another")
    }
  }
  else
  {
    mes <-  paste("The variable",  as.chracter(varCat1), "is not found in", as.chracter(df_m),"...try another variable" )
    stop(mes)
  }
  for (j in 1:ncol(df_S))
  {
    df_S[,j]<-as.numeric(df_S[,j])
  }
  dfs<-scale(df_S)
  nb <- NbClust(dfs, distance = dist, min.nc = minc,  max.nc = maxc, method = meth)
  fviz_nbclust(nb) + theme_gray()
  listResults <- list()
  for(i in 1:3) {
    listResults <- nb$All.index
    listResults <- nb$Best.nc
    listResults <-nb$Best.partition
   }
return(listResults)
}

VisualOptClusters<- function(df_m, meth="hclust", meth2="ward.D2", dist="euclidean", varCat1="Genre", value="Lactobacillus", nc=3, graphs=TRUE){
  # dist: "euclidean", "manhattan", "maximum", "canberra", "binary",  "minkowski", "pearson", "spearman", "kendall"
  # meth: "kmeans", "pam", "clara", "fanny", "hclust","agnes", "diana"
  # meth2(hclust): "ward.D","ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", 
  
  ni <- 1
  nf <- 9
  if (varCat1 %in% colnames(df_m)) {
    if (value == "All")
      chosenrows <- df_m$varCat1
    else {
      if (value1 %in% df_m$varCat1) {
        chosenrows <- (df_m$varCat1 == value)
        if (length(which(chosenrows)) > 3)  {
          df_S <- df_m[chosenrows, ]
          df_S <- df_S[, -ni:-nf]
          for (j in 1:ncol(df_S)) {
            df_S[, j] <- as.numeric(df_S[, j])
          }
        }
        else
          stop("too small number of individuals...try another value")
      }
      else
        stop("this value is not found ...try another")
    }
  }
  else
  {
    mes <-  paste("The variable",  as.chracter(varCat1), "is not found in", as.chracter(df_m),"...try another variable" )
    stop(mes)
  }
  for (j in 1:ncol(df_S))
  {
    df_S[,j]<-as.numeric(df_S[,j])
  }
  df_S<-scale(df_S)
 out <- eclust(df_S, FUNcluster=meth, k = nc,hc_metric =dist, hc_method = meth2, graph = FALSE)
if (graphs) { 
  if (meth=="kmeans")  {
     fviz_cluster(out, geom = "point", frame.type = "norm")    }
  else if (meth=="pam") {
     fviz_cluster(out, geom = "point", frame.type = "norm")   }
  else if (meth=="hclust") {
    fviz_dend(out, rect = TRUE, show_labels = TRUE)  }
    fviz_silhouette(out) + scale_fill_brewer(palette = "Set1") +  scale_color_brewer(palette = "Set1") +  theme_gray()+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    fviz_gap_stat(out$gap_stat)   
    }
 print(out$silinfo)
return(out$silinfo)
}



PointClusterVal <- function(df_m, meth="ward.D2", d="euclidean", varCat1="Genre", value="Lactobacillus", VarCat2="Nutrition"){
# https://www.rdocumentation.org/packages/fpc/versions/2.2-8/topics/cluster.stats
# https://www.datanovia.com/en/lessons/cluster-validation-statistics-must-know-methods/
#  External clustering validation

ni <- 1
nf <- 9
if (varCat1 == varCat2)
  stop("both variables cannot be equal")
else
{
  df_A <- df_m[, ni:nf]
  for (j in 1:ncol(df_A)) {
    df_A[, j] <- as.numeric(factor (df_A[, j]))
  }
  if (varCat2 %in% colnames(df_A)) {
    cl <- df_A$varCat2
    Lcl <- length(unique(cl))
    if (Lcl > 50) {
      mes <-
        paste("The variable",  as.chracter(varCat2),"has more than 50 levels!")
      warning(mes)
    }
  }
  else
    if (varCat2 %in% colnames(df_m)) {
      mes <-paste("The variable", as.chracter(varCat2),  "is a peak intensity...try another variable")
      stop(mes)
    }
  else {
    mes <- paste("The variable", as.chracter(varCat2), "is not found in", as.chracter(df_m), "...try another variable" )
    stop(mes)
  }
  if (varCat1 %in% colnames(df_m)) {
    if (value == "All")
      chosenrows <- df_m$varCat1
    else {
      if (value1 %in% df_m$varCat1) {
        chosenrows <- (df_m$varCat1 == value)
        if (length(which(chosenrows)) > 3)  {
          df_S <- df_m[chosenrows, ]
          df_S <- df_S[, -ni:-nf]
          for (j in 1:ncol(df_S)) {
            df_S[, j] <- as.numeric(df_S[, j])
          }
        }
        else
          stop("too small number of individuals...try another value")
      }
      else
        stop("this value is not found ...try another")
    }
  }
  else
  {
    mes <-  paste("The variable", as.chracter(varCat1), "is not found in",as.chracter(df_m), "...try another variable" )
    stop(mes)
  }
}
for (j in 1:ncol(df_S))
{
  df_S[,j]<-as.numeric(df_S[,j])
}
df_S<-scale(df_S)
s <- get_clust_tendency(df_S, 40, graph = TRUE)

MD <- dist(df_S, method =dist)
out <- eclust(df_S, FUNcluster=meth, k = nc,hc_metric =dist, hc_method = meth2, graph = FALSE)
out_stats <- cluster.stats(MD, out$cluster)
table(cl, out$cluster)

out_ext_stats <- cluster.stats(MD, cl, out$cluster)

return(c(s$hopkins_stat, out_stats$dunn,out_stats$dunn2,out_stats$pearsongamma,out_stats$avg.silwidth, out_ext_stats$corrected.rand,out_ext_stats$vi))
}
