
library("NbClust")
library("factoextra")
library("fpc")

# varCat1="Genre", value="Lactobacillus" 
# OptClusters(df_Peaks, varCat1="Taxonomie"	, value="Enterococcus faecalis")
# OptClusters(df_Peaks, varCat1="Taxonomie"	, value="All")
# OptClusters(df_Peaks, varCat1="Taxonomie"	, value="All", ind="gap statistics")
# OptClusters(df_Peaks, meth="hclust", varCat1="Taxonomie"	, value="All", ind="gap statistics")

OptClusters <- function(df_m,  meth="kmeans", dist="euclidean", varCat1, value, minc=2, maxc=10, ind="average silhouette width"){
  # meth: "ward.D","ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid",  "kmeans"
  # d: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"
  # varCat1: "Identifiant_MALDI"	"Taxonomie"	"Genre"	"Date.d.analyse"	"Origine"	"Ruche"	"Nutrition"	"Date.de.récolte"	"Lieu.de.la.ruche"
  ni <- 1
  nf <- 10
  f<-eval(parse(text=(paste0("df_m$", varCat1))))
  if (varCat1 %in% colnames(df_m)) {
    if (value == "All")
      chosenrows <- rep (TRUE,length(f))
    else {
      if (value %in% f) {
        chosenrows <- (f == value)
      }
      else
        stop("this value is not found ...try another")
    }
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
  {
    mes <-  paste("The variable",  varCat1, "is not found in this dataframe ...try another variable" )
    stop(mes)
  }
  
  dfs<-scale(df_S)
  # https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
 if (meth=="kmeans"){
  if (ind=="total within sum of square") {
   
    print(fviz_nbclust(dfs,kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "Elbow method"))
   
                          }
  else if (ind=="average silhouette width") {
     print(fviz_nbclust(dfs, kmeans, method = "silhouette")+
    labs(subtitle = "Silhouette method"))
  }
  else  if (ind=="gap statistics") {
  set.seed(1677)  
  print(fviz_nbclust(dfs, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
    labs(subtitle = "Gap statistic method"))
  }}
  else if (meth=="pam"){
    if (ind=="total within sum of square") {
      
      print(fviz_nbclust(dfs,pam, method = "wss") +
              geom_vline(xintercept = 4, linetype = 2)+
              labs(subtitle = "Elbow method"))
      
    }
    else if (ind=="average silhouette width") {
      print(fviz_nbclust(dfs, pam, method = "silhouette")+
              labs(subtitle = "Silhouette method"))
    }
    else  if (ind=="gap statistics") {
      set.seed(1677)  
      print(fviz_nbclust(dfs, pam, nstart = 25,  method = "gap_stat", nboot = 50)+
              labs(subtitle = "Gap statistic method"))
    }}
  else if (meth=="hclust"){   #"ward.D2", "euclidean"
    if (ind=="total within sum of square") {
      
      print(fviz_nbclust(dfs,hcut, method = "wss") +
              geom_vline(xintercept = 4, linetype = 2)+
              labs(subtitle = "Elbow method"))
      
    }
    else if (ind=="average silhouette width") {
      print(fviz_nbclust(dfs, hcut, method = "silhouette")+
              labs(subtitle = "Silhouette method"))
    }
    else  if (ind=="gap statistics") {
      set.seed(1677)  
      print(fviz_nbclust(dfs, hcut, nstart = 25,  method = "gap_stat", nboot = 50)+
              labs(subtitle = "Gap statistic method"))
    }}
  
  
  # nb <- NbClust(dfs, distance = dist, min.nc = minc,  max.nc = maxc, method = "kmeans", index=c(  "silhouette", "duda", "ratkowsky","mcclain",  "gplus",  "dunn")) #, 
  # print(fviz_nbclust(nb) + theme_gray())
}


VisualOptClusters<- function(df_m, meth="hclust", meth2="ward.D2", dist="euclidean", varCat1="Genre", value="Lactobacillus", nc=3, graphs=TRUE){
  # dist: "euclidean", "manhattan", "maximum", "canberra", "binary",  "minkowski", "pearson", "spearman", "kendall"
  # meth: "kmeans", "pam", "clara", "fanny", "hclust","agnes", "diana"
  # meth2(hclust): "ward.D","ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", 
  
  ni <- 1
  nf <- 10
  f<-eval(parse(text=(paste0("df_m$", varCat1))))
  if (varCat1 %in% colnames(df_m)) {
    if (value == "All")
      chosenrows <- rep (TRUE,length(f))
    else {
      if (value %in% f) {
        chosenrows <- (f == value)
      }
      else
        stop("this value is not found ...try another")
    }
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
  {
    mes <-  paste("The variable",  varCat1, "is not found in this dataframe ...try another variable" )
    stop(mes)
  }
  
  df_S<-scale(df_S)

if (graphs) { 
  if (meth=="kmeans")  {
    out <- eclust(df_S, kmeans, k = nc,hc_metric =dist, hc_method = meth2, graph = FALSE)
    print(fviz_cluster(out, geom = "point", frame.type = "norm") )   }
  else if (meth=="pam") {
    out <- eclust(df_S, pam, k = nc,hc_metric =dist, hc_method = meth2, graph = FALSE)
    print(fviz_cluster(out, geom = "point", frame.type = "norm") )}
  else if (meth=="hclust") {
    out <- eclust(df_S, hclust, k = nc,hc_metric =dist, hc_method = meth2, graph = FALSE)
    print(fviz_dend(out, rect = TRUE, show_labels = TRUE))  }
  
   print(fviz_silhouette(out) + scale_fill_brewer(palette = "Set1") +  scale_color_brewer(palette = "Set1") +  theme_gray()+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
   print(fviz_gap_stat(out$gap_stat) )
   
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
