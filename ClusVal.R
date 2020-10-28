
library("NbClust")
library("factoextra")
library("fpc")
library("cluster")


#' @title OptClusters
#' @description OptClusters is a wrapper of several functions to visualize and compute optimal clusters for different clustering and evaluation methods
#' @param       df_m: dataframe containing peaks and metadata
#' @param       meth: clustering algorithms (meth="kmeans", default value), other values: "pam" or "hclust"
#' @param       dist:  distances ("euclidean", default),  "maximum", "manhattan", "canberra", "binary" "minkowski"
#' @param       varCat1: categorical variable for choosing isolates, examples: "Taxonomie"	,"Genre",	"Date.d.analyse"	,"Origine","Ruche", "Nutrition"	, "Date.de.récolte"	, "Lieu.de.la.ruche"
#' @param       value: level of catVar1, examples: "Lactobacillus" ("Genre"), Taxonomie("Pediococcus pentosaceus"), "Erica cinerea" ("Nutrition"),...
#' @param       minc: minimal number of clusters (minc=2, default value)
#' @param       maxc: maximal number of clusters (maxc=10, default value)
#' @param       ind: methods to evaluate clustering algorithms:"total within sum of squares" ,"average silhouette width"  and "gap statistics"
#' @param        nb: number of bootstrap samples (nb=100, default value)
#' @return       figures and statistics
#' @examples      OptClusters(df_Peaks, varCat1="Taxonomie"	, value="Enterococcus faecalis"),
#'                OptClusters(df_Peaks, varCat1="Taxonomie"	, value="All")
#'                OptClusters(df_Peaks, meth="pam", varCat1="Taxonomie"	, value="All", ind="gap statistics"),
#'                OptClusters(df_Peaks, meth="hclust", varCat1="Taxonomie"	, value="All", ind="gap statistics")
#'        source: https://rpubs.com/pg2000in/OptimumClusters
#'                http://rstudio-pubs-static.s3.amazonaws.com/265632_3ad9e0b981244e15887677f8dffb39a0.html#using-30-different-indices
#'                https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
#'                https://www.rdocumentation.org/packages/factoextra/versions/1.0.7/topics/fviz_nbclust
#'@export

OptClusters <- function(df_m,  meth="kmeans", dist="euclidean", varCat1, value, minc=2, maxc=10, ind="average silhouette width", nb=100){

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
      stop("too small number of isolates...try another value")
  }
  else
  {
    mes <-  paste("The variable",  varCat1, "is not found in this dataframe ...try another variable" )
    stop(mes)
  }
  dfs<-scale(df_S)
  if (meth=="kmeans") {
    if (ind == "total within sum of squares") {
      print(
        fviz_nbclust(dfs, kmeans, method = "wss") +
          geom_vline(xintercept = 4, linetype = 2) +
          labs(subtitle = "Elbow method")
      )
    }
    else if (ind == "average silhouette width") {
      print(
        fviz_nbclust(dfs, kmeans, method = "silhouette") +
          labs(subtitle = "Silhouette method")
      )
    }
    else  if (ind == "gap statistics") {
      set.seed(1677)
      print(
        fviz_nbclust(
          dfs,
          kmeans,
          nstart = 25,
          method = "gap_stat",
          nboot = nb
        ) +
          labs(subtitle = "Gap statistic method")
      )
    }
  }
  else if (meth == "pam") {
    if (ind == "total within sum of squares") {
      print(
        fviz_nbclust(dfs, pam, method = "wss") +
          geom_vline(xintercept = 4, linetype = 2) +
          labs(subtitle = "Elbow method")
      )
    }
    else if (ind == "average silhouette width") {
      print(
        fviz_nbclust(dfs, pam, method = "silhouette") +
          labs(subtitle = "Silhouette method")
      )
    }
    else  if (ind == "gap statistics") {
      set.seed(1677)
      print(
        fviz_nbclust(
          dfs,
          pam,
          # nstart = 25,
          method = "gap_stat",
          nboot = nb
        ) +
          labs(subtitle = "Gap statistic method")
      )
    }
  }
  else if (meth == "hclust") {
    #"ward.D2", "euclidean"
    if (ind == "total within sum of square") {
      print(
        fviz_nbclust(dfs, hcut, method = "wss") +
          geom_vline(xintercept = 4, linetype = 2) +
          labs(subtitle = "Elbow method")
      )

    }
    else if (ind == "average silhouette width") {
      print(
        fviz_nbclust(dfs, hcut, method = "silhouette") +
          labs(subtitle = "Silhouette method")
      )
    }
    else  if (ind == "gap statistics") {
      set.seed(1677)
      print(
        fviz_nbclust(
          dfs,
          hcut,
          nstart = 25,
          method = "gap_stat",
          nboot = nb
        ) +
          labs(subtitle = "Gap statistic method")
      )
    }
  }
  # nb <- NbClust(dfs, distance = dist, min.nc = minc,  max.nc = maxc, method = "kmeans", index=c(  "silhouette", "duda", "ratkowsky","mcclain",  "gplus",  "dunn")) #,
  # print(fviz_nbclust(nb) + theme_gray())
}


#' @title VisualOptClusters
#' @description visualization and statistics for clustering validation
#' @param       df_m: dataframe containing peaks and metadata
#' @param       meth: clustering algorithms ("hclust", default value), other values: "kmeans", "pam", "clara", "fanny", "hclust","agnes", "diana"
#' @param       dist:  distances ("euclidean", default value),  "maximum", "manhattan", "canberra", "binary" "minkowski"
#' @param       meth2: hc methods ("ward.D2", default value), average", "ward.D", "single", "complete", "mcquitty", "median" or "centroid"
#' @param       varCat1: categorical variable for choosing isolates, examples: "Taxonomie"	,"Genre",	"Date.d.analyse"	,"Origine","Ruche", "Nutrition"	, "Date.de.récolte"	, "Lieu.de.la.ruche"
#' @param       value: level of catVar1, examples: "Lactobacillus" ("Genre"), Taxonomie("Pediococcus pentosaceus"), "Erica cinerea" ("Nutrition"),...
#' @param       nc: number of clusters (nc=3, default value)
#' @return          figures and statistics
#' @examples        s<-VisualOptClusters(df_Peaks, varCat1="Genre", value="Enterococcus"),
#'                  s<-VisualOptClusters(df_Peaks, meth="pam", varCat1="Genre", value="Lactobacillus", nc=2)
#'                 s<-VisualOptClusters(df_Peaks, meth="hclust", meth2="average", dist="pearson",varCat1="Genre", value="All", nc=5),
#'                 s<-VisualOptClusters(df_Peaks, meth="kmeans", dist="euclidean",varCat1="Genre", value="All", nc=3)
#'        source:  https://www.rdocumentation.org/packages/factoextra/versions/1.0.7/topics/fviz_cluster
#'                 https://afit-r.github.io/kmeans_clustering
#'                 http://rstudio-pubs-static.s3.amazonaws.com/265632_3ad9e0b981244e15887677f8dffb39a0.html#using-30-different-indices
#'@export

  VisualOptClusters<- function(df_m, meth="hclust", meth2="ward.D2", dist="euclidean", varCat1, value, nc=3){
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
      stop("too small number of isolates...try another value")
  }
  else
  {
    mes <-  paste("The variable",  varCat1, "is not found in this dataframe ...try another variable" )
    stop(mes)
  }
  df_S<-scale(df_S)
  if (meth=="kmeans")  {
    out <- eclust(df_S, meth, k = nc,hc_metric =dist, graph = FALSE)
    print(fviz_silhouette(out) + scale_fill_brewer(palette = "Set1") +  scale_color_brewer(palette = "Set1") +  theme_gray()+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
    print(fviz_cluster(out, geom = "point", ellipse.type = "norm") )   }
  else if (meth=="pam") {
    out <- eclust(df_S, meth, k = nc,hc_metric =dist, hc_method = meth2, graph = FALSE)
    print(fviz_silhouette(out) + scale_fill_brewer(palette = "Set1") +  scale_color_brewer(palette = "Set1") +  theme_gray()+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
    print(fviz_cluster(out, geom = "point", ellipse.type = "norm") )}
  else if (meth=="hclust") {
    out <- eclust(df_S, meth, k = nc,hc_metric =dist, hc_method = meth2, graph = FALSE)
    print(fviz_silhouette(out) + scale_fill_brewer(palette = "Set1") +  scale_color_brewer(palette = "Set1") +  theme_gray()+theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()))
    print(fviz_dend(out, cex = 0.4, rect = FALSE, show_labels = FALSE))  }

 print(out$silinfo)
return(out$silinfo)
}

#' @title PointClusterVal
#' @description wrapper of functions for clustering tendency and validation statistics
#' @param       df_m: dataframe containing peaks and metadata
#' @param       meth: clustering algorithms ("hclust", default value), other values: "kmeans", "pam", "clara", "fanny", "hclust","agnes", "diana"
#' @param       dist:  distances ("euclidean", default value),  "maximum", "manhattan", "canberra", "binary" "minkowski"
#' @param       meth2: hc methods ("ward.D2", default value), average", "ward.D", "single", "complete", "mcquitty", "median" or "centroid"
#' @param       varCat1: categorical variable for choosing isolates, examples: "Taxonomie","Genre","Date.d.analyse","Origine","Ruche", "Nutrition", "Date.de.récolte", "Lieu.de.la.ruche"
#' @param       value: level of catVar1, examples: "Lactobacillus" ("Genre"), Taxonomie("Pediococcus pentosaceus"), "Erica cinerea" ("Nutrition"),...
#' @param       varCat2: categorical variable for partitioning the set of isolates chosen by using varCat1
#' @return      figures and statistics
#' @examples    ff<-PointClusterVal(df_Peaks, varCat1="Genre", value="Lactobacillus", varCat2="Nutrition")
#'        source:https://www.rdocumentation.org/packages/fpc/versions/2.2-8/topics/cluster.stats
#'               https://www.datanovia.com/en/lessons/cluster-validation-statistics-must-know-methods/
#'@export

PointClusterVal <- function(df_m, meth="hclust",  dist="euclidean",  meth2 = "ward.D2", varCat1, value, varCat2){

  ni <- 1
  nf <- 10
  f <- eval(parse(text = (paste0("df_m$", varCat1))))

  if (varCat1 == varCat2)
    stop("both variables cannot be equal")
  else
  {
    if (varCat1 %in% colnames(df_m)) {
      if (value == "All")
        chosenrows <- rep (TRUE, length(f))
      else {
        if (value %in% f) {
          chosenrows <- (f == value)
        }
        else
          stop("this value is not found ...try another")
      }
      if (length(which(chosenrows)) > 3)  {
        df_SS <- df_m[chosenrows,]
        df_S <- df_SS[,-ni:-nf]
        for (j in 1:ncol(df_S)) {
          df_S[, j] <- as.numeric(df_S[, j])
        }
      }
      else
        stop("too small number of isolates...try another value")
    }
    else
    {
      mes <-
        paste("The variable",
              varCat1,
              "is not found in this dataframe ...try another variable")
      stop(mes)
    }
    df_A <- df_SS[, ni:nf]
    if (varCat2 %in% colnames(df_A)) {
      cl <- eval(parse(text = (paste0(
        "df_A$", varCat2
      ))))
      cl <- as.numeric(factor (cl))
      Lcl <- length(unique(cl))
      if (Lcl > 50) {
        mes <- paste("The variable", varCat2, "has more than 50 levels!")
        warning(mes)
      }
    }
    else
      if (varCat2 %in% colnames(df_m)) {
        mes <-
          paste("The variable",
                varCat2,
                "is a peak intensity...try another variable")
        stop(mes)
      }
    else {
      mes <-
        paste("The variable",
              varCat2,
              "is not found in this dataframe ...try another variable")
      stop(mes)
    }
  }
  df_S <- scale(df_S)
  s <- get_clust_tendency(df_S, 40, graph = TRUE)
  print("Clustering Tendency (Hopkins statistic)")
  print(s)

  MD <- dist(df_S, method = dist)
  out <-
    eclust(
      df_S,
      FUNcluster = meth,
      k =  Lcl,
      hc_metric = dist,
      hc_method = meth2,
      graph = FALSE
    )
  out_stats <- cluster.stats(MD, out$cluster)
  print(out_stats)
  out_ext_stats <- cluster.stats(MD, cl, out$cluster)
  print("external clustering validation (corrected.rand, vi)")
  print("Corrected Rand index")
  print(out_ext_stats$corrected.rand)
  print("Meila variation of information (VI) index")
  print(out_ext_stats$vi)

return(c(s$hopkins_stat, out_stats$dunn,out_stats$dunn2,out_stats$pearsongamma,out_stats$avg.silwidth, out_ext_stats$corrected.rand,out_ext_stats$vi))

}
