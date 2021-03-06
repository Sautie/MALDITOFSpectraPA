% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClusVal.R
\name{OptClusters}
\alias{OptClusters}
\title{OptClusters}
\usage{
OptClusters(
  df_m,
  meth = "kmeans",
  dist = "euclidean",
  varCat1,
  value,
  minc = 2,
  maxc = 10,
  ind = "average silhouette width",
  nb = 100
)
}
\arguments{
\item{df_m:}{dataframe containing peaks and metadata}

\item{meth:}{clustering algorithms (meth="kmeans", default value), other values: "pam" or "hclust"}

\item{dist:}{distances ("euclidean", default),  "maximum", "manhattan", "canberra", "binary" "minkowski"}

\item{varCat1:}{categorical variable for choosing isolates, examples: "Taxonomie"    ,"Genre",    "Date.d.analyse"    ,"Origine","Ruche", "Nutrition"    , "Date.de.récolte"    , "Lieu.de.la.ruche"}

\item{value:}{level of catVar1, examples: "Lactobacillus" ("Genre"), Taxonomie("Pediococcus pentosaceus"), "Erica cinerea" ("Nutrition"),...}

\item{minc:}{minimal number of clusters (minc=2, default value)}

\item{maxc:}{maximal number of clusters (maxc=10, default value)}

\item{ind:}{methods to evaluate clustering algorithms:"total within sum of squares" ,"average silhouette width"  and "gap statistics"}

\item{nb:}{number of bootstrap samples (nb=100, default value)}
}
\value{
figures and statistics
}
\description{
OptClusters is a wrapper of several functions to visualize and compute optimal clusters for different clustering and evaluation methods
}
\examples{
     OptClusters(df_Peaks, varCat1="Taxonomie"	, value="Enterococcus faecalis"),
               OptClusters(df_Peaks, varCat1="Taxonomie"	, value="All")
               OptClusters(df_Peaks, meth="pam", varCat1="Taxonomie"	, value="All", ind="gap statistics"),
               OptClusters(df_Peaks, meth="hclust", varCat1="Taxonomie"	, value="All", ind="gap statistics")
       source: https://rpubs.com/pg2000in/OptimumClusters
               http://rstudio-pubs-static.s3.amazonaws.com/265632_3ad9e0b981244e15887677f8dffb39a0.html#using-30-different-indices
               https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/
               https://www.rdocumentation.org/packages/factoextra/versions/1.0.7/topics/fviz_nbclust
}
