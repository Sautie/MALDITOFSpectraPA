library("ggplot2")
library("factoextra")
library("FactoMineR")
library("corrplot")


#' @title PCA_Clus
#' @description Principal Component Analysis and clustering of Maldi_Tof spectra
#' @param       df_m: dataframe containing peaks and metadata
#' @param       dist:  distances, "euclidean", (default value),  "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" or "kendall".
#' @param       varCat1: categorical variable for choosing isolates, examples: "Taxonomie"	,"Genre",	"Date.d.analyse"	,"Origine","Ruche", "Nutrition"	, "Date.de.récolte"	, "Lieu.de.la.ruche"
#' @param       value: level of catVar1, examples: "Lactobacillus" ("Genre"), Taxonomie("Pediococcus pentosaceus"), "Erica cinerea" ("Nutrition"),...
#' @param       meth: clustering method, ward(default), "average" ,"single","complete‘
#' @param      graph: visual analysis, "dendf", "dendh", (dendograms) "factorMapf", "factorMaph", "factorMapClus",(factor maps) (default value, graph="factorMapClus" )
#' @param      pc: number of principal components (pc=3, default value)
#' @param      ni,nf:  first and last columns corresponding to categorical variables  (default values, ni=1, nf=10)
#' @return     figures and statistics
#' @examples  pc<-PCA_Clus(df_Peaks, varCat1="Genre", value="Lactobacillus"),
#'            pc<-PCA_Clus(df_Peaks, varCat1="Genre", value="Lactobacillus", graph="dendh")
#'            pc<-PCA_Clus(df_Peaks, varCat1="Genre", value="Lactobacillus", graph="dendf"),
#'            pc<-PCA_Clus(df_Peaks, varCat1="Genre", value="Lactobacillus", graph="factorMapf")
#'     source: https://www.rdocumentation.org/packages/FactoMineR/versions/2.2/topics/PCA
#'            https://www.rdocumentation.org/packages/FactoMineR/versions/2.2/topics/HCPC
#'            http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/117-hcpc-hierarchical-clustering-on-principal-components-essentials/
#'            https://rpkgs.datanovia.com/factoextra/
#'            http://factominer.free.fr/factomethods/hierarchical-clustering-on-principal-components.html
#'@export

PCA_Clus<-function(df_m, varCat1, value, meth="ward", dist="euclidean", graph="factorMapClus", pc=3, ni=1, nf=10) {

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


pca <- PCA(df_S, scale.unit = TRUE, ncp = pc, graph = FALSE)
hc <- HCPC(pca, nb.clust=-1, metric=dist, method=meth, graph=FALSE)
print("contribution of variables for each cluster")
print(hc$desc.var$quanti)  # contribution of variables for each cluster
print("contribution of principal components for each cluster")
print(hc$desc.axes$quanti) # contribution of principal components for each cluster

if (graph=="dendf"){
 print( fviz_dend(hc, show_labels = FALSE))
  print(fviz_dend(hc,  cex = 0.6,  palette =  "Dark2", rect = TRUE, rect_fill = TRUE, rect_border = "Dark2"))
}
else if (graph=="dendh")
  plot(hc, choice="tree")
else if (graph=="factorMapf") {
  print(fviz_cluster(hc, geom = "point",  labelsize = 1, main = "Factor map"))
  print(fviz_cluster(hc, repel = TRUE,  show.clust.cent = TRUE, palette = "Dark2",  show_labels = FALSE, ggtheme = theme_gray(),main = "Factor map"))
}
else if (graph=="factorMaph")
  plot(hc, choice="map", draw.tree=FALSE, ind.names=FALSE)
else if (graph=="factorMapClus") {
  plot(hc, choice="3D.map", angle=80, ind.names=TRUE)
  plot(hc, choice="3D.map", angle=60, ind.names=FALSE)
                   }

return (hc$desc.axes$quanti)
}

#' @title SPCA
#' @description function for Principal Component Analysis of Maldi_Tof spectra with external clusters based on metadata information
#' @param       df_m: dataframe containing peaks and metadata
#' @param       varCat1: categorical variable for choosing isolates, examples: "Taxonomie"	,"Genre",	"Date.d.analyse"	,"Origine","Ruche", "Nutrition"	, "Date.de.récolte"	, "Lieu.de.la.ruche"
#' @param       value: level of catVar1, examples: "Lactobacillus" ("Genre"), Taxonomie("Pediococcus pentosaceus"), "Erica cinerea" ("Nutrition"),...
#' @param       varCat2: categorical variable for partitioning the chosen isolates, examples: Taxonomie,	Genre,	Date.d.analyse,	Origine,	Ruche, Nutrition,	Date.de.récolte, Lieu.de.la.ruche
#' @param       contDim:  graph and statistics for PC contributions (default, contDim=TRUE)
#' @param       contVar:  graph and statistics for variable contributions(default,contVar= FALSE)
#' @param       contInd:   graph and statistics for isolate contributions (default, contInd=TRUE)
#' @examples       hg<-SPCA(df_Peaks, varCat1="Genre", value="Lactobacillus", varCat2="Nutrition")
#'@export
#'
SPCA<- function(df_m,  varCat1, value, varCat2, contDim=TRUE, contVar= FALSE, contInd=FALSE){
  # ni, nf: columns corresponding to categorical variables
  ni <- 1
  nf <- 10
  f<-eval(parse(text=(paste0("df_m$", varCat1))))

  if (varCat1 == varCat2)
    stop("both variables cannot be equal")
  else
  {
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
        df_SS <- df_m[chosenrows, ]
        df_S <- df_SS[, -ni:-nf]
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
    df_A <- df_SS[, ni:nf]
    if (varCat2 %in% colnames(df_A)) {
      cl <- eval(parse(text=(paste0("df_A$", varCat2))))
      cln <- as.numeric(factor (cl))
      Lcl <- length(unique(cln))
      if (Lcl > 50) {
        mes <-paste("The variable", varCat2,"has more than 50 levels!")
        warning(mes)
      }
    }
    else
      if (varCat2 %in% colnames(df_m)) {
        mes <-paste("The variable", varCat2,  "is a peak intensity...try another variable")
        stop(mes)
      }
    else {
      mes <- paste("The variable", varCat2, "is not found in this dataframe ...try another variable" )
      stop(mes)
    }
  }
   ncp<-5
    pca <- PCA(df_S, scale.unit = TRUE, ncp, graph = FALSE)
    if (contDim)
      print(fviz_eig(pca, addlabels = TRUE, ylim = c(0, 50)))
    if (contVar) {
       v <- get_pca_var(pca)
      # print(corrplot(v$contrib, is.corr=FALSE) )
      print(fviz_pca_var(pca, col.var = "black"))
    }
    if (contInd){
      print(fviz_contrib(pca, choice = "ind", axes = 1:ncp))
      print(fviz_pca_ind(pca, pointsize=4,col.ind="contrib", geom = "point") +
      scale_color_gradient2(low="lightsalmon", mid="maroon",
                            high=,"indianred1", midpoint=4) ) }

    colors = c("gray48", "darkseagreen1", "rosybrown", "firebrick", "olivedrab", "red", "blue", "green", "tan", "black", "orange", "yellow", "purple", "pink", "brown","sandybrown", "violet","cyan", "maroon1","lightsalmon","yellow4", "darkgreen","powderblue","gray","orangered","aquamarine", "steelblue","wheat4","indianred1","palevioletred1", "turquoise1","maroon","goldenrod3","yellowgreen","orchid3","ivory2","lightskyblue","navajowhite","royalblue1","darkslategray4","purple4","rosybrown1","khaki3","deeppink","rosybrown3","burlywood4","hotpink1", "gray29","mediumspringgreen","mediumpurple1" )
    p<-fviz_pca_ind(pca,  geom.ind = "point", label="none",  habillage=factor(cl),  addEllipses = TRUE,  legend.title = "Groups")
    print(p + scale_color_brewer(palette="Dark2") + theme_minimal())

    dimVar<- dimdesc(pca, axes = 1:ncp, proba = 0.05)  # contribution of variables for ncp dimensions
    print(dimVar)
    return(dimVar)

  }




