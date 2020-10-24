library("ggplot2")
library("factoextra")
library("FactoMineR")
library("corrplot")


PCA_Clus<-function(df_m, varCat1, varCat1="Genre", value="Lactobacillus", meth="ward", dist="euclidean", graph="factorMapClus") {
# varCat1: categorical variable, value: level of the chosen categorical variable
# ni, nf: columns corresponding to categorical variables 
 #  dist: "euclidean" and "manhattan"
 # meth: "average" ([unweighted pair-]group [arithMetic] average method, aka ‘UPGMA’), "single" (single linkage), "complete" (complete linkage), "ward" (Ward's method), "weighted" (weighted average linkage, aka ‘WPGMA’), its generalization "flexible" which uses (a constant version of) the Lance-Williams formula
# graph: "dendf", "dendh", "factorMapf", "factorMaph", "factorMapClus", 
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

pca <- PCA(df_S, scale.unit = TRUE, ncp = 5, graph = FALSE)
hc <- HCPC(pca, nb.clust=-1, metric=dist, method=meth, graph=FALSE)
hc$desc.var$quanti  # contribution of variables for each cluster
hc$desc.axes$quanti # contribution of principal components for each cluster

if (graph=="dendf"){
  fviz_dend(hc, show_labels = FALSE)
  fviz_dend(hc,  cex = 0.6,  palette =  "Dark2", rect = TRUE, rect_fill = TRUE, rect_border = "Dark2")    
}
else if (graph=="dendh")
  plot(hc, choice="tree")
else if (graph=="factorMapf") {
  fviz_cluster(hc, geom = "point", main = "Factor map")
  fviz_cluster(hc, repel = TRUE,  show.clust.cent = TRUE, palette = "Dark2",  show_labels = FALSE, ggtheme = theme_gray(),main = "Factor map")
}  
else if (graph=="factorMaph") 
  plot(hc, choice="map", draw.tree=FALSE, ind.names=FALSE)
else if (graph=="factorMapClus") {
  plot(hc, choice="3D.map", angle=60, ind.names=FALSE)
  plot(hc, choice="3D.map", angle=80, ind.names=TRUE)
                   } 
 
return (hc$desc.axes$quanti)
}


SPCA<- function(df_m,  varCat1="Genre", value="Lactobacillus", VarCat2="Nutrition", contDim=TRUE, contVar= FALSE, contInd=FALSE){
  # varCat1, varCat2: categorical variable, value: level of the chosen categorical variable
  # varCat1, varCat2: Taxonomie,	Genre,	Date.d.analyse,	Origine,	Ruche, Nutrition,	Date.de.récolte, Lieu.de.la.ruche
  # ni, nf: columns corresponding to categorical variables 
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
    pca <- PCA(df_S, scale.unit = TRUE, ncp = 5, graph = FALSE)
    if (contDim)
      fviz_eig(pca, addlabels = TRUE, ylim = c(0, 50))
    if (contVar) {
       v <- get_pca_var(pca)
       corrplot(v$contrib, is.corr=FALSE) 
       fviz_pca_var(pca, col.var = "black")
    } 
    if (contInd)
      fviz_contrib(pca, choice = "ind", axes = 1:ncp)
    
    colors = c("gray48", "rosybrown", "firebrick", "olivedrab", "red", "blue", "green", "tan", "black", "orange", "yellow", "purple", "pink", "brown","sandybrown", "violet","cyan", "maroon1","lightsalmon","yellow4", "darkgreen","powderblue","gray","orangered","aquamarine", "steelblue","wheat4", "darkseagreen1","indianred1","palevioletred1", "turquoise1","maroon","goldenrod3","yellowgreen","orchid3","ivory2","lightskyblue","navajowhite","royalblue1","darkslategray4","purple4","rosybrown1","khaki3","deeppink","rosybrown3","burlywood4","hotpink1", "gray29","mediumspringgreen","mediumpurple1" )
    fviz_pca_ind(pca,  geom.ind = "point",   col.ind = cl,  palette = colors,  addEllipses = TRUE,  legend.title = "Groups" )
    fviz_pca_ind(pca, geom.ind = "point",   col.ind = cl,  palette = colors,  addEllipses = TRUE, ellipse.type = "convex",  legend.title = "Groups")
    
    
    dimVar<- dimdesc(pca, axes = 1:ncp, proba = 0.05)  # contribution of variables for ncp dimensions
    
    return(dimVar)
    
  }
  
  
  
  
  