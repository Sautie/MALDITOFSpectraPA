library(ggpubr)

MDS_Clus<-function(df_m, varCat1="Genre", value="Lactobacillus", dist="euclidean", nc=3, graph="mdsclus") {
  # varCat1: categorical variable, value: level of the chosen categorical variable
  # ni, nf: columns corresponding to categorical variables 
# dist: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
  # graph: mdsclus, lab_mdsclus, mds, lb_mds

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

d <- dist(df_S, method = dist, p=2)  # compute distance matrix
mds <- cmdscale(d)
colnames(mds) <- c("Dim.1", "Dim.2")

clu <- kmeans(mds, centers=nc)
groups<-as.factor(clu$cluster)
mds <- data.frame(mds, groups)

if (graph=="lab_mdsclus")
ggscatter(mds, x = "Dim.1", y = "Dim.2", label = rownames(df_S), color = "groups", palette = "Dark2",   size = 1, ellipse = TRUE,  ellipse.type = "convex",  repel = TRUE)
else if (graph=="mdsclus")
ggscatter(mds, x = "Dim11", y = "Dim.2",  color = "groups",  palette = "Dark2",   size = 1,  ellipse = TRUE, ellipse.type = "convex", repel = TRUE)
else if (graph=="lab_mds")
ggscatter(mds, x = "Dim.1", y = "Dim.2", label = rownames(df_S),  size = 1,  repel = TRUE)
else if (graph=="mds")
ggscatter(mds, x = "Dim.1", y = "Dim.2",  size = 1,  repel = TRUE)

}


SMDS<- function(df_m,  varCat1="Genre", value="Lactobacillus", VarCat2="Nutrition", contDim=TRUE, contVar= FALSE, contInd=FALSE){
  # varCat1, varCat2: categorical variable, value: level of the chosen categorical variable
  # varCat1, varCat2: Taxonomie,	Genre,	Date.d.analyse,	Origine,	Ruche, Nutrition,	Date.de.récolte, Lieu.de.la.ruche
  # ni, nf: columns corresponding to categorical variables 
  # dist: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
  ni <- 1
  nf <- 9
  if (varCat1 == varCat2)
    stop("both variables cannot be equal")
  else
  {
    df_A <- df_m[, ni:nf]
    if (varCat2 %in% colnames(df_A)) {
      cl <-  as.numeric(factor(df_A$varCat2))
      Lcl <- length(unique(cl))
      if (Lcl > 50) {
        mes <-paste("The variable",  as.chracter(varCat2),"has more than 50 levels!")
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

  d <- dist(df_S, method = dist, p=2)  # compute distance matrix
  mds <- cmdscale(d)
  colnames(mds) <- c("Dim.1", "Dim.2")
  groups<-as.factor(df_A$varCat2)
  mds <- data.frame(mds, groups)
  
  if (graph=="lab_mdsGroups")
    ggscatter(mds, x = "Dim.1", y = "Dim.2", label = rownames(df_m), color = "groups", palette = "Dark2",   size = 1, ellipse = TRUE,  ellipse.type = "convex",  repel = TRUE)
  else if (graph=="mdsGroups")
    ggscatter(mds, x = "Dim11", y = "Dim.2",  color = "groups",  palette = "Dark2",   size = 1,  ellipse = TRUE, ellipse.type = "convex", repel = TRUE)
  else if (graph=="lab_mds")
    ggscatter(mds, x = "Dim.1", y = "Dim.2", label = rownames(df_m),  size = 1,  repel = TRUE)
  else if (graph=="mds")
    ggscatter(mds, x = "Dim.1", y = "Dim.2",  size = 1,  repel = TRUE)
  
  }



