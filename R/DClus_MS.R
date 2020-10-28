library("ggpubr")
library("ggplot2")


#'@title MDS_Clus,
#'@description function for Multidimensional scaling and kmeans-based analysis of Maldi_Tof spectra
#'@param       df_m: dataframe containing peaks and metadata
#'@param       dist:  distances dist: "euclidean" (default value), "maximum", "manhattan", "canberra", "binary" or "minkowski"
#'@param       varCat1: categorical variable for choosing isolates, examples: "Taxonomie"	,"Genre",	"Date.d.analyse"	,"Origine","Ruche", "Nutrition"	, "Date.de.récolte"	, "Lieu.de.la.ruche"
#'@param       value: level of the chosen categorical variable catVar1, examples: "Lactobacillus" ("Genre"), Taxonomie("Pediococcus pentosaceus"), "Erica cinerea" ("Nutrition"),...
#'@param      grah: graphs: lab_mdsclus,(default value) lb_mds, mdsclaing
#'@return     figures and statistics
#'@examples    g<-MDS_Clus(df_Peaks, varCat1="Genre", value="Lactobacillus")
#'      source: https://rstudio-pubs-static.s3.amazonaws.com/274936_050c742fb3514bbaa87ce6ee2686af8c.html
#'               http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/122-multidimensional-scaling-essentials-algorithms-and-r-code/
#'               http://ugrad.stat.ubc.ca/R/library/mva/html/cmdscale.html
#'@export

MDS_Clus<-function(df_m, varCat1, value, dist="euclidean", nc=3, grah="lab_mdsclus") {
  # ni, nf: columns corresponding to categorical variables
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

d <- dist(df_S, method = dist)  # compute distance matrix, p:power for minkowski distance
mds <- cmdscale(d)
colnames(mds) <- c("Dim.1", "Dim.2")

clu <- kmeans(mds, centers=nc)
groups<-as.factor(clu$cluster)
mds <-data.frame(mds, groups)
print(mds)
print((mds[,1]))
colors = c("gray48", "darkseagreen1", "rosybrown", "firebrick", "olivedrab", "red", "blue", "green", "tan", "black", "orange", "yellow", "purple", "pink", "brown","sandybrown", "violet","cyan", "maroon1","lightsalmon","yellow4", "darkgreen","powderblue","gray","orangered","aquamarine", "steelblue","wheat4","indianred1","palevioletred1", "turquoise1","maroon","goldenrod3","yellowgreen","orchid3","ivory2","lightskyblue","navajowhite","royalblue1","darkslategray4","purple4","rosybrown1","khaki3","deeppink","rosybrown3","burlywood4","hotpink1", "gray29","mediumspringgreen","mediumpurple1" )
if (grah=="lab_mdsclus")
print(ggscatter(mds, x = "Dim.1", y = "Dim.2", label = rownames(mds), color = "groups", palette = "colors",   size = 1, ellipse = TRUE,  ellipse.type = "convex",  repel = TRUE))
else if (grah=="mdscaling"){
  print(ggplot(mds, aes(x=Dim.1, y=Dim.2)) +
    geom_text(label=rownames(mds)))
                               }
else if (grah=="lab_mds")
 print(ggscatter(mds, x = "Dim.1", y = "Dim.2", label = rownames(mds),  color = "groups", palette = "colors",  size = 1,  repel = TRUE))
}

#'@title SMDS
#'@description function for Multidimensional scaling and  external cluster-based analysis of Maldi_Tof spectra
#'@param       df_m: dataframe containing peaks and metadata
#'@param       dist:  distances: "euclidean" (default value), "maximum", "manhattan", "canberra", "binary" or "minkowski"
#'@param       varCat1: categorical variable for choosing isolates, examples: "Taxonomie"	,"Genre",	"Date.d.analyse"	,"Origine","Ruche", "Nutrition"	, "Date.de.récolte"	, "Lieu.de.la.ruche"
#'@param       value: level of the chosen categorical variable catVar1, examples: "Lactobacillus" ("Genre"), Taxonomie("Pediococcus pentosaceus"), "Erica cinerea" ("Nutrition"),...
#'@param      grah: graphs: "lab_mdsGroups" (default value), "mdsGroups"
#'@return     figures and statistics
#'@examples SMDS(df_Peaks, varCat1="Genre", value="Lactobacillus",  varCat2="Ruche",dist="euclidean",  grah="lab_mdsGroups")
#'           SMDS(df_Peaks, varCat1="Genre", value="Lactobacillus" varCat2="Ruche",)
#'           SMDS(df_Peaks, varCat1="Genre", value="Lactobacillus", varCat2="Nutrition")
#'           SMDS(df_Peaks, varCat1="Genre", value="All",varCat2="Taxonomie")
#'@export
#'
SMDS<- function(df_m, varCat1,  value, varCat2, dist="euclidean",  grah="lab_mdsGroups"){
  # ni, nf: columns corresponding to categorical variables
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
      cln <- eval(parse(text = (paste0(
        "df_A$", varCat2
      ))))
      cl <- as.numeric(factor (cln))
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

  d <- dist(df_S, method = dist, p = 2)  # compute distance matrix
  mds <- cmdscale(d)
  colnames(mds) <- c("Dim.1", "Dim.2")
  groups <- as.factor(cln)
  mds <- data.frame(mds, groups)
  colors = c("gray48", "darkseagreen1", "rosybrown", "firebrick", "olivedrab", "red", "blue", "green", "tan", "black", "orange", "yellow", "purple", "pink", "brown","sandybrown", "violet","cyan", "maroon1","lightsalmon","yellow4", "darkgreen","powderblue","gray","orangered","aquamarine", "steelblue","wheat4","indianred1","palevioletred1", "turquoise1","maroon","goldenrod3","yellowgreen","orchid3","ivory2","lightskyblue","navajowhite","royalblue1","darkslategray4","purple4","rosybrown1","khaki3","deeppink","rosybrown3","burlywood4","hotpink1", "gray29","mediumspringgreen","mediumpurple1" )
  if (grah=="lab_mdsGroups")
    print(
      ggscatter(
        mds,
        x = "Dim.1",
        y = "Dim.2",
        label = rownames(mds),
        color = "groups",
        palette = "colors",
        size = 1,
        ellipse = TRUE,
        ellipse.type = "convex",
        repel = TRUE
      )
    )
  else if (grah == "mdsGroups")
    print(
      ggscatter(
        mds,
        x = "Dim.1",
        y = "Dim.2",
        label = rownames(mds),
        color = "groups",
        palette = "colors",
        size = 1,
        repel = TRUE
      )
    )

  }



