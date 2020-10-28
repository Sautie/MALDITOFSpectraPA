library("pvclust")
library("ape")
library("factoextra")
library("cluster")
library("dendextend")
library("tidyverse")



#' @title fBHclus
#' @description function wrapper of pvclust for building dendrograms labeled with bootstrap probabilities for isolates chosen according to varCat1 levels
#' @param       df_m: dataframe containing peaks and metadata
#' @param       meth: hierarchical clustering algorithm ("ward2", default value), other values: "average", "ward.D", "single", "complete", "mcquitty", "median" or "centroid"
#' @param       dist:  distance ("euclidean", default value), "maximum", "manhattan", "canberra", "binary" "correlation", "uncentered", "abscor"
#' @param       varCat1: categorical variable, example: "Genre","Taxonomie", " Nutrition", "Ruche"...
#' @param       value: value of catVar1 "Lactobacillus" ("Genre"), Taxonomie("Pediococcus pentosaceus"), "Erica cinerea" ("Nutrition"),...
#' @param       nb: n of bootstrap iterations (nb=100, default)
#' @param       fig: boolean variable to indicate output figure (fig=TRUE, default value)
#' @param       ni,nf:  first and last columns corresponding to categorical variables  (default values, ni=1, nf=10)
#' @return       output cluster and figure
#' @examples     dft<-BHclus(df_Peaks, varCat1="Genre", value="All", nb=500),
#'               dft<-BHclus(df_Peaks, varCat1="Genre", value="Lactobacillus", nb=500)
#'               dft<-BHclus(df_Peaks, varCat1="Taxonomie", value="Pediococcus pentosaceus"),
#'               dft<-BHclus(df_Peaks, meth="complete", dist="canberra",varCat1="Taxonomie", value="Pediococcus pentosaceus")
#'       source: https://academic.oup.com/bioinformatics/article/22/12/1540/207339
#'               https://www.rdocumentation.org/packages/pvclust/versions/2.2-0/topics/pvclust
#'               https://www.rdocumentation.org/packages/stats/versions/3.2.1/topics/dist
#'@export

BHclus <- function(df_m, meth="ward.D2", dist="euclidean", varCat1, value, nb=100,  fig=TRUE, ni=1, nf= 10){

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

  set.seed(873)
  pv <- pvclust(t(df_S),
                method.hclust=meth,
                method.dist=dist, nboot=nb)
  if (fig) {
    plot(pv, hang = -1, cex = 0.5)
    pvrect(pv)
  }
  clus<- pvpick(pv)
  clus<-as.matrix(clus)
  return (clus)
}


#' @title Phyclus
#' @description builds different types of dendrograms for isolates chosen according to varCat1 levels
#' @param       df_m: dataframe containing peaks and metadata
#' @param       meth: hierarchical clustering algorithm ("ward2", default value), other values: "average", "ward.D", "single", "complete", "mcquitty", "median" or "centroid"
#' @param       dist:  distance ("euclidean", defaul), "euclidean", "maximum", "manhattan", "canberra", "binary" "minkowski"
#' @param       varCat1: categorical variable, example: "Genre","Taxonomie", " Nutrition", "Ruche"...
#' @param       value: value of catVar1 "Lactobacillus" ("Genre"), Taxonomie("Pediococcus pentosaceus"), "Erica cinerea" ("Nutrition"),...
#' @param       nc: n of clusters (nc=4, default)
#' @param       dendrogram: dendrogram type ("phylogram", default value), "cladogram", "unrooted", "fan" and "radial"
#' @param       ni,nf:  first and last columns corresponding to categorical variables  (default values, ni=1, nf=10)
#' @return       output cluster and figure
#' @examples    c<-Phyclus(df_Peaks, varCat1="Taxonomie",value="Pediococcus pentosaceus"),
#'              c<-Phyclus(df_Peaks, varCat1="Nutrition",value="Erica cinerea")
#'              c<-Phyclus(df_Peaks, varCat1="Taxonomie",value="Pediococcus pentosaceus", dendogram="cladogram")
#'              source: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust
#'                      https://www.rdocumentation.org/packages/ape/versions/5.4-1
#'@export

Phyclus <- function(df_m, meth="ward.D2", d="euclidean", varCat1, value,  dendrogram="phylogram", nc=4,  ni=1, nf=10){

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
  dd <- dist(scale(df_S), method = d)
  hc <- hclust(dd, method = meth)
  cl <- cutree(hc, nc)
  colors = c("red", "blue", "yellow4", "black", "orange", "yellow", "green","purple", "pink", "brown","sandybrown", "violet","cyan", "maroon1","lightsalmon", "darkgreen","powderblue","firebrick","gray","orangered","aquamarine", "steelblue","wheat4","olivedrab", "darkseagreen1","indianred1","palevioletred1", "turquoise1","maroon","goldenrod3","yellowgreen","orchid3","ivory2","lightskyblue","gray48", "rosybrown","tan","navajowhite","royalblue1","darkslategray4","purple4","rosybrown1","khaki3","deeppink","rosybrown3","burlywood4","hotpink1", "gray29","mediumspringgreen","mediumpurple1" )
  if (dendrogram=="phylogram")
    plot(as.phylo(hc), edge.color = "black", tip.color = colors[cl], cex = 0.6, label.offset = 1)
  else if (dendrogram=="cladogram")
    plot(as.phylo(hc), edge.color = "black",  type = "cladogram", tip.color = colors[cl], cex = 0.6, label.offset = 1)
  else if (dendrogram=="unrooted")
    plot(as.phylo(hc), edge.color = "black",  type = "unrooted", tip.color = colors[cl], cex = 0.6,   no.margin = TRUE)
  else if (dendrogram=="fan")
    plot(as.phylo(hc), edge.color = "black",  type = "fan", tip.color = colors[cl], cex=0.6)
  else if (dendrogram=="radial")
    plot(as.phylo(hc), type = "radial", tip.color = colors[cl], cex=0.6)

  coph <- cophenetic(hc)
  ccoph<-cor(dd, coph)
  print("Cophenetic coefficient")
  print(ccoph)

  return(ccoph)
}

#' @title PhyclusVar
#' @description builds different types of dendrograms for isolates selected and categorized based on varCat2 and varCat1 levels
#' @param       df_m: dataframe containing peaks and metadata
#' @param       meth: hierarchical clustering algorithm ("ward2", default value), other values: "average", "ward.D", "single", "complete", "mcquitty", "median" or "centroid"
#' @param       dist:  distance ("euclidean", default value), "euclidean", "maximum", "manhattan", "canberra", "binary" "minkowski"
#' @param       varCat1: categorical variable, example: "Taxonomie"	,"Genre",	"Date.d.analyse"	,"Origine","Ruche", "Nutrition"	, "Date.de.récolte"	, "Lieu.de.la.ruche"
#' @param       value: value of catVar1 "Lactobacillus" ("Genre"), Taxonomie("Pediococcus pentosaceus"), "Erica cinerea" ("Nutrition"),...
#' @param       varCat2: categorical variable
#' @param       nc: n of clusters (nc=4, default value)
#' @param       dendrogram: dendrogram type ("phylogram", default value), "cladogram", "unrooted", "fan" and "radial"
#' @param       ni,nf:  first and last columns corresponding to categorical variables  (default values, ni=1, nf=10)
#' @return       output cophenetic coefficient and figure
#'  @examples    dft<-PhyclusVar(df_Peaks, varCat1="Genre", value="Lactobacillus", varCat2="Nutrition"),
#'               dft<-PhyclusVar(df_Peaks, varCat1="Taxonomie", value="Lactobacillus plantarum", varCat2="Ruche")
#'               dft<-PhyclusVar(df_Peaks, varCat1="Taxonomie", value="Pediococcus pentosaceus", varCat2="Nutrition", dendrogram="cladogram"),
#'               dft<-PhyclusVar(df_Peaks, varCat1="Taxonomie", value="All", varCat2="Date.de.récolte", dendrogram="fan")
#'               dft<-PhyclusVar(df_Peaks, varCat1="Taxonomie", value="Lactobacillus plantarum", varCat2="Ruche", dendrogram="fan"),
#'               dft<-PhyclusVar(df_Peaks, varCat1="Taxonomie", value="Lactobacillus plantarum", varCat2="Ruche", dendrogram="unrooted")
#'              source: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust
#'                      https://www.rdocumentation.org/packages/ape/versions/5.4-1
#'@export

PhyclusVar <- function(df_m, meth="ward.D2", d="euclidean", varCat1, value, varCat2,  dendrogram="phylogram", ni=1, nf=10){

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
      cl <- as.numeric(factor (cl))
      Lcl <- length(unique(cl))
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
  dd <- dist(scale(df_S), method = d)
  hc <- hclust(dd, method = meth)
  colors = c("red", "blue", "maroon", "black", "orange",  "purple",  "brown", "firebrick", "yellowgreen", "violet","maroon1","darkgreen",  "sandybrown", "cyan",  "lightsalmon","yellow4", "powderblue", "yellow","green", "pink","gray","orangered","aquamarine", "steelblue","wheat4","olivedrab", "darkseagreen1","indianred1","palevioletred1", "turquoise1","goldenrod3","orchid3","ivory2","lightskyblue","gray48", "rosybrown","tan","navajowhite","royalblue1","darkslategray4","purple4","rosybrown1","khaki3","deeppink","rosybrown3","burlywood4","hotpink1", "gray29","mediumspringgreen","mediumpurple1" )
  if (dendrogram=="phylogram")
    plot(as.phylo(hc), edge.color = "black", tip.color = colors[cl], cex = 0.8, label.offset = 1)
  else if (dendrogram=="cladogram")
    plot(as.phylo(hc), edge.color = "black",  type = "cladogram", tip.color = colors[cl], cex = 0.6, label.offset = 1)
  else if (dendrogram=="unrooted")
    plot(as.phylo(hc), edge.color = "black",  type = "unrooted", tip.color = colors[cl], cex = 0.6,  no.margin = TRUE) #
  else if (dendrogram=="fan")
    plot(as.phylo(hc), edge.color = "black",  type = "fan", tip.color = colors[cl], cex = 0.6)
  else if (dendrogram=="radial")
    plot(as.phylo(hc), type = "radial", tip.color = colors[cl], cex=0.6)

  coph <- cophenetic(hc)
  ccoph<-cor(dd, coph)
  print("Cophenetic coefficient")
  print(ccoph)
  return(ccoph)
}

#' @title Vizclus
#' @description clusters and builds dendrograms for isolates selected according to varCat1 levels
#' @param       df_m: dataframe containing peaks and metadata
#' @param       meth: hierarchical clustering algorithm ("ward2", default value), other values: "average", "ward.D", "single", "complete", "mcquitty", "median" or "centroid"
#' @param       dist:  distance ("euclidean", default value), "euclidean", "maximum", "manhattan", "canberra", "binary" "minkowski"
#' @param       varCat1: categorical variable for choosing isolates, examples: "Taxonomie"	,"Genre",	"Date.d.analyse"	,"Origine","Ruche", "Nutrition"	, "Date.de.récolte"	, "Lieu.de.la.ruche"
#' @param       value: level of catVar1, examples, "Lactobacillus" ("Genre"), Taxonomie("Pediococcus pentosaceus"), "Erica cinerea" ("Nutrition"),...
#' @param       nc: n of clusters (nc=4, default value)
#' @param       dendrogram: dendrogram and factor map
#' @param       ni,nf:  first and last columns corresponding to categorical variables  (default values, ni=1, nf=10)
#' @return       output figures and statistics
#'  @examples    Lt<-Vizclus(df_Peaks,varCat1="Taxonomie", value="Pediococcus pentosaceus"),
#'               Lt<-Vizclus(df_Peaks,varCat1="Ruche", value="VI")
#'               Lt<-Vizclus(df_Peaks,varCat1="Nutrition", value="Caraway", nc=2),
#'               Lt<-Vizclus(df_Peaks,varCat1="Genre", value="Lactobacillus", nc=4)
#'               Lt<-Vizclus(df_Peaks,varCat1="Nutrition", value="All"),
#'               Lt<-Vizclus(df_Peaks,varCat1="Nutrition", value="Caraway", nc=3, graph="fm")
#'      source: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust
#'              https://www.rdocumentation.org/packages/factoextra/versions/1.0.7/topics/fviz
#'              http://www.sthda.com/english/wiki/factoextra-r-package-easy-multivariate-data-analyses-and-elegant-visualization#:~:text=factoextra%20is%20an%20R%20package,exploratory%20multivariate%20data%20analyses%2C%20including%3A&text=Factor%20Analysis%20of%20Mixed%20Data,both%20quantitative%20and%20qualitative%20variables.
#'              http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
#'@export

Vizclus <- function(df_m, meth="ward.D2", dist="euclidean", varCat1, value, nc=4, graph="phylogram", ni=1, nf=10)
{
  f<-eval(parse(tex=(paste0("df_m$", varCat1))))
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

dd <- dist(scale(df_S), method = dist)
hc <- hclust(dd, method = meth)
# print(hc$labels)

cl <- cutree(hc, k = nc)
print("clustering vector: cluster assignment to each isolate")
print(cl)

tclust<-table(cl)
print("The size of each cluster")
print(tclust)

coph <- cophenetic(hc)
print("Cophenetic coeficient")
cc<-cor(dd, coph)
print(cc)

indc<-1:nc
colors = c("red", "blue", "green", "black", "orange",  "purple", "yellow", "pink", "brown","sandybrown", "violet","cyan", "maroon1","lightsalmon","yellow4", "darkgreen","powderblue","firebrick","gray","orangered","aquamarine", "steelblue","wheat4","olivedrab", "darkseagreen1","indianred1","palevioletred1", "turquoise1","maroon","goldenrod3","yellowgreen","orchid3","ivory2","lightskyblue","gray48", "rosybrown","tan","navajowhite","royalblue1","darkslategray4","purple4","rosybrown1","khaki3","deeppink","rosybrown3","burlywood4","hotpink1", "gray29","mediumspringgreen","mediumpurple1" )
if (graph=="phylogram")
   print(fviz_dend(hc, k = nc,  cex = 0.5, k_colors = colors[1:nc], color_labels_by_k = TRUE,  rect = TRUE,  horiz = TRUE ))
else
   print(fviz_cluster(list(data = df_S, cluster = cl),  palette = "jco", ellipse.type = "convex",  repel = TRUE, show.clust.cent = FALSE, labelsize = 9, ggtheme = theme_minimal()))
listresults <- list()
listresults[[1]] <- hc
listresults[[2]] <- nc
listresults[[3]] <- cl
return(listresults)
}


#' @title Dendrogram_pairComp
#' @description  Dendrogram_pairComp makes pairwise dendrogram comparisons and correlations
#' @param       df_m: dataframe containing peaks and metadata
#' @param       methCorr: dendrogram correlation method ("cophenetic", default value), others: "baker", "common_nodes", "FM_index"
#' @param       meth1,meth2...meth6: hierarchical clustering algorithms "ward2" ,"average", "ward.D", "single", "complete", "mcquitty", "median" or "centroid"
#' @param       dist1,dist2,...dist6:  distances "euclidean", "maximum", "manhattan", "canberra", "binary" "minkowski"
#' @param       varCat1: categorical variable for choosing isolates, examples: "Taxonomie"	,"Genre",	"Date.d.analyse"	,"Origine","Ruche", "Nutrition"	, "Date.de.récolte"	, "Lieu.de.la.ruche"
#' @param       value: level of catVar1, examples: "Lactobacillus" ("Genre"), Taxonomie("Pediococcus pentosaceus"), "Erica cinerea" ("Nutrition"),...
#' @return       untangled and tangled dendrograms, dendrogram correlation matrix and statistics
#'  @examples    mc<-Dendrogram_pairComp(df_Peaks, varCat1="Taxonomie", value="Pediococcus pentosaceus"),
#'               mc<-Dendrogram_pairComp(df_Peaks, varCat1="Nutrition", value="Marrubium vulgare")
#'               source: https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html#:~:text=The%20dendextend%20package%20offers%20a,its%20branches%2C%20nodes%20and%20labels.
#'                      https://www.rdocumentation.org/packages/dendextend/versions/1.14.0
#'                      https://www.r-graph-gallery.com/340-custom-your-dendrogram-with-dendextend.html
#'                      https://academic.oup.com/bioinformatics/article/31/22/3718/240978
#'                      https://cran.rstudio.com/web/packages/dendextend/vignettes/Cluster_Analysis.html
#'                      https://www.datanovia.com/en/lessons/comparing-cluster-dendrograms-in-r/
#'                      https://rdrr.io/cran/dendextend/man/untangle.html
#'@export

Dendrogram_pairComp<- function(df_m, varCat1, value, methCorr="cophenetic",meth1="ward.D2", meth2="mcquitty",
                    meth3="complete", meth4="centroid", meth5="single", meth6="average",
                    dist1="euclidean", dist2="euclidean", dist3="euclidean", dist4="euclidean", dist5="euclidean", dist6="euclidean", graph="Complex", nc=3)
{
  ni <- 1
  nf <- 10
  f <- eval(parse(tex = (paste0("df_m$", varCat1))))
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
    mes <-
      paste("The variable",
            varCat1,
            "is not found in this dataframe ...try another variable")
    stop(mes)
  }

  hcA <- hclust((dist(df_S, method = dist1)), method = meth1)
  hcB <- hclust((dist(df_S, method = dist2)), method = meth2)
  dend1 <- as.dendrogram((hcA))
  dend2 <- as.dendrogram((hcB))
  dends <- dendlist(dend1, dend2)
  if (graph == "Complex")
    tanglegram(dend1, dend2)
  else {
    de <- dendlist(
      dend1 %>%
        set(
          "labels_col",
          value = c("skyblue", "orange", "grey"),
          k = nc
        ) %>%
        set("branches_lty", 1) %>%
        set(
          "branches_k_color",
          value = c("skyblue", "orange", "grey"),
          k = nc
        ),
      dend2 %>%
        set(
          "labels_col",
          value = c("skyblue", "orange", "grey"),
          k = nc
        ) %>%
        set("branches_lty", 1) %>%
        set(
          "branches_k_color",
          value = c("skyblue", "orange", "grey"),
          k = nc
        )

    )

    tanglegram(
      de,
      common_subtrees_color_lines = FALSE,
      highlight_distinct_edges  = TRUE,
      highlight_branches_lwd = FALSE,
      margin_inner = 7,
      lwd = 2
    )
  }

  if (graph == "Complex") {
    tg <- untangle(dends, method = "step1side")
  }
  else{
    tg <- untangle(de, method = "step1side")
  }
  tanglegram(tg, main_left = "untangled dendrograms")
  print("Alignment score")
  print(entanglement(tg))
  print("Cophenetic correlation coefficient")
  print(cor_cophenetic(dend1, dend2))
  print("Baker correlation coefficient")
  print(cor_bakers_gamma(dend1, dend2))


  hcD <- hclust((dist(df_S, method = dist3)), method = meth3)
  hcE <- hclust((dist(df_S, method = dist4)), method = meth4)
  hcF <- hclust((dist(df_S, method = dist5)), method = meth5)
  hcG <- hclust((dist(df_S, method = dist6)), method = meth6)
  dend3 <- as.dendrogram(hcD)
  dend4 <- as.dendrogram(hcE)
  dend5 <- as.dendrogram(hcF)
  dend6 <- as.dendrogram(hcG)
  dendMultiples <-
    dendlist(
      n1 = dend1,
      n2 = dend2,
      n3 = dend3,
      n4 = dend4,
      n5 = dend5,
      n6 = dend6
    )
  print(paste("n1:", meth1, dist1))
  print(paste("n2:", meth2, dist2))
  print(paste("n3:", meth3, dist3))
  print(paste("n4:", meth4, dist4))
  print(paste("n5:", meth5, dist5))
  print(paste("n6:", meth6, dist6))
  print("Tree correlation matrix")
  cors <- cor.dendlist(dendMultiples, method = methCorr)
  print(round(cors, 2))
  return (cors)
}



