library("pvclust")
library("ape")
library("factoextra")
library("cluster")


# default values dft<-BHclus(df)
# dft<-BHclus(df_Peaks, varCat1="Genre", value="All")
# dft<-BHclus(df_Peaks, varCat1="Genre", value="Lactobacillus")
# dft<-BHclus(df_Peaks, varCat1="Taxonomie", value="Pediococcus pentosaceus")
# dft<-BHclus(df_Peaks, meth="complete", dist="canberra",varCat1="Taxonomie", value="Pediococcus pentosaceus")
# meth: "average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median" or "centroid"
# dist:"euclidean", "maximum", "manhattan", "canberra", "binary" "minkowski" "correlation", "uncentered", "abscor"
# Genre Taxonomie Nutrition Ruche   varCat1="Genre", value="Lactobacillus" pspe="Pediococcus pentosaceus", fl="Erica cinerea", ruch="VI",

BHclus <- function(df_m, meth="ward.D2", dist="euclidean", varCat1, value, nb=100,  fig=TRUE){
  
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
namefile<-paste(varCat1, "txt", sep = ".")
write.table(clus, namefile)
return (clus)
}

# dist:"euclidean", "maximum", "manhattan", "canberra", "binary" "minkowski" 
# meth:"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
# Genre Taxonomie Nutrition Ruche   varCat1="Genre", value="Lactobacillus" pspe="Pediococcus pentosaceus", fl="Erica cinerea", ruch="VI"
# c<-Phyclus(df_Peaks, varCat1="Taxonomie",value="Pediococcus pentosaceus")
# c<-Phyclus(df_Peaks, varCat1="Nutrition",value="Erica cinerea") 
Phyclus <- function(df_m, meth="ward.D2", d="euclidean", varCat1, value,  dendrogram="phylogram", nc=4){
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

  dd <- dist(scale(df_S), method = d)
  hc <- hclust(dd, method = meth)
  cl <- cutree(hc, nc)
  colors = c("red", "blue", "green", "black", "orange", "yellow", "purple", "pink", "brown","sandybrown", "violet","cyan", "maroon1","lightsalmon","yellow4", "darkgreen","powderblue","firebrick","gray","orangered","aquamarine", "steelblue","wheat4","olivedrab", "darkseagreen1","indianred1","palevioletred1", "turquoise1","maroon","goldenrod3","yellowgreen","orchid3","ivory2","lightskyblue","gray48", "rosybrown","tan","navajowhite","royalblue1","darkslategray4","purple4","rosybrown1","khaki3","deeppink","rosybrown3","burlywood4","hotpink1", "gray29","mediumspringgreen","mediumpurple1" )
  if (dendrogram=="phylogram")
    plot(as.phylo(hc), edge.color = "black", tip.color = colors[cl], cex = 0.5, label.offset = 1)  
  else if (dendrogram=="cladogram")
    plot(as.phylo(hc), edge.color = "black",  type = "cladogram", tip.color = colors[cl], cex = 0.5, label.offset = 1)
  else if (dendrogram=="unrooted")
    plot(as.phylo(hc), edge.color = "black",  type = "unrooted", tip.color = colors[cl], cex = 0.6,   no.margin = TRUE)
  else if (dendrogram=="fan")
    plot(as.phylo(hc), edge.color = "black",  type = "fan", tip.color = colors[cl], cex=0.5)
  else if (dendrogram=="radial")
    plot(as.phylo(hc), type = "radial", tip.color = colors[cl])
  
  coph <- cophenetic(hc)
  ccoph<-cor(dd, coph)
  return(ccoph)
}


# dft<-PhyclusVar(df_Peaks, varCat1="Genre", value="Lactobacillus", varCat2="Nutrition")
# dft<-PhyclusVar(df_Peaks, varCat1="Taxonomie", value="Lactobacillus plantarum", varCat2="Ruche")

PhyclusVar <- function(df_m, meth="ward.D2", d="euclidean", varCat1, value, varCat2,  dendrogram="phylogram"){
                    
  # varCat1, varCat2: Taxonomie	Genre	Date.d.analyse	Origine	Ruche	Nutrition	Date.de.récolte	Lieu.de.la.ruche
  # value:"Pediococcus pentosaceus",.. Nutrition: "Erica cinerea", Ruche:"VI", 
  ni <- 1
  nf <- 10
  f<-eval(parse(text=(paste0("df_m$", varCat1))))
  
  if (varCat1 == varCat2)
    stop("both variables cannot be equal")
  else
  {
    df_A <- df_m[, ni:nf]
    # for (j in 1:ncol(df_A)) {
    #   df_A[, j] <- as.numeric(factor (df_A[, j]))
    # }
    if (varCat2 %in% colnames(df_A)) {
      cl <- eval(parse(text=(paste0("df_A$", varCat2)))) 
      cl <- as.numeric(factor (cl))
      Lcl <- length(unique(cl))
      if (Lcl > 50) {
        mes <-
          paste("The variable", varCat2,"has more than 50 levels!")
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
    
  dd <- dist(scale(df_S), method = d)
  hc <- hclust(dd, method = meth)
  colors = c("red", "blue", "green", "black", "orange", "yellow", "purple", "pink", "brown","sandybrown", "violet","cyan", "maroon1","lightsalmon","yellow4", "darkgreen","powderblue","firebrick","gray","orangered","aquamarine", "steelblue","wheat4","olivedrab", "darkseagreen1","indianred1","palevioletred1", "turquoise1","maroon","goldenrod3","yellowgreen","orchid3","ivory2","lightskyblue","gray48", "rosybrown","tan","navajowhite","royalblue1","darkslategray4","purple4","rosybrown1","khaki3","deeppink","rosybrown3","burlywood4","hotpink1", "gray29","mediumspringgreen","mediumpurple1" )
  if (dendrogram=="phylogram")
    plot(as.phylo(hc), edge.color = "black", tip.color = colors[cl], cex = 0.6, label.offset = 1)  
  else if (dendrogram=="cladogram")
    plot(as.phylo(hc), edge.color = "black",  type = "cladogram", tip.color = colors[cl], cex = 0.6, label.offset = 1)
  else if (dendrogram=="unrooted")
    plot(as.phylo(hc), edge.color = "black",  type = "unrooted", tip.color = colors[cl], cex = 0.6,   no.margin = TRUE)
  else if (dendrogram=="fan")
    plot(as.phylo(hc), edge.color = "black",  type = "fan", tip.color = colors[cl])
  else if (dendrogram=="radial")
    plot(as.phylo(hc), type = "radial", tip.color = colors[cl])
  
  coph <- cophenetic(hc)
  ccoph<-cor(dd, coph)
  }
 
  return(ccoph)
}
# pspe="Pediococcus pentosaceus", fl="Erica cinerea", ruch="VI",vcol="genre", pgen="Lactobacillus"
# (df_m, meth="ward.D2", d="euclidean", varCat1, value, varCat2,  dendrogram="phylogram")
# Lt<-Vizclus(df_Peaks,varCat1="Taxonomie", value="Pediococcus pentosaceus") 
# Lt<-Vizclus(df_Peaks,varCat1="Ruche", value="VI")
# Lt<-Vizclus(df_Peaks,varCat1="Nutrition", value="Caraway", nc=2)
# Lt<-Vizclus(df_Peaks,varCat1="Nutrition", value="All")
Vizclus <- function(df_m, meth="ward.D2", dist="euclidean", varCat1, value, nc=4, graph="dendrogram") 
{
  ni <- 1
  nf <- 10
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
      stop("too small number of individuals...try another value")
  }
  else
  {
    mes <-  paste("The variable",  varCat1, "is not found in this dataframe ...try another variable" )
    stop(mes)
  }

dd <- dist(scale(df_S), method = dist)
hc <- hclust(dd, method = meth)
cl <- cutree(hc, k = nc)
print("Individuals of each cluster")
print(cl)

tclust<-table(cl)
print("The size of each cluster")
print(tclust)

coph <- cophenetic(hc)
print("Cophenetic coeficient")
cc<-cor(dd, coph)
print(cc)

indc<-1:nc
colors = c("red", "blue", "green", "black", "orange", "yellow", "purple", "pink", "brown","sandybrown", "violet","cyan", "maroon1","lightsalmon","yellow4", "darkgreen","powderblue","firebrick","gray","orangered","aquamarine", "steelblue","wheat4","olivedrab", "darkseagreen1","indianred1","palevioletred1", "turquoise1","maroon","goldenrod3","yellowgreen","orchid3","ivory2","lightskyblue","gray48", "rosybrown","tan","navajowhite","royalblue1","darkslategray4","purple4","rosybrown1","khaki3","deeppink","rosybrown3","burlywood4","hotpink1", "gray29","mediumspringgreen","mediumpurple1" )
if (graph=="dendrogram")
print(fviz_dend(hc, k = nc,  cex = 0.5, k_colors = colors[1:nc], color_labels_by_k = TRUE,  rect = TRUE,  horiz = TRUE ))
else
print(fviz_cluster(list(data = df_S, cluster = cl),  palette = "jco", ellipse.type = "convex",  repel = TRUE, show.clust.cent = FALSE, labelsize = 7, ggtheme = theme_minimal()))
listresults <- list()
   listresults[[1]] <-hc
   listresults[[2]] <-nc
   listresults[[3]] <-cl
return(listresults)
}






