library("pvclust")
library("ape")
library("factoextra")
library("cluster")
library("dendextend")


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
# dft<-PhyclusVar(df_Peaks, varCat1="Taxonomie", value="Pediococcus pentosaceus", varCat2="Nutrition", dendrogram="cladogram")
# dft<-PhyclusVar(df_Peaks, varCat1="Taxonomie", value="All", varCat2="Date.de.récolte", dendrogram="fan")
# dft<-PhyclusVar(df_Peaks, varCat1="Taxonomie", value="Lactobacillus plantarum", varCat2="Ruche", dendrogram="fan")
# dft<-PhyclusVar(df_Peaks, varCat1="Taxonomie", value="Lactobacillus plantarum", varCat2="Ruche", dendrogram="unrooted")

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
        stop("too small number of individuals...try another value")
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
    plot(as.phylo(hc), type = "radial", tip.color = colors[cl])
  
  coph <- cophenetic(hc)
  ccoph<-cor(dd, coph)
  

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
#methCorr="cophenetic", "baker", "common_nodes", "FM_index"
# Dendrogram_pairComp(df_Peaks, varCat1="Taxonomie", value="Pediococcus pentosaceus")
# Dendrogram_pairComp(df_Peaks, varCat1="Nutrition", value="Marrubium vulgare")
Dendrogram_pairComp<- function(df_m, varCat1, value, methCorr="cophenetic",meth1="ward.D2", meth2="mcquitty",
                    meth3="complete", meth4="centroid", meth5="single", meth6="average", 
                    dist1="euclidean", dist2="euclidean", dist3="euclidean", dist4="euclidean", dist5="euclidean", dist6="euclidean")
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
  
  hcA <- hclust((dist(df_S, method=dist1)), method=meth1)
  hcB <- hclust((dist(df_S, method=dist2)), method=meth2)
  dend1 <- as.dendrogram((hcA))
  dend2 <- as.dendrogram((hcB))
  dends <- dendlist(dend1, dend2)
  tanglegram(dend1, dend2)
  
  tg<-untangle(dends, method = "step1side")
  tanglegram(tg)
  print("Alignment score")
  print(entanglement(tg))
  print("Cophenetic correlation coefficient")
  print(cor_cophenetic(dend1, dend2))
  print("Baker correlation coefficient")
  print(cor_bakers_gamma(dend1, dend2))
  
  
  hcD <- hclust((dist(df_S, method=dist3)), method=meth3)
  hcE<- hclust((dist(df_S, method=dist4)), method=meth4)
  hcF<- hclust((dist(df_S, method=dist5)), method=meth5)
  hcG<- hclust((dist(df_S, method=dist6)), method=meth6)
  dend3<-as.dendrogram(hcD)
  dend4<-as.dendrogram(hcE)
  dend5<-as.dendrogram(hcF)
  dend6<-as.dendrogram(hcG)
  dendMultiples <- dendlist(n1= dend1, n2= dend2, n3 = dend3, n4 = dend4, n5= dend5,n6 = dend6 )
  print(paste("n1:", meth1, dist1))
  print(paste("n2:", meth2, dist2))
  print(paste("n3:", meth3, dist3))
  print(paste("n4:", meth4, dist4))
  print(paste("n5:", meth5, dist5))
  print(paste("n6:", meth6, dist6))
  print("Tree correlation matrix")
  cors <- cor.dendlist( dendMultiples, method=methCorr)
  print(round(cors, 2))
  
}




