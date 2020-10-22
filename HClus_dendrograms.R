library("pvclust")
library("ape")
library("factoextra")
library("cluster")


# default values dft<-BHclus(df)
# dft<-BHclus(df, vcol="all")
BHclus <- function(df_m, meth="ward.D2", d="euclidean", nb=100, vcol="genre", pgen="Lactobacillus", 
                   pspe="Pediococcus pentosaceus", fl="Erica cinerea", ruch="VI",  fig=TRUE){
 
   if (vcol=="genre") {
     df_m<-df_m[df_m$Genre==pgen, ]
     namefile<-paste(pgen, "txt", sep = ".")
   }
  else if (vcol=="species") {
    df_m<-df_m[df_m$Taxonomie==pspe, ] 
    namefile<-paste(pspe, "txt", sep = ".")
    
  }
  else if (vcol=="flower") {
    df_m<-df_m[df_m$Nutrition==fl, ]
    namefile<-paste(fl, "txt", sep = ".")
  }
  else if (vcol=="Ruche"){
    df_m<-df_m[df_m$Ruche==ruch, ] 
    namefile<-paste(ruch, "txt", sep = ".")
  }

set.seed(873) 
pv <- pvclust(t(df_m[,-1:-9 ]),
              method.hclust=meth,
              method.dist=d, nboot=nb)
if (fig) {
  plot(pv, hang = -1, cex = 0.5)
  pvrect(pv)
}
clus<- pvpick(pv)
clus<-as.matrix(clus)
write.table(clus, namefile)
return (clus)
}


Phyclus <- function(df_m, meth="ward.D2", d="euclidean", varCat1="Genre", value="Lactobacillus", 
                     dendrogram="phylogram", nc=4){
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
  dd <- dist(scale(df_S), method = d)
  hc <- hclust(dd, method = meth)
  cl <- cutree(hc, nc)
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
  return(ccoph)
}


PhyclusVar <- function(df_m, meth="ward.D2", d="euclidean", varCat1="Genre", value="Lactobacillus", VarCat2="Nutrition",  dendrogram="phylogram"){
                    
  # varCat1, varCat2: Taxonomie	Genre	Date.d.analyse	Origine	Ruche	Nutrition	Date.de.récolte	Lieu.de.la.ruche
  # value:"Pediococcus pentosaceus",.. Nutrition: "Erica cinerea", Ruche:"VI", 
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


Vizclus <- function(df_m, meth="ward.D2", d="euclidean", vcol="genre", pgen="Lactobacillus", 
                    pspe="Pediococcus pentosaceus", fl="Erica cinerea", ruch="VI", dendrogram="phylogram", nc=4){

if (vcol=="genre") {
  df_m<-df_m[df_m$Genre==pgen, ]
  namefile<-paste(pgen, "txt", sep = ".")
}
else if (vcol=="species") {
  df_m<-df_m[df_m$Taxonomie==pspe, ] 
  namefile<-paste(pspe, "txt", sep = ".")
}
else if (vcol=="flower") {
  df_m<-df_m[df_m$Nutrition==fl, ]
  namefile<-paste(fl, "txt", sep = ".")
}
else if (vcol=="Ruche"){
  df_m<-df_m[df_m$Ruche==ruch, ] 
  namefile<-paste(ruch, "txt", sep = ".")
}

df_A<-df_m[,1:9 ]
df_S<-df_m[,-1:-9 ]
for (j in 1:ncol(df_S)) { df_S[,j]<-as.numeric(df_S[,j])  }
for (j in 1:ncol(df_A)) { df_A[,j]<-as.numeric( factor (df_A[,j]) )  }

dd <- dist(scale(df_S), method = d)
hc <- hclust(dd, method = meth)
cl <- cutree(hc, k = nc)

tclust<-table(cl)

coph <- cophenetic(hc)
ccoph<-cor(dd, coph)

indC<-1:nc
colors = c("red", "blue", "green", "black", "orange", "yellow", "purple", "pink", "brown", "violet","cyan", "lightsalmon", "darkgreen","firebrick","gray","aquamarine", "steelblue","wheat4","olivedrab" )

fviz_dend(hc, k = nc, 
          cex = 0.6, 
          k_colors = colors[indc],
          color_labels_by_k = TRUE, 
          rect = TRUE, 
          horiz = TRUE
     )

fviz_cluster(list(data = df_S, cluster = cl),
             palette = "jco", 
             ellipse.type = "convex", 
             repel = TRUE, 
             show.clust.cent = FALSE, labelsize = 6,
             ggtheme = theme_minimal()
)

listresults <- list()
for(i in 1:(nc+2)){
  if (i<=nc)  listresults[[i]] <- rownames(df_A)[cl == i]
  else if (i==(nc+1)) listresults[[i]] <-tclust
  else if (i==(nc+2)) listresults[[i]] <-ccoph
}
listresults<-as.matrix(listresults)
write.table(listresults, namefile)
return(listresults)
}







