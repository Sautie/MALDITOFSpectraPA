library(factoextra)
library(corrplot)


Visual_CorrDistM <- function(df_m, CorrDist="Corr", dist="euclidean", Ls=6, CorrFig="color", layout="upper", ord="FPC", pv=TRUE, sig=0.05){
  # dist: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" or "kendall".
  # CorrDist:"Corr" or "Dist"
  # CorrFig "circle" , "square", "ellipse", "number", "pie", "shade" and "color"
  # layout:"full" "upper" "lower"
  # ord: "AOE" "FPC" "hclust" 
  # https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
  # https://www.rdocumentation.org/packages/corrplot/versions/0.84/topics/corrplot
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
  
  if (CorrDist=="Corr")
       { MC<-cor(df_S) 
          if (pv) {
             r1 <- cor.mtest(df_S, conf.level = .95)
             corrplot(MC, method = CorrFig, type=layout, order=ord, p.mat = r1$p, sig.level = sig, insig = "blank")
          }
          else
            corrplot(MC, method = CorrFig, type=layout, order=ord)  } 
 else {
     df_S <- scale(df_S)
     dM <- get_dist(df_S, method = dist)
     fviz_dist(dM, lab_size = Ls) }

}

