library(factoextra)
library(corrplot)

# Visual_CorrDistM (df_Peaks, varCat1="Taxonomie", value="Pediococcus pentosaceus")
Visual_CorrDistM <- function(df_m,varCat1, value, CorrDist="Corr", dist="euclidean", Ls=6, CorrFig="color", layout="upper", ord="FPC", pv=TRUE, sig=0.05)
  {
  # dist: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" or "kendall".
  # CorrDist:"Corr" or "Dist"
  # CorrFig "circle" , "square", "ellipse", "number", "pie", "shade" and "color"
  # layout:"full" "upper" "lower"
  # ord: "AOE" "FPC" "hclust" 
  # https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
  # https://www.rdocumentation.org/packages/corrplot/versions/0.84/topics/corrplot
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
  
  if (CorrDist=="Corr")
       { MC<-cor(df_S) 
       print(round(MC, 2))
          if (pv) {
             r1 <- cor.mtest(df_S, conf.level = .95)
             corrplot(MC, method = CorrFig, type=layout, order=ord, tl.cex=0.4,  p.mat = r1$p, sig.level = sig, insig = "blank")
          }
          else
            corrplot(MC, method = CorrFig, type=layout, order=ord, tl.cex=0.5)
       write.csv(as.data.frame(MC), "DMatrix.csv")
       } 
 else {
     df_S <- scale(df_S)
     dM <- get_dist(df_S, method = dist)
     fviz_dist(dM, lab_size = Ls) }

}

