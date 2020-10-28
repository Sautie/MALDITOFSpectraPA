library("factoextra")
library("corrplot")


#' @title Visual_CorrDistM
#' @description Visual_CorrDistM is based on functions for computing and visualizing distance matrix and correlograms
#' @param       df_m: dataframe containing peaks and metadata
#' @param       dist:  distances "euclidean", (default value),  "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" or "kendall".
#' @param       varCat1: categorical variable for choosing isolates, examples: "Taxonomie"	,"Genre",	"Date.d.analyse"	,"Origine","Ruche", "Nutrition"	, "Date.de.récolte"	, "Lieu.de.la.ruche"
#' @param       value: level of catVar1, examples: "Lactobacillus" ("Genre"), Taxonomie("Pediococcus pentosaceus"), "Erica cinerea" ("Nutrition"),...
#' @param       Visual: correlogram "Corr", (default value) or distance matrix figure ("Dist") (this second option is more general because it includes correlation, see dist argument)
#' @param       CorrFig: correlogram type, "circle" , "square", "ellipse", "number", "pie", "shade" and "color" (default value)
#' @param       Ord: correlogram arrangement methods  to find patterns (ord="FPC", default value), "AOE" ,"hclust"
#' @param       layout: correlogram  layout: "full" "upper" "lower", (layout="upper", default value)
#' @param       pv: boolean variable (including or not probability values) (pv=TRUE, default value)
#' @param       sig: significance level (sig="0.05", default value)
#' @return        figures and statistics
#' @examples      Visual_CorrDistM (df_Peaks, varCat1="Taxonomie", value="Pediococcus pentosaceus"),
#'                Visual_CorrDistM (df_Peaks, varCat1="Taxonomie", value="Pediococcus pentosaceus", VisualM="dist")
#'                Visual_CorrDistM (df_Peaks, varCat1="Taxonomie", value="Pediococcus pentosaceus", VisualM="dist", dist="pearson")
#'         source: https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
#'                 https://www.rdocumentation.org/packages/corrplot/versions/0.84/topics/corrplot
#'@export
Visual_CorrDistM <- function(df_m,varCat1, value, VisualM="Corr", dist="euclidean", Ls=6, CorrFig="color", layout="upper", ord="FPC", pv=TRUE, sig=0.05)
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
      df_S <- df_m[chosenrows,]
      df_S <- df_S[,-ni:-nf]
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

  if (VisualM == "Corr")
  {
    MC <- cor(df_S)
    print(round(MC, 2))
    print("running...")
    if (pv) {
      r1 <- cor.mtest(df_S, conf.level = .95)
      corrplot(
        MC,
        method = CorrFig,
        type = layout,
        order = ord,
        tl.cex = 0.4,
        p.mat = r1$p,
        sig.level = sig,
        insig = "blank"
      )
    }
    else
      corrplot(
        MC,
        method = CorrFig,
        type = layout,
        order = ord,
        tl.cex = 0.5
      )
    write.csv(as.data.frame(MC), "CMatrix.csv")
  }
  else {
    df_S <- scale(df_S)
    dM <- get_dist(df_S, method = dist)
    print(fviz_dist(dM, lab_size = Ls))
  }

}

