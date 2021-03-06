% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VisualDM_Correlograms.R
\name{Visual_CorrDistM}
\alias{Visual_CorrDistM}
\title{Visual_CorrDistM}
\usage{
Visual_CorrDistM(
  df_m,
  varCat1,
  value,
  VisualM = "Corr",
  dist = "euclidean",
  Ls = 6,
  CorrFig = "color",
  layout = "upper",
  ord = "FPC",
  pv = TRUE,
  sig = 0.05
)
}
\arguments{
\item{df_m:}{dataframe containing peaks and metadata}

\item{dist:}{distances "euclidean", (default value),  "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" or "kendall".}

\item{varCat1:}{categorical variable for choosing isolates, examples: "Taxonomie"    ,"Genre",    "Date.d.analyse"    ,"Origine","Ruche", "Nutrition"    , "Date.de.récolte"    , "Lieu.de.la.ruche"}

\item{value:}{level of catVar1, examples: "Lactobacillus" ("Genre"), Taxonomie("Pediococcus pentosaceus"), "Erica cinerea" ("Nutrition"),...}

\item{Visual:}{correlogram "Corr", (default value) or distance matrix figure ("Dist") (this second option is more general because it includes correlation, see dist argument)}

\item{CorrFig:}{correlogram type, "circle" , "square", "ellipse", "number", "pie", "shade" and "color" (default value)}

\item{Ord:}{correlogram arrangement methods  to find patterns (ord="FPC", default value), "AOE" ,"hclust"}

\item{layout:}{correlogram  layout: "full" "upper" "lower", (layout="upper", default value)}

\item{pv:}{boolean variable (including or not probability values) (pv=TRUE, default value)}

\item{sig:}{significance level (sig="0.05", default value)}
}
\value{
figures and statistics
}
\description{
Visual_CorrDistM is based on functions for computing and visualizing distance matrix and correlograms
}
\examples{
     Visual_CorrDistM (df_Peaks, varCat1="Taxonomie", value="Pediococcus pentosaceus"),
               Visual_CorrDistM (df_Peaks, varCat1="Taxonomie", value="Pediococcus pentosaceus", VisualM="dist")
               Visual_CorrDistM (df_Peaks, varCat1="Taxonomie", value="Pediococcus pentosaceus", VisualM="dist", dist="pearson")
        source: https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html
                https://www.rdocumentation.org/packages/corrplot/versions/0.84/topics/corrplot
}
