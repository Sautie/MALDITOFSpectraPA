% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DProcessing.R
\name{SlowXml26Cleaning}
\alias{SlowXml26Cleaning}
\title{SlowXml26Cleaning}
\usage{
SlowXml26Cleaning(Dirp, newdir, op = "delete")
}
\arguments{
\item{Dirp:}{folder containing the Maldi_Tof spectrum files,}

\item{newdir:}{folder where the chosen files are copied,}

\item{op:}{"delete"(default value) or "copy", metadata xlsx file}
}
\description{
copies all mzXML files with 26 lines to a new folder or delete those with less than 26 lines, slow and old version of XmlNCleaning
}
\examples{
 SlowXml26Cleaning("newdata5", "newdata"),  SlowXml26Cleaning("newdata5", "newdata", op="copy")
}
