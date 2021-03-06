% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.R
\name{getCellCycleGenes}
\alias{getCellCycleGenes}
\title{If files are passed both s and g2m phase lists must be supplied.
 - text files with a single headerless column of symbols are
   expected.
If files are not passed, the lists from Seurat will be used.
If species is "mm" and lists from Seurat are being used, the
symbols will be converted with bioMart.}
\usage{
getCellCycleGenes(
  sgenes_file = NULL,
  g2mgenes_file = NULL,
  species = c("mm", "hs")
)
}
\arguments{
\item{sgenes_file}{A text file with a single column of gene symbols}

\item{g2mgenes_file}{A text file with a single column of gene symbols}

\item{species}{Either "hs" or "mm" - ignored if files are specified}
}
\description{
If files are passed both s and g2m phase lists must be supplied.
 - text files with a single headerless column of symbols are
   expected.
If files are not passed, the lists from Seurat will be used.
If species is "mm" and lists from Seurat are being used, the
symbols will be converted with bioMart.
}
