% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plot.R
\name{violinPlotSection}
\alias{violinPlotSection}
\title{Make a latex section containing a set of violin plots}
\usage{
violinPlotSection(
  data,
  seurat_object,
  cluster_ids,
  type = "positive",
  group.by = opt$testfactor,
  ident.include = opt$identinclude,
  vncol = 4,
  vnrow = 3,
  pt_size = 0.1,
  outdir = opt$outdir,
  analysis_title = "violin plots",
  fc_type = "fold change",
  plot_dir_var = "plotsDir",
  to_pdf = TRUE,
  use.minfc = FALSE
)
}
\arguments{
\item{data}{The data}

\item{seurat_object}{The seurat object}

\item{cluster_ids}{The cluster ids}

\item{type}{Either "positive" or "negative"}

\item{group.by}{A factor to group by}

\item{ident.include}{Identity to include}

\item{outdir}{The directory for the output}

\item{analysis_title}{Title for this section}

\item{fc_type}{The metric by which the violin plots are ordered}

\item{use.minfc}{Use minimum foldchange}

\item{ncol}{Number of columns in the figure}
}
\description{
Make a latex section containing a set of violin plots
}
