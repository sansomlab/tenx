## Generic helper functions

#' Tidy up numbers in results tables
#' @param results_table A dataframe or matrix of results contaning numeric data.
#' @param nsignif Number of significant figures (used for numbers > 1)
#' @param nround Number of digits (used for numbers < 1)
tidyNumbers <- function(results_table, nsignif=3, nround=2)
{
  x <- as.data.frame(results_table)
  for(col in colnames(x))
  {
    if(is.numeric(x[[col]]))
    {
      a <- abs(x[[col]])< 1 & !is.na(x[[col]])
      b <- abs(x[[col]])> 1 & !is.na(x[[col]])
      x[[col]][a] <- signif(x[[col]][a],nsignif)
      x[[col]][b] <- round(x[[col]][b],nround)
    }
  }
  x
}

#' format numbers in results tables as character strings
#'
#' @param results_table A dataframe or matrix of results that contains some numeric columns.
#' @param number_fmt Sprintf format for numeric columns
sprintfResults <- function(results_table,
                           number_fmt="%0.3g")
{
    for(col in colnames(results_table))
    {
        if(is.numeric(results_table[[col]]))
        {
            xx <- results_table[[col]]

            nas <- is.na(xx)

            xx[!nas] <- sapply(xx[!nas],
                               sprintf,
                               fmt=number_fmt)

            results_table[[col]] <- xx
        }
    }

    results_table
}


#' Function to add "top" indicator column to degenes
#' @param data The data
#' @param m_col The column containing the log2 ratio
#' @param id_col A column containing unique identifiers
#' @param ngenes The number of genes to demarcate
topGenes <- function(data, m_col = "avg_logFC",
                     use_fc = TRUE,
                     id_col="gene", ngenes=7)
{

  tmp <- data[order(data$p.adj),id_col][1:ngenes]
  if(use_fc)
  {
    # ensure non-redunant
    data <- data[!data[[id_col]] %in% tmp,]
    tmp <- c(tmp, data[rev(order(abs(data[[m_col]]))),id_col][1:ngenes])
  }
  tmp <- unique(tmp)
  tmp
}

#' Function to add "top" and "sig" indicator columns to degenes list
#' @param data The data
#' @param m_col The column containing the log2 ratio
#' @param p_col The column containing the p-value
#' @param p_threshold The p-threshold below which genes are considered significant
#' @param m_col The column containing the log2 ratio
#' @param ngenes The number of genes to demarcate
#' @param id_col A column containing a unique identifier
categoriseGenes <- function(data,m_col="avg_logFC", use_fc=TRUE,
                            p_col="p.adj", p_threshold=0.05,
                            ngenes=7,
                            id_col="gene")
{

  tmp <- topGenes(data[data[[m_col]] > 0,],
                  m_col=m_col, use_fc=use_fc,
                  ngenes=ngenes,id_col=id_col)
  tmp2 <- topGenes(data[data[[m_col]] < 0,],
                   m_col=m_col, use_fc=use_fc,
                   ngenes=ngenes,id_col=id_col)
  data$top <- FALSE
  data$top[data[[id_col]] %in% unique(c(tmp,tmp2))] <- TRUE


  data$sig <- FALSE
  data$sig[data[[p_col]] < p_threshold] <- TRUE
  data
}
