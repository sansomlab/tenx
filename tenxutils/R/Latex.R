## Helper functions for writing latex output from Rscripts

## Function definitions
#' Write latex for subsection title
#' @param subsection_name The subsection title.
getSubsectionTex <- function(subsection_name)
{
paste("\\subsection{",gsub("_","\\\\_",subsection_name),"}
      ")
}

#' Write latex for subsubsection title
#' @param subsubsection_name The subsubsection title.
getSubsubsectionTex <- function(subsubsection_name)
{
paste("\\subsubsection{",gsub("_","\\\\_",subsubsection_name),"}
      ")
}

#' Write latex for a figure
#' @param image_file The image file to show in the figure.
#' @param caption The figure caption.
#' @param plot_dir_var The directory containing the image_file
#' @param height The height of the figure as a fraction of the page height.
getFigureTex <- function(image_file, caption, plot_dir_var="plotDir", height=0.9)
{
    caption <- gsub("_","\\\\_", caption)
    ## image_file  <- gsub(".pdf","", pdf_file)

    paste0("\\begin{figure}[H]
            \\includegraphics[width=1.0\\textwidth,height=",height,"\\textheight,keepaspectratio]{{{\\",
           plot_dir_var,
           "/",
           image_file,
           "}}}
               \\caption{",
           caption,
           "}
           \\end{figure}")
}

#' Write latex for a subfigure
#' @param image_file The image file to show in the figure.
#' @param caption The figure caption.
#' @param plot_dir_var The directory containing the image_file
getSubFigureTex <- function(image_file, caption, plot_dir_var="plotDir")
{
    caption <- gsub("_","\\\\_", caption)
    ## pdf_file  <- gsub(".pdf","", pdf_file)

    paste0("\\begin{subfigure}[b]{1.0\\textwidth}
            \\includegraphics[width=1.0\\textwidth]{{{\\",
           plot_dir_var,
           "/",
           image_file,
           "}}}
               \\caption{",
           caption,
           "}
           \\end{subfigure}")
}

#' Function for writing a tex snippet to a file
#' @param file_path The file to write the latex snippet to.
#' @param tex The tex snippet to write (a vector of lines).
writeTex <- function(file_path, tex)
    {
        ## write out latex snippet
        ## (for \input{} in report)
        con <- file(file_path, "w")
        writeLines(tex, con = con)
        close(con)
    }
