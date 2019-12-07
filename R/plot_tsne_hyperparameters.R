## Title ----
##
## Visualisation of effect of tSNE hyperparameters
##
## Description ----
##
## Make plots e.g. of tSNE runs with different perplexity values.
##
## Details ----
##
##
## Usage ----
##
## See options.

stopifnot(
  require(ggplot2),
  require(reshape2),
  require(Matrix),
  require(optparse),
  require(tenxutils)
)

# Options ----

option_list <- list(
    make_option(
      c("--table"),
      default="none",
      help="A table containing the tsne coordinates and phenotype information"
      ),
    make_option(
      c("--shapefactor"),
      default="none",
      help="A column in the cell metadata to use for deterimining the shape of points on the tSNE plots"
      ),
    make_option(
      c("--colorfactor"),
      default="none",
      help="Column in the cell metadata to use for deterimining the color of points on the tSNE. One plot will be made per color factor."
      ),
    make_option(
      c("--hyperparameter"),
      default="perplexity",
      help="The column containing the hyperparameter to facet on"
    ),
    make_option(
      c("--plotdirvar"),
      default="tsneDir",
      help="latex var for plot location"
      ),

    make_option(
      c("--pointsize"),
      default=0.5,
      help="The point size for the tSNE plots"),
    make_option(
      c("--pointalpha"),
      default=0.8,
      help="The alpha setting for the points on the tSNE plots"
      ),
    make_option(
      c("--outdir"),
      default="seurat.out.dir",
      help="outdir")
    )

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)


##
message("Reading in rdims table")
plot_data <- read.table(opt$table, sep="\t", header=TRUE)

tex = ""

print("Making tSNE plots for the grouping variables")

color_var <- opt$colorfactor
plot_data[,color_var] <- factor(plot_data[,color_var], 
				levels=sort(as.numeric(unique(plot_data[,color_var]))))


print(paste("Making",color_var,"tSNE plot"))

## convert the hyperparameter into a factor with correct
## level ordering

param_values <- as.numeric(plot_data[[opt$hyperparameter]])
param_levels <- sort(unique(param_values))
param_factor <- factor(as.character(param_values),levels=as.character(param_levels))

plot_data[[opt$hyperparameter]]  <- param_factor

if(opt$shapefactor=="none" | !(opt$shapefactor %in% colnames(plot_data)))
    {
        gp <- ggplot(plot_data, aes_string("tSNE_1", "tSNE_2", color=color_var))
    } else {
        gp <- ggplot(plot_data, aes_string("tSNE_1", "tSNE_2",
                                           color=color_var, shape=opt$shapefactor))
    }

gp <- gp + geom_point(size=opt$pointsize)
gp <- gp + facet_wrap(as.formula(paste("~", opt$hyperparameter)),
                      scales="free",
                      ncol=3)

# calculate the figure height.
h = ceiling(length(param_levels)/3) * 2.5

plotfilename = paste("tSNE", opt$hyperparameter,sep=".")

save_ggplots(file.path(opt$outdir, plotfilename),
             gp,
             width=8,
             height=h)


texCaption <- paste("Effect of",opt$hyperparameter,"on the tSNE projection (plots are colored by cluster)")

tex <- paste(tex,
             getFigureTex(plotfilename,texCaption,
             plot_dir_var=opt$plotdirvar),
             sep="\n")

print("Writing out latex snippet")
## write out latex snippet
## (for \input{} in report)
tex_file <- file.path(opt$outdir,
                      paste0("tSNE.", opt$hyperparameter,".tex"))

writeTex(tex_file, tex)
