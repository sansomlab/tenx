## Run SingleR in a given Seurat object

# Libraries ----

stopifnot(
  require(optparse),
  require(Seurat),
  require(dplyr),
  require(tenxutils),
  require(SingleR),
  require(scater),
  require(tidyr),
  require(ggforce),
  require(MASS),
  require(viridis)
)

# Options ----

option_list <- list(
  make_option(c("--seuratobject"), default="begin.Robj",
              help="A seurat object after PCA"),
  make_option(c("--reference"),
              help="reference dataset for celltype assignment; see https://bioconductor.org/packages/3.11/bioc/vignettes/SingleR/inst/doc/SingleR.html#5_available_references" ),
  make_option(c("--predictions"), default = NULL,
              help="rds object containig the result of the call to singleR"),
  make_option(c("--show_annotation_in_plots"), default=NULL,
              help="Column names from the metadata slot to show in the annotation of the output plots"),
  make_option(c("--outdir"), default="seurat.out.dir",
              help="outdir")
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

s <- loadSeurat(path=opt$seuratobject)

a <- ifelse("SCT" %in% names(s), yes = "SCT", no = "RNA")

sce <- as.SingleCellExperiment(s, assay = a) # Seems to be using SCT assay (so does it take the defualt?)

# read in the predctions
pred <- readRDS(opt$predictions)

common <- intersect(rownames(sce), rownames(pred))
sce <- sce[common,]


# Plot label assignment scores
cat("Plotting prediction results \n")
annotation_frame <- as.data.frame(colData(sce))
qc_cols <- c("nCount_RNA", "nFeature_RNA", "percent.mito")

if (!is.null(opt$show_annotation_in_plots)) {
    meta_cols <- unlist(strsplit(opt$show_annotation_in_plots, split = ","))

  keep_annotation <- c(qc_cols, meta_cols)
} else {
  keep_annotation <- qc_cols
}

keep_annotation <- keep_annotation[keep_annotation %in% colnames(annotation_frame)]

annotation_frame <- annotation_frame[,keep_annotation]

# Explore predicted scores ----------------
scores <- as.data.frame(pred$scores)
scores$barcode <- rownames(pred)
scores <- pivot_longer(data = scores, cols = -barcode, names_to = "label", values_to = "scores")
pruned_labels <- data.frame(barcode=rownames(pred), pruned=is.na(pred$pruned.labels))
scores <- left_join(scores, pruned_labels, by="barcode")
n_pruned <- scores %>% group_by(barcode) %>% summarise(p=unique(pruned)) %>% pull(p) %>%table
cat("\nNumber of cells with ambiguous labels:", format(n_pruned["TRUE"]), " \n")

# "delta.med", the difference between the score of the assigned label and the median of all scores for each cell.
png(paste0(opt$outdir, "/singler_delta_med.png"), height=480, width=540)
plotScoreDistribution(pred, show = "delta.med", ncol = 4, show.nmads = 3)
dev.off()

# Calculate median values per label
scores.summary <- scores %>% group_by(label) %>%
  summarise(mean=mean(scores)) %>% ungroup() %>%
  arrange(desc(mean))
scores$label <- factor(scores$label, levels=scores.summary$label)

# Plot score distribution per label
p <- ggplot(scores, aes(x=label, y=scores)) +
  geom_sina(alpha=0.25) +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

save_ggplots(paste0(opt$outdir, "/score_distribution_per_label.png"),
             p)

# Explore first and second scores
top_scores <- scores %>% group_by(barcode) %>% mutate(mean_scores=mean(scores)) %>%
  top_n(2, scores) %>% arrange(desc(scores), .by_group = TRUE) %>%
  summarise(first = label[1], second=label[2],
            first_score=scores[1], second_score=scores[2],
            diff_score=scores[1]-scores[2], mean_score=unique(mean_scores), pruned=unique(pruned))

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
top_scores$density <- get_density(top_scores$first_score, top_scores$second_score, n = 100)

p <- ggplot(top_scores, aes(x=first_score, y=second_score, col=density)) +
  geom_point(alpha=0.20, size=2) +
  geom_abline(intercept = 0, slope = 1, col="blue", lty="dashed") +
  geom_vline(xintercept = top_scores[top_scores$density==max(top_scores$density),]$first_score[1],
             lty="dashed", col="grey") +
  geom_hline(yintercept = top_scores[top_scores$density==max(top_scores$density),]$second_score[1],
             lty="dashed", col="grey") +
  ylim(0, max(c(top_scores$first_score, top_scores$second_score))) +
  xlim(0, max(c(top_scores$first_score, top_scores$second_score))) +
  xlab("First score") +
  ylab("Second score") +
    scale_color_viridis()

save_ggplots(paste0(opt$outdir, "/first_vs_second_density_scatter.pdf"),
             p,
             height=5,
             width=5.5)

p <- ggplot(top_scores, aes(x=first_score, y=second_score, col=pruned)) +
  geom_point(alpha=0.20, size=2) +
  geom_abline(intercept = 0, slope = 1, col="blue", lty="dashed") +
  geom_vline(xintercept = top_scores[top_scores$density==max(top_scores$density),]$first_score[1],
             lty="dashed", col="grey") +
  geom_hline(yintercept = top_scores[top_scores$density==max(top_scores$density),]$second_score[1],
             lty="dashed", col="grey") +
  ylim(0, max(c(top_scores$first_score, top_scores$second_score))) +
  xlim(0, max(c(top_scores$first_score, top_scores$second_score))) +
  xlab("First score") +
  ylab("Second score") +
  guides(col = guide_legend(title = "Pruned"))

save_ggplots(paste0(opt$outdir, "/first_vs_second_prunned_scatter.pdf"),
             p,
             width=5,
             height=5.5)


# Density plots
p <- ggplot(top_scores, aes(x= diff_score, col=first)) +
  geom_density() +
  guides(col = guide_legend(ncol = 2, title = "First label"))

save_ggplots(paste0(opt$outdir, "/diff_density.pdf"),
             p,
             height=4.5,
             width=10)


# Heatmap of first and second after prunning
df <- top_scores %>% filter(pruned==FALSE) %>% dplyr::select(first, second) %>% table %>% as.data.frame()
levels <- df %>% arrange(desc(Freq)) %>% pull(first) %>% unique()
df$first <- factor(df$first, levels = rev(levels))
df$Freq[df$Freq ==0 ] <- NA

p <- ggplot(df, aes(x=second, y=first, fill=Freq)) +
  geom_tile() +
  scale_fill_gradient(low="grey90", high="red") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle=45, hjust=0))

save_ggplots(paste0(opt$outdir, "/first_vs_second_heatmap.pdf"),
             p,
             height=5,
             width=6)

# Plot markers used to assign labels
all.markers <- metadata(pred)$de.genes
sce$labels <- pred$labels


do_plot <- function() {
plotHeatmap(sce, order_columns_by=list(I(pred$labels)),
            features=unique(unlist(all.markers[[lab]])))
}

for (lab in unique(pred$labels)) {
    save_plots(paste0(opt$outdir, "/markers_heatmap_", lab,".png"),
               plot_fn=do_plot,
               height = 3000, width=800,
               to_pdf = FALSE)
}


# PCA
pruned_labels <- data.frame(barcode=rownames(pred), pruned_label=as.character(pred$pruned.labels))
top_scores <- left_join(top_scores, pruned_labels, by="barcode")
top_scores$pruned_label <- as.character(top_scores$pruned_label)
top_scores$pruned_label[is.na(top_scores$pruned_label)]="Pruned"
top_scores$density <- NULL
scores <- as.data.frame(pred$scores)
rownames(scores) <- rownames(pred)
pca <- prcomp(scores)
df <- as.data.frame(pca$x)
df$barcode <- rownames(df)
df <- left_join(df, top_scores, by="barcode")
df_bg <- df[,c("PC1", "PC2")]
eigs <- pca$sdev^2

p <-ggplot(df, aes(x=PC1, y=PC2, color=pruned_label)) +
  geom_point(data=df_bg, color="grey", alpha=0.55) +
  geom_point(alpha=0.5) +
  facet_wrap(~pruned_label) +
  guides(col=FALSE) +
  xlab(paste0("PC1 ", round( (eigs[1] / sum(eigs) * 100), digits = 2), "%")) +
    ylab(paste0("PC2 ", round( (eigs[2] / sum(eigs) * 100), digits = 2), "%"))

save_ggplots(paste0(opt$outdir, "/pca_scores.png"),
             p,
             height=6,
             width=6)

cat("Saving output tables \n")
# Output table
pc_out <- df %>% dplyr::select(PC1, PC2, barcode)
colnames(pc_out) <- c("singler_scores_pc1", "singler_scores_pc2", "barcode")
top_scores <- left_join(top_scores, pc_out, by="barcode")
write.table(top_scores, gzfile(paste0(opt$outdir, "/singler_predictions.tsv.gz")),
            sep="\t", quote=FALSE, row.names=FALSE)

# Output summary table
t <- table(top_scores$pruned_label) %>% as.data.frame()
colnames(t) <- c("label", "n")
t <- t[order(t$n, decreasing = TRUE),]
write.table(t,paste0(opt$outdir, "/singler_predictions_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

message("Completed\n\n")
sessionInfo()


# Correlation with QC scores
# meta <- colData(sce) %>% as.data.frame() %>% dplyr::select(nCount_RNA, nFeature_RNA, percent.mito)
# meta$barcode <- rownames(meta)
# top_scores <- left_join(top_scores, meta, by="barcode")

# ggplot(top_scores, aes(x=first_score, y=nFeature_RNA, col=percent.mito)) +
#  geom_point(alpha=0.5) +
#  facet_wrap(~first)
