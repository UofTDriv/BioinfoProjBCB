
#example script for zscore and hclust 
#may 13 2025

library(dplyr)
library(ggplot2)
library(openxlsx)
library(tidyr)

sigsub <- read.csv("outputs/DESeq2_results/Severe_acute_respiratory_syndrome_coronavirus_2_DESeq2_results.csv", header = TRUE)
# rename gene to Geneid
sigsub <- rename(sigsub, Geneid = gene)
counttable <- read.csv("data/processed/filtered_gene_counts.csv", row.names = 1, header = TRUE)
# remove geneid col
counttable <- counttable[, !colnames(counttable) %in% c("Geneid")]
head(counttable)

######do zscore of log10counts+1#######
sig_counts_log10<-subset(counttable, rownames(counttable) %in% sigsub$Geneid)
head(sig_counts_log10)
#do log10 of counts+1
sig_counts_log10<-log10(sig_counts_log10+1)
sig_counts_log10_m<-as.matrix(sig_counts_log10)

sdlist <- apply(sig_counts_log10_m, 1, sd)
meanlist <- rowMeans(sig_counts_log10_m)

#get sd and mean 
sig_counts_log10$sd<-sdlist
sig_counts_log10$mean<-meanlist

# find number of samples
nsamples <- ncol(sig_counts_log10) - 2

# calculate z-scores
zscores_log10 <- (sig_counts_log10[,1:nsamples] - sig_counts_log10$mean) / sig_counts_log10$sd

##then quick plots showing neg vs pos
#merge the zscores and kept id list and keep only those that merge 
zscores_log10<-as.matrix(zscores_log10)

#cluster the x and y axes using hierarchical clustering to prepare the ordering and then force the order 
d_zscores_log10 <- dist(zscores_log10, method = "euclidean")
hc_genes_log10 <- hclust(d_zscores_log10, method = "average"  )
hc_genes_log10$labels
hc_genes_log10$order
geneorder_log10<-hc_genes_log10$labels[c(hc_genes_log10$order)]

zscores_log10_t<-t(zscores_log10)
d_zscores_log10_1 <- dist(zscores_log10_t, method = "euclidean")
hc_genes_log10_1 <- hclust(d_zscores_log10_1, method = "average"  )
hc_genes_log10_1$labels
hc_genes_log10_1$order
sampleorder_log10<-hc_genes_log10_1$labels[c(hc_genes_log10_1$order)]

###plot hm
zscores_log10_df<-as.data.frame(zscores_log10)
zscores_log10_df$gene<-rownames(zscores_log10_df)

zscores_log10_df <- zscores_log10_df %>% tidyr::pivot_longer(cols = -gene, names_to = "sample", values_to = "zscore")

colnames(zscores_log10_df)<-c("gene", "sample", "zscore")

zscores_log10_df$gene<-factor(zscores_log10_df$gene, levels = geneorder_log10 )
zscores_log10_df$sample<-factor(zscores_log10_df$sample, levels = sampleorder_log10 )

summary(zscores_log10_df$zscore)

#if there are big values they can be set down to a max value 
zscores_log10_df$mod_zscore<-zscores_log10_df$zscore
zscores_log10_df[zscores_log10_df$zscore > 3,]$mod_zscore <-3

# select top 100 most significant genes (padj preferred), keep only those present in the zscore df
sig_col <- if ("padj" %in% colnames(sigsub)) "padj" else if ("pvalue" %in% colnames(sigsub)) "pvalue" else stop("No p-value column (padj or pvalue) in sigsub")
top_genes <- sigsub %>% filter(!is.na(.data[[sig_col]])) %>% arrange(.data[[sig_col]]) %>% slice_head(n = 100) %>% pull(Geneid)
top_genes <- intersect(top_genes, unique(zscores_log10_df$gene))
if (length(top_genes) == 0) top_genes <- unique(zscores_log10_df$gene)

# subset long df to top genes and preserve ordering from clustering if possible
zscores_plot_df <- zscores_log10_df %>% filter(gene %in% top_genes)
plot_gene_levels <- geneorder_log10[geneorder_log10 %in% top_genes]
if (length(plot_gene_levels) != length(unique(zscores_plot_df$gene))) {
  plot_gene_levels <- sort(unique(zscores_plot_df$gene))
}
zscores_plot_df$gene <- factor(zscores_plot_df$gene, levels = plot_gene_levels)

# fixed sample order provided by user (will keep only samples present in the data; append any extras)
sample_levels <- c(
  "X30306045","W90303867","W90303891","W81803328","W72103811","W81802779","W81703168","W61006922",
  "W71701367","X30101334","X22302231","W71607287","X30504387","W71700992","W71700985","W81802080",
  "W71701045","W81704242","W61401321","W61402112","W90303291","W90406146","W90504037","W71607804",
  "X22404828","X01703004","W90502794","W61402109","W81703961","W73109775","W90404502","Neg_Control_W",
  "W90504181","W71606182","W81703049","W61303224","W73007934","W71607867","W81801994","X30306010",
  "X30501984","W81803327","W72803457","W60805434","W72206124","W71606944","W61303468","W61402317",
  "W60805435","W61402333","W61303382","W61402326","W61303414","W73103274","W61402335","W61402295",
  "W61402322","W61402303","W73007312","W61303509","W72406818","W61401221","W71700722","X30203412",
  "W71702778","W61303457","W71701863","W71607726","W90202950","W90504879","W72900530","W72504187",
  "W72407764","W73008317","W73103566","W61400968","W71607235","W61501463","W73007677","X11007358",
  "W61501903","W71607203","W61006900","X22106916","X30108590","W61501250","W73103801","W71607585",
  "W73008045","W61304557","W71605335","W71607748","W72902029","X30108596","W71701812","X30203379",
  "W50504632","W73006205","W72102578","W72102619","W71608660","W61303057","W71607236","W73008772",
  "X30108583","W61304603","W61401551","W73007482","W90407426","W71703213","W71607201","W61304616",
  "W61304621","W61402710","W61303838","W71607588","W71703204","W61303884","W71701651","W72003672",
  "W72803699","W60804985","W61402281"
)

# keep only samples present and preserve given order; append any other samples not listed
present_samples <- unique(zscores_plot_df$sample)
sample_levels <- sample_levels[sample_levels %in% present_samples]
extra_samples <- setdiff(present_samples, sample_levels)
if (length(extra_samples) > 0) sample_levels <- c(sample_levels, sort(extra_samples))

zscores_plot_df$sample <- factor(zscores_plot_df$sample, levels = sample_levels)
sampleorder_log10 <- sample_levels

#create heatmap and save to pdf (top 100)
pdf("outputs/plots/initial_heatmap_DEGs_human_log10plus1_top100.pdf", width = 17, height = 17)
ggplot(zscores_plot_df, aes(x = sample, y = gene, fill = mod_zscore)) + geom_tile() + 
  scale_fill_gradientn(colors=c(low = 'darkblue', mid = 'white', high = 'red4'), 
                        limits=c(-3, 3), breaks=c(-3, -1.5, 0, 1.5, 3) ) +
  xlab("Sample") + ylab("Gene") + labs(fill="Zscore of log10(Counts+1)") + theme_bw()  +
  theme(axis.title = element_text(size=16), legend.position = "bottom", plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, size=1), 
      axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1), axis.text.y = element_text(hjust = 1, size = 12)  ) 
dev.off()