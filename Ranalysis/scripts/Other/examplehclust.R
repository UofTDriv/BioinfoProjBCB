
#example script for zscore and hclust 
#may 13 2025

library(dplyr)
library(ggplot2)
library(openxlsx)

sigsub <- read.csv("outputs/DESeq2_results/Severe_acute_respiratory_syndrome-related_coronavirus_DESeq2_results.csv")
head(sigsub)
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
d_zscores_log10_1 <- dist(zscorest, method = "euclidean")
hc_genes_log10_1 <- hclust(d_zscores_log10_1, method = "average"  )
hc_genes_log10_1$labels
hc_genes_log10_1$order
sampleorder_log10<-hc_genes_log10_1$labels[c(hc_genes_log10_1$order)]

###plot hm
zscores_log10_df<-as.data.frame(zscores_log10)
zscores_log10_df$gene<-rownames(zscores_log10_df)

zscores_log10_df<-melt(zscores_log10_df, id.vars = c("gene"))
colnames(zscores_log10_df)<-c("gene", "sample", "zscore")

zscores_log10_df$gene<-factor(zscores_log10_df$gene, levels = geneorder_log10 )
zscores_log10_df$sample<-factor(zscores_log10_df$sample, levels = sampleorder_log10 )

summary(zscores_log10_df$zscore)

#if there are big values they can be set down to a max value 
#zscores_log10_df$mod_zscore<-zscores_log10_df$zscore
#zscores_log10_df[zscores_log10_df$zscore > 3,]$mod_zscore <-3

#create heatmap and save to pdf 
pdf("initial_heatmap_DEGs_human_log10plus1.pdf", width = 17, height = 6)
ggplot(zscores_log10_df, aes(x = sample, y = gene, fill = mod_zscore)) + geom_tile() + 
  scale_fill_gradientn(colors=c(low = 'darkblue', mid = 'white', high = 'red4'), 
                        limits=c(-3, 3), breaks=c(-3, -1.5, 0, 1.5, 3) ) + #may need to adjust scale bounds 
  xlab("Sample") + ylab("Gene") + labs(fill="Zscore of log10(Counts+1)") + theme_bw()  +
  theme(axis.title = element_text(size=16), legend.position = c("bottom"), plot.background = element_blank(),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill=NA, size=1), 
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1), axis.text.y = element_text(hjust = 1, size = 12)  ) 

dev.off()