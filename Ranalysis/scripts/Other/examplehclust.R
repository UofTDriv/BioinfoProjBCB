library(dplyr)
library(ggplot2)
library(openxlsx)
library(tidyr)

# Updated Gene list in the EXACT desired order
# Note: REMOVED Normalized "HLA-C" and "HLA-E" to standard hyphens for matching

target_genes <- c( "SRGN", "PLAUR", "SOD2", "BCL2A1", "CXCL8", "RGS2", "CCL4L2", "PDZK1IP1", "C15orf48", "WFDC21P", "DUOXA2", "CDKN1A", "ARF6", "HES4", "S100P", "HSPB1", "SDCBP2", "SAT1", "MYL12A", "ACTG1", "ACTB", "GLUL", "FTH1", "NFKBIA", "IER2", "GADD45B", "PPP1R15A", "DUSP1", "PMAIP1", "ARL4C", "KRT17", "RND3", "ERO1A", "DUSP5", "KRT4", "CSTA", "METRNL", "C1orf116", "PITX1", "TP53INP2", "IER5", "TPM4", "PPIF", "MXD1", "BIRC3", "TNFAIP3", "CCNL1", "ACSL1", "FOS", "CEACAM5", "TMPRSS11D", "CCN1", "APOBEC3A", "OPTN", "DNAJA1", "IFI44", "MCL1", "MARCKS", "MT2A", "HLA-C", "B2M", "TXNIP", "RPL36AL", "CFL1", "TMSB10", "TXN", "TMSB4X", "FTL", "NOP10", "ARPC5", "IFITM3", "IFITM1", "UBE2L6", "BST2", "SOX4", "PLSCR1", "KLF6", "SELENOK", "DRAP1", "SYF2", "GRB2", "PSME2", "LAP3", "CD74", "C6orf62", "SAMD9L", "OAS1", "PARP14", "MX1", "MX2", "TNFSF10", "IFI6", "IFI35", "EPSTI1", "STAT1", "IFIT3", "IFIT1", "ISG20", "IFIT2", "OASL", "GBP1", "RSAD2", "ISG15", "IRF7", "SAMD9", "YPEL5", "CXCL1", "GBP5", "CARD16", "CASP1", "NMI", "ARRDC3", "ACTR3", "SAA1", "IDO1", "UPK1B", "CXCL11", "CXCL10", "CCL8", "DCP1A", "SPRR1A", "CNFN", "RHCG", "MAL", "KRT13", "EMP1", "KRT16", "RNASE7", "IL23A", "SPINK7" ) 

sigsub <- read.csv("outputs/DESeq2_results/Severe_acute_respiratory_syndrome_coronavirus_2_DESeq2_results.csv", header = TRUE)
sigsub <- rename(sigsub, Geneid = gene)
counttable <- read.csv("data/processed/filtered_gene_counts.csv", row.names = 1, header = TRUE)
counttable <- counttable[, !colnames(counttable) %in% c("Geneid")]

###### Z-score calculation #######
sig_counts_log10 <- subset(counttable, rownames(counttable) %in% sigsub$Geneid)
sig_counts_log10 <- log10(sig_counts_log10 + 1)
sig_counts_log10_m <- as.matrix(sig_counts_log10)

sdlist <- apply(sig_counts_log10_m, 1, sd)
meanlist <- rowMeans(sig_counts_log10_m)

sig_counts_log10$sd <- sdlist
sig_counts_log10$mean <- meanlist
nsamples <- ncol(sig_counts_log10) - 2

zscores_log10 <- (sig_counts_log10[,1:nsamples] - sig_counts_log10$mean) / sig_counts_log10$sd
zscores_log10 <- as.matrix(zscores_log10)

# (Clustering still runs for data prep, but we will override the gene order below)
zscores_log10_df <- as.data.frame(zscores_log10)
zscores_log10_df$gene <- rownames(zscores_log10_df)
zscores_log10_df <- zscores_log10_df %>% tidyr::pivot_longer(cols = -gene, names_to = "sample", values_to = "zscore")

# FILTER AND ORDER BY YOUR LIST
# We use indexing with target_genes to ensure the order is preserved during the intersect
top_genes <- target_genes[target_genes %in% unique(zscores_log10_df$gene)]

if (length(top_genes) == 0) {
  stop("None of the target genes were found in the dataset.")
}

zscores_plot_df <- zscores_log10_df %>% filter(gene %in% top_genes)

# Use rev() so the first gene in your list appears at the TOP of the plot
zscores_plot_df$gene <- factor(zscores_plot_df$gene, levels = rev(top_genes))

# Sample ordering (User defined)
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

present_samples <- unique(zscores_plot_df$sample)
sample_levels <- sample_levels[sample_levels %in% present_samples]
zscores_plot_df$sample <- factor(zscores_plot_df$sample, levels = sample_levels)

# Plotting
zscores_plot_df$mod_zscore <- pmin(pmax(zscores_plot_df$zscore, -3), 3)

pdf("outputs/plots/ordered_heatmap_DEGs.pdf", width = 17, height = 17)
ggplot(zscores_plot_df, aes(x = sample, y = gene, fill = mod_zscore)) +
  geom_tile() +
  scale_fill_gradientn(colors = c('darkblue', 'white', 'red4'), limits = c(-3, 3),
                       guide = guide_colorbar(title = "Zscore of log10(Counts+1)",
                                              title.position = "top",
                                              direction = "horizontal",
                                              barwidth = unit(10, "cm"),
                                              barheight = unit(0.5, "cm"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
        axis.text.y = element_text(size = 8),
        legend.position = "bottom")
dev.off()