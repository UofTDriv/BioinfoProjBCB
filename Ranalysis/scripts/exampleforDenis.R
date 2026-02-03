#code excerpt for denis


rnacnt<-read.table("....CountMat.txt", sep = "\t",header = T)
head(rnacnt)
nrow(rnacnt)
#43388


##rrna
rrnalist<-read.csv("hgnc_rRNAgenelist.csv", header = T)
rrnalist
#use Approved.symbol column to match my counts

rrna_marked<-rnacnt
rrna_marked$group<-"othergenes"
rrna_marked[rrna_marked$Geneid %in% rrnalist$Approved.symbol,]$group<-"rrna_genes"
rrna_marked_agg<- rrna_marked %>% 
  group_by(group) %>% 
  summarise(across(2:13, sum))
rrna_marked_agg_t<-as.data.frame(t(rrna_marked_agg))
colnames(rrna_marked_agg_t)<-rrna_marked_agg_t[1,]; rrna_marked_agg_t<-rrna_marked_agg_t[2:nrow(rrna_marked_agg_t),]
rrna_marked_agg_t$othergenes<-as.numeric(rrna_marked_agg_t$othergenes); rrna_marked_agg_t$rrna_genes<-as.numeric(rrna_marked_agg_t$rrna_genes)

rrna_marked_agg_t$sumgenes<-rrna_marked_agg_t$othergenes + rrna_marked_agg_t$rrna_genes
rrna_marked_agg_t$percentrrna<-rrna_marked_agg_t$rrna_genes / rrna_marked_agg_t$sumgenes * 100
write.table(rrna_marked_agg_t, file = "rrna_rates_datasetRNAseq.txt", sep = "\t", quote = F)
#save rrna table


###remove rrna genes from the table and move on
counts_filt<-subset(rnacnt, !(rnacnt$Geneid %in% rrnalist$Approved.symbol))
nrow(counts_filt) #43349

#drop rows which are 0
counts_filt<-counts_filt %>% rowwise() %>% mutate(total = sum(c_across(where(is.numeric))))
head(counts_filt)
counts_filt<-counts_filt[counts_filt$total > 0,]
nrow(counts_filt) #28617
counts_filt$total<-NULL


##alt splice genes
listofASgenes<-counts_filt[grepl("^.*\\-AS[0-9]$", counts_filt$Geneid),]
listofASgenes$Geneid

counts_filt[grepl("^ACTA", counts_filt$Geneid),]

#trim the AS1/AS2/AS3/AS4/AS5
counts_filt$Geneid.trim<-gsub("-AS[0-9]", "", counts_filt$Geneid,  )
#Example - Acta2
counts_filt[counts_filt$Geneid.trim == "ACTA2",]
counts_filt$Geneid<-NULL
head(counts_filt)

rawtablefiltsum<-counts_filt %>%
  group_by(Geneid.trim) %>%
  summarize_all(sum)
head(rawtablefiltsum)
rawtablefiltsum[grepl("^ACTA", rawtablefiltsum$Geneid.trim),]
#check example again

head(rawtablefiltsum)
rawtablefiltsum<-as.data.frame(rawtablefiltsum)
rownames(rawtablefiltsum)<-rawtablefiltsum$Geneid.trim; rawtablefiltsum$Geneid.trim<-NULL
head(rawtablefiltsum)
nrow(rawtablefiltsum)#27753

#filter again on merged samples
#skip now since going into DESEq2
#rawtablefiltsum<-rawtablefiltsum[rowSums(rawtablefiltsum)>5,]

##filter for the GENE GENE_1 GENE_2 situations 
counts_filt2<-rawtablefiltsum
counts_filt2$Geneid.trim<-gsub("_[0-9]", "", rownames(counts_filt2)  )
#counts_filt2$Geneid<-NULL
head(counts_filt2)
counts_filt2$Geneid.trim
genenamecounts<-as.data.frame(table(counts_filt2$Geneid.trim))
genenamecounts<-genenamecounts %>% arrange(desc(Freq))
head(genenamecounts)
counts_filt2[grepl("^HLA-C", counts_filt2$Geneid.trim),]



rawtablefiltsum2<-counts_filt2 %>%
  group_by(Geneid.trim) %>%
  summarize_all(sum)
head(rawtablefiltsum2)
rawtablefiltsum2[grepl("^HLA-C", rawtablefiltsum2$Geneid.trim),]


head(rawtablefiltsum2)
rawtablefiltsum2<-as.data.frame(rawtablefiltsum2)
rownames(rawtablefiltsum2)<-rawtablefiltsum2$Geneid.trim; rawtablefiltsum2$Geneid.trim<-NULL
head(rawtablefiltsum2)
nrow(rawtablefiltsum2)#27505  so removed ~250 genes

#also drop all the LOC105371466 types
listoflocs<-rawtablefiltsum2[grepl("^LOC[0-9]*", rownames(rawtablefiltsum2)),]
#remove
rawtablefiltsum2_filt2<-rawtablefiltsum2[!(rownames(rawtablefiltsum2) %in% rownames(listoflocs)),]
nrow(rawtablefiltsum2_filt2)
#20369


##also the LINC genes
#also drop all the LOC105371466 types
listoflincs<-rawtablefiltsum2_filt2[grepl("^LINC[0-9]*", rownames(rawtablefiltsum2_filt2)),]
#remove
rawtablefiltsum2_filt2<-rawtablefiltsum2_filt2[!(rownames(rawtablefiltsum2_filt2) %in% rownames(listoflincs)),]
nrow(rawtablefiltsum2_filt2) #19327

#write.table(rawtablefiltsum2_filt2, file = "dataset_filtered_counts.txt", sep = "\t", row.names = T, col.names = T, quote = F)
rawtablefiltsum2_filt2<-read.table("dataset_filtered_counts.txt", header = T)

###Do PCA of these samples.
rawtablefiltsum2_filt2_t<-as.data.frame(t(rawtablefiltsum2_filt2))
rownames(rawtablefiltsum2_filt2_t)

#run pca on genes - double check length of genes remaining and use that value here 
pca_res <- prcomp(rawtablefiltsum2_filt2_t[,1:19327], scale. = TRUE)
autoplot(pca_res, data = rawtablefiltsum2_filt2_t)
#plot other pcs

pca1<-autoplot(pca_res, x = 1, y = 2, data = rawtablefiltsum2_filt2_t, colour = 'group1') + ggtitle("PCA of RNAseq Samples") + theme_bw() + coord_equal() +
  xlim(-1, 1) + ylim(-1, 1)

pdf("dataset_pca_v2.pdf", width = 5, height = 5)
pca1
dev.off()

#continue with Deseq2
