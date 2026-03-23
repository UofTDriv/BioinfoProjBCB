ml <- read.csv("outputs/FinalRanalysis/outputs/unaligned_merged_long260320op.csv")


colnames(ml)
"outputs/FinalRanalysis/outputs/unaligned_merged_long260320op.csv"
# remove all rows where "bracken_reads"       "cladeReads" are both 0
ml <- ml[!(ml$bracken_reads == 0 & ml$cladeReads == 0), ]

write.csv(ml, "outputs/FinalRanalysis/outputs/unaligned_merged_long260320op_filtered.csv", row.names = FALSE)