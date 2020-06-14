
BiocManager::install(version = "3.11")
BiocManager::install("DESeq2")
BiocManager::install("plyr")
BiocManager::install("dplyr")
BiocManager::install("splitstackshape")
BiocManager::install("ggplot2")
BiocManager::install("RColorBrewer")




library("DESeq2")


directory <- getwd()
sampleFiles <- grep ("counts" ,list.files(directory),value=TRUE)
sampleCondition <- sub("(.*counts).*","\\1" ,sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~condition)
dds <- DESeq(ddsHTSeq)
res <- results(dds)


#pca_plot
rld = rlogTransformation(dds)

svg("./figures/pcaplot.svg")

plotPCA(rld, intgroup=c("condition"))

dev.off()


res$id <- rownames(res)


rownames(res) <- NULL


write.csv(res, file = "./rnaseq_analysis/raw_deseqtable.csv")

deseqtable <- read.csv("./rnaseq_analysis/raw_deseqtable.csv")

deseqtable$fold_change <- 2^deseqtable[,3]

library(splitstackshape)

tab <- read.table("hpara.gff", sep="\t", header=FALSE, comment.char="#",
                  na.strings=".", stringsAsFactors=FALSE,
                  quote="", fill=FALSE)

results2 <- concat.split(data = tab, split.col = "V9", sep = "=", drop = TRUE)
gff_file <- concat.split(data = results2, split.col = "V9_2", sep = ";", drop = TRUE)

names(gff_file)[12] <- "id"
names(gff_file)[10] <- "gff_annotation"
filt_gff <- gff_file[,c(10, 12)]

library("dplyr")

merged_gff_deseq<- inner_join(deseqtable, filt_gff, by="id")


filtered_merged_gff_deseq <- merged_gff_deseq[,c(9, 8, 6, 7, 3, 10)]


write.csv(filtered_merged_gff_deseq, file = "./rnaseq_analysis/raw_condition2_vs_condition1.csv")


sig <- filter(filtered_merged_gff_deseq, padj<0.05)
sigup <- filter(sig, fold_change>2)
sigdown <- filter(sig, fold_change<0.5)
colnames(sigdown)
sigdown$fold_change <- -1/sigdown[,1]


sorted_sigup<- arrange(sigup, desc(fold_change))
sorted_sigdown<- arrange(sigdown, fold_change)

write.csv(sig, file = "./rnaseq_analysis/allsig_conditon2_vs_condition1.csv")
write.csv(sorted_sigup, file = "./rnaseq_analysis/upsig_conditon2_vs_condition1.csv")
write.csv(sorted_sigdown, file = "./rnaseq_analysis/downsig_conditon2_vs_condition1.csv")


ko <- dir(pattern="\\_ko.txt$")

ko_numbers <- read.table(ko, header = FALSE, fill = TRUE)


names(ko_numbers)[1] <- "id"
names(ko_numbers)[2] <- "ko"

kegg <- dir(pattern="\\_kegg.csv$")

kegg_data<- read.csv(kegg, header = FALSE, fill = TRUE)

names(kegg_data)[1] <- "ko"
names(kegg_data)[3] <- "kegg_annotation"
names(kegg_data)[4] <- "lowest_category"
names(kegg_data)[5] <- "highest_category"

merged_kegg_ko <- inner_join(kegg_data, ko_numbers, by="ko")

merged_final_kegg <- inner_join(merged_kegg_ko, filtered_merged_gff_deseq, by="id")

filtered_finalkegg <- merged_final_kegg[,c(8, 7, 11, 9, 10, 1, 12, 2, 3, 4, 5)]


write.csv(filtered_finalkegg, file = "./kegg_analysis/raw_kegg_table.csv")

sigkegg <- filter(filtered_finalkegg, padj<0.05)
sigupkegg <- filter(sigkegg, fold_change>2)
sigdownkegg <- filter(sigkegg, fold_change<0.5)

#calculate foldchange and add column
sigdownkegg$fold_change <- -1/sigdownkegg[,1]


sorted_sigupkegg<- arrange(sigupkegg, desc(fold_change))
sorted_sigdownkegg<- arrange(sigdownkegg, fold_change)

write.csv(sigkegg, file = "./kegg_analysis/allsig_condition2_vs_condition1_kegg.csv")
write.csv(sorted_sigupkegg, file = "./kegg_analysis/upsig_condition2_vs_condition1_kegg.csv")
write.csv(sorted_sigdownkegg, file = "./kegg_analysis/downsig_condition2_vs_condition1_kegg.csv")


library(ggplot2)

library(RColorBrewer)

#up_lowest_cat

piechartdata_up <- as.data.frame(sigupkegg)

colourCount_up = length(unique(piechartdata_up$lowest_category))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

p<- ggplot(piechartdata_up, aes(x="", fill=lowest_category)) +
  geom_bar(width = 1,) + coord_polar("y") +
  scale_fill_manual(values = getPalette(colourCount_up))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


svg("./figures/overexpressed_lowest_category.svg", width = 20, height = 10)

p + blank_theme +
  theme(axis.text.x=element_blank()) 

dev.off()


#up_highest_cat


colourCount_up_high = length(unique(piechartdata_up$highest_category))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

q<- ggplot(piechartdata_up, aes(x="", fill=highest_category)) +
  geom_bar(width = 1,) + coord_polar("y") +
  scale_fill_manual(values = getPalette(colourCount_up_high))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


svg("./figures/overexpressed_highest_category.svg", width = 20, height = 10)

q + blank_theme +
  theme(axis.text.x=element_blank()) 

dev.off()

#down_lowest_cat


piechartdata_down <- as.data.frame(sigdownkegg)




colourCount_down = length(unique(piechartdata_down$lowest_category))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

r<- ggplot(piechartdata_down, aes(x="", fill=lowest_category)) +
  geom_bar(width = 1,) + coord_polar("y") +
  scale_fill_manual(values = getPalette(colourCount_down))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


svg("./figures/repressed_lowest_category.svg", width = 20, height = 10)

r + blank_theme +
  theme(axis.text.x=element_blank()) 

dev.off()

#down_highest_cat

colourCount_down_high = length(unique(piechartdata_down$highest_category))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

s<- ggplot(piechartdata_down, aes(x="", fill=highest_category)) +
  geom_bar(width = 1,) + coord_polar("y") +
  scale_fill_manual(values = getPalette(colourCount_down_high))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


svg("./figures/repressed_highest_category.svg", width = 20, height = 10)

s + blank_theme +
  theme(axis.text.x=element_blank()) 

dev.off()

#values in legend
detach("package:dplyr", unload=TRUE)

library("plyr")

summary_up_low<- count(sigupkegg, 'lowest_category')
sorted_sum_up_low <- arrange(summary_up_low, desc(freq))

summary_up_high<- count(sigupkegg, 'highest_category')
sorted_sum_up_high <- arrange(summary_up_high, desc(freq))

summary_down_low<- count(sigdownkegg, 'lowest_category')
sorted_sum_down_low <- arrange(summary_down_low, desc(freq))

summary_down_high<- count(sigdownkegg, 'highest_category')
sorted_sum_down_high <- arrange(summary_down_high, desc(freq))

write.csv(sorted_sum_up_low, file = "./kegg_analysis/overexpressed_kegg_summary_catlow.csv")
write.csv(sorted_sum_up_high, file = "./kegg_analysis/overexpressed_kegg_summary_cathigh.csv")
write.csv(sorted_sum_down_low, file = "./kegg_analysis/repressed_kegg_summary_catlow.csv")
write.csv(sorted_sum_down_high, file = "./kegg_analysis/repressed_kegg_summary_cathigh.csv")





