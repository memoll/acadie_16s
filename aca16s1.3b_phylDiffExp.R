#######################################################################################
# Differential expression analysis of the phyllosphere samples                        #
# Studying the effects of neonicotinoids on the phyllosphere bacterial communities    #
# Data: Miseq-16S - L'Acadie (ACA)                                                    #
# Mona Parizadeh - 2019-2020                                                          #
#######################################################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq") #‘1.27.6’
library(vegan); packageVersion("vegan") #‘2.5.6’
library(ggplot2); packageVersion("ggplot2")
library(DESeq2); packageVersion("DESeq2") #‘1.29.4’
library(dplyr); packageVersion("dplyr") #‘0.8.5’

# Import data (non-rarefied) #### 
setwd("../mp/aca_16s/files/")
ps = readRDS("16S_aca_phyl.rds")
ps

#Phyloseq to deseq2 conversion ####
phTOds = phyloseq_to_deseq2(ps, design = ~ neonic) 
is(phTOds); isS4(phTOds)
#contents
slotNames(phTOds) 
#estimate size factors 
fcs = estimateSizeFactors(phTOds) #no need to calculate geometric means
#Bayesian estimation of dispersion
dsp = estimateDispersions(fcs)
plotDispEsts(dsp)

#DESeq ####
dds = DESeq(dsp, test = "Wald", fitType="local") #DESeq(phTOds, test = "Wald", fitType="local")

#Investigate test results table ####
resultsNames(dds)
res = results(dds) #extracts a table from a DESeq analysis
#contains base means across samples, log2 fold changes, standard errors, test statistics, p-values and adjusted p-values
#baseMean = the average of the normalized counts taken over all samples
#log2FoldChange = log2 fold change between the groups. E.g. value 2 means that the expression has increased 4-fold
#Fold change is a measure describing how much a quantity changes between control and treatment.
#lfcSE = standard error of the log2FoldChange estimate; stat = Wald statistic; pvalue = Wald test p-value; padj = Benjamini-Hochberg adjusted p-value
#rownames = ASVs
res = res[order(res$padj, na.last=NA), ] #remove padj NAs
mcols(res, use.names=TRUE) #or: colnames(aca.neo.res)
class(res); is(res)
slotNames(res)
summary(res)
hist(res$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 P-value", ylab="Frequency")
#How many adjusted p-values were less than 0.1/0.5/0.01?
sum(res$padj < 0.1, na.rm=TRUE)
sum(res$padj < 0.05, na.rm=TRUE)
sum(res$padj < 0.01, na.rm=TRUE)
#MA plot ####
plotMA(res) 
legend("bottomright", legend="differentially abundant", lty=-1, pch=1, col="red", bty='n')
#red points: adjusted p value less than 0.1 (default threshold)

#Set padj the threshold ####
alpha = 0.05 # Threshold on the adjusted p-value
sigtab = res[(res$padj < alpha), ]
#Combine tax with results ####
sigtab = cbind(as(sigtab, "data.frame"), 
               as(tax_table(ps)[rownames(sigtab), ], "matrix"))
dim(sigtab)

#plot
ggplot(sigtab, aes(x=rownames(sigtab), y=log2FoldChange, color=Phylum)) + 
  geom_hline(yintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size = 5))

#genera in phyla ####
genPhlm = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
genPhlm = sort(genPhlm, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(genPhlm))
# Order genera based on their log2fold
gen = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
gen = sort(gen, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(gen))
#order total number of ASVs associated w/ each genus
sort(table(sigtab$Genus), decreasing = TRUE)

#subset ASVs associated w/ neonic (>0)
sigtab_neo = sigtab[which(sigtab$log2FoldChange>0),];dim(sigtab_neo)
#order based on total number of ASVs
sort(table(sigtab_neo$Genus), decreasing = TRUE)
sort(table(sigtab_neo$Phylum), decreasing = TRUE)

             #subset ASVs associated w/ control (<0)
sigtab_ctl = sigtab[which(sigtab$log2FoldChange<0),];dim(sigtab_ctl)
#order based on total number of ASVs
sort(table(sigtab_ctl$Genus), decreasing = TRUE)
sort(table(sigtab_ctl$Phylum), decreasing = TRUE)

#plot
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  theme_classic() +
  geom_hline(yintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=2) + 
  scale_color_manual(values=c("cornflowerblue","indianred1","mediumvioletred","darkolivegreen4")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, hjust = 0, vjust=0.5),
        axis.title = element_text( size = 12, face = "bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text=element_text(size=12),
        legend.position = "left") +
  annotate("text", x = 12.3, y = 5, label = 'atop(bold("Neonicotinoid-treated"))', parse = TRUE, 
           size = 6, colour = "azure4") +
  annotate("text", x = 13.5, y = -2, label = 'atop(bold("Control"))', parse = TRUE, 
           size = 6, colour = "azure4")
