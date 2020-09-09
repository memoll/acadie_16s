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

# Import data #### 
setwd("/data/users/mona/miseq_16S/Mona_16S_all/")
ps = readRDS("article1/16S_aca_phyl.rds")
ps

#Subset hosts ####
ps.sy = subset_samples(ps, sample_data(ps)$host == "soy")
ps.sy = prune_taxa(taxa_sums(ps.sy)>0, ps.sy)
ps.sy 
ps.cr = subset_samples(ps, sample_data(ps)$host == "corn")
ps.cr = prune_taxa(taxa_sums(ps.cr)>0, ps.cr)
ps.cr
#Subset treatments ####
ps.ctl = subset_samples(ps, sample_data(ps)$neonic == "N")
ps.ctl = prune_taxa(taxa_sums(ps.ctl)>0, ps.ctl)
ps.ctl 
ps.neo = subset_samples(ps, sample_data(ps)$neonic == "Y")
ps.neo = prune_taxa(taxa_sums(ps.neo)>0, ps.neo)
ps.neo

#Phyloseq to deseq2 conversion ####
phTOds = phyloseq_to_deseq2(ps, design = ~ neonic) #dds file
is(phTOds); isS4(phTOds)
#contents
slotNames(phTOds) 
#estimate size factors 
fcs = estimateSizeFactors(phTOds) #no need to calculate geometric means
#Bayesian estimation of dispersion
dsp = estimateDispersions(fcs)
plotDispEsts(dsp)

#DESeq ####
dds = DESeq(phTOds, test = "Wald", fitType="local")
kable(head(colData(dds))) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
#size factor: median ratio of a sample over a pseudosample
#sort(sizeFactors(dds))

# investigate test results table ####
#results(): extracts a table from a DESeq analysis
#contains base means across samples, log2 fold changes, standard errors, test statistics, p-values and adjusted p-values
#baseMean = the average of the normalized counts taken over all samples
#log2FoldChange = log2 fold change between the groups. E.g. value 2 means that the expression has increased 4-fold
#Fold change is a measure describing how much a quantity changes between control and treatment.
#lfcSE = standard error of the log2FoldChange estimate
#stat = Wald statistic
#pvalue = Wald test p-value
#padj = Benjamini-Hochberg adjusted p-value
#rownames = ASVs
resultsNames(dds)
res = results(dds)
res = res[order(res$padj, na.last=NA), ] #remove padj NAs
kable(head(res)) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
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
kable(sigtab) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
length(which(sigtab$log2FoldChange>0)) #Y
length(which(sigtab$log2FoldChange<0)) #N
sigtab$Phylum
sigtab$Family
sigtab$Genus  
#Replace the long name of a genus
sigtab$Genus = as.character(sigtab$Genus)
sigtab$Genus[sigtab$Genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] = "Rhizobium"

#plot
#fold change is for measuring the changes in gene expression
ggplot(sigtab, aes(x=rownames(sigtab), y=log2FoldChange, color=Phylum)) + 
  geom_hline(yintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size = 5))

#genera in phyla ####
#remove NAs genera  #no, it causes bias
#sigtabGen = subset(sigtab, !is.na(Genus))
#dim(sigtabGen)
#length(which(sigtabGen$log2FoldChange>0)) #Y
#sigtabGen[which(sigtabGen$log2FoldChange>0),]["Phylum"]
#length(which(sigtabGen$log2FoldChange<0)) #N
#sigtabGen[which(sigtabGen$log2FoldChange<0),]["Phylum"]

# Phylum order based on log2Foldchange (it gives only one value to phylum, but there are common phyla in both treatments that we have to count them as well)
genPhlm = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
length(which(genPhlm>0))
which(genPhlm>0)
length(which(genPhlm<0))
which(genPhlm<0)
genPhlm = sort(genPhlm, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(genPhlm))

# Order genera based on their log2fold
gen = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
gen = sort(gen, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(gen))

#order total number of ASVs associated to each genus
sort(table(sigtab$Genus), decreasing = TRUE)
#subset ASVs associated to neonic (>0)
sigtab_neo = sigtab[which(sigtab$log2FoldChange>0),];dim(sigtab_neo)
#order based on total number of ASVs
sort(table(sigtab_neo$Genus), decreasing = TRUE)
sort(table(sigtab_neo$Phylum), decreasing = TRUE)

#subset ASVs associated to control (<0)
sigtab_ctl = sigtab[which(sigtab$log2FoldChange<0),];dim(sigtab_ctl)
#order based on total number of ASVs
sort(table(sigtab_ctl$Genus), decreasing = TRUE)
sort(table(sigtab_ctl$Phylum), decreasing = TRUE)

#check if the most diff abundant ASVs in neo were present in ctl
any(sigtabGen_ctl$Genus == "Hymenobacter")
#mean(sigtabGen[which(sigtabGen$Genus == "Hymenobacter"),]$padj)
any(sigtabGen_ctl$Genus == "Pseudomonas")
#Methylobacterium: the most in ctl, but also one in ctl
any(sigtabGen_neo$Genus == "Arsenophonus")

#Fig 4A ####
#plot
fig4a = ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
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
           size = 6, colour = "azure4") +
  labs(tag = "A)") + theme(plot.tag = element_text(size = 12, face = "bold"))
library(cowplot)
save_plot("/data/users/mona/miseq_16S/Mona_16S_all/article1/graphs/ParizadehM_Fig4A.pdf", fig4a, ncol = 2, nrow = 1)

#family in phyla ####
# Order genera based on their log2fold
faml = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
faml = sort(faml, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(faml))

#order total number of ASVs associated to each genus
sort(table(sigtab$Family), decreasing = TRUE)
#order neonic based on total number of ASVs
sort(table(sigtab_neo$Family), decreasing = TRUE)
#order control based on total number of ASVs
sort(table(sigtab_ctl$Family), decreasing = TRUE)
#plot
ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + 
  theme_classic() +
  geom_hline(yintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=3) + 
  scale_color_manual(values=c("cornflowerblue","indianred1","mediumvioletred","darkolivegreen4")) +
  theme(axis.text.x = element_text(size = 16, angle = -90, hjust = 0, vjust=0.5),
        axis.title = element_text( size = 16, face = "bold"),
        legend.title = element_text(size=16, face="bold"),
        legend.text=element_text(size=16),
        legend.position = "left") +
  annotate("text", x = 12.8, y = 5, label = 'atop(bold("Neonicotinoid-treated"))', parse = TRUE, 
           size = 8, colour = "azure4") +
  annotate("text", x = 13.5, y = -1.5, label = 'atop(bold("Control"))', parse = TRUE, 
           size = 8, colour = "azure4") +
  labs(tag = "A)") + theme(plot.tag = element_text(size = 18, face = "bold"))

#% Subset 2016 ####
ps.16 = subset_samples(ps, sample_data(ps)$year == "2016")
ps.16 = prune_taxa(taxa_sums(ps.16)>0, ps.16)
#Phyloseq to deseq2 conversion ####
#neonic
neo.phTOds.16 = phyloseq_to_deseq2(ps.16, design = ~ neonic) #dds file
#estimate size factors
neo.fcs.16 = estimateSizeFactors(neo.phTOds.16)
#no need for calculating the geometric means
#Bayesian estimation of dispersion
neo.dsp.16 = estimateDispersions(neo.fcs.16)
plotDispEsts(neo.dsp.16)
#DSEeq ####
neo.dds.16 = DESeq(neo.phTOds.16, test = "Wald", fitType="local")
#Investigate test results table ####
resultsNames(neo.dds.16)
neo.res.16 = results(neo.dds.16)
neo.res.16 = neo.res.16[order(neo.res.16$padj, na.last=NA), ] #remove padj NAs
kable(head(neo.res.16)) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
hist(neo.res.16$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 P-value", ylab="Frequency")
#MA plot ####
plotMA(neo.res.16) 
legend("bottomright", legend="differentially abundant", lty=-1, pch=1, col="red", bty='n')
#red points: adjusted p value less than 0.1 (threshold)
#Set padj the threshold ####
neo.sigtab.16 = neo.res.16[(neo.res.16$padj < alpha), ]
#Combine tax with results ####
neo.sigtab.16 = cbind(as(neo.sigtab.16, "data.frame"), as(tax_table(ps.16)[rownames(neo.sigtab.16), ], "matrix"))
kable(head(neo.sigtab.16)) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
dim(neo.sigtab.16)
dim(neo.sigtab.16[which(neo.sigtab.16$log2FoldChange>0),]) #Y
dim(neo.sigtab.16[which(neo.sigtab.16$log2FoldChange<0),]) #N

#% Subset 2017 ####
ps.17 = subset_samples(ps, sample_data(ps)$year == "2017")
ps.17 = prune_taxa(taxa_sums(ps.17)>0, ps.17)
#Phyloseq to deseq2 conversion ####
#neonic
neo.phTOds.17 = phyloseq_to_deseq2(ps.17, design = ~ neonic) #dds file
#estimate size factors
neo.fcs.17 = estimateSizeFactors(neo.phTOds.17)
#no need for calculating the geometric means
#Bayesian estimation of dispersion
neo.dsp.17 = estimateDispersions(neo.fcs.17)
plotDispEsts(neo.dsp.17)
#DSEeq ####
neo.dds.17 = DESeq(neo.phTOds.17, test = "Wald", fitType="local")
#Investigate test results table ####
resultsNames(neo.dds.17)
neo.res.17 = results(neo.dds.17)
neo.res.17 = neo.res.17[order(neo.res.17$padj, na.last=NA), ] #remove padj NAs
kable(head(neo.res.17)) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
hist(neo.res.17$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 P-value", ylab="Frequency")
#MA plot ####
plotMA(neo.res.17) 
legend("bottomright", legend="differentially abundant", lty=-1, pch=1, col="red", bty='n')
#red points: adjusted p value less than 0.1 (threshold)
#Set padj the threshold ####
neo.sigtab.17 = neo.res.17[(neo.res.17$padj < alpha), ]
#Combine tax with results ####
neo.sigtab.17 = cbind(as(neo.sigtab.17, "data.frame"), as(tax_table(ps.17)[rownames(neo.sigtab.17), ], "matrix"))
kable(neo.sigtab.17) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
dim(neo.sigtab.17)
dim(neo.sigtab.17[which(neo.sigtab.17$log2FoldChange>0),]) #Y
dim(neo.sigtab.17[which(neo.sigtab.17$log2FoldChange<0),]) #N

#% Subset 2018 ####
ps.18 = subset_samples(ps, sample_data(ps)$year == "2018")
ps.18 = prune_taxa(taxa_sums(ps.18)>0, ps.18)
#Phyloseq to deseq2 conversion ####
#neonic
neo.phTOds.18 = phyloseq_to_deseq2(ps.18, design = ~ neonic) #dds file
#estimate size factors
neo.fcs.18 = estimateSizeFactors(neo.phTOds.18) 
#Bayesian estimation of dispersion
neo.dsp.18 = estimateDispersions(neo.fcs.18)
plotDispEsts(neo.dsp.18)
#DSEeq ####
neo.dds.18 = DESeq(neo.dsp.18, test = "Wald", fitType="local")
#Investigate test results table ####
resultsNames(neo.dds.18)
neo.res.18 = results(neo.dds.18)
neo.res.18 = neo.res.18[order(neo.res.18$padj, na.last=NA), ] #remove padj NAs
kable(head(neo.res.18)) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
hist(neo.res.18$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 P-value", ylab="Frequency")
#MA plot ####
plotMA(neo.res.18) 
legend("bottomright", legend="differentially abundant", lty=-1, pch=1, col="red", bty='n')
#red points: adjusted p value less than 0.1 (threshold)
#Set padj the threshold ####
neo.sigtab.18 = neo.res.18[(neo.res.18$padj < alpha), ]
#Combine tax with results ####
neo.sigtab.18 = cbind(as(neo.sigtab.18, "data.frame"), as(tax_table(ps.18)[rownames(neo.sigtab.18), ], "matrix"))
kable(neo.sigtab.18) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
dim(neo.sigtab.18)
dim(neo.sigtab.18[which(neo.sigtab.18$log2FoldChange>0),]) #Y
dim(neo.sigtab.18[which(neo.sigtab.18$log2FoldChange<0),]) #N

save.image("/data/users/mona/miseq_16S/Mona_16S_all/article1/a1_3_p2_deseq_aca_16S_phyl.RData")

