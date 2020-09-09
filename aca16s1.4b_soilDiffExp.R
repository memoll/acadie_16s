# Explanatory analysis of soil 16S - DESeq2 - article
# Studying the effects of neonic on soil bacteria
# Data: Miseq-16S - Acadie
# Mona Parizadeh - June 2020

library(phyloseq); packageVersion("phyloseq") #‘1.27.6’
library(vegan); packageVersion("vegan") #‘2.5.6’
library(ggplot2); packageVersion("ggplot2")
library(DESeq2); packageVersion("DESeq2") #‘1.29.4’
library(dplyr); packageVersion("dplyr") #‘0.8.5’
library(kableExtra); packageVersion("kableExtra") #‘1.1.0.9000’

# Import data #### 
setwd("/data/users/mona/miseq_16S/Mona_16S_all/")
ps = readRDS("article1/16S_aca_soil.rds")
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

#plot
#fold change is for measuring the changes in gene expression
ggplot(sigtab, aes(x=rownames(sigtab), y=log2FoldChange, color=Phylum)) + 
  geom_hline(yintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size = 5))

#genus in phylum ####
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
length(gen)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(gen))

#order total number of ASVs associated to each genus
sort(table(sigtab$Genus), decreasing = TRUE)
#subset ASVs associated to neonic (>0)
sigtab_neo = sigtab[which(sigtab$log2FoldChange>0),];dim(sigtab_neo)
#order based on total number of ASVs
sort(table(sigtab_neo$Genus), decreasing = TRUE)
kable(sort(table(sigtab_neo$Genus), decreasing = TRUE)) %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = F) %>%
  add_header_above(c("Neonicotinoid-treated","")) %>%
  row_spec(1:2, bold = T, color = "white", background = "#D7261E") 

sort(table(sigtab_neo$Phylum), decreasing = TRUE)


#subset ASVs associated to control (<0)
sigtab_ctl = sigtab[which(sigtab$log2FoldChange<0),];dim(sigtab_ctl)
#order based on total number of ASVs
sort(table(sigtab_ctl$Genus), decreasing = TRUE)
kable(sort(table(sigtab_ctl$Genus), decreasing = TRUE)) %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = F) %>%
  add_header_above(c("Control","")) %>%
  row_spec(1:2, bold = T, color = "white", background = "#D7261E") 

sort(table(sigtab_ctl$Phylum), decreasing = TRUE)

#check for beneficial bacteria (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5428506/pdf/41598_2017_Article_472.pdf) repression
any(sigtab_ctl$Genus == "Bacillus")
any(sigtab_neo$Genus == "Bacillus")
#any(sigtab_ctl$Genus == "Agromyces")
#any(sigtab_neo$Genus == "Agromyces")
#any(sigtab_ctl$Genus == "Micromonospora")
#any(sigtab_neo$Genus == "Micromonospora")
#any(sigtab_ctl$Genus == "Pseudonocardia") #T
#any(sigtab_neo$Genus == "Pseudonocardia") #T
#any(sigtab_ctl$Genus == "Acremonium")
#any(sigtab_neo$Genus == "Acremonium")
#any(sigtab_ctl$Genus == "Lysobacter")
#any(sigtab_neo$Genus == "Lysobacter")
any(sigtab_ctl$Genus == "Mesorhizobium")
any(sigtab_neo$Genus == "Mesorhizobium")
any(sigtab_ctl$Genus == "Microvirga")
any(sigtab_neo$Genus == "Microvirga")
any(sigtab_ctl$Genus == "Bradyrhizobium")
any(sigtab_neo$Genus == "Bradyrhizobium")
#any(sigtab_ctl$Genus == "Acremonium")
#any(sigtab_neo$Genus == "Acremonium")
#any(sigtab_ctl$Genus == "Chaetomium")
#any(sigtab_neo$Genus == "Chaetomium")
#any(sigtab_ctl$Genus == "Pseudomonas")
#any(sigtab_neo$Genus == "Pseudomonas")
#any(sigtab_ctl$Genus == "Nocardioides")
#any(sigtab_neo$Genus == "Nocardioides")
#any(sigtab_ctl$Genus == "Micromonospora")
#any(sigtab_neo$Genus == "Micromonospora")
#Benjamin's comment
any(sigtab_ctl$Genus == "Bosea")
any(sigtab_neo$Genus == "Bosea")
any(sigtab_ctl$Genus == "Nitrospira")
any(sigtab_neo$Genus == "Nitrospira")
any(sigtab_ctl$Genus == "Nitrosospira")
any(sigtab_neo$Genus == "Nitrosospira")
any(sigtab_ctl$Genus == "Ammoniphilus")
any(sigtab_neo$Genus == "Ammoniphilus")
any(sigtab_ctl$Genus == "Hyphomicrobium")
any(sigtab_neo$Genus == "Hyphomicrobium")
any(sigtab_ctl$Genus == "Rhodanobacter")
any(sigtab_neo$Genus == "Rhodanobacter")
any(sigtab_ctl$Genus == "Rhizobacter")
any(sigtab_neo$Genus == "Rhizobacter")
#supressed
any(sigtab_ctl$Genus == "Mycobacterium")
any(sigtab_neo$Genus == "Mycobacterium")
any(sigtab_ctl$Genus == "Arthrobacter")
any(sigtab_neo$Genus == "Arthrobacter")
#Filimon 2014:
any(sigtab_ctl$Genus == "Azotobacter")
any(sigtab_neo$Genus == "Azotobacter")
#Pang 2020 (biodegraders):
#any(sigtab_ctl$Genus == "Pseudoxanthomonas")
#any(sigtab_neo$Genus == "Pseudoxanthomonas")
#any(sigtab_ctl$Genus == "Rhizobium")
#any(sigtab_neo$Genus == "Rhizobium")
#any(sigtab_ctl$Genus == "Rhodococcus")
#any(sigtab_neo$Genus == "Rhodococcus")
#any(sigtab_ctl$Genus == "Actinomycetes")
#any(sigtab_neo$Genus == "Actinomycetes")
#any(sigtab_ctl$Genus == "Stenotrophomonas")
#any(sigtab_neo$Genus == "Stenotrophomonas")
#any(sigtab_ctl$Genus == "Sinorhizobium")
#any(sigtab_neo$Genus == "Sinorhizobium")
#any(sigtab_ctl$Genus == "Pseudomonas")
#any(sigtab_neo$Genus == "Pseudomonas")
#any(sigtab_ctl$Genus == "Sphingomonas")
#any(sigtab_neo$Genus == "Sphingomonas")
#any(sigtab_ctl$Genus == "Acinetobacter")
#any(sigtab_neo$Genus == "Acinetobacter")
#any(sigtab_ctl$Genus == "Aerophilus")
#any(sigtab_neo$Genus == "Aerophilus")
any(sigtab_ctl$Genus == "Pigmentiphaga")
any(sigtab_neo$Genus == "Pigmentiphaga")

#Benjamin:
any(sigtab_ctl$Genus == "Mycobacterium") #& Pang
any(sigtab_neo$Genus == "Mycobacterium")
any(sigtab_ctl$Genus == "Arthrobacter")
any(sigtab_neo$Genus == "Arthrobacter")
any(sigtab_ctl$Genus == "Streptomyces")
any(sigtab_neo$Genus == "Streptomyces")
#Fig 4B ####
#Plot 
fig4b = ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  theme_classic() +
  geom_hline(yintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=3) +
  scale_color_manual(values = c("#7CE3D8","darkgoldenrod1","mediumvioletred","indianred1","tan3",
                                "cornflowerblue","seagreen","red2","tan4","yellowgreen","tomato3")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, hjust = 0, vjust=0.5),
        axis.title = element_text( size = 12, face = "bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text=element_text(size=12),
        legend.position = "left") +
  annotate("text", x = 53.5, y = 1.5, label = 'atop(bold("Neonicotinoid-treated"))', parse = TRUE, 
           size = 6, colour = "azure4") +
  annotate("text", x = 58, y = -3, label = 'atop(bold("Control"))', parse = TRUE, 
           size = 6, colour = "azure4") +
  labs(tag = "B)") + theme(plot.tag = element_text(size = 12, face = "bold"))
library(cowplot)
save_plot("/data/users/mona/miseq_16S/Mona_16S_all/article1/graphs/ParizadehM_Fig4B.pdf", fig4b, ncol = 2, nrow = 2)

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
  scale_color_manual(values = c("#7CE3D8","darkgoldenrod1","mediumvioletred","indianred1","tan3",
                                "cornflowerblue","seagreen","red2","tan4","yellowgreen","tomato3")) +
  theme(axis.text.x = element_text(size = 10, angle = -90, hjust = 0, vjust=0.5),
        axis.title = element_text( size = 12, face = "bold"),
        legend.title = element_text(size=12, face="bold"),
        legend.text=element_text(size=12),
        legend.position = "left") +
  annotate("text", x = 50.5, y = 1.5, label = 'atop(bold("Neonicotinoid-treated"))', parse = TRUE, 
           size = 6, colour = "azure4") +
  annotate("text", x = 55.5, y = -0.5, label = 'atop(bold("Control"))', parse = TRUE, 
           size = 6, colour = "azure4") +
  labs(tag = "B)") + theme(plot.tag = element_text(size = 12, face = "bold"))

#% Extra#### 
#not in the article
#Relative Abundance ####
#add a new column including both months and years 
sample_data(ps)$mnt_yr = as.factor(paste(sample_data(ps)$month, sample_data(ps)$year, sep="_"))
sample_data(ps)$mnt_yr = factor(sample_data(ps)$mnt_yr,levels = c("July_2016","Aug_2016","Sep_2016",
                                                                  "July_2017","Aug_2017","Sep_2017",
                                                                  "July_2018","Aug_2018","Sep_2018"))
#change it to numeric
sample_data(ps)$mnt_yr = as.numeric(sample_data(ps)$mnt_yr)
#make a new phyloseq object containing only the ASVs in sigtab
ps.sigtab = ps %>%
  subset_taxa(rownames(tax_table(ps)) %in% rownames(sigtab))
#make dataframe
melt.sigtab = ps.sigtab %>%
  transform_sample_counts(function(otu) otu/sum(otu)) %>%
  psmelt()
summary(melt.sigtab$Abundance)
#choose a threshod of relative abundance
melt.sigtab1 = subset(melt.sigtab, Abundance > 0.02)
kable(melt.sigtab1) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
#number of NA genera
subset_taxa(ps.sigtab,is.na(Genus)) #145/294
#or
dim(subset(sigtab, is.na(Genus)));dim(sigtab)

#ASV Frequency ####
#make dataframe
tab.sigtab = as.data.frame(table(sigtab$Genus))
filt = tab.sigtab %>% 
  filter(Freq>=3) #filter based on the frequency of each ASV
melt.filt = melt.sigtab %>%
  subset(melt.sigtab$Genus %in% filt$Var1)
#filter(Genus == "Ellin6067" | Genus == "Gemmatimonas" | Genus == "Gaiella" | Genus == "Mycobacterium")

#boxplot
melt.filt$month = factor(melt.filt$month,levels = c("July","Aug","Sep"))
ggplot(data = melt.filt, mapping = aes_string(x = "neonic",y = "Abundance",
                                              color = "Genus", fill = "Genus")) +
  theme_bw() +
  ylab("Relative Abundance - log10") +
  geom_boxplot() +
  #geom_jitter(alpha=0.1) +
  facet_grid(month~year) + 
  scale_y_log10() +
  scale_color_manual(values = c("#7CE3D8","darkgoldenrod1","mediumvioletred","seagreen","tan3",
                                "cornflowerblue","indianred1","red2","tan4","yellowgreen",
                                "mediumorchid1","tomato3","lightsalmon","wheat4","darkred",
                                "darkgray")) +
  scale_fill_manual(values = c("#7CE3D8","darkgoldenrod1","mediumvioletred","seagreen","tan3",
                               "cornflowerblue","indianred1","red2","tan4","yellowgreen",
                               "mediumorchid1","tomato3","lightsalmon","wheat4","darkred",
                               "darkgray")) +
  scale_x_discrete(name = "Treatment",
                   labels = c(N="Control",Y="Neonicotinoid-treated")) +
  theme(legend.position="right")

#loess plot
ggplot(melt.filt, aes(x = mnt_yr, y = Abundance, group = Genus, color = Genus)) + 
  theme_bw() +
  #geom_jitter() +
  geom_smooth(alpha=0.05) +
  ylab("Relative Abundance (log 10)") +
  scale_y_log10() +
  scale_x_continuous(name="\nTime", breaks = 1:9,
                     labels= c("July 2016","August 2016","September 2016",
                               "July 2017","August 2017","September 2017",
                               "July 2018","August 2018","September 2018")) +
  scale_color_manual(values = c("#7CE3D8","darkgoldenrod1","mediumvioletred","seagreen","tan3",
                                "cornflowerblue","indianred1","red2","tan4","yellowgreen",
                                "mediumorchid1","tomato3","lightsalmon","wheat4","darkred",
                                "darkgray")) + 
  facet_wrap(~neonic, ncol = 2,
             labeller=labeller(neonic = c(N="Control",Y="Neonicotinoid-treated"))) +
  theme(legend.title = element_text(size=16, face="bold"),
        legend.text=element_text(size=16),
        strip.text.x = element_text(size=16, face="bold"),
        axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=14, hjust = 0, angle = -45),
        axis.title = element_text(size = 16, face = "bold"))

#% Subset 2016 ####
ps.16 = subset_samples(ps, sample_data(ps)$year == "2016")
ps.16 = prune_taxa(taxa_sums(ps.16)>0, ps.16)
#Phyloseq to deseq2 conversion 
#neonic
neo.phTOds.16 = phyloseq_to_deseq2(ps.16, design = ~ neonic) #dds file
#estimate size factors
neo.fcs.16 = estimateSizeFactors(neo.phTOds.16)
#no need for calculating the geometric means
#Bayesian estimation of dispersion
neo.dsp.16 = estimateDispersions(neo.fcs.16)
plotDispEsts(neo.dsp.16)
#DSEeq 
neo.dds.16 = DESeq(neo.phTOds.16, test = "Wald", fitType="local")
#Investigate test results table 
resultsNames(neo.dds.16)
neo.res.16 = results(neo.dds.16)
neo.res.16 = neo.res.16[order(neo.res.16$padj, na.last=NA), ] #remove padj NAs
kable(head(neo.res.16)) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
hist(neo.res.16$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 P-value", ylab="Frequency")
#MA plot 
plotMA(neo.res.16) 
legend("bottomright", legend="differentially abundant", lty=-1, pch=1, col="red", bty='n')
#red points: adjusted p value less than 0.1 (threshold)
#Set padj the threshold 
neo.sigtab.16 = neo.res.16[(neo.res.16$padj < alpha), ]
#Combine tax with results 
neo.sigtab.16 = cbind(as(neo.sigtab.16, "data.frame"), as(tax_table(ps.16)[rownames(neo.sigtab.16), ], "matrix"))
kable(head(neo.sigtab.16)) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
dim(neo.sigtab.16)
dim(neo.sigtab.16[which(neo.sigtab.16$log2FoldChange>0),]) #Y
dim(neo.sigtab.16[which(neo.sigtab.16$log2FoldChange<0),]) #N
#subset ASVs associated to neonic (>0)
neo.sigtab.neo16 = sigtab[which(neo.sigtab.16$log2FoldChange>0),];dim(neo.sigtab.neo16)
#order based on total number of ASVs
sort(table(neo.sigtab.neo16$Genus), decreasing = TRUE)
kable(sort(table(neo.sigtab.neo16$Genus), decreasing = TRUE)) %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = F) %>%
  add_header_above(c("Neonicotinoid-treated 2016","")) %>%
  row_spec(1:2, bold = T, color = "white", background = "#D7261E") 
#subset ASVs associated to control (<0)
neo.sigtab.ctl16 = sigtab[which(neo.sigtab.16$log2FoldChange<0),];dim(neo.sigtab.ctl16)
#order based on total number of ASVs
sort(table(neo.sigtab.ctl16$Genus), decreasing = TRUE)
kable(sort(table(neo.sigtab.ctl16$Genus), decreasing = TRUE)) %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = F) %>%
  add_header_above(c("Control 2016","")) %>%
  row_spec(1:2, bold = T, color = "white", background = "#D7261E") 

genPhlm16 = tapply(neo.sigtab.16$log2FoldChange, neo.sigtab.16$Phylum, function(x) max(x))
sort(table(neo.sigtab.neo16$Phylum), decreasing = TRUE)
sort(table(neo.sigtab.ctl16$Phylum), decreasing = TRUE)

#% Subset 2017 ####
ps.17 = subset_samples(ps, sample_data(ps)$year == "2017")
ps.17 = prune_taxa(taxa_sums(ps.17)>0, ps.17)
#Phyloseq to deseq2 conversion 
#neonic
neo.phTOds.17 = phyloseq_to_deseq2(ps.17, design = ~ neonic) #dds file
#estimate size factors
neo.fcs.17 = estimateSizeFactors(neo.phTOds.17)
#no need for calculating the geometric means
#Bayesian estimation of dispersion
neo.dsp.17 = estimateDispersions(neo.fcs.17)
plotDispEsts(neo.dsp.17)
#DSEeq 
neo.dds.17 = DESeq(neo.phTOds.17, test = "Wald", fitType="local")
#Investigate test results table
resultsNames(neo.dds.17)
neo.res.17 = results(neo.dds.17)
neo.res.17 = neo.res.17[order(neo.res.17$padj, na.last=NA), ] #remove padj NAs
kable(head(neo.res.17)) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
hist(neo.res.17$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 P-value", ylab="Frequency")
#MA plot 
plotMA(neo.res.17) 
legend("bottomright", legend="differentially abundant", lty=-1, pch=1, col="red", bty='n')
#red points: adjusted p value less than 0.1 (threshold)
#Set padj the threshold 
neo.sigtab.17 = neo.res.17[(neo.res.17$padj < alpha), ]
#Combine tax with results 
neo.sigtab.17 = cbind(as(neo.sigtab.17, "data.frame"), as(tax_table(ps.17)[rownames(neo.sigtab.17), ], "matrix"))
kable(neo.sigtab.17) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
dim(neo.sigtab.17)
dim(neo.sigtab.17[which(neo.sigtab.17$log2FoldChange>0),]) #Y
dim(neo.sigtab.17[which(neo.sigtab.17$log2FoldChange<0),]) #N
#subset ASVs associated to neonic (>0)
neo.sigtab.neo17 = sigtab[which(neo.sigtab.17$log2FoldChange>0),];dim(neo.sigtab.neo17)
#order based on total number of ASVs
sort(table(neo.sigtab.neo17$Genus), decreasing = TRUE)
kable(sort(table(neo.sigtab.neo17$Genus), decreasing = TRUE)) %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = F) %>%
  add_header_above(c("Neonicotinoid-treated 2017","")) %>%
  row_spec(1:2, bold = T, color = "white", background = "seagreen") 
#subset ASVs associated to control (<0)
neo.sigtab.ctl17 = sigtab[which(neo.sigtab.17$log2FoldChange<0),];dim(neo.sigtab.ctl17)
#order based on total number of ASVs
sort(table(neo.sigtab.ctl17$Genus), decreasing = TRUE)
kable(sort(table(neo.sigtab.ctl17$Genus), decreasing = TRUE)) %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = F) %>%
  add_header_above(c("Control 2017","")) 

genPhlm17 = tapply(neo.sigtab.17$log2FoldChange, neo.sigtab.17$Phylum, function(x) max(x))
sort(table(neo.sigtab.neo17$Phylum), decreasing = TRUE)
sort(table(neo.sigtab.ctl17$Phylum), decreasing = TRUE)

#% Subset 2018 ####
ps.18 = subset_samples(ps, sample_data(ps)$year == "2018")
ps.18 = prune_taxa(taxa_sums(ps.18)>0, ps.18)
#Phyloseq to deseq2 conversion 
#neonic
neo.phTOds.18 = phyloseq_to_deseq2(ps.18, design = ~ neonic) #dds file
#estimate size factors
neo.fcs.18 = estimateSizeFactors(neo.phTOds.18) 
#Bayesian estimation of dispersion
neo.dsp.18 = estimateDispersions(neo.fcs.18)
plotDispEsts(neo.dsp.18)
#DSEeq 
neo.dds.18 = DESeq(neo.dsp.18, test = "Wald", fitType="local")
#Investigate test results table 
resultsNames(neo.dds.18)
neo.res.18 = results(neo.dds.18)
neo.res.18 = neo.res.18[order(neo.res.18$padj, na.last=NA), ] #remove padj NAs
kable(head(neo.res.18)) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
hist(neo.res.18$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 P-value", ylab="Frequency")
#MA plot 
plotMA(neo.res.18) 
legend("bottomright", legend="differentially abundant", lty=-1, pch=1, col="red", bty='n')
#red points: adjusted p value less than 0.1 (threshold)
#Set padj the threshold 
neo.sigtab.18 = neo.res.18[(neo.res.18$padj < alpha), ]
#Combine tax with results 
neo.sigtab.18 = cbind(as(neo.sigtab.18, "data.frame"), as(tax_table(ps.18)[rownames(neo.sigtab.18), ], "matrix"))
kable(head(neo.sigtab.18)) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)
dim(neo.sigtab.18)
dim(neo.sigtab.18[which(neo.sigtab.18$log2FoldChange>0),]) #Y
dim(neo.sigtab.18[which(neo.sigtab.18$log2FoldChange<0),]) #N
#subset ASVs associated to neonic (>0)
neo.sigtab.neo18 = sigtab[which(neo.sigtab.18$log2FoldChange>0),];dim(neo.sigtab.neo18)
#order based on total number of ASVs
sort(table(neo.sigtab.neo18$Genus), decreasing = TRUE)
kable(sort(table(neo.sigtab.neo18$Genus), decreasing = TRUE)) %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = F) %>%
  add_header_above(c("Neonicotinoid-treated 2018","")) %>%
  row_spec(1:2, bold = T, color = "white", background = "cornflowerblue") 
#subset ASVs associated to control (<0)
neo.sigtab.ctl18 = sigtab[which(neo.sigtab.18$log2FoldChange<0),];dim(neo.sigtab.ctl18)
#order based on total number of ASVs
sort(table(neo.sigtab.ctl18$Genus), decreasing = TRUE)
kable(sort(table(neo.sigtab.ctl18$Genus), decreasing = TRUE)) %>%
  kable_styling(bootstrap_options = c("striped","hover"), full_width = F) %>%
  add_header_above(c("Control 2018","")) %>%
  row_spec(1:2, bold = T, color = "white", background = "cornflowerblue")

genPhlm18 = tapply(neo.sigtab.18$log2FoldChange, neo.sigtab.18$Phylum, function(x) max(x))
sort(table(neo.sigtab.neo18$Phylum), decreasing = TRUE)
sort(table(neo.sigtab.ctl18$Phylum), decreasing = TRUE)




