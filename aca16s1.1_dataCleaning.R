###############################################################
# Cleaning and denoising 16S data                             #
# Data: Miseq-16S - all Runs - Subset L'Acadie (ACA)          #
# Mona Parizadeh - 2019-2020                                  #
###############################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq") #‘1.27.6’
library(vegan); packageVersion("vegan") #‘2.5.6’
library(ggplot2); packageVersion("ggplot2") #‘3.3.0’
library(tidyverse); packageVersion("tidyverse") #‘1.3.0.9000’

# Import data #### 
setwd("../mp/aca_16s/files/")
ps = readRDS("ps.rds") 
# Explore data 
nsamples(ps)
ntaxa(ps)
sample_variables(ps)
rank_names(ps)
subset_samples(ps, sample_data(ps)$habitat == "leaf")
subset_samples(ps, sample_data(ps)$habitat == "soil")
subset_samples(ps, sample_data(ps)$neonic == "Y") #neonic-treated
subset_samples(ps, sample_data(ps)$neonic == "N") #control (non-treated)
subset_samples(ps, sample_data(ps)$sample_or_control == "control") #negative controls
#number of seq per sample
summary(sample_sums(ps))  #or: summary(apply(comm,1,sum))
sd(sample_sums(ps), na.rm=TRUE)/sqrt(length(sample_sums(ps)[!is.na(sample_sums(ps))])) #SE
head(sort(sample_sums(ps),TRUE))
hist(sample_sums(ps))
#ASV richness 
summary(estimate_richness(ps, measures = "Observed"))
#distribution of ASVs
hist(log10(taxa_sums(ps))) 
# ASV (OTU) table 
taxa_names(ps) = paste0("ASV", seq(ntaxa(ps))) #replace sequence w/ ASV
otu_mat = function(ps) as(otu_table(ps), "matrix")
otu_mat(ps)[1:5,1:5] 
# taxonomic table 
tax_mat = function(ps) as(tax_table(ps), "matrix")
tax_mat(ps)[1:5,]

#Explore controls ####
#%Positive controls ####
ct.posID = sample_names(ps)[grep("CTL..00",sample_names(ps))] 
ps.pos.ctl = subset_samples(ps, sample_data(ps)$sampleid %in% ct.posID)
ps.pos.ctl = prune_taxa(taxa_sums(ps.pos.ctl)>0, ps.pos.ctl)
#plot 
ps.pos.ctl %>% 
  tax_glom(taxrank = "Genus", NArm=TRUE) %>% #remove NAs
  #transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to relative abundance
  plot_bar(fill="Genus") +
  labs(title = "Bacteria genera present in the positive controls") +
  xlab("Positive Control IDs") + ylab("Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "left")
#identify
mdt.pos.ctl = ps.pos.ctl %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt() 
#check if the bacteria in the positive controls have been sequenced and identiifed 
mdt.pos.ctl$sampleid[mdt.pos.ctl$Genus == "Clavibacter"] 
mdt.pos.ctl$sampleid[mdt.pos.ctl$Genus == "Pectobacterium"] 
mdt.pos.ctl$sampleid[mdt.pos.ctl$Genus == "Escherichia/Shigella"] 
mdt.pos.ctl$sampleid[mdt.pos.ctl$Genus == "Pantoea"]
mdt.pos.ctl$sampleid[mdt.pos.ctl$Genus == "Xanthomonas"] 

#%Negative controls ####
ps.neg_ctl = subset_samples(ps, sample_data(ps)$sample_or_control == "control")
ps.neg_ctl = prune_taxa(taxa_sums(ps.neg_ctl)>0, ps.neg_ctl)
#plot
ps.neg_ctl %>% 
  tax_glom(taxrank = "Genus", NArm=TRUE) %>% #keep NAs, glom based on genus
  plot_bar(fill="Genus") +
  labs(title = "Bacteria families present in the negative controls \nAbsolute abundance") +
  xlab("Negative Control IDs") + ylab("Abundance") +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1))

#Keep only L'Acadie & negative CTRLs (or remove positive CTRLs) ####
#L'Acadie data contains the 3rd phyllosphere replicates (r3), which need to be removed
ps.aca.neg_ctl = subset_samples(ps, sample_data(ps)$site == "ACA" | sample_data(ps)$sampleid %in% sample_data(ps.neg_ctl)$sampleid)
ps.aca.neg_ctl = prune_taxa(taxa_sums(ps.aca.neg_ctl)>0, ps.aca.neg_ctl)
ps.aca.neg_ctl

#Explore contaminants ####
#Find Archaea 
Archaea = subset_taxa(ps.aca.neg_ctl, Kingdom == "Archaea")
(ntaxa(Archaea)/ntaxa(ps.aca.neg_ctl))*100 
# Find Cyanobacteria 
Cyanobacteria = subset_taxa(ps.aca.neg_ctl, Phylum=="Cyanobacteria")
(ntaxa(Cyanobacteria)/ntaxa(ps.aca.neg_ctl))*100 
# Find Chloroplast 
Chloroplast = subset_taxa(ps.aca.neg_ctl, Order=="Chloroplast")  ####keep it for later
(ntaxa(Chloroplast)/ntaxa(ps.aca.neg_ctl))*100 
# Find Mitochondria 
Mitochondria = subset_taxa(ps.aca.neg_ctl, Family=="Mitochondria")
(ntaxa(Mitochondria)/ntaxa(ps.aca.neg_ctl))*100
# Note: we still keep chloroplast and mitochondria for decontam purposes 
((ntaxa(Chloroplast)+ntaxa(Mitochondria))/ntaxa(ps.aca.neg_ctl))*100 

#Denoising & Decontaminating data ####
#1. Remove undefined phyla ####
ps.aca.neg_ctl.noNa = subset_taxa(ps.aca.neg_ctl, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
ps.aca.neg_ctl.noNa = prune_samples(sample_sums(ps.aca.neg_ctl.noNa)>0, ps.aca.neg_ctl.noNa)
ps.aca.neg_ctl.noNa
100-(ntaxa(ps.aca.neg_ctl.noNa)/ntaxa(ps.aca.neg_ctl))*100 

#2. Remove outliers - NMDS####
#relative abundance
ps.ra = transform_sample_counts(ps.aca.neg_ctl.noNa, function(otu) otu/sum(otu)) 
#ordinate
nmds = ordinate(ps.aca.neg_ctl.noNa, method = "NMDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps.ra, nmds, color = "habitat", shape = "sample_or_control") + 
  theme_bw() + geom_point() + ggtitle("nMDS") +
  geom_text(aes(label = sampleid), check_overlap = FALSE, size = 5) + 
  geom_point(size = 1) + scale_shape_manual(values = c(19, 1))
#Remove outliers 1
out = c("CTL0004","CTL3001")
ps.aca.neg_ctl.noNa1 = prune_samples(!sample_data(ps.aca.neg_ctl.noNa)$sampleid %in% out, ps.aca.neg_ctl.noNa)
ps.aca.neg_ctl.noNa1 = prune_taxa(taxa_sums(ps.aca.neg_ctl.noNa1)>0,ps.aca.neg_ctl.noNa1)
ps1.ra = transform_sample_counts(ps.aca.neg_ctl.noNa1, function(otu) otu/sum(otu)) 
nmds1 = ordinate(ps1.ra, method = "NMDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps1.ra, nmds1, color = "habitat", shape = "sample_or_control") + 
  theme_bw() +
  geom_point() + ggtitle("nMDS") +
  geom_text(aes(label = sampleid), check_overlap = FALSE, size = 3, nudge_y = -0.1) + 
  geom_point(size = 5) + scale_shape_manual(values = c(19, 1))

#3. Keep samples with at least 1000 reads ####
ps.aca.neg_ctl.1krds=prune_samples(sample_sums(ps.aca.neg_ctl.noNa1)>=1000, ps.aca.neg_ctl.noNa1)
ps.aca.neg_ctl.1krds=prune_taxa(taxa_sums(ps.aca.neg_ctl.1krds)>0, ps.aca.neg_ctl.1krds)
ps.aca.neg_ctl.1krds
#lost samples
ps.lost_1kreads = prune_samples(sample_sums(ps.aca.neg_ctl.noNa1)<1000, ps.aca.neg_ctl.noNa1)
ps.lost_1kreads
sample_names(ps.lost_1kreads)[grep("CTL..0.",sample_names(ps.lost_1kreads))] #lost
sample_names(ps.aca.neg_ctl.1krds)[grep("CTL..0.",sample_names(ps.aca.neg_ctl.1krds))] #left

#4. Remove contaminants ####
library(decontam); packageVersion("decontam") #‘1.1.2’
# 4.1. Inspect library sizes ####
# Assign the full sample_data
sample_df = function(ps) as(sample_data(ps), "data.frame")
sam.df = sample_df(ps.aca.neg_ctl.1krds)  # Put sample_data into a ggplot-friendly data.frame
sam.df$LibrarySize = sample_sums(ps.aca.neg_ctl.1krds)
sample_data(ps.aca.neg_ctl.1krds) = sample_data(sam.df)
# The shorthand way, assign one column
sample_data(ps.aca.neg_ctl.1krds)$LibrarySize = sample_sums(ps.aca.neg_ctl.1krds)
df = sam.df[order(sam.df$LibrarySize),]
df$Index = seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample_or_control)) + geom_point()
# The library sizes of the samples primarily goes from 25 000 to 100 000 reads, but there are some high-read outliers.
# As expected, the negative control samples have very fewer reads.

# 4.2. Prevalence - Identify Contaminants ####
sample_data(ps.aca.neg_ctl.1krds)$is.neg = sample_data(ps.aca.neg_ctl.1krds)$sample_or_control == "control"
contamdf.prev = isContaminant(ps.aca.neg_ctl.1krds, method="prevalence", neg="is.neg")
# the default threshold for a contaminant is that it reaches a probability of 0.1
table(contamdf.prev$contaminant) #ASVs being recognized as contaminant (TRUE) or not (FALSE)
which(contamdf.prev$contaminant) 
# getting vector holding the identified contaminant IDs
contam_asvs = row.names(contamdf.prev[contamdf.prev$contaminant == TRUE, ])
length(contam_asvs)
contam_taxa = as.data.frame(tax_mat(ps.aca.neg_ctl.1krds))
#In the prevalence test there is a special value worth knowing,  threshold=0.5
#This will identify as contaminants all sequences that are more prevalent in negative controls than in samples.
contamdf.prev05 = isContaminant(ps.aca.neg_ctl.1krds, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
which(contamdf.prev05$contaminant)
# The prevalence can find more number of contaminants than the frequency method
# but it might miss the very frequent contaminants is all the samples (even the negative ones).

# getting vector holding the identified contaminant IDs
contam05_asvs = row.names(contamdf.prev05[contamdf.prev05$contaminant == TRUE, ])
length(contam05_asvs)
contam05_taxa = as.data.frame(tax_mat(ps.aca.neg_ctl.1krds))

# Visualize abundance of potential contaminant ASVs - prevalence cutoff 0.1
ps.pa = transform_sample_counts(ps.aca.neg_ctl.1krds, function(abund) 1*(abund>0))
ps.pa_ctl = prune_samples(sample_data(ps.pa)$sample_or_control == "control", ps.pa) # negative samples
sample_data(ps.pa_ctl)$sampleid
ps.pa_smp = prune_samples(sample_data(ps.pa)$sample_or_control == "sample", ps.pa)  # the real samples
nsamples(ps.pa_smp)
nsamples(ps.pa_ctl)
# Make data.frame of prevalence in positive and negative samples
df.pa = data.frame(pa.smp=taxa_sums(ps.pa_smp), pa.ctl=taxa_sums(ps.pa_ctl),
                   contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.ctl, y=pa.smp, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") +
  ggtitle("Abundance of potential contaminant ASVs - prevalence cutoff 0.1")

# Visualize abundance of potential contaminant ASVs - prevalence cutoff 0.5
#threshold (cutoff): The probability threshold below which the null-hypothesis (not a contaminant) should be rejected in favor of the alternate hypothesis (contaminant)
#Make data.frame of prevalence in positive and negative samples
df.pa05 = data.frame(pa.smp=taxa_sums(ps.pa_smp), pa.ctl=taxa_sums(ps.pa_ctl),
                     contaminant=contamdf.prev05$contaminant)
ggplot(data=df.pa05, aes(x=pa.ctl, y=pa.smp, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") + 
  ggtitle("Abundance of potential contaminant ASVs - prevalence cutoff 0.5") 
#Samples (red points) seem to split pretty cleanly into two branches that show up mostly in positive samples and another that shows up mostly in negative controls, and the contaminant assignment (at default probability threshold) has done a good job of identifying those mostly in negative controls.

# taxonomic identity and abundance of the potential contaminant ASVs - prevalence cutoff 0.1
subset_contam = subset_taxa(ps.aca.neg_ctl.1krds, contamdf.prev$contaminant==TRUE)
plot_bar(subset_contam, fill="Genus") + ggtitle("Taxonomic identity and abundance of the potential contaminant ASVs \nat Genus level - prevalence cutoff 0.1")
# Abundance of the contaminant ASVs by taxonomy
taxo.contam = as.data.frame(as(tax_table(subset_contam),"matrix"))
taxo.contam$abund = apply(otu_table(subset_contam),2,sum)
rownames(taxo.contam) = NULL
head(taxo.contam)
                                
# taxonomic identity and abundance of the potential contaminant ASVs - prevalence cutoff 0.5
subset_contam05 = subset_taxa(ps.aca.neg_ctl.1krds, contamdf.prev05$contaminant==TRUE)
plot_bar(subset_contam05, fill="Genus") + ggtitle("Taxonomic identity and abundance of the potential contaminant ASVs \nat Genus level - prevalence cutoff 0.5")
# Abundance of the contaminant ASVs by taxonomy
taxo.contam05 = as.data.frame(as(tax_table(subset_contam05),"matrix"))
taxo.contam05$abund = apply(otu_table(subset_contam05),2,sum)
rownames(taxo.contam05) = NULL
head(taxo.contam05)

# 4.2.a. Prevalence 0.1 ####
ps.notcontam = subset_taxa(ps.aca.neg_ctl.1krds, contamdf.prev$contaminant==FALSE)
ps.notcontam = prune_samples(sample_sums(ps.notcontam)>0, ps.notcontam)
ps.notcontam
# 4.2.b. Agressive Prevalence 0.5 ####
ps.notcontam05 = subset_taxa(ps.aca.neg_ctl.1krds, contamdf.prev05$contaminant==FALSE)
ps.notcontam05 = prune_samples(sample_sums(ps.notcontam05)>0, ps.notcontam05)
ps.notcontam05
#lost ASVs
ntaxa(ps.aca.neg_ctl.1krds) - ntaxa(ps.notcontam05)
# 4.2.c. Remove Chloroplast ####
ps.notcontamChlp = subset_taxa(ps.notcontam05, (Order!="Chloroplast" | is.na(Order)))
# 4.2.d. Remove Mitochondria ####
ps.notcontamChlpMitc = subset_taxa(ps.notcontamChlp, (Family!="Mitochondria") | is.na(Family))
ps.notcontamChlpMitc = prune_samples(sample_sums(ps.notcontamChlpMitc) > 0, ps.notcontamChlpMitc)
ps.notcontamChlpMitc

#%Remaining negative controls ####
ps.negctl.rem = subset_samples(ps.notcontamChlpMitc, sample_data(ps.notcontamChlpMitc)$sampleid %in% sample_data(ps.neg_ctl)$sampleid)
ps.negctl.rem = prune_taxa(taxa_sums(ps.negctl.rem)>0, ps.negctl.rem)
ps.negctl.rem
ps.negctl.rem %>% 
  tax_glom(taxrank = "Genus", NArm=TRUE) %>% #keep NAs
  plot_bar(fill="Genus") +
  labs(title = "Bacteria genera present in the remaining negative control \nAfter decontamination") +
  xlab("Negative Control") + ylab("Abundance") +
  theme(legend.position = "bottom")
ctl0001 = subset_samples(ps.negctl.rem, sample_data(ps.negctl.rem)$sampleid == "CTL0001")
ctl0001 = prune_taxa(taxa_sums(ctl0001)>0, ctl0001)
sum(taxa_sums(ctl0001)) #no. of sequences (reads)
ctl0002 = subset_samples(ps.negctl.rem, sample_data(ps.negctl.rem)$sampleid == "CTL0002")
ctl0002 = prune_taxa(taxa_sums(ctl0002)>0, ctl0002)
sum(taxa_sums(ctl0002)) #no. of sequences (reads)

#5. Remove outliers - alpha diversity ####
aca.ctl.notcontam_asv.shn = estimate_richness(ps.notcontamChlpMitc, split=TRUE, measures="Shannon") 
plot_richness(ps.notcontamChlpMitc, "month","neonic", measures = "Shannon")
#samples w/ very low alpha diversity
div.out = rownames(aca.ctl.notcontam_asv.shn)[(which(aca.ctl.notcontam_asv.shn$Shannon<2))]
div.out
#Remove outliers 
ps.aca.ctl.div2 = prune_samples(!sample_data(ps.notcontamChlpMitc)$sampleid %in% div.out, ps.notcontamChlpMitc)
ps.aca.ctl.div2 = prune_taxa(taxa_sums(ps.aca.ctl.div2)>0,ps.aca.ctl.div2)
ps.aca.ctl.div2

#6. Filter ASVs w/ less than 10 reads ####
ps.aca.ctl_seq10 = prune_taxa(taxa_sums(ps.aca.ctl.div2) > 10, ps.aca.ctl.div2)
ps.aca.ctl_seq10 = prune_samples(sample_sums(ps.aca.ctl_seq10)>0,ps.aca.ctl_seq10) 
ps.aca.ctl_seq10
100-ntaxa(ps.aca.ctl_seq10)/ntaxa(ps.aca.ctl.div2)*100
#%Remaining negative controls ####
ps.negctl.rem1 = subset_samples(ps.aca.ctl_seq10, sample_data(ps.aca.ctl_seq10)$sampleid %in% sample_data(ps.neg_ctl)$sampleid)
ps.negctl.rem1 = prune_taxa(taxa_sums(ps.negctl.rem1)>0, ps.negctl.rem1)
ps.negctl.rem1
sum(taxa_sums(ps.negctl.rem1)) #no. of sequences (reads)

#Clean metadata ####
#look for outliers
#nmds
ps.aca.ra = transform_sample_counts(ps.aca.ctl_seq10, function(otu) otu/sum(otu)) 
nmds.aca = ordinate(ps.aca.ctl_seq10, method = "NMDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps.aca.ra, nmds.aca, color = "habitat", shape = "host") + 
  theme_bw() + geom_point() + ggtitle("nMDS") +
  geom_text(aes(label = sampleid), check_overlap = TRUE, size = 5)  

#shannon richness
plot_richness(ps.aca.ctl_seq10,"habitat","host", measures = "Shannon") +
  geom_text(aes(label = sampleid), check_overlap = FALSE, size = 3)

#remove outlier
out1 = c("CTL0001", "SLA7632", "SDA7831", "SDA7852", "SDA8841")
ps.aca = prune_samples(!sample_data(ps.aca.ctl_seq10)$sampleid %in% out1, ps.aca.ctl_seq10)
ps.aca = prune_taxa(taxa_sums(ps.aca)>0,ps.aca)
ps.aca
#nmds
ps.aca.ra = transform_sample_counts(ps.aca, function(otu) otu/sum(otu)) 
nmds1.aca = ordinate(ps.aca, method = "NMDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps.aca.ra, nmds1.aca, color = "habitat", shape = "host") + 
  theme_bw() + geom_point() + ggtitle("nMDS") +
  geom_text(aes(label = sampleid), check_overlap = FALSE, size = 5) 

#clean the metadata table (remove the unnecessary variables) 
sample_variables(ps.aca)
sample_data(ps.aca)$site = NULL
sample_data(ps.aca)$sample_or_control = NULL
ps.aca

#Rarefaction ####
#rarefaction curve
source("../mp/scripts/my_functions.R")
rare_curve = calculate_rarefaction_curves(ps.aca,c('Observed', 'Shannon','Simpson'), c(1000,2000,5000,10000,20000))
summary(rare_curve)
rare_curve_summary = ddply(rare_curve, 
                               c('Depth', 'Sample', 'Measure'), summarise, 
                               Alpha_diversity_mean = mean(Alpha_diversity))
obs = rare_curve_summary[rare_curve_summary$Measure=="Observed",]
# Plot
ggplot(data = obs,mapping = aes(x = Depth,y = Alpha_diversity_mean,colour=Sample,group=Sample)) + 
  theme_bw() +
  scale_x_continuous(breaks = seq(0,20000,1000))+ geom_line() + theme_bw() +theme(legend.position = "none") + 
  #labs(title = "Rarefaction curve") +
  labs(x = "\nNumber of sequences", y = "Observed number of ASVs\n") + geom_point(size = 1) +
  theme(axis.text.x = element_text(size = 14, angle = 90),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 18, face = "bold"),
        title = element_text(size = 14, face = "bold")) 

#% 5,000 cutoff #### 
#we rarefy at the lowest cutoff which is the one from phyllosphere samples 
ps.rare5 = rarefy_even_depth(ps.aca, sample.size = 5000, rngseed = 9306, trimOTUs = TRUE, replace = TRUE)
ps.rare5
rarecurve(otu_table(ps.rare5), step=100,label=FALSE, col = "darkred",
          main = "Rarefy to 5,000 reads per sample")
nsamples(ps.aca) - nsamples(ps.rare5)
ntaxa(ps.aca) - ntaxa(ps.rare5)
#lost samples and ASVs after data denoising and rarefaction at 5,000
100-(nsamples(ps.rare5)/nsamples(ps.aca.neg_ctl))*100
100-(ntaxa(ps.rare5)/ntaxa(ps.aca.neg_ctl))*100 

#Subset samples 
#control (not treated) ####
#Subset control (non-treated) samples 
ps.ctl = subset_samples(ps.rare5, sample_data(ps.rare5)$neonic == "N") 
ps.ctl = prune_taxa(taxa_sums(ps.ctl)>0,ps.ctl)
sample_data(ps.ctl)$neonic = NULL
ps.ctl
saveRDS(ps.ctl, "16S_aca_ctl5000.rds")

#phyllosphere ####
#subset phyllosphere
#not rarefied
ps.phyl.aca = subset_samples(ps.aca, sample_data(ps.aca)$habitat == "leaf") 
ps.phyl.aca = prune_taxa(taxa_sums(ps.phyl.aca)>0,ps.phyl.aca)
sample_data(ps.phyl.aca)$habitat = NULL
ps.phyl.aca
saveRDS(ps.phyl.aca, "16S_aca_phyl.rds")
#rarefied at 5,000
ps.phyl = subset_samples(ps.rare5, sample_data(ps.rare5)$habitat == "leaf") 
ps.phyl = prune_taxa(taxa_sums(ps.phyl)>0,ps.phyl)
sample_data(ps.phyl)$habitat = NULL
ps.phyl
saveRDS(ps.phyl, "16S_aca_phyl5000.rds")

#% 10,000 cutoff #### 
#we rarefy at the lowest cutoff which is the one from soil samples 
ps.rare10 = rarefy_even_depth(ps.aca, sample.size = 10000, rngseed = 8306, trimOTUs = TRUE, replace = TRUE)
ps.rare10
rarecurve(otu_table(ps.rare10), step=100,label=FALSE, col = "darkred",
          main = "Rarefy to 10,000 reads per sample")
nsamples(ps.aca) - nsamples(ps.rare10)
ntaxa(ps.aca) - ntaxa(ps.rare10)
#lost samples and ASVs after data denoising and rarefaction at 10,000
100-(nsamples(ps.rare10)/nsamples(ps.aca.neg_ctl))*100
100-(ntaxa(ps.rare10)/ntaxa(ps.aca.neg_ctl))*100 
#soil ####
#subset soil
#not rarefied
ps.soil.aca = subset_samples(ps.aca, sample_data(ps.aca)$habitat == "soil") 
ps.soil.aca = prune_taxa(taxa_sums(ps.soil.aca)>0,ps.soil.aca)
sample_data(ps.soil.aca)$habitat = NULL
ps.soil.aca
saveRDS(ps.soil.aca, "16S_aca_soil.rds")
#rarefied at 10,000
ps.soil = subset_samples(ps.rare10, sample_data(ps.rare10)$habitat == "soil") 
ps.soil = prune_taxa(taxa_sums(ps.soil)>0,ps.soil)
sample_data(ps.soil)$habitat = NULL
ps.soil
saveRDS(ps.soil, "16S_aca_soil10000.rds")
