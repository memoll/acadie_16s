#######################################################################################
# Explanatory analysis of the phyllosphere and soil non-treated (control) samples     #
# Studying the bacteria communities differences between phyllosphere and soil         #
# Data: Miseq-16S - L'Acadie (ACA)                                                    #
# Mona Parizadeh - 2019-2020                                                          #
#######################################################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq") #‘1.27.6’
library(vegan); packageVersion("vegan") #‘2.5.6’
library(ggplot2); packageVersion("ggplot2") #‘3.3.0’

# Import data (rarefied) #### 
setwd("../mp/aca_16s/files/")
ps = readRDS("16S_aca_ctl5000.rds") 
ps

#Explore data
subset_samples(ps, sample_data(ps)$habitat == "leaf" & sample_data(ps)$host == "soy")
subset_samples(ps, sample_data(ps)$habitat == "leaf" & sample_data(ps)$host == "corn")
subset_samples(ps, sample_data(ps)$habitat == "soil" & sample_data(ps)$host == "soy")
subset_samples(ps, sample_data(ps)$habitat == "soil" & sample_data(ps)$host == "corn")
#stats
asv.rich = estimate_richness(ps, measures = "Observed") #ASV per sample (richness)
summary(asv.rich)
sd(asv.rich$Observed, na.rm=TRUE) /  
  sqrt(length(asv.rich$Observed[!is.na(asv.rich$Observed)])) #SE 

#PERMANOVA ####
#make dataframe
df = as(sample_data(ps), "data.frame")
#bray-curtis distance
dis = phyloseq::distance(ps,  method = "bray")
#permanova
set.seed(411166)
adns.t = adonis2(dis ~ habitat*host*month/year, df) 
adns.t 

# Homogeneity of dispersion test ####
#habitat
beta.hbt = betadisper(dis, df$habitat)
beta.hbt
#ANOVA
set.seed(1233)
permutest(beta.hbt)
#host
beta.hst = betadisper(dis, df$host)
beta.hst
set.seed(1234)
permutest(beta.hst)
#phyl-host
beta.phyl.hst = betadisper(dis.phyl, df.phyl$host)
beta.phyl.hst
set.seed(12341)
permutest(beta.phyl.hst)
#soil-host
beta.soil.hst = betadisper(dis.soil, df.soil$host)
beta.soil.hst
set.seed(12342)
permutest(beta.soil.hst)

#Richness ####
shn.rich = cbind(estimate_richness(ps,measures = 'shannon'),
                 sample_data(ps))
ggplot(shn.rich, aes(x = habitat, y = Shannon, color=habitat)) +  
  geom_boxplot()
ggplot(shn.rich, aes(x = habitat, y = Shannon, color=host)) +  
  geom_boxplot()
#Wilcoxon rank-sum test (Compare the diversity between habitats)
library(ggpubr); packageVersion("ggpubr") #‘0.3.0’
compare_means(Shannon ~ habitat, data = shn.rich, method = "wilcox.test", p.adjust.method = "holm")
compare_means(Shannon ~ host, data = shn.rich, method = "wilcox.test", p.adjust.method = "holm")

#Ordination ####
pcoa = ordinate(ps, method = "MDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps, pcoa, color = "habitat", shape = "host") + geom_point() + ggtitle("PCoA") +
  geom_text(aes(label = sampleid), check_overlap = TRUE, size = 3, nudge_y = -0.01) +
  scale_color_manual(name = "Habitat", values = c("chartreuse4", "darkred"))
pcoa1 = plot_ordination(ps,pcoa)
#define variables as factors
hst = get_variable(ps, "host")
sample_data(ps)$host = factor(hst)
hbt = get_variable(ps, "habitat")
sample_data(ps)$habitat = factor(hbt)
#% variables
group.hbt.hst = paste(hbt,hst, sep = "")
#define labels and colors
hbt.hst.colors = c("darkgoldenrod2","darkgreen", "darkgoldenrod2","darkgreen")

#Empty points of the PCoA replace them with the desired shapes 
pcoa1$layers
pcoa1$layers = pcoa1$layers[-1] #remove the original points to add the desired colors and shapes
pcoa1 +
  theme_bw() +
  #group by habitat and host
  stat_ellipse(aes(fill = group.hbt.hst, group = group.hbt.hst), geom = "polygon", 
               level = 0.99, linetype = 0, alpha = 0.2, show.legend=FALSE) +
  scale_fill_manual(name="Host & Habitat ellipses",values = hbt.hst.colors) +
  geom_point(size = 3, alpha = 0.5, aes(color = group.hbt.hst, shape = group.hbt.hst)) +
  scale_shape_manual(name = "Host & Habitat",
                     labels = c("Corn Phyllosphere", "Soybean Phyllosphere","Corn Soil", "Soybean soil"),
                     values = c(19, 17, 19, 17)) +
  scale_colour_manual(name = "Host & Habitat",
                      labels = c("Corn Phyllosphere", "Soybean Phyllosphere","Corn Soil", "Soybean soil"),
                      values = c("chartreuse4","chartreuse4", "darkred", "darkred")) +   
  theme(axis.text.x = element_text(size = 12), #axis 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold"), #main axis label
        legend.title = element_text(size=12, face="bold"),legend.text=element_text(size=12),
        legend.position = "right") +
  guides(fill = guide_legend(override.aes=list(shape=15))) 

# Fitness ####
comm = otu_table(ps)
taxo = tax_table(ps)
library(seqtools); packageVersion("seqtools") #‘0.1.0’
#family
comm.fam = taxocomm(comm, taxo, "Family")
comm.fam.ra = decostand(comm.fam, method="total") #relative abundance
fam01 = sort(apply(comm.fam.ra[,apply(comm.fam.ra,2,mean)>0.01], 2, mean),TRUE)
fam01.ra = comm.fam.ra[,which(colnames(comm.fam.ra) %in% names(fam01))]
#fit vectors on pcoa
env.fam = envfit(pcoa$vectors, fam01.ra)
env.fam
sort(env.fam$vectors$r, decreasing = TRUE)
#dataframe
env.fam.df = as.data.frame(env.fam$vectors$arrows*sqrt(env.fam$vectors$r))
env.fam.df$family = rownames(env.fam.df)

#plot (bacterial families on PCoA)
p.env.fam = plot_ordination(ps, pcoa) 
#Empty points of the PCoA replace them with the desired shapes 
p.env.fam$layers
p.env.fam$layers = p.env.fam$layers[-1] #remove the original points to add the desired colors and shapes
library(ggrepel); packageVersion("ggrepel") #‘0.9.0’
p.env.fam + 
  #ggtitle("PCoA - Bacterial families who are driving the diversity pattern") +
  theme_bw() +
  geom_point(size = 3, alpha = 0.5, aes(color = group.hbt.hst, shape = group.hbt.hst)) +
  scale_shape_manual(name = "Host & Habitat",
                     labels = c("Corn Phyllosphere", "Soybean Phyllosphere","Corn Soil", "Soybean soil"),
                     values = c(19, 17, 19, 17)) +
  scale_colour_manual(name = "Host & Habitat",
                      labels = c("Corn Phyllosphere", "Soybean Phyllosphere","Corn Soil", "Soybean soil"),
                      values = c("chartreuse4","chartreuse4", "darkred", "darkred")) +   
  scale_fill_manual(name = "Host", values=c("darkgoldenrod2","darkgreen"),
                    labels = c(corn="Corn", soy="Soybean")) +
  geom_segment(data=env.fam.df,aes(x=0,xend=0.4*Axis.1,y=0,yend=0.4*Axis.2), size=0.2,
               arrow = arrow(length = unit(0.01, "npc")),colour="black",inherit.aes=FALSE) +
  geom_label_repel(data=env.fam.df,aes(x=0.4*Axis.1,y=0.4*Axis.2,label=family),size=4, fontface="bold",
                   color = "antiquewhite4", label.padding = unit(0.2, "lines"), label.size = 0.2, show.legend = FALSE) +
  theme(axis.text.x = element_text(size = 10), #axis 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10, face = "bold"), #main axis label
        legend.title = element_text(size=12, face="bold"),legend.text=element_text(size=12),
        legend.position = "none") +
  guides(fill = guide_legend(override.aes=list(shape=15))) 

#Phyllosphere #### 
ps.phyl = subset_samples(ps, sample_data(ps)$habitat == "leaf") 
ps.phyl = prune_taxa(taxa_sums(ps.phyl)>0, ps.phyl)
ps.phyl
#make dataframe
df.phyl = as(sample_data(ps.phyl), "data.frame")  
#bray-curtis distance
dis.phyl = phyloseq::distance(ps.phyl,  method = "bray")
#PERMANOVA 
set.seed(4112)
adns.phyl = adonis2(dis.phyl ~ host*month/year, df.phyl)
adns.phyl
#Shannon richness
shn.rich.phyl = cbind(estimate_richness(ps.phyl,measures = 'shannon'),
                 sample_data(ps.phyl))
summary(shn.rich.phyl)
sd(shn.rich.phyl$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.phyl$Shannon[!is.na(shn.rich.phyl$Shannon)])) #SE
#Wilcoxon test
compare_means(Shannon ~ year, data = shn.rich.phyl, method = "wilcox.test", p.adjust.method = "holm")
#%Soybean phyllosphere #### 
ps.phyl.sy = subset_samples(ps.phyl, sample_data(ps.phyl)$host == "soy") 
ps.phyl.sy = prune_taxa(taxa_sums(ps.phyl.sy)>0, ps.phyl.sy)
ps.phyl.sy
#shannon richness
shn.rich.phyl.sy = estimate_richness(ps.phyl.sy, measures = 'shannon')
shn.rich.phyl.sy = cbind(estimate_richness(ps.phyl.sy,measures = 'shannon'),
                         sample_data(ps.phyl.sy))
summary(shn.rich.phyl.sy)
sd(shn.rich.phyl.sy$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.phyl.sy$Shannon[!is.na(shn.rich.phyl.sy$Shannon)])) #SE
#%Corn phyllosphere #### 
ps.phyl.cr = subset_samples(ps.phyl, sample_data(ps.phyl)$host == "corn") 
ps.phyl.cr = prune_taxa(taxa_sums(ps.phyl.cr)>0, ps.phyl.cr)
ps.phyl.cr
#shannon richness
shn.rich.phyl.cr = estimate_richness(ps.phyl.cr, measures = 'shannon')
shn.rich.phyl.cr = cbind(estimate_richness(ps.phyl.cr,measures = 'shannon'),
                         sample_data(ps.phyl.cr))
summary(shn.rich.phyl.cr)
sd(shn.rich.phyl.cr$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.phyl.cr$Shannon[!is.na(shn.rich.phyl.cr$Shannon)])) #SE

#Soil ####
ps.soil = subset_samples(ps, sample_data(ps)$habitat == "soil") 
ps.soil = prune_taxa(taxa_sums(ps.soil)>0, ps.soil)
ps.soil
#make dataframe
df.soil = as(sample_data(ps.soil), "data.frame")
#bray-curtis distance
dis.soil = phyloseq::distance(ps.soil,  method = "bray")
#PERMANOVA 
set.seed(4112)
adns.soil = adonis2(dis.soil ~ host*month/year, df.soil)
adns.soil
#Shannon richness
shn.rich.soil = cbind(estimate_richness(ps.soil,measures = 'shannon'),
                      sample_data(ps.soil))
summary(shn.rich.soil)
sd(shn.rich.soil$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.soil$Shannon[!is.na(shn.rich.soil$Shannon)])) #SE
#Wilcoxon test
compare_means(Shannon ~ year, data = shn.rich.soil, method = "wilcox.test", p.adjust.method = "holm")
#%Soybean soil #### 
ps.soil.sy = subset_samples(ps.soil, sample_data(ps.soil)$host == "soy") 
ps.soil.sy = prune_taxa(taxa_sums(ps.soil.sy)>0, ps.soil.sy)
ps.soil.sy
#shannon richness
shn.rich.soil.sy = estimate_richness(ps.soil.sy, measures = 'shannon')
shn.rich.soil.sy = cbind(estimate_richness(ps.soil.sy,measures = 'shannon'),
                         sample_data(ps.soil.sy))
summary(shn.rich.soil.sy)
sd(shn.rich.soil.sy$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.soil.sy$Shannon[!is.na(shn.rich.soil.sy$Shannon)])) #SE
#%Corn Soil #### 
ps.soil.cr = subset_samples(ps.soil, sample_data(ps.soil)$host == "corn") 
ps.soil.cr = prune_taxa(taxa_sums(ps.soil.cr)>0, ps.soil.cr)
ps.soil.cr
#shannon richness
shn.rich.soil.cr = estimate_richness(ps.soil.cr, measures = 'shannon')
shn.rich.soil.cr = cbind(estimate_richness(ps.soil.cr,measures = 'shannon'),
                         sample_data(ps.soil.cr))
summary(shn.rich.soil.cr)
sd(shn.rich.soil.cr$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.soil.cr$Shannon[!is.na(shn.rich.soil.cr$Shannon)])) #SE

#Phyllosphere-time ####
#% year interactions w/ phyllosphere - shannon richness ####
#2016
ps.phyl16 = subset_samples(ps.phyl, sample_data(ps.phyl)$year == "2016") 
ps.phyl16 = prune_taxa(taxa_sums(ps.phyl16)>0, ps.phyl16)
ps.phyl16
shn.rich.phyl16 = cbind(estimate_richness(ps.phyl16,measures = 'shannon'),
                      sample_data(ps.phyl16))
summary(shn.rich.phyl16)
sd(shn.rich.phyl16$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.phyl16$Shannon[!is.na(shn.rich.phyl16$Shannon)])) #SE
#2017
ps.phyl17 = subset_samples(ps.phyl, sample_data(ps.phyl)$year == "2017") 
ps.phyl17 = prune_taxa(taxa_sums(ps.phyl17)>0, ps.phyl17)
ps.phyl17
shn.rich.phyl17 = cbind(estimate_richness(ps.phyl17,measures = 'shannon'),
                        sample_data(ps.phyl17))
summary(shn.rich.phyl17)
sd(shn.rich.phyl17$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.phyl17$Shannon[!is.na(shn.rich.phyl17$Shannon)])) #SE
#2018
ps.phyl18 = subset_samples(ps.phyl, sample_data(ps.phyl)$year == "2018") 
ps.phyl18 = prune_taxa(taxa_sums(ps.phyl18)>0, ps.phyl18)
ps.phyl18
shn.rich.phyl18 = cbind(estimate_richness(ps.phyl18,measures = 'shannon'),
                        sample_data(ps.phyl18))
summary(shn.rich.phyl18)
sd(shn.rich.phyl18$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.phyl18$Shannon[!is.na(shn.rich.phyl18$Shannon)])) #SE
#Wilcoxon test
compare_means(Shannon ~ year, data = shn.rich.phyl,method = "wilcox.test", p.adjust.method = "holm")

#% month interactions w/ phyllosphere - shannon richness ####
#July
ps.phyl.j = subset_samples(ps.phyl, sample_data(ps.phyl)$month == "July") 
ps.phyl.j = prune_taxa(taxa_sums(ps.phyl.j)>0, ps.phyl.j)
ps.phyl.j
shn.rich.phyl.j = cbind(estimate_richness(ps.phyl.j,measures = 'shannon'),
                        sample_data(ps.phyl.j))
summary(shn.rich.phyl.j)
sd(shn.rich.phyl.j$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.phyl.j$Shannon[!is.na(shn.rich.phyl.j$Shannon)])) #SE
#August
ps.phyl.a = subset_samples(ps.phyl, sample_data(ps.phyl)$month == "Aug") 
ps.phyl.a = prune_taxa(taxa_sums(ps.phyl.a)>0, ps.phyl.a)
ps.phyl.a
shn.rich.phyl.a = cbind(estimate_richness(ps.phyl.a,measures = 'shannon'),
                        sample_data(ps.phyl.a))
summary(shn.rich.phyl.a)
sd(shn.rich.phyl.a$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.phyl.a$Shannon[!is.na(shn.rich.phyl.a$Shannon)])) #SE
#September
ps.phyl.s = subset_samples(ps.phyl, sample_data(ps.phyl)$month == "Sep") 
ps.phyl.s = prune_taxa(taxa_sums(ps.phyl.s)>0, ps.phyl.s)
ps.phyl.s
shn.rich.phyl.s = cbind(estimate_richness(ps.phyl.s,measures = 'shannon'),
                        sample_data(ps.phyl.s))
summary(shn.rich.phyl.s)
sd(shn.rich.phyl.s$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.phyl.s$Shannon[!is.na(shn.rich.phyl.s$Shannon)])) #SE
#Wilcoxon test
compare_means(Shannon ~ month, data = shn.rich.phyl,method = "wilcox.test", p.adjust.method = "holm")

#Soil-time ####
#% year interactions w/ soil - shannon richness ####
#2016
ps.soil16 = subset_samples(ps.soil, sample_data(ps.soil)$year == "2016") 
ps.soil16 = prune_taxa(taxa_sums(ps.soil16)>0, ps.soil16)
ps.soil16
shn.rich.soil16 = cbind(estimate_richness(ps.soil16,measures = 'shannon'),
                        sample_data(ps.soil16))
summary(shn.rich.soil16)
sd(shn.rich.soil16$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.soil16$Shannon[!is.na(shn.rich.soil16$Shannon)])) #SE
#2017
ps.soil17 = subset_samples(ps.soil, sample_data(ps.soil)$year == "2017") 
ps.soil17 = prune_taxa(taxa_sums(ps.soil17)>0, ps.soil17)
ps.soil17
shn.rich.soil17 = cbind(estimate_richness(ps.soil17,measures = 'shannon'),
                        sample_data(ps.soil17))
summary(shn.rich.soil17)
sd(shn.rich.soil17$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.soil17$Shannon[!is.na(shn.rich.soil17$Shannon)])) #SE
#2018
ps.soil18 = subset_samples(ps.soil, sample_data(ps.soil)$year == "2018") 
ps.soil18 = prune_taxa(taxa_sums(ps.soil18)>0, ps.soil18)
ps.soil18
shn.rich.soil18 = cbind(estimate_richness(ps.soil18,measures = 'shannon'),
                        sample_data(ps.soil18))
summary(shn.rich.soil18)
sd(shn.rich.soil18$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.soil18$Shannon[!is.na(shn.rich.soil18$Shannon)])) #SE
#Wilcoxon test
compare_means(Shannon ~ year, data = shn.rich.soil,method = "wilcox.test", p.adjust.method = "holm")

#% month interactions w/ soil - shannon richness ####
#July
ps.soil.j = subset_samples(ps.soil, sample_data(ps.soil)$month == "July") 
ps.soil.j = prune_taxa(taxa_sums(ps.soil.j)>0, ps.soil.j)
ps.soil.j
shn.rich.soil.j = cbind(estimate_richness(ps.soil.j,measures = 'shannon'),
                        sample_data(ps.soil.j))
summary(shn.rich.soil.j)
sd(shn.rich.soil.j$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.soil.j$Shannon[!is.na(shn.rich.soil.j$Shannon)])) #SE
#August
ps.soil.a = subset_samples(ps.soil, sample_data(ps.soil)$month == "Aug") 
ps.soil.a = prune_taxa(taxa_sums(ps.soil.a)>0, ps.soil.a)
ps.soil.a
shn.rich.soil.a = cbind(estimate_richness(ps.soil.a,measures = 'shannon'),
                        sample_data(ps.soil.a))
summary(shn.rich.soil.a)
sd(shn.rich.soil.a$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.soil.a$Shannon[!is.na(shn.rich.soil.a$Shannon)])) #SE
#September
ps.soil.s = subset_samples(ps.soil, sample_data(ps.soil)$month == "Sep") 
ps.soil.s = prune_taxa(taxa_sums(ps.soil.s)>0, ps.soil.s)
ps.soil.s
shn.rich.soil.s = cbind(estimate_richness(ps.soil.s,measures = 'shannon'),
                        sample_data(ps.soil.s))
summary(shn.rich.soil.s)
sd(shn.rich.soil.s$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.soil.s$Shannon[!is.na(shn.rich.soil.s$Shannon)])) #SE
#Wilcoxon test
compare_means(Shannon ~ month, data = shn.rich.soil,method = "wilcox.test", p.adjust.method = "holm")
