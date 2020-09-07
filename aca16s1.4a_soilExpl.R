############################################################################
# Explanatory analysis of soil samples                                     #
# Studying the effects of neonicotinoids on soil bacterial communities     #
# Data: Miseq-16S - L'Acadie (ACA)                                         #
# Mona Parizadeh - 2019-2020                                               #
############################################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq") #‘1.27.6’
library(vegan); packageVersion("vegan") #‘2.5.6’
library(ggplot2); packageVersion("ggplot2") #‘3.3.0’

# Import data #### 
setwd("../mp/aca_16s/files/")
ps = readRDS("16S_aca_soil10000.rds")
ps

#stats
asv.rich = estimate_richness(ps, measures = "Observed") #ASV per sample (richness)
summary(asv.rich)
sd(asv.rich$Observed, na.rm=TRUE) /  
  sqrt(length(asv.rich$Observed[!is.na(asv.rich$Observed)])) #SE
#Subset hosts ####
ps.sy = subset_samples(ps, sample_data(ps)$host == "soy")
ps.sy = prune_taxa(taxa_sums(ps.sy)>0, ps.sy)
ps.sy 
ps.cr = subset_samples(ps, sample_data(ps)$host == "corn")
ps.cr = prune_taxa(taxa_sums(ps.cr)>0, ps.cr)
ps.cr

#PERMANOVA ####
#make dataframe
df = as(sample_data(ps), "data.frame")
#bray-curtis distance
dis = phyloseq::distance(ps,  method = "bray")
set.seed(31108)
adns = adonis2(dis ~ host*year*month*neonic, df) #distance = bray
adns

#Ordination ####
pcoa = ordinate(ps, method = "MDS", k = 2, try = 100, distance = "bray")

#define variables as factors
hst = get_variable(ps, "host")
sample_data(ps)$host = factor(hst)
neo = get_variable(ps, "neonic")
sample_data(ps)$neonic = factor(neo)
yr = get_variable(ps, "year")
sample_data(ps)$year = factor(yr)
#group
group.hst.neo.yr = paste(hst, neo, yr, sep = "")
#order as factor
sample_data(ps)$month = factor(sample_data(ps)$month,levels=c("July","Aug","Sep"))
group.hst.neo.yr.fac = factor(group.hst.neo.yr,levels=c("soyN2016","soyY2016","cornN2017",
                                                        "cornY2017","soyN2018","soyY2018"))
#define labels and colors 
hst.neo.labels = c(soyN2016="Control Soybean 2016", soyY2016="Neonicotinoid-treated Soybean 2016",
                   cornN2017="Control Corn 2017", cornY2017="Neonicotinoid-treated Corn 2017",
                   soyN2018="Control Soybean 2018", soyY2018="Neonicotinoid-treated Soybean 2018")
soil.hst.neo.colors = c("chartreuse4", "darkred", "chartreuse4", "darkred","chartreuse4", "darkred")
  
#ordinate
pcoa1 = plot_ordination(ps, pcoa) 
pcoa1$layers
pcoa1$layers = pcoa1$layers[-1] #remove the original points to add the desired colors and shapes
#plot - months seperated
pcoa1 + 
  theme_bw() +
  xlim(-0.6, 0.6) +
  ylim(-0.5, 0.5) +
  #group by neonic
  stat_ellipse(aes(fill = neonic, group=neonic), geom = "polygon", 
               level = 0.95, linetype = 0, alpha = 0.4, show.legend=TRUE) +
  scale_fill_manual(name="Treatement ellipses",values = c("cornflowerblue","darkgoldenrod2"),
                    labels = c(N="Control", Y="Neonicotinoid-treated")) + 
  geom_point(aes(color = group.hst.neo.yr.fac, shape = group.hst.neo.yr.fac), size=3, alpha=0.75) + 
  #facet_grid(year~month) +
  facet_wrap(~ month, ncol = 2,
             labeller=labeller(month = c(May="May",July="July",
                                         Aug="August",Sep="September"))) +
  #change facet font size
  theme(strip.text.x = element_text(size=16, face="bold")) +
  #change legend labels 
  scale_color_manual(name = "Soil: Treatment, Host & Year",
                     values = soil.hst.neo.colors, 
                     labels = hst.neo.labels) +
  scale_shape_manual(name = "Soil: Treatment, Host & Year", 
                     values = c(19,19,17,17,15,15) ,
                     labels = hst.neo.labels) +
  theme(legend.title = element_text(size=16, face="bold"),
        legend.text=element_text(size=16),
        legend.position = "right",
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text( size = 16, face = "bold")) +
  labs(tag = "C)") + theme(plot.tag = element_text(size = 18, face = "bold")) 

#HOSTS ####

#% Soybean 
#PERMANOVA
#make dataframe
df.sy = as(sample_data(ps.sy), "data.frame")
#bray-curtis distance
dis.sy = phyloseq::distance(ps.sy,  method = "bray")
set.seed(311081)
adns.sy = adonis2(dis.sy ~ year*month*neonic, df.sy) #distance = bray
adns.sy

#% Corn ####
#PERMANOVA ####
#make dataframe
df.cr = as(sample_data(ps.cr), "data.frame")
#bray-curtis distance
dis.cr = phyloseq::distance(ps.cr,  method = "bray")
set.seed(311082)
adns.cr = adonis2(dis.cr ~ month*neonic, df.cr) #distance = bray
adns.cr

#Richness ####
shn.rich = cbind(estimate_richness(ps,measures = 'shannon'),
                 sample_data(ps))
summary(shn.rich)
sd(shn.rich$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich$Shannon[!is.na(shn.rich$Shannon)])) #SE
# Wilcoxon rank-sum test (Mann-Whitney)
library(ggpubr); packageVersion("ggpubr") #‘0.3.0’
compare_means(Shannon ~ neonic, shn.rich, method = "wilcox.test", paired = FALSE, p.adjust.method = "holm") 

#Shannon richness - treatment ####
#control
ps.ctl = subset_samples(ps, sample_data(ps)$neonic == "N")
ps.ctl = prune_taxa(taxa_sums(ps.ctl)>0, ps.ctl)
shn.rich.ctl = cbind(estimate_richness(ps.ctl,measures = 'shannon'),
                     sample_data(ps.ctl))
summary(shn.rich.ctl)
sd(shn.rich.ctl$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.ctl$Shannon[!is.na(shn.rich.ctl$Shannon)])) #SE
#neonic-treated
ps.neo = subset_samples(ps, sample_data(ps)$neonic == "Y")
ps.neo = prune_taxa(taxa_sums(ps.neo)>0, ps.neo)
shn.rich.neo = cbind(estimate_richness(ps.neo,measures = 'shannon'),
                     sample_data(ps.neo))
summary(shn.rich.neo)
sd(shn.rich.neo$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.neo$Shannon[!is.na(shn.rich.neo$Shannon)])) #SE

#Shannon richness - hosts individually in interaction w/ treatment ####
#soybean
shn.rich.sy = cbind(estimate_richness(ps.sy,measures = 'shannon'),
                    sample_data(ps.sy))
compare_means(Shannon ~ neonic, shn.rich.sy, method = "wilcox.test", paired = FALSE, p.adjust.method = "holm") 
##soybean control
ps.sy.ctl = subset_samples(ps.sy, sample_data(ps.sy)$neonic == "N")
ps.sy.ctl = prune_taxa(taxa_sums(ps.sy.ctl)>0, ps.sy.ctl)
shn.rich.sy.ctl = cbind(estimate_richness(ps.sy.ctl,measures = 'shannon'),
                        sample_data(ps.sy.ctl))
summary(shn.rich.sy.ctl)
sd(shn.rich.sy.ctl$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.sy.ctl$Shannon[!is.na(shn.rich.sy.ctl$Shannon)])) #SE
##soybean neonic-treated
ps.sy.neo = subset_samples(ps.sy, sample_data(ps.sy)$neonic == "Y")
ps.sy.neo = prune_taxa(taxa_sums(ps.sy.neo)>0, ps.sy.neo)
shn.rich.sy.neo = cbind(estimate_richness(ps.sy.neo,measures = 'shannon'),
                        sample_data(ps.sy.neo))
summary(shn.rich.sy.neo)
sd(shn.rich.sy.neo$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.sy.neo$Shannon[!is.na(shn.rich.sy.neo$Shannon)])) #SE

#corn
shn.rich.cr = cbind(estimate_richness(ps.cr,measures = 'shannon'),
                    sample_data(ps.cr))
compare_means(Shannon ~ neonic, shn.rich.cr, method = "wilcox.test", paired = FALSE, p.adjust.method = "holm") 
#corn control
ps.cr.ctl = subset_samples(ps.cr, sample_data(ps.cr)$neonic == "N")
ps.cr.ctl = prune_taxa(taxa_sums(ps.cr.ctl)>0, ps.cr.ctl)
shn.rich.cr.ctl = cbind(estimate_richness(ps.cr.ctl,measures = 'shannon'),
                        sample_data(ps.cr.ctl))
summary(shn.rich.cr.ctl)
sd(shn.rich.cr.ctl$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.cr.ctl$Shannon[!is.na(shn.rich.cr.ctl$Shannon)])) #SE
#corn neonic-treated
ps.cr.neo = subset_samples(ps.cr, sample_data(ps.cr)$neonic == "Y")
ps.cr.neo = prune_taxa(taxa_sums(ps.cr.neo)>0, ps.cr.neo)
shn.rich.cr.neo = cbind(estimate_richness(ps.cr.neo,measures = 'shannon'),
                        sample_data(ps.cr.neo))
summary(shn.rich.cr.neo)
sd(shn.rich.cr.neo$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.cr.neo$Shannon[!is.na(shn.rich.cr.neo$Shannon)])) #SE

#Shannon richness - months in interaction w/ treatment ####
compare_means(Shannon ~ neonic, shn.rich, method = "wilcox.test", paired = FALSE, p.adjust.method = "holm",
              group.by = "month") 
#July - control
shn.j.n = shn.rich %>% 
  filter(month == "July" & neonic == "N") 
summary(shn.j.n)
sd(shn.j.n$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.j.n$Shannon[!is.na(shn.j.n$Shannon)])) #SE
#July - neonic-treated
shn.j.y = shn.rich %>% 
  filter(month == "July" & neonic == "Y") 
summary(shn.j.y)
sd(shn.j.y$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.j.y$Shannon[!is.na(shn.j.y$Shannon)])) #SE
#August - control
shn.a.n = shn.rich %>% 
  filter(month == "Aug" & neonic == "N") 
summary(shn.a.n)
sd(shn.a.n$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.a.n$Shannon[!is.na(shn.a.n$Shannon)])) #SE
#August - neonic-treated
shn.a.y = shn.rich %>% 
  filter(month == "Aug" & neonic == "Y") 
summary(shn.a.y)
sd(shn.a.y$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.a.y$Shannon[!is.na(shn.a.y$Shannon)])) #SE
#September - control
shn.s.n = shn.rich %>% 
  filter(month == "Sep" & neonic == "N") 
summary(shn.s.n)
sd(shn.s.n$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.s.n$Shannon[!is.na(shn.s.n$Shannon)])) #SE
#September - neonic-treated
shn.s.y = shn.rich %>% 
  filter(month == "Sep" & neonic == "Y") 
summary(shn.s.y)
sd(shn.s.y$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.s.y$Shannon[!is.na(shn.s.y$Shannon)])) #SE

#Shannon richness - years in interaction w/ treatment ####
compare_means(Shannon ~ neonic, shn.rich, method = "wilcox.test", paired = FALSE, p.adjust.method = "holm",
              group.by = "year") 
#2016 - control
shn.16.n = shn.rich %>% 
  filter(year == "2016" & neonic == "N") 
summary(shn.16.n)
sd(shn.16.n$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.16.n$Shannon[!is.na(shn.16.n$Shannon)])) #SE
#2016 - neonic-treated
shn.16.y = shn.rich %>% 
  filter(year == "2016" & neonic == "Y") 
summary(shn.16.y)
sd(shn.16.y$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.16.y$Shannon[!is.na(shn.16.y$Shannon)])) #SE
#2017 - control
shn.17.n = shn.rich %>% 
  filter(year == "2017" & neonic == "N") 
summary(shn.17.n)
sd(shn.17.n$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.17.n$Shannon[!is.na(shn.17.n$Shannon)])) #SE
#2017 - neonic-treated
shn.17.y = shn.rich %>% 
  filter(year == "2017" & neonic == "Y") 
summary(shn.17.y)
sd(shn.17.y$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.17.y$Shannon[!is.na(shn.17.y$Shannon)])) #SE
#2018 - control
shn.18.n = shn.rich %>% 
  filter(year == "2018" & neonic == "N") 
summary(shn.18.n)
sd(shn.18.n$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.18.n$Shannon[!is.na(shn.18.n$Shannon)])) #SE
#2018 - neonic-treated
shn.18.y = shn.rich %>% 
  filter(year == "2018" & neonic == "Y") 
summary(shn.18.y)
sd(shn.18.y$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.18.y$Shannon[!is.na(shn.18.y$Shannon)])) #SE

#Linear model ####
shn.lm = lm(Shannon ~ neonic * month * year, data=shn.rich)
summary(shn.lm)
anova(shn.lm)
