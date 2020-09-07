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

#unifrac:
#ordu = ordinate(ps, "PCoA", "unifrac", weighted=TRUE)
#plot_ordination(ps, ordu, color="neonic", shape="neonic") +
#  facet_grid(~month)+
#  stat_ellipse(aes(fill = neonic, group=neonic), geom = "polygon", 
#               level = 0.95, linetype = 0, alpha = 0.4, show.legend=TRUE)

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
#Fig. 2C ####
p.fig2c = pcoa1 + 
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
        legend.position = "none",
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text( size = 16, face = "bold")) +
  #guides(col = guide_legend(nrow = 6, title.position = "top")) +
  labs(tag = "C)") + theme(plot.tag = element_text(size = 18, face = "bold")) 

#with legend (to be cropped)
p.fig2c.leg = pcoa1 + 
  theme_bw() +
  #group by neonic
  stat_ellipse(aes(fill = neonic, group=neonic), geom = "polygon", 
               level = 0.95, linetype = 0, alpha = 0.4, show.legend=TRUE) +
  scale_fill_manual(name="Treatement ellipses",values = c("cornflowerblue","darkgoldenrod2"),
                    labels = c(N="Control", Y="Neonicotinoid-treated")) + 
  geom_point(aes(color = group.hst.neo.yr.fac, shape = group.hst.neo.yr.fac), size=3, alpha=0.75) + 
  #facet_grid(year~month) +
  facet_wrap(~ month, ncol = 3,
             labeller=labeller(month = c(May="May",July="July",
                                         Aug="August",Sep="September"))) +
  #change facet font size
  theme(strip.text.x = element_text(size=10, face="bold"),
        strip.text.y = element_text(size=10, face="bold")) +
  #change legend labels 
  scale_color_manual(name = "Soil: Treatment, Host & Year",
                     values = soil.hst.neo.colors, 
                     labels = hst.neo.labels) +
  scale_shape_manual(name = "Soil: Treatment, Host & Year", 
                     values = c(19,19,17,17,15,15) ,
                     labels = hst.neo.labels) +
  theme(legend.title = element_text(size=12, face="bold"),
        legend.text=element_text(size=12),
        legend.position = "right",
        #legend.box = "horizontal",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text( size = 10, face = "bold")) +
  guides(col = guide_legend(nrow = 6, title.position = "top"))

#ggsave("/data/users/mona/miseq_16S/Mona_16S_all/article1/graphs/ParizadehM_Fig2C.pdf", p.fig2c,dpi = 300, width = 340, units = "mm")
library(cowplot)
save_plot("/data/users/mona/miseq_16S/Mona_16S_all/article1/graphs/ParizadehM_Fig2C.pdf", p.fig2c, ncol = 2, nrow = 2)
save_plot("/data/users/mona/miseq_16S/Mona_16S_all/article1/graphs/ParizadehM_Fig2Cleg.pdf", p.fig2c.leg, ncol = 3, nrow = 2)

#Group
#% ellipses
group.hst.neo = paste(hst, neo, sep = "")
pcoa1 + 
  theme_bw() +
  #group by neonic
  stat_ellipse(aes(fill = group.hst.neo, group=group.hst.neo), geom = "polygon", 
               level = 0.95, linetype = 0, alpha = 0.2, show.legend=FALSE) +
  scale_fill_manual(name="neonic & host",values = c("chartreuse4", "darkred","chartreuse4", "darkred")) +
  geom_point(aes(color = group.hst.neo.yr.fac, shape = group.hst.neo.yr.fac), size=4, alpha=0.75) + 
  #change facet font size
  theme(strip.text.x = element_text(size=14, face="bold"),
        strip.text.y = element_text(size=14, face="bold")) +
  #change legend labels 
  scale_color_manual(name = "Soil: Treatment, Host & Year",
                     values = soil.hst.neo.colors, 
                     labels = hst.neo.labels) +
  scale_shape_manual(name = "Soil: Treatment, Host & Year", 
                     values = c(19,19,17,17,15,15) ,
                     labels = hst.neo.labels) +
  theme(legend.title = element_text(size=16, face="bold"),
        legend.text=element_text(size=16),
        axis.title = element_text( size = 14, face = "bold"),
        legend.position = "bottom") +
  guides(col = guide_legend(nrow = 6, title.position = "top")) 
  #labs(tag = "B)") + theme(plot.tag = element_text(size = 18, face = "bold"))

#hosts - PERMANOVA####
#% Soybean ####
#PERMANOVA ####
#make dataframe
df.sy = as(sample_data(ps.sy), "data.frame")
#bray-curtis distance
dis.sy = phyloseq::distance(ps.sy,  method = "bray")
set.seed(311081)
adns.sy = adonis2(dis.sy ~ year*month*neonic, df.sy) #distance = bray
adns.sy
#Ordination ####
pcoa.sy = ordinate(ps.sy, method = "MDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps.sy, pcoa.sy, color = "month", shape = "neonic") + geom_point() + ggtitle("PCoA-soybean") +
  geom_text(aes(label = sampleid), check_overlap = FALSE, size = 4, nudge_y = -0.05) + #nudge_x to seperate id from point
  geom_point(size=4)
#Group and define variables as factors
mnt.sy = get_variable(ps.sy, "month")
sample_data(ps.sy)$month = factor(mnt.sy)
neo.sy = get_variable(ps.sy, "neonic")
sample_data(ps.sy)$neonic = factor(neo.sy)
#group
group.mnt.neo.sy = paste(mnt.sy, neo.sy, sep = "")
#define as factor
group.mnt.neo.sy.fac = factor(group.mnt.neo.sy, levels = c("JulyN","JulyY","AugN","AugY","SepN","SepY"))

#define labels and colors
mnt.neo.labels = c(JulyN="July - Control", JulyY="July - Neonicotinoid-treated",
                   AugN="August - Control", AugY="August - Neonicotinoid-treated",
                   SepN="September - Control", SepY="September - Neonicotinoid-treated")
mnt.neo.colors = c("chartreuse4", "darkred","chartreuse4", "darkred","chartreuse4", "darkred")
#final ordination
p.sy = plot_ordination(ps.sy, pcoa.sy) 
#Empty points of the PCoA replace them with the desired shapes 
p.sy$layers
p.sy$layers = p.sy$layers[-1] #remove the original points to add the desired colors and shapes
p.syA = p.sy +
  theme_bw() +
  #xlim(-0.7, 0.7) +
  ylim(-0.7, 0.7) +
  stat_ellipse(aes(fill = neonic, group=neonic), geom = "polygon", 
               level = 0.95, linetype = 0, alpha = 0.2, show.legend=FALSE) +
  scale_fill_manual(name="Treatment ellipses",labels = c(N="Control",Y="Neonicotinoid-treated"),
                    values = mnt.neo.colors) +
  geom_point(aes(color = group.mnt.neo.sy.fac, shape = group.mnt.neo.sy.fac), size=4, alpha=0.75) + 
  facet_wrap(~ year, labeller=labeller(year = c("2016"="Soybean 2016", "2018"="Soybean 2018"))) +
  #change facet font size
  theme(strip.text.x = element_text(size=14, face="bold")) +
  #change legend labels 
  scale_shape_manual(name = "Month & Treatement", 
                     labels = mnt.neo.labels,
                     values = c(19, 19, 17, 17, 15, 15)) +
  scale_color_manual(name = "Month & Treatement", 
                     labels = mnt.neo.labels,
                     values = mnt.neo.colors) +
  theme(legend.position = "none",
        #legend.title = element_text(size=16, face="bold"),
        #legend.text=element_text(size=14),
        axis.title = element_text(size = 14, face = "bold"),
        legend.justification = "left") 


#% Corn ####
#PERMANOVA ####
#make dataframe
df.cr = as(sample_data(ps.cr), "data.frame")
#bray-curtis distance
dis.cr = phyloseq::distance(ps.cr,  method = "bray")
set.seed(311082)
adns.cr = adonis2(dis.cr ~ month*neonic, df.cr) #distance = bray
adns.cr

#Ordination ####
pcoa.cr = ordinate(ps.cr, method = "MDS", k = 2, try = 100, distance = "bray")
plot_ordination(ps.cr, pcoa.cr, color = "month", shape = "neonic") + geom_point() + ggtitle("PCoA-corn") +
  geom_text(aes(label = sampleid), check_overlap = FALSE, size = 4, nudge_y = -0.05) + #nudge_x to seperate id from point
  geom_point(size=4)
#Group and define variables as factors
mnt.cr = get_variable(ps.cr, "month")
sample_data(ps.cr)$month = factor(mnt.cr)
neo.cr = get_variable(ps.cr, "neonic")
sample_data(ps.cr)$neonic = factor(neo.cr)
#group
group.mnt.neo.cr = paste(mnt.cr, neo.cr, sep = "")
#define as factor
group.mnt.neo.cr.fac = factor(group.mnt.neo.cr, levels = c("JulyN","JulyY","AugN","AugY","SepN","SepY"))

#define labels and colors: same as for sybean (above)
#final ordination
p.cr = plot_ordination(ps.cr, pcoa.cr) 
#Empty points of the PCoA replace them with the desired shapes 
p.cr$layers
p.cr$layers = p.cr$layers[-1] #remove the original points to add the desired colors and shapes
p.crB = p.cr +
  theme_bw() +
  #xlim(-0.7, 0.7) +
  ylim(-0.7, 0.7) +
  stat_ellipse(aes(fill = neonic, group=neonic), geom = "polygon", 
               level = 0.95, linetype = 0, alpha = 0.2, show.legend=TRUE) +
  scale_fill_manual(name="Treatment ellipses",labels = c(N="Control",Y="Neonicotinoid-treated"),
                    values = mnt.neo.colors) +
  geom_point(aes(color = group.mnt.neo.cr.fac, shape = group.mnt.neo.cr.fac), size=4, alpha=0.75) + 
  #change facet font size
  facet_wrap(~ year, labeller=labeller(year = c("2017"="Corn 2017"))) +
  theme(strip.text.x = element_text(size=14, face="bold")) +
  #change legend labels 
  scale_shape_manual(name = "Month & Treatement", 
                     labels = mnt.neo.labels,
                     values = c(19, 19, 17, 17, 15, 15)) +
  scale_color_manual(name = "Month & Treatement", 
                     labels = mnt.neo.labels,
                     values = mnt.neo.colors) +
  theme(legend.title = element_text(size=16, face="bold"),
        legend.text=element_text(size=14),
        axis.title = element_text( size = 14, face = "bold"),
        #legend.justification = "left"
        legend.box = "horizontal") 
#patch ordinations
library(patchwork);packageVersion("patchwork") #‘1.0.0’
p.syA / p.crB + 
  plot_layout(nrow = 2, height = c(1,1), width = c(1, 1), guides = 'collect') +
  plot_annotation(tag_levels = 'A', tag_suffix = ')') &
  theme(plot.tag = element_text(size = 18, face="bold")) 

#%%%%%%%% ####
#Richness ####
shn.rich = cbind(estimate_richness(ps,measures = 'shannon'),
                 sample_data(ps))
summary(shn.rich)
sd(shn.rich$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich$Shannon[!is.na(shn.rich$Shannon)])) #SE
# Wilcoxon rank-sum test (Mann-Whitney)
library(ggpubr); packageVersion("ggpubr") #‘0.3.0’
compare_means(Shannon ~ neonic, shn.rich, method = "wilcox.test", paired = FALSE, p.adjust.method = "holm") 
ggplot(shn.rich, aes(x = neonic, y = Shannon, color=neonic)) +  
  geom_boxplot() +
  facet_wrap(~month)

#control
ps.ctl = subset_samples(ps, sample_data(ps)$neonic == "N")
ps.ctl = prune_taxa(taxa_sums(ps.ctl)>0, ps.ctl)
shn.rich.ctl = cbind(estimate_richness(ps.ctl,measures = 'shannon'),
                     sample_data(ps.ctl))
summary(shn.rich.ctl)
sd(shn.rich.ctl$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.ctl$Shannon[!is.na(shn.rich.ctl$Shannon)])) #SE
#neonic
ps.neo = subset_samples(ps, sample_data(ps)$neonic == "Y")
ps.neo = prune_taxa(taxa_sums(ps.neo)>0, ps.neo)
shn.rich.neo = cbind(estimate_richness(ps.neo,measures = 'shannon'),
                     sample_data(ps.neo))
summary(shn.rich.neo)
sd(shn.rich.neo$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.neo$Shannon[!is.na(shn.rich.neo$Shannon)])) #SE
#or (Tukey post-hoc)
library(agricolae); packageVersion("agricolae")
aov = aov(Shannon ~ neonic, shn.rich)
hsd = HSD.test(aov, "neonic", group = TRUE)
print(hsd)

#soybean
shn.rich.sy = cbind(estimate_richness(ps.sy,measures = 'shannon'),
                    sample_data(ps.sy))
compare_means(Shannon ~ neonic, shn.rich.sy, method = "wilcox.test", paired = FALSE, p.adjust.method = "holm") 
#control
ps.sy.ctl = subset_samples(ps.sy, sample_data(ps.sy)$neonic == "N")
ps.sy.ctl = prune_taxa(taxa_sums(ps.sy.ctl)>0, ps.sy.ctl)
shn.rich.sy.ctl = cbind(estimate_richness(ps.sy.ctl,measures = 'shannon'),
                        sample_data(ps.sy.ctl))
summary(shn.rich.sy.ctl)
sd(shn.rich.sy.ctl$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.sy.ctl$Shannon[!is.na(shn.rich.sy.ctl$Shannon)])) #SE
#neonic
ps.sy.neo = subset_samples(ps.sy, sample_data(ps.sy)$neonic == "Y")
ps.sy.neo = prune_taxa(taxa_sums(ps.sy.neo)>0, ps.sy.neo)
shn.rich.sy.neo = cbind(estimate_richness(ps.sy.neo,measures = 'shannon'),
                        sample_data(ps.sy.neo))
summary(shn.rich.sy.neo)
sd(shn.rich.sy.neo$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.sy.neo$Shannon[!is.na(shn.rich.sy.neo$Shannon)])) #SE
#or
aov.sy = aov(Shannon ~ neonic, shn.rich.sy)
hsd.sy = HSD.test(aov.sy, "neonic", group = TRUE)
print(hsd.sy)

#corn
shn.rich.cr = cbind(estimate_richness(ps.cr,measures = 'shannon'),
                    sample_data(ps.cr))
compare_means(Shannon ~ neonic, shn.rich.cr, method = "wilcox.test", paired = FALSE, p.adjust.method = "holm") 
#control
ps.cr.ctl = subset_samples(ps.cr, sample_data(ps.cr)$neonic == "N")
ps.cr.ctl = prune_taxa(taxa_sums(ps.cr.ctl)>0, ps.cr.ctl)
shn.rich.cr.ctl = cbind(estimate_richness(ps.cr.ctl,measures = 'shannon'),
                        sample_data(ps.cr.ctl))
summary(shn.rich.cr.ctl)
sd(shn.rich.cr.ctl$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.cr.ctl$Shannon[!is.na(shn.rich.cr.ctl$Shannon)])) #SE
#neonic
ps.cr.neo = subset_samples(ps.cr, sample_data(ps.cr)$neonic == "Y")
ps.cr.neo = prune_taxa(taxa_sums(ps.cr.neo)>0, ps.cr.neo)
shn.rich.cr.neo = cbind(estimate_richness(ps.cr.neo,measures = 'shannon'),
                        sample_data(ps.cr.neo))
summary(shn.rich.cr.neo)
sd(shn.rich.cr.neo$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.rich.cr.neo$Shannon[!is.na(shn.rich.cr.neo$Shannon)])) #SE
#or
aov.cr = aov(Shannon ~ neonic, shn.rich.cr)
hsd.cr = HSD.test(aov.cr, "neonic", group = TRUE)
print(hsd.cr)



#########################????
ps.ra = transform_sample_counts(ps, function(otu) otu/sum(otu)) 
ps.ra %>% 
  subset_taxa(Genus=="Gaiella") %>%
  plot_bar(fill="month") +
  facet_grid(~year+neonic, scales = "free")
ps.ra %>% 
  subset_taxa(Genus=="Mycobacterium") %>%
  plot_bar(fill="month") +
  facet_grid(~year+neonic, scales = "free")
ps.ra %>% 
  subset_taxa(Genus=="Gemmatimonas") %>%
  plot_bar(fill="month") +
  facet_grid(~year+neonic, scales = "free")
ps.ra %>% 
  subset_taxa(Genus=="Ellin6067") %>%
  plot_bar(fill="month") +
  facet_grid(~year+neonic, scales = "free")

#Linear model ####
shn.lm = lm(Shannon ~ neonic * month * year, data=shn.rich)
summary(shn.lm)
anova(shn.lm)
#quick code to find the mean for the comparisons (instead of checking the shn.rich matrix)
(HSD.test(aov(shn.lm), c("neonic","month","year"), group = TRUE))$groups 

#Richness-interaction####
#test
compare_means(Shannon ~ neonic, shn.rich, method = "wilcox.test", paired = FALSE, p.adjust.method = "holm",
              group.by = "month") 
#subset Month ####
#July
shn.j.n = shn.rich %>% 
  filter(month == "July" & neonic == "N") 
summary(shn.j.n)
sd(shn.j.n$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.j.n$Shannon[!is.na(shn.j.n$Shannon)])) #SE
shn.j.y = shn.rich %>% 
  filter(month == "July" & neonic == "Y") 
summary(shn.j.y)
sd(shn.j.y$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.j.y$Shannon[!is.na(shn.j.y$Shannon)])) #SE
#Aug
shn.a.n = shn.rich %>% 
  filter(month == "Aug" & neonic == "N") 
summary(shn.a.n)
sd(shn.a.n$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.a.n$Shannon[!is.na(shn.a.n$Shannon)])) #SE
shn.a.y = shn.rich %>% 
  filter(month == "Aug" & neonic == "Y") 
summary(shn.a.y)
sd(shn.a.y$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.a.y$Shannon[!is.na(shn.a.y$Shannon)])) #SE
#Sep
shn.s.n = shn.rich %>% 
  filter(month == "Sep" & neonic == "N") 
summary(shn.s.n)
sd(shn.s.n$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.s.n$Shannon[!is.na(shn.s.n$Shannon)])) #SE
shn.s.y = shn.rich %>% 
  filter(month == "Sep" & neonic == "Y") 
summary(shn.s.y)
sd(shn.s.y$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.s.y$Shannon[!is.na(shn.s.y$Shannon)])) #SE

#test
compare_means(Shannon ~ neonic, shn.rich, method = "wilcox.test", paired = FALSE, p.adjust.method = "holm",
              group.by = "year") 
#Subset year ####
#2016
shn.16.n = shn.rich %>% 
  filter(year == "2016" & neonic == "N") 
summary(shn.16.n)
sd(shn.16.n$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.16.n$Shannon[!is.na(shn.16.n$Shannon)])) #SE
shn.16.y = shn.rich %>% 
  filter(year == "2016" & neonic == "Y") 
summary(shn.16.y)
sd(shn.16.y$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.16.y$Shannon[!is.na(shn.16.y$Shannon)])) #SE
#2017
shn.17.n = shn.rich %>% 
  filter(year == "2017" & neonic == "N") 
summary(shn.17.n)
sd(shn.17.n$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.17.n$Shannon[!is.na(shn.17.n$Shannon)])) #SE
shn.17.y = shn.rich %>% 
  filter(year == "2017" & neonic == "Y") 
summary(shn.17.y)
sd(shn.17.y$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.17.y$Shannon[!is.na(shn.17.y$Shannon)])) #SE
#2018
shn.18.n = shn.rich %>% 
  filter(year == "2018" & neonic == "N") 
summary(shn.18.n)
sd(shn.18.n$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.18.n$Shannon[!is.na(shn.18.n$Shannon)])) #SE
shn.18.y = shn.rich %>% 
  filter(year == "2018" & neonic == "Y") 
summary(shn.18.y)
sd(shn.18.y$Shannon, na.rm=TRUE) /  
  sqrt(length(shn.18.y$Shannon[!is.na(shn.18.y$Shannon)])) #SE

#tukey
#shn.lm.hsd = lm(Shannon ~ neonic * month, data=shn.rich) %>%
#  tukey_hsd()
#shn.lm.hsd = shn.lm.hsd[shn.lm.hsd$p.adj<0.05,]
#kable(shn.lm.hsd) %>%
 # kable_styling(bootstrap_options = "striped", full_width = F)
#neonic
compare_means(Shannon ~ neonic, shn.rich, method = "wilcox.test", paired = FALSE, p.adjust.method = "holm") 


#neonic-month ####
stat.test = compare_means(Shannon ~ neonic, data = shn.rich,method = "wilcox.test", 
                          p.adjust.method = "holm", group.by = "month")
stat.test
#neonic-year ####
compare_means(Shannon ~ neonic, data = shn.rich,method = "wilcox.test", 
              p.adjust.method = "holm", group.by = "year")
#(Shannon ~ year, data = shn.rich,method = "wilcox.test", 
 #             p.adjust.method = "holm", group.by = "neonic")
#neonic-month&year
compare_means(Shannon ~ neonic, data = shn.rich,method = "wilcox.test", 
                          p.adjust.method = "holm", group.by = c("month","year"))
#plot neonic-month
ggplot(shn.rich, aes(x = neonic, y = Shannon, color=neonic, alpha = 0.1)) +  
  theme_bw() +
  geom_boxplot(aes(fill=neonic)) + 
  scale_x_discrete(labels = c(N="Control",Y="Neonicotinoid-treated")) +
  scale_color_manual(values = c("cornflowerblue","mediumvioletred")) +
  scale_fill_manual(values = c("cornflowerblue","mediumvioletred")) +
  facet_wrap(~month, labeller=labeller(month = c(July="July",Aug="August",Sep="September"))) +
  stat_pvalue_manual(data = stat.test, label = "p.adj", xmin = "group1", xmax = "group2", size = 5,
                     y.position = c(7.5,7.5), hide.ns = TRUE, bracket.size = 0) +
  stat_compare_means(comparisons = list(c("N","Y")), size = 5, 
                     method = "wilcox.test", label = "p.signif", label.y.npc = 0.5, label.y = c(7.55,7.55,7.55)) +
  geom_jitter(alpha = 0.5) +
  ylab("Soil Shannon Richness") +
  theme(legend.position = 'none',
        strip.text.x = element_text(size=14, face="bold"),
        axis.text.x = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(face = "bold", size=12),
        axis.title.x = element_blank()) 

#neonic-time####
compare_means(Shannon ~ neonic, shn.rich, method = "wilcox.test", paired = FALSE, 
                            p.adjust.method = "holm", group.by = c("month","year")) %>%
  filter(p.format<0.05)



########################
#sample_data(ps)$neo_mnt = as.factor(paste(sample_data(ps)$neonic, sample_data(ps)$month, sep="_"))

ps.ra = transform_sample_counts(ps, function(otu) otu/sum(otu)) 
ps.fam = tax_glom(ps.ra, "Family")
melt.fam = psmelt(ps.fam)
melt.fam03 = subset(melt.fam, Abundance > 0.03)
length(melt.fam03$Family)
ggplot(data = melt.fam03, mapping = aes_string(x = "neonic",y = "Abundance",
                                            color = "Family", fill = "Family")) +
  geom_boxplot() +
  geom_point(size = 1, alpha = 0.3,
             position = position_jitter(width = 0.3)) +
  facet_grid(month~year) # scale_y_log10()
  #theme(legend.position="none")
  
  
ps.gen = tax_glom(ps.ra, "Genus")
melt.gen = psmelt(ps.gen)
mean(melt.gen$Abundance)
melt.gen01 = subset(melt.gen, Abundance > 0.01)
length(melt.gen01$Genus)
ggplot(data = melt.gen01, mapping = aes_string(x = "neonic",y = "Abundance",
                                               color = "Genus", fill = "Genus")) +
  ylab("Relative Abundance - log10") +
  geom_boxplot() +
  geom_point(size = 1, alpha = 0.3,
             position = position_jitter(width = 0.3)) +
  facet_grid(month~year) + 
  scale_y_log10()+
  theme(legend.position="right")
  
ggplot(melt.gen01, aes(x=year, y=Abundance, group=Family, color=Family)) +
  geom_smooth(method = "lm")+
  facet_wrap(~neonic)


ggplot(melt.gen01, aes(x=year, y=Abundance, group=Family, color=Family)) +
  geom_point(size = 2) +                                               # Changing point size
  geom_smooth(method = "lm", aes(fill = Family))

save.image("/data/users/mona/miseq_16S/Mona_16S_all/article1/a1_4_s1_exp_aca_16S_soil.RData")
