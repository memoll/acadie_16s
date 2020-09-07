#######################################################################################
# Explanatory analysis of the phyllosphere and samples                                #
# Studying the effects of neonicotinoids on the phyllosphere bacterial communities    #
# Data: Miseq-16S - L'Acadie (ACA)                                                    #
# Mona Parizadeh - 2019-2020                                                          #
#######################################################################################

# Load libraries
library(phyloseq); packageVersion("phyloseq") #‘1.27.6’
library(vegan); packageVersion("vegan") #‘2.5.6’
library(ggplot2); packageVersion("ggplot2") #‘3.3.0’

# Import data #### 
setwd("../mp/aca_16s/files/")
ps = readRDS("16S_aca_phyl5000.rds")
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
set.seed(95018)
adns = adonis2(dis ~ host*year*month*neonic, df) #distance = bray
adns

#Ordination #### 
pcoa = ordinate(ps, method = "MDS", k = 2, try = 100, distance = "bray")

#unifrac:
#ordu = ordinate(ps, "PCoA", "unifrac", weighted=TRUE)
#plot_ordination(ps, ordu, color="neonic", shape="neonic") +
#  geom_text(aes(label = sampleid), check_overlap = FALSE, size = 5) +
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
group.hst.neo.yr.fac = factor(group.hst.neo.yr,levels=c("soyN2016","soyY2016","cornN2017",
                                                        "cornY2017","soyN2018","soyY2018"))
#define labels and colors
hst.neo.labels = c(soyN2016="Control Soybean 2016", soyY2016="Neonicotinoid-treated Soybean 2016",
                   cornN2017="Control Corn 2017", cornY2017="Neonicotinoid-treated Corn 2017",
                   soyN2018="Control Soybean 2018", soyY2018="Neonicotinoid-treated Soybean 2018")
hst.neo.colors = c("cornflowerblue","mediumvioletred","cornflowerblue","mediumvioletred","cornflowerblue", "mediumvioletred")

#ordinate
pcoa = ordinate(ps, method = "MDS", k = 2, try = 100, distance = "bray")
#define variables as factors
hst = get_variable(ps, "host")
sample_data(ps)$host = factor(hst)
neo = get_variable(ps, "neonic")
sample_data(ps)$neonic = factor(neo)
yr = get_variable(ps, "year")
sample_data(ps)$year = factor(yr)
#% variables
group.hst.neo.yr = paste(hst, neo, yr, sep = "")
#order as factor
group.hst.neo.yr.fac = factor(group.hst.neo.yr,levels=c("soyN2016","soyY2016","cornN2017",
                                                        "cornY2017","soyN2018","soyY2018"))
#order as factor
sample_data(ps)$month = factor(sample_data(ps)$month,levels=c("July","Aug","Sep"))
#define labels and colors
hst.neo.colors = c("cornflowerblue","mediumvioletred","cornflowerblue","mediumvioletred","cornflowerblue", "mediumvioletred")

#ordinate
pcoa1 = plot_ordination(ps, pcoa) 
pcoa1$layers
pcoa1$layers = pcoa1$layers[-1] #remove the original points to add the desired colors and shapes

#group
#% ellipses
group.hst.neo = paste(hst, neo, sep = "")

p.phyl = pcoa1.phyl + 
  theme_bw() +
  #group by neonic
  stat_ellipse(aes(fill = group.hst.neo, group=group.hst.neo), geom = "polygon", 
               level = 0.95, linetype = 0, alpha = 0.4, show.legend=FALSE) +
  scale_fill_manual(name="neonic & host",values = c("cornflowerblue","darkgoldenrod2","cornflowerblue","darkgoldenrod2")) +
  geom_point(aes(color = group.hst.neo.yr.fac, shape = group.hst.neo.yr.fac), size=3, alpha=0.75) + 
  #change facet font size
  theme(strip.text.x = element_text(size=10, face="bold"),
        strip.text.y = element_text(size=10, face="bold")) +
  #change legend labels 
  scale_color_manual(name = "Phyllosphere: Treatment, Host & Year",
                     values = phyl.hst.neo.colors, 
                     labels = hst.neo.labels) +
  scale_shape_manual(name = "Phyllosphere: Treatment, Host & Year", 
                     values = c(19,19,17,17,15,15) ,
                     labels = hst.neo.labels) +
  annotate("text", x = 0.3, y = 0.3, label = 'atop(bold("Soybean"))', parse = TRUE, 
           size = 5, colour = "#283747") +
  annotate("text", x = -0.3, y = -0.4, label = 'atop(bold("Corn"))', parse = TRUE, 
           size = 5, colour = "#283747") +
  theme(legend.title = element_text(size=12, face="bold"),
        legend.text=element_text(size=12),
        legend.position = "none",
        axis.title = element_text( size = 12, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  guides(col = guide_legend(nrow = 6, title.position = "top")) +
  labs(tag = "A)") + theme(plot.tag = element_text(size = 12, face = "bold")) 

p.phyl.mnt = pcoa1 + 
  theme_bw() +
  #group by neonic
  stat_ellipse(aes(fill = neonic, group=neonic), geom = "polygon", 
               level = 0.95, linetype = 0, alpha = 0.4, show.legend=TRUE) +
  scale_fill_manual(name="Treatment ellipses",values = c("cornflowerblue","darkgoldenrod2"),
                    labels = c(N="Control",Y="Neonicotinoid-treated")) +
  geom_point(aes(color = group.hst.neo.yr.phyl.fac, shape = group.hst.neo.yr.phyl.fac), size=3, alpha=0.75) + 
  facet_wrap(~ month,
             labeller=labeller(month = c(May="May",July="July",
                                         Aug="August",Sep="September"))) +
  #change facet font size
  theme(strip.text.x = element_text(size=14, face="bold")) +
  #change legend labels 
  scale_color_manual(name = "Phyllosphere: Treatment, Host & Year",
                     values = phyl.hst.neo.colors, 
                     labels = hst.neo.labels) +
  scale_shape_manual(name = "Phyllosphere: Treatment, Host & Year", 
                     values = c(19,19,17,17,15,15) ,
                     labels = hst.neo.labels) +
  theme(legend.title = element_text(size=14, face="bold"),
        legend.text=element_text(size=14),
        legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text( size = 14, face = "bold")) +
  guides(col = guide_legend(nrow = 6, title.position = "top")) +
  labs(tag = "B)") + theme(plot.tag = element_text(size = 16, face = "bold")) 

#with legend (to be cropped)
p.phyl.leg = pcoa1 + 
  theme_bw() +
  #group by neonic
  stat_ellipse(aes(fill = neonic, group=neonic), geom = "polygon", 
               level = 0.95, linetype = 0, alpha = 0.4, show.legend=TRUE) +
  scale_fill_manual(name="Treatment ellipses",values = c("cornflowerblue","darkgoldenrod2"),
                    labels = c(N="Control",Y="Neonicotinoid-treated")) +
  geom_point(aes(color = group.hst.neo.yr.phyl.fac, shape = group.hst.neo.yr.phyl.fac), size=3, alpha=0.75) + 
  facet_wrap(~ month,
             labeller=labeller(month = c(May="May",July="July",
                                         Aug="August",Sep="September"))) +
  #change facet font size
  theme(strip.text.x = element_text(size=10, face="bold"),
        strip.text.y = element_text(size=10, face="bold")) +
  #change legend labels 
  scale_color_manual(name = "Phyllosphere: Treatment, Host & Year",
                     values = phyl.hst.neo.colors, 
                     labels = hst.neo.labels) +
  scale_shape_manual(name = "Phyllosphere: Treatment, Host & Year", 
                     values = c(19,19,17,17,15,15) ,
                     labels = hst.neo.labels) +
  theme(legend.title = element_text(size=12, face="bold"),
        legend.text=element_text(size=12),
        legend.position = "right",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text( size = 10, face = "bold")) +
  guides(col = guide_legend(nrow = 6, title.position = "top"))  

#Fig. 2A&B ####
#library(patchwork); packageVersion("patchwork")
#p.phyl + p.phyl.mnt + plot_layout(widths = c(1, 1.5), guides = 'collect') +
 # plot_annotation(tag_levels = 'A', tag_suffix = ')') &
 # theme(plot.tag = element_text(size = 12, face="bold"))
library(cowplot)
save_plot("/data/users/mona/miseq_16S/Mona_16S_all/article1/graphs/ParizadehM_Fig2A.pdf", p.phyl, ncol = 1, nrow = 1)
save_plot("/data/users/mona/miseq_16S/Mona_16S_all/article1/graphs/ParizadehM_Fig2B.pdf", p.phyl.mnt, ncol = 2, nrow = 2)
save_plot("/data/users/mona/miseq_16S/Mona_16S_all/article1/graphs/ParizadehM_Fig2Bleg.pdf", p.phyl.leg, ncol = 2, nrow = 1)

#HOSTS ####
#% Soybean ####
#PERMANOVA ####
#make dataframe
df.sy = as(sample_data(ps.sy), "data.frame")
#bray-curtis distance
dis.sy = phyloseq::distance(ps.sy,  method = "bray")
set.seed(950181)
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
mnt.neo.colors = c("cornflowerblue","mediumvioletred","cornflowerblue","mediumvioletred","cornflowerblue", "mediumvioletred")
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
  geom_point(aes(color = group.mnt.neo.sy.fac, shape = group.mnt.neo.sy.fac), size=2, alpha=0.75) + 
  facet_wrap(~ year, labeller=labeller(year = c("2016"="Soybean 2016", "2018"="Soybean 2018"))) +
  #change facet font size
  theme(strip.text.x = element_text(size=10, face="bold")) +
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
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        legend.justification = "left") +
  labs(tag = "A)") + theme(plot.tag = element_text(size = 12, face = "bold"))

#% Corn ####
#PERMANOVA ####
#make dataframe
df.cr = as(sample_data(ps.cr), "data.frame")
#bray-curtis distance
dis.cr = phyloseq::distance(ps.cr,  method = "bray")
set.seed(950182)
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
  geom_point(aes(color = group.mnt.neo.cr.fac, shape = group.mnt.neo.cr.fac), size=2, alpha=0.75) + 
  #change facet font size
  facet_wrap(~ year, labeller=labeller(year = c("2017"="Corn 2017"))) +
  theme(strip.text.x = element_text(size=10, face="bold")) +
  #change legend labels 
  scale_shape_manual(name = "Month & Treatement", 
                     labels = mnt.neo.labels,
                     values = c(19, 19, 17, 17, 15, 15)) +
  scale_color_manual(name = "Month & Treatement", 
                     labels = mnt.neo.labels,
                     values = mnt.neo.colors) +
  theme(legend.title = element_text(size=12, face="bold"),
        legend.text=element_text(size=12),
        axis.text = element_text(size = 10),
        axis.title = element_text( size = 12, face = "bold"),
        #legend.justification = "left"
        legend.box = "horizontal") +
  labs(tag = "B)") + theme(plot.tag = element_text(size = 12, face = "bold"))

#Fig. 3 ####
#patch ordinations
library(patchwork);packageVersion("patchwork") #‘1.0.0’
p.syA / p.crB + 
  plot_layout(nrow = 2, height = c(1,1), width = c(1, 1), guides = 'collect') +
  plot_annotation(tag_levels = 'A', tag_suffix = ')') &
  theme(plot.tag = element_text(size = 9, face="bold")) 
ggsave("/data/users/mona/miseq_16S/Mona_16S_all/article1/graphs/ParizadehM_Fig3.pdf", fig3, dpi = 300, width = 170, units = "mm")
#ggsave("/data/users/mona/miseq_16S/Mona_16S_all/article1/graphs/ParizadehM_Fig3.jpeg", fig3, scale = 2, dpi = 1200, width = 6.6)
save_plot("/data/users/mona/miseq_16S/Mona_16S_all/article1/graphs/ParizadehM_Fig3A.pdf", p.syA, ncol = 2, nrow = 1)
save_plot("/data/users/mona/miseq_16S/Mona_16S_all/article1/graphs/ParizadehM_Fig3B.pdf", p.crB, ncol = 2, nrow = 1)

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
compare_means(Shannon ~ neonic, shn.rich, method = "wilcox.test", paired = FALSE, p.adjust.method = "holm",
              group.by = "month") 
compare_means(Shannon ~ neonic, shn.rich, method = "wilcox.test", paired = FALSE, p.adjust.method = "holm",
              group.by = "year") 


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

#Linear model ####
shn.lm = lm(Shannon ~ neonic * month, data=shn.rich)
summary(shn.lm)
anova(shn.lm)
shn.lm.hsd = lm(Shannon ~ neonic * month, data=shn.rich) %>%
  tukey_hsd()
shn.lm.hsd = shn.lm.hsd[shn.lm.hsd$p.adj<0.05,]
kable(shn.lm.hsd) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)

ggplot(shn.rich, aes(x = neonic, y = Shannon, color=month, alpha = 0.1)) + #alpha: color intensity)) 
  theme_bw() +
  geom_boxplot(aes(fill=month))

#Temporal variation-composition ####
set.seed(950188)
adns.t = adonis2(dis ~ year*month*neonic, df, strata = host)
adns.t

save.image("/data/users/mona/miseq_16S/Mona_16S_all/article1/a1_3_p1_exp_aca_16S_phyl.RData")
