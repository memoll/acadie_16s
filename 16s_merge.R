#####################################################################
# Script to merge multiple sequencing runs results of DADA2         #
# Based on DADA2 workflow for Big Data by Benjamin Callahan         #
# Data: Miseq-16S-seqtab                                            #
# Mona Parizadeh - 2018-2019                                        #
#####################################################################

# Load libarary ####
library(dada2); packageVersion("dada2") #‘1.12.1’
library(phyloseq); packageVersion("phyloseq") #‘1.27.6’

# merge multiple runs ####
st_A = readRDS("../mp/aca_16s/files/seqtab_A.rds")
st_B = readRDS("../mp/aca_16s/files/seqtab_B.rds")
st_C = readRDS("../mp/aca_16s/files/seqtab_C.rds")
st_all = mergeSequenceTables(st_A, st_B, st_C)

# Assign taxonomy using Silva 132 ####
#Download the Silva taxonomic training data formatted for DADA2, derived from the Silva Project's version 132 release httfps://zenodo.org/record/1172783#.XHlu0-JKjfY
taxaSilva <- assignTaxonomy(st_all, "../mp/aca_16s/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
unname(head(taxaSilva))
dim(taxaSilva) 

# Add species:
system.time({taxaSilva.spe = addSpecies(taxaSilva, "../mp/aca_16s/tax/silva_species_assignment_v132.fa.gz", verbose = TRUE)
colnames(taxaSilva.spe) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
unname(taxaSilva.spe)})

# Inspect the taxonomic assignments 
taxaSilva.spe.print <- taxaSilva.spe # Removing sequence rownames for display only
rownames(taxaSilva.spe.print) <- NULL
head(taxaSilva.spe.print)

# Import the biosample (metadata) file ####
#It contains informations about samples and treatments
map = import_qiime_sample_data("../mp/aca_16s/files/aca_biosample.csv")

# Make phyloseq object ####
ps = phyloseq(otu_table(seqtab, taxa_are_rows = FALSE), sample_data(map), tax_table(taxaSilva.spe))
#save
saveRDS(ps, file="../mp/aca_16s/files/ps.rds")
