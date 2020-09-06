#################################################
# Script to analyze 16S data with DADA2         #
# Based on DADA2 tutorial by Benjamin Callahan  #
# Data: Miseq-16S-demultiplexed                 #
# Mona Parizadeh - 2018-2019                    #
#################################################

# I have used the same script for all the three sequencing runs; however, I had to set few parameters differently among the runs. 
# You can find the parameter settings for the 2nd and the 3rd runs (if they were different from the 1st one) in the comments below the related function.

# Load libarary ####
library(dada2); packageVersion("dada2") #‘1.12.1’
require(parallel)

# Define path ####
path <- "../mp/aca_16s/dada2"
list.files(path)

# Filter and trim ####
#Lists forward and reverse fastq files and sorts the files to make sure they are all in the same direction (R1 and R2)
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE)) # F 
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE)) # R 
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)
head(sample.names)

# Examine quality profiles of forward and reverse reads
plotQualityProfile(fnFs[1:2]) # forward reads 
plotQualityProfile(fnRs[1:2]) # reverse reads 

# Perform filtering and trimming 
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

detectCores() #in order to choose multithread

# Filter the forward and reverse reads
#takes time
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,270), 
                     maxN=0, maxEE=c(2,2), minLen = 50, truncQ=2, rm.phix=FALSE, compress=TRUE, 
                     verbose = TRUE, matchIDs = TRUE , multithread=TRUE) 

#run2: truncLen=c(220,260); run3: truncLen=c(230,270)

#maxN : max number of ambiguous bases allowed, since dada2 cannot handle ambiguous bases keep it at 0.
#maxEE : max number of estimated errors allowed by an individual read; increase it only if you have low quality reads
#truncQ : truncates reads at the first instance of a Q score less than or equal to the value specified
#fastqPairedFilter(): funtion to check each pair (one forward and the related reverse) seperately to troubleshoot filterAndTrim() 

head(out)

# Calculate error rates ####
#takes time
errF <- learnErrors(filtFs, multithread=TRUE)  
errR <- learnErrors(filtRs, multithread=TRUE) 
#play with the arguments to improve the convergence, if needed; example: learnErrors(..., nbases=1e8, MAX_CONSIST = 20, multithread = TRUE)

#In the error plot, the red line shows the error rates expected under the nominal definition of the Q-value.
#As error rate decrease, the Q (quality) score increases.

# Dereplicate the filtered fastq files ####
#This combines all identical sequencing reads into “unique sequences” with a corresponding “abundance”: the number of reads with that unique sequence, to reduce computation time by eliminating redundant comparisons
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Infer the sequence variants in each sample ####
#This prodeuces an object that describes dada2 denoising results (removes all sequencing errors)
#Using the developed error model, it calculates abundance p-values for each unique sequence and checks if a given error rate is too abundant in a sample to be explained by sequencing errors
#takes time
dadaFs <- dada(derepFs, err=errF, pool="pseudo", multithread=TRUE)   
dadaRs <- dada(derepRs, err=errR, pool="pseudo", multithread=TRUE)

#if p-value is low: there are more reads of sequence than can be explained by sequencing errors
#if p-value is high: read is caused by error
#pool = pseudo is recommended when the diversity is high in samples; it approximates pooling in linear time and so takes less time.
#pool = TRUE, the algorithm will pool together all samples prior to sample inference.
#Pooling can increase sensitivity to rare per-sample variants, but it takes time.
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16
dadaFs[[1]] 
dadaRs[[1]] 
dadaFs[[3]] 
dadaRs[[3]] 

# Merge the denoised forward and reverse reads ####
#This merges error-free forward and reverse reads (after being denoised) and removes paired reads that did not exactly overlapped
#takes time
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 10, maxMismatch = 0, returnRejects = TRUE, verbose=TRUE) 
#may need to comprimise the previous trimming to make reads overlap
#maxMismatch = 0 (since we have already removed errors)
#minOverlap = 10; by defaulyt it's 20 nt to overlap, which might be a lot for high variable regions
#verbose = TRUE: to print a summary of the function results
#Inspect the merger data.frame from the first sample
#returnRejects = TRUE : the pairs that were rejected based on mismatches in the overlap region are retained in the return data.frame.
#if the primers overhang (hang or extend outward over), use overhang = TRUE; see ITS workflow.
#trimOverhang = TRUE; in case any of our reads go passed their opposite primers (not necessary)
head(mergers[[1]])

# Construct sequence table ####
seqtab <- makeSequenceTable(mergers) 
#The sequences being tabled vary in length; analogous to an OTU table
dim(seqtab) 
#seqtab is a matrix w/ rows corresponding to (and named by) the samples, and columns corresponding to (and named by) the sequence variants. 
#The lengths of the merged sequences should all fall within the expected range for this V4 amplicon.

# Inspect distribution of sequence lengths (length of ASVs):
table(nchar(getSequences(seqtab)))
#top: size (nt), below: frequency

#In order to merge this dataset with the other sequencing runs, for each run we need to save the seqtab at this moment:
# Run 1:
saveRDS(seqtab, file = "../mp/aca_16s/files/seqtab_A.rds")
# Run 2: saveRDS(seqtab, file = "../mp/aca_16s/files/seqtab_B.rds")
# Run 3: saveRDS(seqtab, file = "../mp/aca_16s/files/seqtab_C.rds")




