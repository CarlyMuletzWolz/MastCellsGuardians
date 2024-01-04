### By: Carly Muletz Wolz

### DADA2 pipeline and creating feature table, taxonomy table and sequence file for later analyses


#############  INITIAL processing of files in terminal to get this into dada2 format #############

## From basespace, the files are downloaded with each sample having a folder 
# and within that folder are the forward and reverse reads

## To make it compatible with dada2 and R we need to copy all the fastq files to one folder 
# Then you can delete the old folder that used to hold the reverse and forward reads. 

## You need to navigate to the project folder (most likely called FASTQ_Generation...) in terminal 
# cd /Users/lkgentry/Library/CloudStorage/OneDrive-SmithsonianInstitution/Kelsey\'s\ Project/Kelsey-2-371840470/FASTQ_Generation_2022-11-15_19_49_06Z-629625305
# Move all the files within each folders up one folder to the FASTQ files, which will then hold all
# of the fastq files in one main folder as opposed to per sample folders
# then remove the folders with nothing in them now, and unzip

# mv -v *_L001*/* .     ### move all files within those folders that have L001 in them up one folder (forward and reverse reads) 

# rm -r *_L001*/        ### remove all the empty folders

### YOU MUST UNZIP all of the files

# gunzip *_L001*

###############  INSTALL DADA2 if you don't already have it.  ################
### Follow tutorial on how to install (follow 1. and 2.) https://benjjneb.github.io/dada2/dada-installation.html

library(dada2)
packageVersion("dada2")

library(devtools)
devtools::install_github("benjjneb/dada2")

## received warning:no function found corresponding to methods exports from ‘GenomicAlignments’ for: ‘concatenateObjects’


########## DADA2 tutorial is very helpful: https://benjjneb.github.io/dada2/tutorial.html ##########


## Run data is here. On home desktop. 

setwd("/Users/lkgentry/Library/CloudStorage/OneDrive-SmithsonianInstitution/Kelsey\'s\ Project/Kelsey-2-371840470/FASTQ_Generation_2022-11-15_19_49_06Z-629625305")


path <- "/Users/lkgentry/Library/CloudStorage/OneDrive-SmithsonianInstitution/Kelsey\'s\ Project/Kelsey-2-371840470/FASTQ_Generation_2022-11-15_19_49_06Z-629625305"
list.files(path)



#FILTER AND TRIM 

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#INSPECTION OF QUALITY PROFILES

plotQualityProfile(fnFs[8:16])


plotQualityProfile(fnRs[8:16])


## Quality looks great for forward, starts dropping on reverse around 175
## Since need to trim forward primer off 19 bp, will do 275. We want 450 bp merged

#FILTERING AND TRIMMING 

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))


## parameters for filtering data

## Not uncommon to lose 6 to 10k sequences between reads.in and reads.out
##maxee to 2,5 

#out3 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,190),
                   #   maxN=0, maxEE=c(2,2), trimLeft = 19, trimRight = 23, 
                    #  truncQ=2, rm.phix=TRUE,
                    #  compress=TRUE, multithread=TRUE) 
out2 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,190),
                      maxN=0, maxEE=c(2,5), trimLeft = 19, trimRight = 23, 
                      truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE)
##out1 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,190),
                     # maxN=0, maxEE=c(2,5), trimLeft = 19, trimRight = 23, 
                      #truncQ=2, rm.phix=TRUE,
                      #compress=TRUE, multithread=TRUE)

#out4 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,190),
                     # maxN=0, maxEE=c(2,5), trimLeft = 19, trimRight = 23, 
                     # truncQ=2, rm.phix=TRUE,
                     # compress=TRUE, multithread=TRUE)


head(out2)
str(out2)
mean(out2[,2])
mean(out2[,1])

## on average 80% with this. Seems good. 
mean(out2[,2])/mean(out2[,1])


## NEED to run whichever parameter last so that filtFs and filtRs are from this
# out3 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,180), 
                      #maxN=0, maxEE=c(2,2), trimLeft = 19, trimRight = 23,
                      #truncQ=2, rm.phix=TRUE,
                      #compress=TRUE, multithread=TRUE)
#head(out3)
#mean(out3[,2])

## on average 80% with this, sticking with first one. they are the same... 
#mean(out3[,2])/mean(out3[,1])

## GOING with out3

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


plotErrors(errF, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]


mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

setwd("/Users/lkgentry/OneDrive - Smithsonian Institution/GrayferLabExperiments/KelseyAnalysisNov2022/")
      
## Used out 1 parameters
saveRDS(seqtab, "XengK2022exp2_seqtab.rds")

##Had to install packages and updates below

#if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.16")

#BiocManager::install("DECIPHER")

#Do this code later and see what it looks like

## If want to read back in 
seqtab <- readRDS("XengK2022exp2_seqtab.rds")

dim(seqtab)

## Distribution of amplicon sizes in bp
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)

## Still retain 94.5% of sequences after chimera removal
sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))

## Make sure you change 'out1' to whatever output you end up selecting
track <- cbind(out2, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(track) <- sample.names

head(track)


write.csv(track, "dada2_output.csv")


## https://benjjneb.github.io/dada2/training.html
## Assign taxonomy is up three directories so that I can use these files for multiple projects
taxa <- assignTaxonomy(seqtab.nochim, "/Users/lkgentry/Documents/Microbiome Analysis/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa <- addSpecies(taxa, "/Users/lkgentry/Documents/Microbiome Analysis/silva_species_assignment_v138.1 (1).fa.gz")



#  inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


library(phyloseq); packageVersion("phyloseq")

## combine feature table and taxonomy table in same order
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))
ps

## rename ASVs to numbers
new.names <- paste0("ASV", seq(ntaxa(ps))) # Define new names ASV1, ASV2, ...
seqs <- taxa_names(ps) # Store sequences
names(seqs) <- new.names # Make map from ASV1 to full sequence
taxa_names(ps) <- new.names # Rename to human-friendly format


## convert feature table to matrix
site_species <-as(otu_table(ps), "matrix")

## need to change this to match mapping file later
rownames(site_species)

# not sure what this is doing. needed it for another project. Might be useful for underscores?
samples.out <- rownames(site_species)
rownames(site_species) <- sapply(strsplit(samples.out, "f"), `[`, 1)
rownames(site_species)

## transpose to make a species by site matrix

species_site <- t(site_species)

# taxon table 
tax <- as(tax_table(ps), "matrix")



getwd()
setwd("/Users/lkgentry/OneDrive - Smithsonian Institution/GrayferLabExperiments/KelseyAnalysisNov2022/")

## Write this file out and look at it. Determine what to do with ASVs in Neg Controls
## select whole worksheet and filter by reads in neg controls, 
## MAKE SURE you have selected all columns/rows before you filter

## Remove ASVs if they are in all negative controls or in one negative control, but almost all samples
## pretty subjective, no guidelines really on what to do here
## Once done making neg control decisions, remove negative controls from file 
## and rename _final.csv
## Can also use 'decontam' package. That is what we are doing now in most cases.




write.csv(species_site, "XENGK2022EXP2_feature_table.csv")

write.csv(tax, "XENGK2022EXP2_taxonomy.csv")
write.csv(seqs, 'XENGK2022EXP2_feature_DNAsequences.csv')




writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))  }
  
fileConn<-file(filename)
writeLines(fastaLines, fileConn)
close(fileConn)
}

seqsDf <- as.data.frame(seqs)
seqsDf$ASV <- row.names(seqsDf)

seqsDf <- seqsDf[ ,c(2,1)]

colnames(seqsDf) <- c('name','seq')

writeFasta(seqsDf, "XENGK2022EXP2_DNAsequences.fasta")


##library(seqRFLP) No longer available
## you want to generate a fasta file of sequences for reading back in 
seq_data <- dataframe2fas(seqs, file = "XENGK2022EXP2_feature_DNAsequence.fasta")


####### FILTERING #########

## For this case with fewer samples and fewer neg controls, just going to manually assess.
## Removed ASV270 was only ASV in first negPCR control and not in any samples
## No ASVs in second negPCR control
## Removed this sample, T21SXENG4, only had 1 sequence, all other samples great coverage


###called XENGK2022_feature_table2.csv

## still going to take over to decontam for the neg extraction control, which has several ASVs in it 








