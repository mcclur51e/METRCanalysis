##########################################################################################################
##### DESCRIPTION                                                                                    #####
##### Code for preliminary processing of .fastq files from Illumina sequencing. Processing includes: #####
##### 1. standard DADA2 pipeline                                                                     #####
##### 2. assigning taxonomy using a modified silva_nr_v132 database                                  #####
##### 3. passing data into phyloseq object                                                           #####
##### 4. correcting a few mis-identified important taxa                                              #####
##### 5. adding columns to mapping file based on calculations from processed sequencing data         #####
##### 6. final output = physeqR                                                                      #####
##### Use this file to prepare initial raw phyloseq object for all future downstream analysis.       #####
###########################################################################################################################################################
##### This code requires as input:                                                                                                                    #####
##### 1. A working repository (default "~/Desktop/Wannigan_METRC/")                                                                                   #####
##### 2. Folder within working repository containing all .fastq files (default "~/Desktop/Wannigan_METRC/raw/16S/")                                   #####
##### 3. .csv file containing sample metadata listed with unique sampleIDs in first column (default "~/Desktop/Wannigan_METRC/raw/table_map.csv")     #####
##### 4. .fa file containing reference taxonomy database (default "~/Masters/silva_nr_v132_train_set_RossMod.fa")                                     #####                                   
##### 5. .csv file containing a list of sequences to be assigned to Homo sapiens (default "~/Desktop/Wannigan_METRC/raw/list_homoASVs.csv")           #####
###########################################################################################################################################################

#####################################################################
########## Processing Set −up #######################################
#####################################################################
# R version 4.5.2 (2025-10-31) -- "[Not] Part in a Rumble"
# All code written by Emily Ann McClure. No AI was used at any stage in writing or editing this code. 
# github.com/mcclure51e/
########## Call libraries for use ##########
# Packages from CRAN
library("ggplot2") # version 4.0.0
library("data.table") # version 1.17.8
library("dplyr") # version 1.1.4

# Packages from Bioconductor
library("phyloseq") # version 1.52.0
library("decontam") # version 1.28.0 # identify contaminant ASVs 
library("dada2") # version 1.36.0

##################################################################################
########## If running from scrum (or similar), use the following block ###########
##################################################################################
#input <- commandArgs(trailingOnly = TRUE) # Read arguments when beginning script in scrum
#path <- as.character(input[1])
#setwd(path) # set working directory
#dir.create("output") # create directory for output files to go

################################################################################
########## If running by line (i.e. RStudio) use the following block ###########
################################################################################
path <- "~/Desktop/Wannigan_METRC/" # change this to directory where you will be working
setwd(path) # set working directory
#dir.create("outputPrelim") # create directory for output files to go
#dir.create("plots") # create directory for plots to go

################################################################################
##### ASV assignment with DADA2 ################################################
################################################################################
fnFs <- sort(list.files(paste0(path,"raw/16S/"), pattern="R1_001", full.names = TRUE))
fnRs <- sort(list.files(paste0(path,"raw/16S/"), pattern="R2_001", full.names = TRUE))

if(length(fnFs) != length(fnRs)) stop("At least one sample is unpaired. Please check forward and reverse reads are present for all samples")

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#sample.namesR <- sapply(strsplit(basename(fnRs), "_"), `[`, 1) # to check if reverse reads have same sample names as forward
#write.csv(sample.names,"outputPrelim/sampleNames.csv") # print list of sample names to check if they match map file
table_map <- data.frame(read.csv("raw/table_mapComplete.csv", header = TRUE, row.names = 1, check.names=FALSE)) # reads csv file into data.frame with row names in column 1

save(fnFs,file=("outputPrelim/output_fnFs.RData")) # Save in a .RData file 
save(fnRs,file=("outputPrelim/output_fnRs.RData")) # Save in a .RData file 

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
save(filtFs,file=("outputPrelim/output_filtFs.RData")) # Save in a .RData file
save(filtRs,file=("outputPrelim/output_filtRs.RData")) # Save in a .RData file

#27F primer used for v1v3 sequencing=20bp "AGRGTTYGATYMTGGCTCAG"
#515R primer used for v1v3 sequencing=19bp "TBACCGCGGCTGCTGGCAC"
filtered <- filterAndTrim(fnFs, filtFs, 
                          fnRs, filtRs,
                          trimLeft=20, trimRight=5,
                          maxLen=530, minLen = 200,
                          maxN=0, maxEE=c(3,5), truncQ=2, rm.phix=TRUE,
                          compress=TRUE, verbose=TRUE, multithread=TRUE) # filtering values set here are specific to this study. Modify as appropriate when sequencing other regions.
save(filtered,file=("outputPrelim/output_filtered.RData")) # Save .RData file

### Subset filtFs and filtRs to include files with > 0 reads (some files may have been emptied during error analysis) 
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]
save(filtFs,file=("outputPrelim/output_filtFs.RData")) # Save .RData file
save(filtRs,file=("outputPrelim/output_filtRs.RData")) # Save .RData file 

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
save(errF,file=("outputPrelim/output_errF.RData")) # Save .RData file
save(errR,file=("outputPrelim/output_errR.RData")) # Save .RData file
#plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
save(dadaFs,file=("outputPrelim/output_dadaFs.RData")) # Save in a .RData file
save(dadaRs,file=("outputPrelim/output_dadaRs.RData")) # Save in a .RData file
save(mergers,file=("outputPrelim/output_mergers.RData")) # Save in a .RData file

seqtab <- makeSequenceTable(mergers)
# table(nchar(getSequences(seqtab))) # use to check distribution of sequence lengths
# hist(nchar(getSequences(seqtab)), main="ASV lengths in bp") # make histogram of ASV lengths
seqtab.all <- seqtab[,nchar(colnames(seqtab)) %in% 250:510] # these cut-off values have been chosen specific to this study. Modify as appropriate when sequencing other regions.
seqtab.noBim <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE)
save(seqtab.noBim,file=("outputPrelim/output_seqtabNoBim.RData")) # Save in a .RData file


getNreads <- function(x) sum(getUniques(x))
track <- cbind(filtered, sapply(dadaFs, getNreads), sapply(dadaRs, getNreads), sapply(mergers, getNreads), rowSums(seqtab), rowSums(seqtab.noBim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "tabled", "nonBim")
track

#####################################################################
########## Transfer data into Phyloseq ##############################
##### This preliminary Phyloseq object will be used to trim #########
##### singletons from the dataset before assigning taxonomy #########
#####################################################################

OTUraw = otu_table(seqtab.noBim, taxa_are_rows=FALSE) # assigns ASV table from DADA output
MAP = sample_data(table_map) # reads csv file into data.frame with row names in column 1
phyW = phyloseq(OTUraw, MAP) # prepare phyloseq object
save(phyW, file=("outputPrelim/physeq_raw.RData"))
save(OTUraw,file=("outputPrelim/output_otuRaw.RData")) # Save the phyloseq data object in a .RData file 

taxRaw <- assignTaxonomy(OTUraw, "~/Masters/silva_nr_v132_train_set_RossModFresh.fa", multithread=TRUE) # assign taxonomy to trimmed dataset
save(taxRaw,file=("outputPrelim/output_taxRaw.RData")) # Save the phyloseq data object in a .RData file 

#####################################################################
########## Transfer data into Phyloseq ##############################
#####################################################################
#load("outputPrelim/output_otuRaw.RData")
#load("outputPrelim/output_taxRaw.RData")
#table_map <- data.frame(read.csv("raw/table_map.csv", header = TRUE, row.names = 1, check.names=FALSE)) # reads csv file into data.frame with row names in column 1

OTU = otu_table(OTUraw, taxa_are_rows=FALSE) # assigns ASV table from trimmed DADA output
TAX = tax_table(taxRaw) # assigns taxonomy table
MAP = sample_data(table_map) # assigns metadata table
physeq = phyloseq(OTU, TAX, MAP) # prepare phyloseq object

### add data to taxonomy table ###
table_spp <- data.frame(read.csv("raw/table_SppDescription.csv", header = TRUE, check.names=FALSE)) # reads csv file into data.frame with row names in column 1
bind.asv <- data.frame(cbind(tax_table(physeq), paste0("ASV", seq(ntaxa(physeq))))) # add column listing ASVs numerically
bind.seq <- cbind(bind.asv, row.names(tax_table(physeq))) # add column listing ASV sequences
bind.o <- cbind(bind.seq, seq(ntaxa(physeq))) # adding a column for sorting table later
bind.length <- cbind(bind.o, nchar(bind.o[,10])) # add a column listing length of sequences
bind.count <- cbind(bind.length, taxa_sums(physeq)) # add a column listing total reads for ASV
bind.prev <- cbind(bind.count, rowSums(t(otu_table(physeq)) != 0)) # add a column listing total number of samples in which ASV appears
bind.prevSample <- cbind(bind.prev, rowSums(t(otu_table(subset_samples(physeq, Control%in%c("sample","positive")))) != 0)) # add a column listing total number of samples (excluding controls) in which ASV appears
bind.CP <- cbind(bind.prevSample, as.numeric(bind.prevSample[,13]) / as.numeric(bind.prevSample[,14])) # calculate ~average reads/sample
colnames(bind.CP) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain","ASV",
                       "Sequence","Sort","Length","Count","Prevalence","PrevSamples","CP") # rename columns
bind.sppData2 <- merge(bind.CP, table_spp[,c("pathogenStatus","commonSource","commonSourceSpecific","atmosphere","Genus","Species")], by=c("Genus","Species"), all.x=TRUE) # add extra info about species
bind.sppData <- bind.sppData2[order(as.numeric(bind.sppData2$Sort)), ] # put back in original order (merge function above scrambles order)
bind.sppData$atmosphere <- with(bind.sppData, ifelse(Genus=="Escherichia/Shigella" & is.na(Species), atmosphere=="aerobic",
                                                     ifelse(Genus=="Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium" & is.na(Species), atmosphere=="aerobic",
                                                            atmosphere))) # correct a few atmosphere assignments
bind.mod <- bind.sppData[,c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain","ASV","Length","Count",
                            "Prevalence","PrevSamples","CP","pathogenStatus","commonSource","commonSourceSpecific","atmosphere","Sequence")] # rename columns
rownames(bind.mod) <- bind.mod$ASV # reassign rownames to ASV names
TAX2 = tax_table(as.matrix(bind.mod)) # define new taxonomy table
taxa_names(physeq) <- tax_table(TAX2)[,"ASV"]
physeqR = phyloseq(otu_table(physeq), TAX2, MAP) # prepare phyloseq object modified with ASV numbers

#######################################
##### add columns to mapping data #####
#######################################

sample_data(physeqR) <- sample_data(physeqR) %>%
  cbind(sampleID = sample_names(physeqR), # sample names
        LibrarySize = sample_sums(physeqR), # total reads count per sample
        BacteriaReads = sample_sums(subset_taxa(physeqR, Kingdom=="Bacteria")), # total bacteria reads per sample
        AnimalReads = sample_sums(subset_taxa(physeqR, Kingdom=="Animalia")) # total bacteria reads per sample
  )
sample_data(physeqR) <- sample_data(physeqR) %>%
  cbind(RatioReads = sample_data(physeqR)$BacteriaReads / sample_data(physeqR)$AnimalReads, # ratio of bacteria reads:animal reads per sample
        BacteriaPercent = sample_data(physeqR)$BacteriaReads / sample_data(physeqR)$LibrarySize # percent of total reads that are bacteria
) 


########################################################################################
##### save phyloseq object and its components for easy access in future processing #####
########################################################################################
save(physeqR,file=("outputPrelim/physeqR.RData")) # Save the phyloseq data object in a .RData file 
write.csv(tax_table(physeqR),"outputPrelim/table_tax.csv") # Save taxonomy table as .csv
write.csv(otu_table(physeqR),"outputPrelim/table_otu.csv") # Save ASV table as .csv
write.csv(data.frame(sample_data(physeqR)),"outputPrelim/table_mapModified.csv") # Save modified sample data as .csv

##########################################################################################################
##### You have now successfully imported all the raw data into R and performed preliminary clean-up. #####
##### You are ready to proceed with further analysis                                                 #####
##########################################################################################################

##### All code written by Emily Ann McClure. No AI was used at any stage in writing or editing this code. #####
