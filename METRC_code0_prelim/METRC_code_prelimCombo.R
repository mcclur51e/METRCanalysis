##########################################################################################################
##### DESCRIPTION                                                                                    #####
##### Code for preliminary processing of NGS and culture data                                        #####
##########################################################################################################
##### This code requires as input:                                                                   #####
##### 1. "outputPrelim/physeqR.RData : a phyloseq object containing paired and trimmed 16S reads     #####
##### 2. "outputPrelim/table_compileCulture.csv" : a data table containing compiled culture data     #####
##########################################################################################################

#####################################################################
########## Processing Set −up #######################################
#####################################################################
# R version 4.5.2 (2025-10-31) -- "[Not] Part in a Rumble"
# All code written by Emily Ann McClure. No AI was used at any stage in writing or editing this code. 
# github.com/mcclure51e/
########## Call libraries for use ##########
# Packages from CRAN
library("data.table") # version 1.17.4
library("dplyr") # version 1.1.4 , allows use of '%>%' for chained functions
library("reshape") # version 0.8.9 , use for cast functions

# Packages from Bioconductor
library("phyloseq") # version 1.52.0

########## Call functions for use ##########
source("~/Masters/color_palettes.R",local=TRUE) # load common color palettes for use

########## Define master input and output locations ##########    
path <- "~/Desktop/Wannigan_METRC/" # change this to directory where you will be working
setwd(path) # set working directory
#dir.create("output") # create directory for output files
#dir.create("plots") # create directory to collect output plots

### Load preliminary culture data
dt.cultCompile <- data.table(read.csv("outputPrelim/table_compileCulture.csv", header = TRUE, check.names = FALSE)) # reads csv file into data.table with no row names
df.cultCount <- setNames(cast(data.frame(dt.cultCompile), study_id+timepoint~., value="abundance", sum),
                         c("study_id","timepoint","isolatesPerSample")) # count number of cultured isolates per sample

### Load sequencing data
load("outputPrelim/physeqR.RData") # load the phyloseq data object from a .RData file 
phyR <- subset_taxa(physeqR, taxa_sums(physeqR)>1) # keep only taxa with more than one read in the dataset

########## Find study ids by sample groups ##########
#ls.cultIDs <- unique(dt.cultCompile$study_id) # study_ids from culture data
ls.cultFUIDs <- unique(subset(dt.cultCompile, timepoint=="follow-up")$study_id) # study_ids from culture follow-up data
ls.cultBaseIDs <- unique(subset(dt.cultCompile, timepoint=="baseline")$study_id) # study_ids from culture baseline data

ls.seqFUIDs <- unique(sample_data(subset_samples(phyR, Control%in%c("sample","positive") & SampleType=="sequencing" & timepoint=="follow-up" & sample_sums(phyR)>0))$study_id) # study_ids with reads from sequencing at follow-up
ls.seqBaseIDs <- unique(sample_data(subset_samples(phyR, Control%in%c("sample","positive") & SampleType=="sequencing" & timepoint=="baseline" & sample_sums(phyR)>0))$study_id) # study_ids with reads from sequencing at baseline

ls.baseIDs <- unique(union(ls.cultBaseIDs, ls.seqBaseIDs)) # study_ids sampled at baseline
ls.FUIDs <- unique(union(ls.cultFUIDs, ls.seqFUIDs)) # study_ids sampled at follow-up
ls.allIDs <- c(ls.baseIDs, ls.FUIDs)

########## Describe patient outcomes ##########
### Outcomes listed by METRC metadata
ls.outDeepInfection <- unique(sample_data(subset_samples(phyR, outcomeText=="deepInfection" & SampleType!="none"))$study_id) # list study_id resulting in deep infection
ls.outFlapFailure <- unique(sample_data(subset_samples(phyR, outcomeText=="flapFailure" & SampleType!="none"))$study_id) # list study_id resulting in flap failure
ls.outNonUnion <- unique(sample_data(subset_samples(phyR, outcomeText=="nonUnion" & SampleType!="none"))$study_id) # list study_id resulting in non-union
ls.outOther <- unique(sample_data(subset_samples(phyR, outcomeText=="Other" & SampleType!="none"))$study_id) # list study_id resulting in "other"
ls.outListed <- c(ls.outDeepInfection, ls.outFlapFailure, ls.outNonUnion, ls.outOther) # list study_id with METRC listed outcome
ls.healed <- unique(setdiff(ls.baseIDs, union(ls.outListed, ls.FUIDs))) # list study_id resulting in healing (no follow-up samples or outcome listed)
ls.outUnknown <- unique(setdiff(ls.allIDs, union(ls.outListed, ls.healed))) # list study_id with unknown outcome (follow-up samples, but no listed outcome)
sample_data(phyR)$outcomeText <- with(sample_data(phyR), ifelse(study_id%in%ls.healed, "healed", outcomeText)) # Update outcomeText in sample data tables

######################################################################
########## Identify subset of ASVs to represent human reads ##########
######################################################################
### NOTE: this must be done before trimming ASVs for length
df.taxHomo <- data.frame(tax_table(subset_taxa(phyR, Phylum=="Chordata"))) # list taxa that are assigned to primates
df.taxHomo <- df.taxHomo %>% mutate_at(c("Length","Count","Prevalence","PrevSamples","CP"), as.numeric) # convert listed columns to numeric
ls.homoReads <- subset(df.taxHomo, Count>nrow(df.taxHomo) & # must have more reads in the dataset than number of samples in the dataset
                         Prevalence > quantile(df.taxHomo$Prevalence, na.rm = T, probs = c(0.995)))$ASV  # keep ASVs with top .05% prevalence counts (present in the most samples)


####################################################################
########## Add human reads count to sample data #######################
####################################################################
sd.comp <- merge(data.frame(sample_data(phyR)), unique(dt.cultCompile[,c("study_id","timepoint","isolatesPerSample")]), by=c("study_id","timepoint"), all.x=TRUE)
rownames(sd.comp) <- sd.comp$sampleID
phyR.sdComp <- phyloseq(otu_table(phyR), tax_table(phyR), sample_data(sd.comp))
sample_data(phyR.sdComp)$isolatesPerSample <- with(sample_data(phyR.sdComp), ifelse(isolatesPerSample>0, isolatesPerSample, 0))
sample_data(phyR.sdComp) <- sample_data(phyR.sdComp) %>%
  cbind(baseCultureResults = ifelse(sample_data(phyR.sdComp)$study_id%in%ls.cultBaseIDs & sample_data(phyR.sdComp)$timepoint=="baseline", "positive",
                                    ifelse(!sample_data(phyR.sdComp)$study_id%in%ls.cultBaseIDs & sample_data(phyR.sdComp)$timepoint=="baseline","negative", "untested")), # culture results when sampled at baseline
        fuCultureResults = ifelse(sample_data(phyR.sdComp)$study_id%in%ls.cultFUIDs & sample_data(phyR.sdComp)$timepoint=="follow-up", "positive",
                                  ifelse(!sample_data(phyR.sdComp)$study_id%in%ls.cultFUIDs & sample_data(phyR.sdComp)$timepoint=="follow-up","negative", "untested")), # culture results when sampled at follow-up
        outcomeSimple = ifelse(sample_data(phyR.sdComp)$outcomeText%in%c("deepInfection","flapFailure","nonUnion"),"infection","unInfected"), # simplified outcome
        HomoReads = sample_sums(subset_taxa(phyR, ASV%in%ls.homoReads)) # add column listing H.sapiens read count to sample data
  )

###################################################################
########## Reduce dataset based on length of amplicon #############
###################################################################
tax_table(phyR.sdComp)[,"Prevalence"] <- as.numeric(tax_table(phyR.sdComp)[,"Prevalence"])
tax_table(phyR.sdComp)[,"Length"] <- as.numeric(tax_table(phyR.sdComp)[,"Length"])
val.minLength <- as.numeric(min(tax_table(subset_taxa(phyR.sdComp, Kingdom=="Bacteria" & !is.na(Phylum) & Prevalence>1 & Length>400))[,c("Length")])) # find shortest bacterial amplicon present as at least 10 reads in dataset
val.maxLength <- as.numeric(max(tax_table(subset_taxa(phyR.sdComp, Kingdom=="Bacteria" & !is.na(Phylum) & Prevalence>1 & Length<600))[,c("Length")])) # find longest bacterial amplicon present as at least 10 reads in dataset
phyR.red <- subset_taxa(phyR.sdComp, Length>=val.minLength & Length<=val.maxLength) # trim sequences above or below the lengths of bacterial reads

####################################################################
########## Modify sample_data (reduced ASVs) #######################
####################################################################
phyR.sdXtend <- phyR.red
sample_data(phyR.sdXtend) <- sample_data(phyR.sdXtend) %>%
  cbind(LibrarySizeRed = sample_sums(phyR.sdXtend), # add column listing number of reads per sample
        BacteriaReadsRed = sample_sums(subset_taxa(phyR.sdXtend, Kingdom=="Bacteria")), # add column listing number of bacterial reads per sample
        estimate_richness(subset_taxa(phyR.sdXtend, Kingdom=="Bacteria"), split = TRUE, measures = NULL) # include alpha diversity calculations  
        )
sample_data(phyR.sdXtend) <- sample_data(phyR.sdXtend) %>%
  cbind(HomoPercent = ifelse(sample_data(phyR.sdXtend)$LibrarySizeRed>0, sample_data(phyR.sdXtend)$HomoReads / sample_data(phyR.sdXtend)$LibrarySize, 
                             ifelse(sample_data(phyR.sdXtend)$HomoReads>0, 1, 0)),
        BacteriaPercentRed = ifelse(sample_data(phyR.sdXtend)$LibrarySizeRed>0, sample_data(phyR.sdXtend)$BacteriaReadsRed / sample_data(phyR.sdXtend)$LibrarySizeRed, 0), # add column listing percent of total reads that are bacteria
        BHR = sample_data(phyR.sdXtend)$BacteriaReads / sample_data(phyR.sdXtend)$HomoReads, # add column listing ratio of bacterial to human reads
        BHRred = sample_data(phyR.sdXtend)$BacteriaReadsRed / sample_data(phyR.sdXtend)$HomoReads, # ratio of bacterial to human reads
        logBHR = log10(as.numeric(sample_data(phyR.sdXtend)$BacteriaReadsRed / sample_data(phyR.sdXtend)$HomoReads)), # add column listing the log of BHR of reduced dataset   
        outcomeTP = paste0(sample_data(phyR.sdXtend)$outcomeSimple, "\n", sample_data(phyR.sdXtend)$timepoint), # create new groups of result + timepoint
        allCultureResults = ifelse(sample_data(phyR.sdXtend)$baseCultureResults=="positive"|sample_data(phyR.sdXtend)$fuCultureResults=="positive", "positive",
                                   ifelse(sample_data(phyR.sdXtend)$baseCultureResults=="negative"|sample_data(phyR.sdXtend)$fuCultureResults=="negative", "negative", "untested")) # list of culture results
        #logBactHomo = log10(BHR), # add column listing the log of BHR
  )
sample_data(phyR.sdXtend) <- sample_data(phyR.sdXtend) %>%
  cbind(cultureTP = paste0(sample_data(phyR.sdXtend)$allCultureResults, "\n", sample_data(phyR.sdXtend)$timepoint) # create new groups of culture result + timepoint
)
save(phyR.sdXtend,file=("outputPrelim/physeq_reduced.RData")) # Save the phyloseq data object in a .RData file 

##########################################################################################################
##### You have now successfully pre-processed the combine (culture and sequencing) dataset. #####
##### You are ready to proceed with further analysis                                                 #####
##########################################################################################################

##### All code written by Emily Ann McClure. No AI was used at any stage in writing or editing this code. #####
