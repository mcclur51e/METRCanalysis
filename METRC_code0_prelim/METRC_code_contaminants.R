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
##### first run prelimCombo
##### 1. A working repository (default "~/Desktop/Wannigan_METRC/")                                                                                   #####
##### 2. "outputPrelim/table_compileCulture.csv" : a data table containing compiled culture data                                     #####
##### 3. .csv file containing sample metadata listed with unique sampleIDs in first column (default "~/Desktop/Wannigan_METRC/raw/table_map.csv")     #####
##### 4. .fa file containing reference taxonomy database (default "~/Masters/silva_nr_v132_train_set_RossMod.fa")                                     #####                                   
###########################################################################################################################################################

#####################################################################
########## Processing Set −up #######################################
#####################################################################
# R version 4.5.2 (2025-10-31) -- "[Not] Part in a Rumble"
# All code written by Emily Ann McClure. No AI was used at any stage in writing or editing this code. 
# github.com/mcclure51e/
########## Call libraries for use ##########
# Packages from CRAN
library("data.table") # version 1.17.4
library("ggplot2") # version 4.0.0

# Packages from Bioconductor
library("decontam") # version 1.28.0, use for identifying contaminants
library("dplyr") # version 1.1.4 , allows use of '%>%' for chained functions
library("phyloseq") # version 1.52.0

# run code "prelimCombo"

########## Define master input and output locations ##########    
path <- "~/Desktop/Wannigan_METRC/" # change this to directory where you will be working
setwd(path) # set working directory
#dir.create("output") # create directory for output files
#dir.create("plots") # create directory to collect output plots

load("outputPrelim/physeq_reduced.RData") # load prepared phyloseq object, phyR.sdXtend
dt.cultCompile <- data.table(read.csv("outputPrelim/table_compileCulture.csv", header = TRUE, check.names = FALSE)) # Load culture data, dt.cultCompile

#############################################################
########## Identify and remove likely contaminants ##########
#############################################################
### Calculate ASVs to remove based on counts in the dataset
phyR.neg2 <- subset_samples(phyR.sdXtend, Control%in%c("negExtract","negSeq","negative")) # subset data to negative controls
phyR.neg <- subset_taxa(phyR.sdXtend, taxa_sums(phyR.neg2) > 1 & Kingdom=="Bacteria") # subset to taxa present in negative controls
df.neg <- data.frame(otu_table(phyR.neg)) 
phyR.taxTrim <- subset_taxa(phyR.sdXtend, taxa_sums(phyR.sdXtend) > quantile(as.matrix(df.neg[df.neg>0]), na.rm = T, probs = c(0.7))) # keep taxa present at counts above the 70% quantile from negative controls (i.e. have to have more reads than most ASVs from negative control samples)

### identify contaminants using decontam
ps <- subset_taxa(phyR.sdXtend, Kingdom=="Bacteria")
sample_data(ps)$is.neg <- sample_data(ps)$Control == "negative"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)
#table(contamdf.prev$contaminant)
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Control%in%c("negative","negEtract","negSeq"), ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Control=="sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#phyR.noContam <- subset_taxa(phyR.taxTrim, !ASV%in%rownames(subset(df.pa, contaminant=="TRUE")))
#phyR.noContamSamples <- subset_samples(phyR.noContam, sample_sums(phyR.noContam)>0 & Control%in%c("sample","positive"))

phyR.contam <- subset_taxa(phyR.sdXtend, ASV%in%rownames(subset(df.pa, contaminant=="TRUE")))
phyR.blank <- subset_taxa(phyR.contam, !Genus%in%c(dt.cultCompile$genus))

phyR.noContam <- subset_taxa(phyR.sdXtend, !ASV%in%tax_table(phyR.blank)[,c("ASV")])
phyR.noContamSamples <- subset_samples(phyR.noContam, Control%in%c("sample","positive"))

#############################################################
########## Remove replicates ################################
#############################################################
### Find replicates
agg.maxReads <- setNames(aggregate(x=sample_data(phyR.noContamSamples)$LibrarySizeRed, 
                                   by=list(sample_data(phyR.noContamSamples)$study_id, sample_data(phyR.noContamSamples)$timepoint),
                                   FUN=max), 
                         c("study_id","timepoint","MaxLibrarySizeRed"))
### Trim dataset to remove replicates and samples with 0 reads

ps.noContamPos <- psmelt(phyR.noContamSamples)

ps.noContamPosMax <- merge(ps.noContamPos,
      agg.maxReads,
      by=c("study_id","timepoint"),
      all=TRUE)

ps.temp <- subset(ps.noContamPosMax, LibrarySizeRed==MaxLibrarySizeRed)

ls.seqPosBase <- unique(subset(ps.temp, timepoint=="baseline")$sampleID)
ls.seqPosFU <- unique(subset(ps.temp, timepoint=="follow-up")$sampleID)

phyD.seqPosBase <- subset_samples(phyR.noContam, timepoint=="baseline" & sampleID%in%ls.seqPosBase)  
phyD.seqPosFU <- subset_samples(phyR.noContam, timepoint=="follow-up" & sampleID%in%ls.seqPosFU)  
phyD.seqPos <- merge_phyloseq(phyD.seqPosBase, phyD.seqPosFU)
phyD.seqPosBact <- subset_taxa(phyD.seqPos, taxa_sums(phyD.seqPos)>1 & Kingdom=="Bacteria")

###########################################################
########## Rarefaction curve to define failed samples #####
###########################################################
##### prepare data
### separate taxonomy table and convert columns to numbers for taxa selection later
dfD.seqPos <- data.frame(tax_table(phyD.seqPosBact)) # convert taxonomy table to data frame
dfD.seqPos <- dfD.seqPos %>% mutate_at(c("Length","Count","Prevalence","PrevSamples","CP"), as.numeric) # convert listed columns to numeric
sample_data(phyD.seqPosBact)$rareID <- paste0("xx", sample_data(phyD.seqPosBact)$LibrarySize, sample_data(phyD.seqPosBact)$LibrarySizeRed, sample_data(phyD.seqPosBact)$study_id) # assign temporary sample names to begin with alpha character, because R
sample_names(phyD.seqPosBact) <- sample_data(phyD.seqPosBact)$rareID # change sample names to begin with alpha character, because R
set.seed(42) # random number here
psdata <- phyD.seqPosBact

### calculate rarefaction curves
calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- reshape2::melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    molten_alpha_diversity
  }
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  rarefaction_curve_data
}

rarefaction_curve_data <- calculate_rarefaction_curves(psdata, c('Observed', 'Shannon'), rep(c(1, 5, 10, 25, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500, 750, 1000, 2500, 5000, 1:100 * 10000), each = 10))
summary(rarefaction_curve_data)

rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(psdata)), by.x = 'Sample', by.y = 'row.names')

### plot rarefaction curve
pRare <- ggplot(data = subset(rarefaction_curve_data_summary_verbose, Measure=="Observed"),
                mapping = aes(x = Depth, y = Alpha_diversity_mean,
                              ymin = Alpha_diversity_mean - Alpha_diversity_sd, ymax = Alpha_diversity_mean + Alpha_diversity_sd,
                              group = Sample)) +
  geom_line() + 
  geom_pointrange() + 
  facet_wrap(timepoint~outcomeSimple, scales = 'free_y') + 
  scale_x_log10() +
  theme_bw() +
  theme(text=element_text(size=12),axis.text.x = element_text(),legend.position="none",panel.spacing = unit(0, "lines"))
pRare
ggsave(pRare, filename="Bioburden/Plots/METRCburden_rarefaction.pdf", dpi="retina",width=9,height=6,units="in") # Save figure to .pdf file

ls.depth <- unique(rarefaction_curve_data$Depth) # these were set in line 337, when rarefaction_curve_data was first defined
dt.rcData <- data.table(rarefaction_curve_data) # convert to data table, because dcast

df.rare <- data.frame(dcast(subset(dt.rcData, Measure=="Observed"), Sample ~ Depth, value.var = "Alpha_diversity", fun.aggregate=median)) # prepare table with all depths for each sample in different column
ls.depthNames <- paste0("depth", ls.depth) # create list of slope depths
colnames(df.rare) <- c("Sample", paste0("depth",ls.depth))[1:ncol(df.rare)] # rename columns with rarefation depth

df.rareStretch <- merge(df.rare, data.frame(sample_data(phyD.seqPosBact))[,c("BacteriaReads","LibrarySize","LibrarySizeRed","Observed")], by.x = 'Sample', by.y = 'row.names') # add columns to table
ls.slopeNames <- paste0("slope",ls.depth[2:length(ls.depth)]) # create list of slope depths
ls.addNames <- paste0("add",ls.depth[2:length(ls.depth)]) # create list of column names for ASVs added at each depth

# add columns calculating slope between sampling depths
for(d in ls.depthNames[2:length(ls.depthNames)]){
  for(r in rownames(df.rareStretch)){
    s <- which(ls.depthNames==d)-1
    df.rareStretch[r,ls.slopeNames[s]] <- (df.rareStretch[r,d] - df.rareStretch[r,ls.depthNames[s]])/(ls.depth[s+1] - ls.depth[s])
  }
}
# add columns calculating how many ASVs are added between sampling depths
for(d in ls.depthNames[2:length(ls.depthNames)]){
  for(r in rownames(df.rareStretch)){
    s <- which(ls.depthNames==d)-1
    df.rareStretch[r,ls.addNames[s]] <- df.rareStretch[r,ls.slopeNames[s]] * df.rareStretch[r,"Observed"] 
  }
}

df.rareAdd <- data.frame(apply(df.rareStretch[, ls.addNames], 2, FUN=range, na.rm=TRUE)[2,]) # find derivative of rarefaction curves
colnames(df.rareAdd) <- "ASVsAdded" # rename column
df.rareAdd$depth <- ls.depth[2:(nrow(df.rareAdd)+1)] # label depth at which the derivative was calculated


########## DEFINE RAREFACTION LEVEL HERE ##########
### choose rarefaction cut-off as step where 0.5% or fewer taxa are added
df.rareAddLow <- subset(df.rareAdd, round(df.rareAdd[,1],1)<=(ntaxa(phyD.seqPosBact)*.005))
#df.rareAddLow
var.rareMin <- df.rareAddLow[1,2] - 1 # define rarefaction minimum as: (the smallest depth with slope<0.1) - 1

### separate samples by whether they meet rarefaction cut-off
phyD.seqPosBase2 <- subset_samples(phyD.seqPosBase, BacteriaReadsRed>var.rareMin)
phyD.seqPosFU2 <- subset_samples(phyD.seqPosFU, BacteriaReadsRed>var.rareMin)
phyD.posNGS <- merge_phyloseq(phyD.seqPosBase2, phyD.seqPosFU2)
sample_data(phyD.posNGS)$seqOutcome <- "positive" 

phyD.seqNegBase2 <- subset_samples(phyD.seqPosBase, BacteriaReadsRed<=var.rareMin)
phyD.seqNegFU2 <- subset_samples(phyD.seqPosFU, BacteriaReadsRed<=var.rareMin)
phyD.negNGS <- merge_phyloseq(phyD.seqNegBase2, phyD.seqNegFU2)
sample_data(phyD.negNGS)$seqOutcome <- "negative" 

sample_data(phyD.posNGS)$NGSresult <- "positive"
sample_data(phyD.negNGS)$NGSresult <- "negative"

save(phyD.posNGS,file=("outputPrelim/physeq_posNGS.RData")) # save the phyloseq data object in a .RData file 
save(phyD.negNGS,file=("outputPrelim/physeq_negNGS.RData")) # save the phyloseq data object in a .RData file 







