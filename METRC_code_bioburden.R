# R version 4.4.1 (2024-06-14) -- "Race for Your Life"

########## Call libraries for use ##########
# Packages from CRAN
#library("ggplot2") # version 3.4.2

# Packages from Bioconductor
library("phyloseq") # version 1.44.0
library("microViz") # version 0.10.10
library("ggplot2") # version 3.4.2
library("pheatmap") # version 1.0.12 
library("tidyr") # version 1.3.0 
library("reshape2") # version 1.4.4 , use for melt and cast functions
library("data.table") # version 1.14.8

########## Call functions for use ##########
source("~/Masters/color_palettes.R",local=TRUE) # load common color palettes for use

##################################
########## End preamble ##########
##################################                     

########## Define master input and output locations ##########    
path <- "~/Desktop/Wannigan_METRC/" # change this to directory where you will be working
setwd(path) # set working directory
#dir.create("output") # create directory for output files
#dir.create("plots") # create directory to collect output plots



load("outputPrelim/output_mergers.RData")
seqtab <- makeSequenceTable(mergers)

# table(nchar(getSequences(seqtab))) # use to check distribution of sequence lengths
#seqtab.all <- seqtab[,nchar(colnames(seqtab)) %in% 443:529]
seqtab.all <- seqtab
seqtab.noBim <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE)

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

taxRaw <- assignTaxonomy(OTUraw, "~/Masters/silva_nr_v132_train_set_RossMod.fa", multithread=TRUE) # assign taxonomy to trimmed dataset
#save(taxRaw,file=("outputPrelim/output_taxRaw.RData")) # Save the phyloseq data object in a .RData file 


#####################################################################
########## Transfer data into Phyloseq ##############################
#####################################################################
#load("outputPrelim/output_otuRaw.RData")
#load("outputPrelim/output_taxRaw.RData")
#table_map <- data.frame(read.csv("raw/table_map.csv", header = TRUE, row.names = 1, check.names=FALSE)) # reads csv file into data.frame with row names in column 1

OTU = otu_table(OTUraw, taxa_are_rows=FALSE) # assigns ASV table from trimmed DADA output
TAX = tax_table(taxRaw) # assigns taxonomy table
MAP = sample_data(table_map) # reads csv file into data.frame with row names in column 1
physeq = phyloseq(OTU, TAX, MAP) # prepare phyloseq object

### create and assign ASV numbers to consensus sequences for easier reference ###
ASV <- paste0("ASV", seq(ntaxa(physeq))) 
Sequence <- row.names(tax_table(physeq)) 
bind.length <- cbind(tax_table(physeq),nchar(rownames(tax_table(physeq))))
bind.count <- cbind(bind.length,taxa_sums(physeq))
bind.asv <- cbind(bind.count,ASV)
bind.seq <- cbind(bind.asv,Sequence) 

bind.seq <- as.data.frame(bind.seq)

TAX2 = tax_table(as.matrix(bind.seq)) # define new taxonomy table
colnames(TAX2) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species","Strain","Length","Count","ASV","Sequence")
physeqR = phyloseq(OTU, TAX2, MAP) # prepare phyloseq object modified with ASV numbers
taxa_names(physeqR) <- tax_table(physeqR)[,"ASV"]

### Load culture data
df.baseCult <- data.frame(read.csv("outputPrelim/table_cultureBase.csv", header = TRUE, check.names=FALSE))[,-1] # reads csv file into data.frame with row names in column 1
df.fuCult <- data.frame(read.csv("outputPrelim/table_cultureFU.csv", header = TRUE, check.names=FALSE))[,-1] # reads csv file into data.frame with row names in column 1
df.baseCult$timepoint <- "baseline"
df.fuCult$timepoint <- "follow-up"
df.cult <- rbind(df.baseCult, df.fuCult)
df.cult$genus <- with(df.cult, ifelse(genus=="Propionibacterium", "Cutibacterium", genus)) # correct naming convention
write.csv(df.cult, "outputPrelim/table_cultureResultsRaw.csv")

########## Find study ids by sample groups ##########
ls.cultIDs <- unique(df.cult$study_id) # list all study IDs (patient IDs) from culture data
ls.cultFUIDs <- unique(subset(df.cult, timepoint=="follow-up")$study_id) # list all study IDs (patient IDs) from culture follow-up data
ls.cultBaseIDs <- unique(subset(df.cult, timepoint=="baseline")$study_id) # list all study IDs (patient IDs) from culture baseline data

ls.seqIDs <- unique(sample_data(subset_samples(physeqR,Control=="sample"))$study_id) # study_ids with reads from sequencing
ls.seqFUIDs <- unique(sample_data(subset_samples(physeqR,Control=="sample" & SampleType=="sequencing" & timepoint=="follow-up"))$study_id) # study_ids with reads from sequencing at follow-up
ls.seqBaseIDs <- unique(sample_data(subset_samples(physeqR,Control=="sample" & SampleType=="sequencing" & timepoint=="baseline"))$study_id) # study_ids with reads from sequencing at baseline

ls.FUIDs <- unique(union(ls.cultFUIDs, ls.seqFUIDs))
ls.FUoverlapIDs <- unique(intersect(ls.cultFUIDs,ls.seqFUIDs))
ls.baseIDs <- unique(union(ls.cultBaseIDs, ls.seqBaseIDs))
ls.baseOverlapIDs <- unique(intersect(ls.cultBaseIDs,ls.seqBaseIDs))
ls.allFourTests <- intersect(ls.FUIDs,ls.baseIDs)
ls.allIDs <- unique(union(ls.cultIDs, ls.seqIDs))

########## Describe patient outcomes ##########
### Outcomes listed by METRC metadata
ls.outAmputation <- unique(subset(table_map, outcomeText=="Amputation" & SampleType!="none")$study_id) # list study_id resulting in amputation
ls.outBoneGraft <- unique(subset(table_map, outcomeText=="boneGraft" & SampleType!="none")$study_id) # list study_id resulting in bone graft
ls.outDeepInfection <- unique(subset(table_map, outcomeText=="deepInfection" & SampleType!="none")$study_id) # list study_id resulting in deep infection
ls.outFlapFailure <- unique(subset(table_map, outcomeText=="flapFailure" & SampleType!="none")$study_id) # list study_id resulting in flap failure
ls.outNonUnion <- unique(subset(table_map, outcomeText=="nonUnion" & SampleType!="none")$study_id) # list study_id resulting in non-union
ls.outOther <- unique(subset(table_map, outcomeText=="Other" & SampleType!="none")$study_id) # list study_id resulting in "other
ls.outListed <- unique(union(union(union(ls.outAmputation,ls.outBoneGraft),union(ls.outDeepInfection,ls.outFlapFailure)),union(ls.outNonUnion,ls.outOther))) # list study_id with METRC listed outcome

### Outcomes inferred from absence of info in METRC metadata
ls.healed <- unique(setdiff(ls.baseIDs,union(ls.outListed,ls.FUIDs))) # list study_id resulting in healing (no follow-up samples or outcome listed)
ls.outUnknown <- unique(setdiff(ls.allIDs,union(ls.outListed,ls.healed))) # list study_id with unknown outcome (follow-up samples, but no listed outcome)
# NOTE: some samples have an outcome listed without a follow-up sample result

### Table of outcomes listed by study_id
df.outcome <- unique(table_map[,c("study_id","SampleType","outcomeText")]) # create simple table listing outcome by study_id
df.outcome <- subset(df.outcome, !is.na(study_id)) # remove samples with no study_id (controls)
df.outcome$outcomeText <- with(df.outcome, ifelse(study_id%in%ls.healed, "healed",
                                                  ifelse(study_id%in%ls.outUnknown, "unk",
                                                         outcomeText))) # add "healed" value to outcomeText
write.csv(data.frame(df.outcome), "output/table_outcomeText.csv")

### Update outcomeText in sample data tables
df.cleanNew <- data.frame(sample_data(physeqR))
df.cleanNewTrim <- df.cleanNew[ , !(names(df.cleanNew)%in%c("outcomeText"))]
df.cleanNewTrim$sampleID <- rownames(df.cleanNewTrim)
df.cleanNewOutcome <- merge(data.frame(df.cleanNewTrim), df.outcome[,c("study_id","outcomeText")], by="study_id", all.x=TRUE)
rownames(df.cleanNewOutcome) <- df.cleanNewOutcome$sampleID
phyT.cleanNew2 <- phyloseq(tax_table(physeqR), otu_table(physeqR),sample_data(df.cleanNewOutcome))
#phyT.cleanNew <- subset_samples(phyT.cleanNew2, !sample_names(phyT.cleanNew2)%in%c("1388","1388-1","1673","1963","2427")) # remove empty replicates
#save(phyT.cleanNew,file=("outputPrelim/physeq_clean.RData")) # Save the phyloseq data object in a .RData file 

### add columns to mapping data
df.map <- as.data.frame(sample_data(phyT.cleanNew2))
df.map$study_id_tp <- with(df.map, ifelse(timepoint=="baseline",paste(study_id,"base"),
                                          ifelse(timepoint=="follow-up",paste(study_id,"fu"),
                                                 paste(Control,I5_Index_ID, I7_Index_ID))))
df.map$sampleID <- rownames(df.map) # sample names
df.map$LibrarySize <- sample_sums(physeqR) # total reads count per sample
df.map$AnimalReads <- sample_sums(subset_taxa(physeqR, Kingdom=="Animalia")) # total animal reads per sample
df.map$BacteriaReads <- sample_sums(subset_taxa(physeqR, Kingdom=="Bacteria")) # total bacteria reads per sample
df.map$ChloroplastReads <- sample_sums(subset_taxa(physeqR, Class=="Oxyphotobacteria")) # total chloroplast reads per sample
df.map$RatioReads <- df.map$BacteriaReads / df.map$AnimalReads # ratio of bacteria reads:animal reads per sample
df.map$BacteriaPercent <- df.map$BacteriaReads / df.map$LibrarySize # percent of total reads that are bacteria
df.map$AnimalPercent <- df.map$AnimalReads / df.map$LibrarySize # percent of total reads that are animal
df.map$StenotrophomonasReads <- sample_sums(subset_taxa(physeqR, Genus=="Stenotrophomonas")) # total animal reads per sample
df.map$DelftiaReads <- sample_sums(subset_taxa(physeqR, Genus=="Delftia")) # total animal reads per sample
df.map$StaphReads <- sample_sums(subset_taxa(physeqR, Genus=="Staphylococcus"))
df.map$coagNegStaphReads <- sample_sums(subset_taxa(physeqR, Genus=="Staphylococcus" & Species%in%c("capitis", "caprae", "epidermidis","hominis", "haemolyticus", "lugdunensis","pettenkoferi", "saprophyticus", "simulans", "warneri")))
df.map$ASV17 <- as.numeric(sample_sums(subset_taxa(physeqR,ASV=="ASV23")))
df.map$ASV69 <- as.numeric(sample_sums(subset_taxa(physeqR,ASV=="ASV22")))
df.map$ASV85 <- as.numeric(sample_sums(subset_taxa(physeqR,ASV=="ASV17")))
df.map$ASV66 <- as.numeric(sample_sums(subset_taxa(physeqR,ASV=="ASV32")))
df.map$ASV52 <- as.numeric(sample_sums(subset_taxa(physeqR,ASV=="ASV27")))
df.map$ASV84 <- as.numeric(sample_sums(subset_taxa(physeqR,ASV=="ASV84")))
df.map$ASV5 <- as.numeric(sample_sums(subset_taxa(physeqR,ASV=="ASV5")))
df.map$ASV46 <- as.numeric(sample_sums(subset_taxa(physeqR,ASV=="ASV46")))
df.map$ASV62 <- as.numeric(sample_sums(subset_taxa(physeqR,ASV=="ASV62")))
df.map$ASV157 <- as.numeric(sample_sums(subset_taxa(physeqR,ASV=="ASV157")))
df.map$ASV23bact <- as.numeric(df.map$ASV17/df.map$BacteriaReads)
df.map$ASV22bact <- as.numeric(df.map$ASV69/df.map$BacteriaReads)
df.map$ASV17bact <- as.numeric(df.map$ASV85/df.map$BacteriaReads)
df.map$ASV32bact <- as.numeric(df.map$ASV66/df.map$BacteriaReads)
df.map$ASV27bact <- as.numeric(df.map$ASV52/df.map$BacteriaReads)
df.map$ASV84bact <- as.numeric(df.map$ASV84/df.map$BacteriaReads)
df.map$ASV5bact <- as.numeric(df.map$ASV5/df.map$BacteriaReads)
df.map$ASV46bact <- as.numeric(df.map$ASV46/df.map$BacteriaReads)
df.map$ASV62bact <- as.numeric(df.map$ASV62/df.map$BacteriaReads)
df.map$ASV157bact <- as.numeric(df.map$ASV157/df.map$BacteriaReads)
df.map$ASVsum <- df.map$ASV17 + df.map$ASV69 + df.map$ASV85 + df.map$ASV66 + df.map$ASV52 + df.map$ASV84 + df.map$ASV5 + df.map$ASV46 + df.map$ASV62 + df.map$ASV157
df.map$ASVsum <- with(df.map, ifelse(ASVsum==0, 0.001, ASVsum))
df.map$ASVsumBact <- df.map$BacteriaReads/df.map$ASVsum
df.map$outcomeSimple <- with(df.map, ifelse(outcomeText%in%c("deepInfection","flapFailure","nonUnion"),"infection","un-infected"))
df.map$baseCultureResults <- with(df.map, ifelse(study_id%in%ls.cultBaseIDs, "Pos","Neg"))
df.map$FUCultureResults <- with(df.map, ifelse(study_id%in%ls.cultFUIDs, "Pos","Neg"))
df.map$allCultureResults <- with(df.map, ifelse(study_id%in%ls.cultBaseIDs & timepoint=="baseline", "Pos",
                                                ifelse(study_id%in%ls.cultFUIDs & timepoint=="follow-up", "Pos",
                                                       "Neg")))

physeqR = phyloseq(tax_table(physeqR),otu_table(physeqR),sample_data(df.map)) # return map to phyloseq object
save(physeqR,file=("outputPrelim/physeq_initial.RData")) # Save the phyloseq data object in a .RData file 
#load("outputPrelim/physeq_initial.RData")
table_map <- sample_data(physeqR) # 

## Removing samples with fewer than 3 OTUs in the entire sample, removing this from the dataset
physeqR <- subset_samples(physeqR, !sampleID%in%c("M4-1043-RT1","ORL-1043RT1","RIH-2097RDT","UOK-2427","2014",
                                                  "CMC-2134","RYD-2138","1848","1614","AGY-2287",
                                                  "UOK-2427","CMC-2134") & LibrarySize>10)
physeqR <- subset_samples(physeqR, !sample_names(physeqR)%in%c("1388","1388-1","1673","1963","2427")) # remove empty replicates
physeqR <- subset_samples(physeqR, !Control%in%c("negative","positive"))
physeqR <- subset_samples(physeqR, !sample_names(physeqR)%in%c("1043-RT1", "1517-POOL-5UL","1517-POOL-METRC6","CMC-1517RDT")) # remove extra replicates (hand curated via plot_bar)


### Temporary plot_bar


#plot_bar(subset_samples(physeqR, study_id%in%c("1043","1517")), fill="Order")




### Counting values for results section ###
# list samples with positive sequencing results (>10 reads is a VERY low threshold to set here)
ls.seqPosBaseIDs <- unique(sample_data(subset_samples(physeqR, sample_sums(physeqR)>10 & timepoint=="baseline"))$study_id)
ls.seqPosFUIDs <- unique(sample_data(subset_samples(physeqR, sample_sums(physeqR)>10 & timepoint=="follow-up"))$study_id)
# patients with culture + sequencing data at each timepoint
num.basePts <- length(ls.seqPosBaseIDs)
num.fuPts <- length(ls.seqPosFUIDs)
# patients with infections at each timepoint
num.baseInfect <- length(unique(subset(df.outcome, outcomeText%in%c("deepInfection","flapFailure","nonUnion") & study_id%in%ls.seqPosBaseIDs))$study_id)
num.fuInfect <- length(unique(subset(df.outcome, outcomeText%in%c("deepInfection","flapFailure","nonUnion") & study_id%in%ls.seqPosFUIDs))$study_id)
# patients culture positive at each timepoint
num.baseCultPos <- length(intersect(ls.cultBaseIDs,ls.seqPosBaseIDs))
num.fuCultPos <- length(intersect(ls.cultFUIDs, ls.seqPosFUIDs))
# print results
print(c("patients at baseline",num.basePts, 
        "patients at follow-up",num.fuPts,
        "infections after baseline", num.baseInfect,
        "infections at follow-up", num.fuInfect,
        "culture positive at baseline", num.baseCultPos,
        "culture positive at follow-up", num.fuCultPos
        ))

### Calculate stats for bacteria:human read ratios
phyR.cultPosBase <- subset_samples(physeqR, study_id%in%ls.cultBaseIDs & timepoint=="baseline" & sample_sums(physeqR)>10)
phyR.cultPosFU <- subset_samples(physeqR, study_id%in%ls.cultFUIDs & timepoint=="follow-up"& sample_sums(physeqR)>10)

med.baseCultPosRatio <- median(sample_data(phyR.cultPosBase)$ASVsumBact)
med.fuCultPosRatio <- median(sample_data(phyR.cultPosFU)$ASVsumBact)
sd.baseCultPosRatio <- sd(sample_data(subset_samples(phyR.cultPosBase, ASVsum>1))$ASVsumBact)
sd.fuCultPosRatio <- sd(sample_data(subset_samples(phyR.cultPosFU, ASVsum>1))$ASVsumBact)

phyR.cultNegBase <- subset_samples(physeqR, !study_id%in%ls.cultBaseIDs & timepoint=="baseline" & sample_sums(physeqR)>10)
phyR.cultNegFU <- subset_samples(physeqR, !study_id%in%ls.cultFUIDs & timepoint=="follow-up" & sample_sums(physeqR)>10)

med.baseCultNegRatio <- mean(sample_data(phyR.cultNegBase)$ASVsumBact)
med.fuCultNegRatio <- mean(sample_data(phyR.cultNegFU)$ASVsumBact)
sd.baseCultNegRatio <- sd(sample_data(subset_samples(phyR.cultNegBase, ASVsum>1))$ASVsumBact)
sd.fuCultNegRatio <- sd(sample_data(subset_samples(phyR.cultNegFU, ASVsum>1))$ASVsumBact)

med.cultPosRatio <- median(union(sample_data(phyR.cultPosBase)$ASVsumBact, sample_data(phyR.cultPosFU)$ASVsumBact))
med.cultNegRatio <- median(union(sample_data(phyR.cultNegBase)$ASVsumBact, sample_data(phyR.cultNegFU)$ASVsumBact))
sd.cultPosRatio <- sd(union(sample_data(subset_samples(phyR.cultPosBase, ASVsum>1))$ASVsumBact, sample_data(subset_samples(phyR.cultPosFU, ASVsum>1))$ASVsumBact))
sd.cultNegRatio <- sd(union(sample_data(subset_samples(phyR.cultNegBase, ASVsum>1))$ASVsumBact, sample_data(subset_samples(phyR.cultNegFU, ASVsum>1))$ASVsumBact))


print(c("patients at baseline with positive culture", nsamples(phyR.cultPosBase), 
        #"patients removed due to 0 human reads", nsamples(phyR.cultPosBaseRemoved),
        "patients at follow-up with positive culture", nsamples(phyR.cultPosFU),
        "bacteria:human ratio with positive culture at baseline", med.baseCultPosRatio,
        "stDev", sd.baseCultPosRatio,
        "bacteria:human ratio with positive culture at follow-up", med.fuCultPosRatio,
        "stDev", sd.fuCultPosRatio,
        "bacteria:human ratio with positive culture", med.cultPosRatio,
        "stDev", sd.cultPosRatio,
        "patients at baseline with negative culture", nsamples(phyR.cultNegBase), 
        "patients at follow-up with negative culture", nsamples(phyR.cultNegFU),
        "bacteria:human ratio with negative culture at baseline", med.baseCultNegRatio,
        "stDev", sd.baseCultNegRatio,
        "bacteria:human ratio with negative culture at follow-up", med.fuCultNegRatio,
        "stDev", sd.fuCultNegRatio,
        "bacteria:human ratio with negative culture", med.cultNegRatio,
        "stDev", sd.cultNegRatio
        
))


### Mann-Whitney U test
#physeqR2 <- subset_samples(physeqR, ASVsumBact!="Inf" & !is.na(ASVsumBact))
map.R <- data.frame(sample_data(physeqR))
# Is ratio significant between culture results
print("test if baseline culture pos/neg bacteria:human ratio is significant")
wilcox.test(map.R$ASVsumBact ~ map.R$baseCultureResults, data = subset(map.R, timepoint=="baseline"), exact = FALSE)
print("test if follow-up culture pos/neg bacteria:human ratio is significant")
wilcox.test(map.R$ASVsumBact ~ map.R$FUCultureResults, data = subset(map.R, timepoint=="follow-up"), exact = FALSE)
print("test if all culture pos/neg bacteria:human ratio is significant")
wilcox.test(map.R$ASVsumBact ~ map.R$allCultureResults, data = map.R, exact = FALSE)


library(ggpubr)
plot_richness(subset_samples(physeqR, timepoint=="baseline"), x="baseCultureResults", measures=c("Chao1", "Shannon")) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test")

plot_richness(subset_samples(physeqR, timepoint=="follow-up"), x="baseCultureResults", measures=c("Chao1", "Shannon")) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test")

df.rich <- estimate_richness(physeqR, split = TRUE, measures = NULL)
df.mapXtend <- cbind(sample_data(physeqR),df.rich)

pDot.shannon <- ggplot(data=subset(df.mapXtend, Shannon>0 & ASVsumBact>0.01), aes(x=ASVsumBact, y=Shannon, color=allCultureResults)) +
  geom_point() +
  #facet_wrap(sex~indole) +
  #scale_y_log10() +
  scale_x_log10() +
  #ylim(0,100) +
  theme_bw() +
  #ggtitle("Library Size") +
  theme(text=element_text(size=12),axis.text.x = element_text(face="italic", angle=45, hjust=1),panel.spacing = unit(0, "lines")) +
  scale_color_manual(values = pal.CB)
pDot.shannon

df.mapXtend$ASVsumBactLog <- log10(df.mapXtend$ASVsumBact)

lm.shannon = lm(ASVsumBactLog~Shannon, data = subset(df.mapXtend, ASVsumBact>0.01)) #Create the linear regression
summary(lm.shannon) #Review the results

pDot.temp <- ggplot(data=subset(df.mapXtend, Shannon>0 & ASVsumBact>0.01), aes(x=LibrarySize, y=ASVsumBact, color=allCultureResults)) +
  geom_point() +
  #facet_wrap(sex~indole) +
  #scale_y_log10() +
  #scale_x_log10() +
  #ylim(0,100) +
  theme_bw() +
  #ggtitle("Library Size") +
  theme(text=element_text(size=12),axis.text.x = element_text(face="italic", angle=45, hjust=1),panel.spacing = unit(0, "lines")) +
  scale_color_manual(values = pal.CB)
pDot.temp






### culture alpha diversity

alpha.cult <- unique(dt.cultCompile[,c("study_id","isolatesPerSample","source")])
alpha.cult$timepoint <- with(alpha.cult, ifelse(source=="culture","baseline",
                                                ifelse(source=="cultureFU","follow-up","oops")))
df.mapX2 <- merge(df.mapXtend, alpha.cult[,c("study_id","timepoint","isolatesPerSample")], by=c("study_id","timepoint"))
df.mapX2$isolatesPerSample <- as.character(df.mapX2$isolatesPerSample)

pDot.alpha <- ggplot(data=df.mapX2, aes(x=isolatesPerSample, y=ASVsumBact, color=timepoint)) +
  geom_jitter()+
  geom_boxplot() +
  #facet_wrap(sex~indole) +
  #scale_y_log10() +
  #scale_x_log10() +
  #ylim(0,100) +
  theme_bw() +
  #ggtitle("Library Size") +
  theme(text=element_text(size=12),axis.text.x = element_text(face="italic", angle=45, hjust=1),panel.spacing = unit(0, "lines")) +
  scale_color_manual(values = pal.CB)
pDot.alpha

df.mapT <- subset(df.mapX2, timepoint=="baseline")
cor(df.mapT$ASVsumBact, as.numeric(df.mapT$isolatesPerSample))






##############################################
# old scratch area
##############################################

### Trim to followed patients
phy.baseOverlap <- subset_samples(phyT.cleanNew, study_id%in%ls.baseOverlapIDs)

map_baseOverlap <- subset(table_map, study_id%in%ls.baseOverlapIDs)


map_baseOverlapOutcome <- unique(map_baseOverlap[,c("study_id","outcomeText")])


agg.baseOutcome <- setNames(aggregate(x = map_baseOverlapOutcome$study_id,       # Specify data column
                                      by = list(map_baseOverlapOutcome$outcomeText, ),   # Specify group indicator
                                      FUN = length),                          # Specify function (i.e. sum)
                            c("outcome","count"))        


df.outcomeInfect <- subset(df.outcome, outcomeText%in%c("deepInfection","flapFailure","nonUnion"))
df.outcomeHeal <- subset(df.outcome, outcomeText=="healed")


df.outcomeTrim <- rbind(df.outcomeInfect, df.outcomeHeal)

unique(map_baseOverlap$study_id)
