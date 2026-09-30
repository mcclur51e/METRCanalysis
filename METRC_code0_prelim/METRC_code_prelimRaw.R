#####################################################################
########## Processing Set −up #######################################
#####################################################################
# R version 4.5.0 (2025-04-11) -- "How About a Twenty-Six"

########## Call libraries for use ##########
# Packages from CRAN
library("tidyr") # version 1.3.0
library("tidyverse") # version 2.0.0
library("stringr") # version 1.5.0

##################################
########## End preamble ##########
##################################                     

########## Define master input and output locations ##########    
path <- "~/Desktop/Wannigan_METRC/" # change this to directory where you will be working
setwd(path) # set working directory
dir.create("outputPrelim") # create directory for output files
#dir.create("plots") # create directory to collect output plots

########## Import raw data ##########
df.cultBase <- data.frame(read.csv("raw/bioburdenCultureBaseList.csv", header = TRUE, check.names=FALSE)) # reads csv file into data.frame with row names in column 1
df.cultBase$bugName <- str_to_title(df.cultBase$MasterBaseBug) # Duplicate row with capitalized first letters and change column name
df.cultFU <- data.frame(read.csv("raw/bioburdenCultureFollowList.csv", header = TRUE, check.names=FALSE)) # reads csv file into data.frame with row names in column 1
df.cultFU$bugName <- str_to_title(df.cultFU$MasterFollowBug) # Duplicate row with capitalized first letters and change column name

### List the study IDs for datasets and find paired/unique/etc
#ls.studyIDs <- unique(c(df.pcr$study_id,df.cultBase$study_id,df.pcrFU$study_id,df.cultFU$study_id)) # list all study IDs
ls.cultIDs <- unique(df.cultBase$study_id) # list all study IDs (patient IDs) from culture data
ls.cultFUIDs <- unique(df.cultFU$study_id) # list all study IDs (patient IDs) from culture follow data
#ls.cultFUuniqueIDs <- setdiff(ls.cultFUIDs,unique(c(ls.pcrFUIDs,ls.cultIDs,ls.pcrIDs))) # samples that only have follow-up culture data

#ls.healIDs <- setdiff(ls.pcrIDs,union(ls.pcrFUIDs,ls.cultFUIDs)) # list all study IDs (patient IDs) with no follow-up pcr
#ls.baseFUIDs <- intersect(ls.pcrIDs,union(ls.pcrFUIDs,ls.cultFUIDs)) # list all study IDs (patient IDs) with follow-up pcr
#ls.noCultbaseFUIDs <- setdiff(union(ls.pcrFUIDs,ls.cultFUIDs), ls.pcrIDs) # list IDs with no baseline pcr data
#ls.noPCRbaseFUIDs <- setdiff(union(ls.pcrFUIDs,ls.cultFUIDs), ls.cultIDs) # list IDs with no baseline culture data
#ls.nobaseIDs <- setdiff(union(ls.pcrFUIDs,ls.cultFUIDs),union(ls.pcrIDs,ls.cultIDs)) # list IDs present only in follow-up data
#ls.follow4 <- Reduce(intersect, list(ls.pcrIDs,ls.pcrFUIDs,ls.cultIDs,ls.cultFUIDs)) # list IDs present in all 4 datasets

###################################################
########## Compile data into master sheet #########
###################################################


########## culture data ##########
### Add new columns 'genus' and 'species' from split 'bugName' column in culture data
### this is not perfect and only takes 1st ID, does not account for potential multiple matches listed
df.cultBase$genus <- df.cultBase$bugName
df.cultBaseSplit <- df.cultBase %>% separate(genus, c("genus","identifier"))
df.cultBaseSplit$abundance <- 1 # create abundance column indicating the taxon was cultured

### Remove yeast/fungus
df.cultBaseSplit <- subset(df.cultBaseSplit, !genus%in%c("Alternaria","Aspergillus", "Bipolaris", "Botryotrichum", 
                                                 "Candida", "Chaetomium", "Curvularia", "Dactylaria","Fusarium", 
                                                 "Mold", "Mucor", "Penicillium", "Scopulariopsis", "Trichosporon",
                                                 "Yeast", "Diptheroids", "Diphtheroids","Nonfermenter","Fermenters"))

### Shift mrsa/mssa to 'identifier' column
df.cultBaseSplit$identifier <- with(df.cultBaseSplit, ifelse(genus=="Mssa","MSSA",
                                                         ifelse(genus=="Mrsa","MRSA",
                                                                ifelse(identifier=="Coagulase","CoagNeg",identifier))))
### Correct a bunch of mis-typing and label MRSA/MSSA as Staphylococcus
df.cultBaseSplit$genus <- with(df.cultBaseSplit, ifelse(genus=="Bacteriodes","Bacteroides",
                                                        ifelse(genus=="Clostridial","Clostridium",
                                                               ifelse(genus=="Flavimonas","Pseudomonas",
                                                                      ifelse(genus=="Gamella","Gemella",
                                                                             ifelse(genus=="Nocardia","Nocardioides",
                                                                                    ifelse(genus=="Orchrobactrum","Ochrobactrum", 
                                                                                           ifelse(genus=="Peptostreptococci","Peptostreptococcus",
                                                                                                  ifelse(genus=="Propionibacterium", "Cutibacterium",
                                                                                                         ifelse(genus=="Stenotophomonas","Stenotrophomonas",
                                                                                                                ifelse(genus=="Mssa","Staphylococcus",
                                                                                                                       ifelse(genus=="Mrsa","Staphylococcus", genus))))))))))))
### Add identifier for empty cells
df.cultBaseSplit$identifier <- with(df.cultBaseSplit, ifelse(is.na(identifier),"spp.",
                                                             ifelse(identifier=="Species","spp.",
                                                                    ifelse(identifier=="Aureus","aureus",
                                                                           identifier))))

### Calculate taxa sum across dataset (from how many samples a taxon was cultured)
ls.cultBaseSampSum <- setNames(aggregate(x = df.cultBaseSplit$abundance,           # Specify data column
                                     by = list(df.cultBaseSplit$study_id),     # Specify group indicator
                                     FUN = sum),                               # Specify function (i.e. sum)
                           c("study_id","isolatesPerSample"))                          # specify column headers

### Add prevalence column to culture data
df.cultBaseSplit <- merge(df.cultBaseSplit[,c("study_id","genus","identifier","abundance")],
                      ls.cultBaseSampSum,by="study_id") # 
df.cultBaseSplit$source <- "culture" # add column to indicate that taxon came from culture
df.cultBaseSplit$fraction <- df.cultBaseSplit$abundance / df.cultBaseSplit$isolatesPerSample  # add column taxon abundance in the specified sample
df.cultBaseSplit$timepoint <- "baseline" # add column to indicate all samples collected at baseline
write.csv(df.cultBaseSplit,"outputPrelim/table_cultureBase.csv")

########## Follow-up culture data ##########
### Add new columns 'genus' and 'species' from split 'bugName' column in culture data
### this is not perfect and only takes 1st ID, does not account for potential multiple matches listed
df.cultFU$genus <- df.cultFU$bugName
df.cultFUSplit <- df.cultFU %>% separate(genus, c("genus","identifier"))
df.cultFUSplit$abundance <- 1 # create abundance column indicating the taxon was cultured

### Remove yeast/fungus
df.cultFUSplit <- subset(df.cultFUSplit, !genus%in%c("Alternaria","Aspergillus", "Bipolaris", "Botryotrichum", 
                                                     "Candida", "Chaetomium", "Curvularia", "Dactylaria","Fusarium", 
                                                     "Mold", "Mucor", "Penicillium", "Scopulariopsis", "Trichosporon",
                                                     "Yeast", "Diptheroids", "Diphtheroids","Nonfermenter","Fermenters"))
### Shift mrsa/mssa to 'identifier' column
df.cultFUSplit$identifier <- with(df.cultFUSplit, ifelse(genus=="Mssa","MSSA",
                                                         ifelse(genus=="Mrsa","MRSA",
                                                                ifelse(identifier=="Coagulase","CoagNeg",identifier))))
### Correct a bunch of mis-typing and label MRSA/MSSA as Staphylococcus
df.cultFUSplit$genus <- with(df.cultFUSplit, ifelse(genus=="Bacteriodes","Bacteroides",
                                                    ifelse(genus=="Clostridial","Clostridium",
                                                           ifelse(genus=="Flavimonas","Pseudomonas",
                                                                  ifelse(genus=="Gamella","Gemella",
                                                                         ifelse(genus=="Nocardia","Nocardioides",
                                                                                ifelse(genus=="Orchrobactrum","Ochrobactrum", 
                                                                                       ifelse(genus=="Peptostreptococci","Peptostreptococcus",
                                                                                              ifelse(genus=="Propionibacterium", "Cutibacterium",
                                                                                                     ifelse(genus=="Stenotophomonas","Stenotrophomonas",
                                                                                                            ifelse(genus=="Mssa","Staphylococcus",
                                                                                                                   ifelse(genus=="Mrsa","Staphylococcus", genus))))))))))))
### Add identifier for empty cells
df.cultFUSplit$identifier <- with(df.cultFUSplit, ifelse(is.na(identifier),"spp.",
                                                         ifelse(identifier=="Species","spp.",
                                                                ifelse(identifier=="Aureus","aureus",
                                                                       identifier))))
### Calculate taxa sum across dataset (from how many samples a taxon was cultured)
ls.cultFUSampSum <- setNames(aggregate(x = df.cultFUSplit$abundance,           # Specify data column
                                       by = list(df.cultFUSplit$study_id),     # Specify group indicator
                                       FUN = sum),                               # Specify function (i.e. sum)
                             c("study_id","isolatesPerSample"))                 # specify column headers
### Add prevalence column to culture data
df.cultFUSplit <- merge(df.cultFUSplit[,c("study_id","genus","identifier","abundance")],
                        ls.cultFUSampSum,by="study_id") # 
df.cultFUSplit$source <- "culture" # add column to indicate that taxon came from culture in follow-up sample
df.cultFUSplit$timepoint <- "follow-up" # add column to indicate all samples collected at follow-up
df.cultFUSplit$fraction <- df.cultFUSplit$abundance / df.cultFUSplit$isolatesPerSample  # add column calculating abundance of taxon in sample
write.csv(df.cultFUSplit,"outputPrelim/table_cultureFollowUp.csv")

########## Compile raw datasets together ##########
df.cultCompile <- rbind(df.cultBaseSplit, df.cultFUSplit)
write.csv(df.cultCompile,"outputPrelim/table_compileCulture.csv")


####################################################################
########## You have now successfully imported all the ##############
########## raw culture data into R and performed ###################
########## preliminary clean-up so that you are ready ##############
########## to proceed with further analysis ########################
####################################################################

                                          