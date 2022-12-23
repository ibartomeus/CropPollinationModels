# This script reads data from CropPol and cleans it for analysis.
# clean data is stored in /data

#first load needed packages----

library(tidyverse)
library(RCurl)

#now we read data from CropPol----

x <- getURL("https://raw.githubusercontent.com/ibartomeus/OBservData/master/Final_Data/CropPol_field_level_data.csv")
field_data <- read.csv(text = x)

# data processing: study IDs----

# in order to maximize the sample size within each study, 
# we keep all the years of each study together.  
# CropPol currently gives each study_year a separate study_id, 
# so I make a "study_id2" column that lumps years within study, 
# and then fix various inconsistencies
##IB note: I wonder if keeping years separated makes also sense

field_data$study_id2 <- NA
for (i in 1:nrow(field_data)) {
  thisrow = field_data[i,]
  # get rid of the year using strsplit() to separate the differnt parts 
  # of studyID and then subsetting away the year part
  field_data$study_id2[i] <- paste(subset(
    unlist(strsplit(thisrow$study_id, split='_')), 
    unlist(strsplit(thisrow$study_id, split='_')) 
    != thisrow$sampling_year), collapse="_")
}

# fix naming issues that didn't process correctly above
##IB note: it would be great if we can do that in a generic way. Why are badly processes?
field_data[which(field_data$study_id2=="Alexandra_Maria_Klein_Coffea_arabica_Indonesia_2000_2001"),]$study_id2 <- "Alexandra_Maria_Klein_Coffea_arabica_Indonesia"
field_data[which(field_data$study_id2=="Alexandra_Maria_Klein_Coffea_canephora_Indonesia_2000_2001"),]$study_id2 <- "Alexandra_Maria_Klein_Coffea_canephora_Indonesia"
field_data[which(field_data$study_id2=="Alice_Classen_Coffea_arabica_Tanzania_2011_2012"),]$study_id2 <- "Alice_Classen_Coffea_arabica_Tanzania"
field_data[which(field_data$study_id2=="Bryony_Willcox_Mangifera_indica_Australia_2"),]$study_id2 <- "Bryony_Willcox_Mangifera_indica_Australia"
field_data[which(field_data$study_id2=="Davi_L_Ramos_Phaseolus_vulgaris L_Brazil_2015_2016"),]$study_id2 <- "Davi_L_Ramos_Phaseolus_vulgaris L_Brazil"
field_data[which(field_data$study_id2=="Margaret_Mayfield_Actinidia_deliciosa_New_Zealand_NA"),]$study_id2 <- "Margaret_Mayfield_Actinidia_deliciosa_New_Zealand"

# there are some "studies" where the different years come from different papers.
# keep these separate if the methods differed by year, otherwise lump them

# Dara_Stanley_Brassica_napus_Ireland
# 2009 and 2010 data come from different papers using differnt methods, 
# so I will kept separate
field_data[which(field_data$study_id=="Dara_Stanley_Brassica_napus_Ireland_2009"),]$study_id2 <- "Dara_Stanley_Brassica_napus_Ireland_study_A"
field_data[which(field_data$study_id=="Dara_Stanley_Brassica_napus_Ireland_2010"),]$study_id2 <- "Dara_Stanley_Brassica_napus_Ireland_study_B"

# David_Kleijn_Malus_domestica_Netherlands
# methods for 2010-11 appear different from 2013-14, so keep these separate
field_data[which(field_data$study_id=="David_Kleijn_Malus_domestica_Netherlands_2010"),]$study_id2 <- "David_Kleijn_Malus_domestica_Netherlands_study_A"
field_data[which(field_data$study_id=="David_Kleijn_Malus_domestica_Netherlands_2011"),]$study_id2 <- "David_Kleijn_Malus_domestica_Netherlands_study_A"
field_data[which(field_data$study_id=="David_Kleijn_Malus_domestica_Netherlands_2013"),]$study_id2 <- "David_Kleijn_Malus_domestica_Netherlands_study_B"
field_data[which(field_data$study_id=="David_Kleijn_Malus_domestica_Netherlands_2014"),]$study_id2 <- "David_Kleijn_Malus_domestica_Netherlands_study_B"

# Alexandra_Maria_Klein_Prunus_dulcis_USA
# 2008 and 2009 come from different papers, but methods appear identical, so I will keep them lumped

# Rachael_Winfree_Citrullus_lanatus_USA
# 2004-2007 are net data, 2008-12 are the visit data.  
# Note: these won't be used in the main analysis anyway since there is no associated yield data.
field_data[which(field_data$study_id=="Rachael_Winfree_Citrullus_lanatus_USA_2004"),]$study_id2 <- "Rachael_Winfree_Citrullus_lanatus_USA_study_A"
field_data[which(field_data$study_id=="Rachael_Winfree_Citrullus_lanatus_USA_2005"),]$study_id2 <- "Rachael_Winfree_Citrullus_lanatus_USA_study_A"
field_data[which(field_data$study_id=="Rachael_Winfree_Citrullus_lanatus_USA_2007"),]$study_id2 <- "Rachael_Winfree_Citrullus_lanatus_USA_study_A"
field_data[which(field_data$study_id=="Rachael_Winfree_Citrullus_lanatus_USA_2008"),]$study_id2 <- "Rachael_Winfree_Citrullus_lanatus_USA_study_B"
field_data[which(field_data$study_id=="Rachael_Winfree_Citrullus_lanatus_USA_2010"),]$study_id2 <- "Rachael_Winfree_Citrullus_lanatus_USA_study_B"
field_data[which(field_data$study_id=="Rachael_Winfree_Citrullus_lanatus_USA_2011"),]$study_id2 <- "Rachael_Winfree_Citrullus_lanatus_USA_study_B"
field_data[which(field_data$study_id=="Rachael_Winfree_Citrullus_lanatus_USA_2012"),]$study_id2 <- "Rachael_Winfree_Citrullus_lanatus_USA_study_B"

# fix naming issues
field_data[which(field_data$study_id2=="David_Biddinger_Malus_pumila_USA"),]$study_id2 <- "David_Biddinger_Malus_domestica_USA"
field_data[which(field_data$study_id2=="Julianna_Wilson_Malus_pumila_USA"),]$study_id2 <- "Julianna_Wilson_Malus_domestica_USA"
field_data[which(field_data$study_id2=="Elizabeth_Elle_Vaccinium_corymbosum_USA"),]$study_id2 <- "Elizabeth_Elle_Vaccinium_corymbosum_Canada"
field_data[which(field_data$crop=="Malus pumila"),]$crop <- "Malus domestica"

unique(field_data$study_id2)	# 144 studies

# data processing: drop data that we can't use----

# drop pan trap data
field_data2a = subset(field_data, !study_id2 %in% c(
  "Breno_M_Freitas_Gossypium_hirsutum_Brazil",
  "Georg_Andersson_Fragaria_ananassa_Sweden",
  "Heather_Lee_Grab_Fragaria_ananassa_USA",
  "Rachel_Mallinger_Malus_domestica_USA",
  "Yi_Zou_Brassica_napus_China"
)) ##IB note it would be great if we can do that in a generic way.

# drop enclosure data
field_data2a = subset(field_data2a, !study_id2 %in% c(	
  "Willis_Chan_Raine_Cucurbita_pepo_Canada"
))	##IB note: it would be great if we can do that in a generic way.

# drop for missing HB data (HB data not recorded)	
field_data2a = subset(field_data2a, !study_id2 %in% c(	
  "Luisa_G_Carvalheiro_Helianthus_annuus_South_Africa",
  "Rachael_Winfree_Capsicum_annuum_USA",
  "Rachael_Winfree_Cucumis_melo_USA",
  "Rachael_Winfree_Malus_Domestica_USA",
  "Rachael_Winfree_Solanum_lycopersicum_USA",
  "Rachael_Winfree_Vaccinium_corymbosum_USA",
  "Rachael_Winfree_Citrullus_lanatus_USA_study_A",
  "Rachael_Winfree_Vaccinium_macrocarpon_USA",
  "Yael_Mandelik_Citrullus_lanatus_Israel",
  "Yael_Mandelik_Helianthus_annuus_Israel"
)) ##IB note: it would be great if we can do that in a generic way.

# drop studies with unclear methods for now, unless we learn that they were ok--see below
field_data2a = subset(field_data2a, !study_id2 %in% c(	
  "Breno_M_Freitas_Bixa_orellana_Brazil",
  "Johan_Ekroos_Vicia_faba_Sweden",
  "Rebecca_Steward_Fragaria_ananassa_Sweden"
))

# notes on why these 3 studies were dropped:
# Breno_M_Freitas_Bixa_orellana_Brazil
# abundance method listed as "other" in the database, so it might be pan trap data, etc, and it also might be bees only. The paper listed as Dainese et al 2019, but I can't find it in the data provided
# Johan_Ekroos_Vicia_faba_Sweden
# taxa listed as "only bees (bumblebees)"  does this mean that HB were not sampled, or were there just not any?  Were other wild bees than bumble bees sampled?  There is no paper associated with this data
# Rebecca_Steward_Fragaria_ananassa_Sweden
# there is only data for bees in the wild bees column.  The paper says HB were sampled, but most bees were bumblebees.  But there are 0 HB in the data--where are they?  

unique(field_data2a$study_id2)	# 125 studies

# data processing: visit rates and abundance----

# In the CropPol database, insect visitors are recorded either as "abundance" 
# or "visitation_rate".  Abundance is typically = net data, 
# and visitation_rate = timed observations.  
# For this analysis, we will use visitation rate if it is available 
# and the data looks reasonable, and otherwise use the abundance data.

# I wrote a loop to do this (not included in this file), but then I had to change a few studies manually after checking the data, so it is currently being done with a manually created csv.
# In particular:
# I use abundance instead of visits for the Garratt studies because it looks like the visit data might not include all pollinators, whereas the abundance data does
# I use visits instead of abundance for Taylor_Ricketts_Coffea_arabica_Costa_Rica because the visit data has only 4 fewer rows and that is probably worth it
# I use visits instead of abundance for  Thijs_Fijen_Allium_porrum_France, Thijs_Fijen_Allium_porrum_Italy, and Virginie_Boreux_Coffea_canephora_India because of potential data issues in the visit data

# read in a csv that lists whether to use the visitation_rate or abundance column for each study
field_data_summary_edited <- read.csv(file="../Reilly_models/field_data_summary_edited_v3.csv", header=T)
field_data_summary_edited 
##IB note: it would be great if we add this info to CropPol to do that in a generic way.


# After going back to the relevant papers, I also made some manual edits to the richness_restriction column in the CropPol database, which lists which insect groups were recorded for each study. These edits are in the "taxa_recorded_edited" column in this csv.  
# In particular:
# Saul_A_Cunningham_Brassica_napus_Australia
# from the paper: "Pollinators were counted in three categories: feral honeybees, hover flies and native bees."  I changed the category from NA to "bees and syrphids"
# Blande_Viana_Passiflora_edulis_Brazil
# this data comes from a thesis and was also in Ricketts et al 2008.  It is clearly bees only, so I changed the category from NA to "bees"
# Alejandro_Trillo_Fragaria_ananassa_Spain:
# the methods in the paper indicate that they used different methods to sample bees and other insects. This is a problem for comparisons, so I changed the category from "all visitors considered" to "bees", so that only the bees (including honey bee) part will be used. 

# merge this csv with the field data
field_data2b <- merge(field_data2a, field_data_summary_edited, by="study_id2")

# check to make sure we still have the correct number of rows
nrow(field_data2a)	# 2993
nrow(field_data2b)	# 2993

#setdiff(unique(field_data2a$study_id2), unique(field_data2b$study_id2))

# data processing: clean up 0's and NAs----

# there is some inconsistency in the database about whether missing insect data is recorded as 0s or NAs
# This can be a serious issue, especially when using na.rm=T, as we do in some of the functions below.
# To be safe, the following code will go through and manually reset NAs and 0's depending on which taxa 
# were actually recorded based on the "taxa_recorded_edited" column

# these are the CropPol insect categories:  
# honeybee
# bombus
# wildbees # note that this DOES NOT include bombus!!
# syrphids
# humbleflies
# other_flies
# beetles
# lepidoptera
# nonbee_hymenoptera
# others

# ***JR code note: what I do in the loop below is subsetting a single row from the data,
# doing some processing on it, then rebuilding the dataframe from the processed rows
# using rbind.  There are various other ways this can be accomplished, but I 
# find this method more intuitive and easier to troubleshoot.

# set up a new empty dataframe
field_data2c = data.frame()

# loop through each row, replacing NAs and 0's depending on which taxa were listed as being recorded
for (i in 1:nrow(field_data2b)) {
  thisrow = field_data2b[i,]
  # if this is an abundance study, we need to do this for the abundance columns (ab_)
  if (thisrow$abund_or_visits_edited=="use_abundance") {
    if (!is.na(thisrow$abundance)) {
      if (thisrow$taxa_recorded_edited=="bees") {
        if (is.na(thisrow$ab_honeybee)) { thisrow$ab_honeybee <- 0 }
        if (is.na(thisrow$ab_bombus)) { thisrow$ab_bombus <- 0 }
        if (is.na(thisrow$ab_wildbees)) { thisrow$ab_wildbees <- 0 }
        thisrow$ab_syrphids <- NA
        thisrow$ab_humbleflies <- NA
        thisrow$ab_other_flies <- NA
        thisrow$ab_beetles <- NA
        thisrow$ab_lepidoptera <- NA
        thisrow$ab_nonbee_hymenoptera <- NA
        thisrow$ab_others <- NA
      } else if (thisrow$taxa_recorded_edited=="all visitors") {
        if (is.na(thisrow$ab_honeybee)) { thisrow$ab_honeybee <- 0 }
        if (is.na(thisrow$ab_bombus)) { thisrow$ab_bombus <- 0 }
        if (is.na(thisrow$ab_wildbees)) { thisrow$ab_wildbees <- 0 }			
        if (is.na(thisrow$ab_syrphids)) { thisrow$ab_syrphids <- 0 }	
        if (is.na(thisrow$ab_humbleflies)) { thisrow$ab_humbleflies <- 0 }	
        if (is.na(thisrow$ab_other_flies)) { thisrow$ab_other_flies <- 0 }	
        if (is.na(thisrow$ab_beetles)) { thisrow$ab_beetles <- 0 }	
        if (is.na(thisrow$ab_lepidoptera)) { thisrow$ab_lepidoptera <- 0 }	
        if (is.na(thisrow$ab_nonbee_hymenoptera)) { thisrow$ab_nonbee_hymenoptera <- 0 }
        if (is.na(thisrow$ab_others)) { thisrow$ab_others <- 0 }	
      } else if (thisrow$taxa_recorded_edited=="bees and leps") { 
        if (is.na(thisrow$ab_honeybee)) { thisrow$ab_honeybee <- 0 }
        if (is.na(thisrow$ab_bombus)) { thisrow$ab_bombus <- 0 }
        if (is.na(thisrow$ab_wildbees)) { thisrow$ab_wildbees <- 0 }			
        thisrow$ab_syrphids <- NA
        thisrow$ab_humbleflies <- NA
        thisrow$ab_other_flies <- NA
        thisrow$ab_beetles <- NA
        if (is.na(thisrow$ab_lepidoptera)) { thisrow$ab_lepidoptera <- 0 }	
        thisrow$ab_nonbee_hymenoptera <- NA
        thisrow$ab_others <- NA			
      } else if (thisrow$taxa_recorded_edited=="bees and syrphids") { 
        if (is.na(thisrow$ab_honeybee)) { thisrow$ab_honeybee <- 0 }
        if (is.na(thisrow$ab_bombus)) { thisrow$ab_bombus <- 0 }
        if (is.na(thisrow$ab_wildbees)) { thisrow$ab_wildbees <- 0 }			
        if (is.na(thisrow$ab_syrphids)) { thisrow$ab_syrphids <- 0 }	
        thisrow$ab_humbleflies <- NA
        thisrow$ab_other_flies <- NA
        thisrow$ab_beetles <- NA
        thisrow$ab_lepidoptera <- NA
        thisrow$ab_nonbee_hymenoptera <- NA
        thisrow$ab_others <- NA			
      } else if (thisrow$taxa_recorded_edited=="bees, syrphids and wasps") { 
        if (is.na(thisrow$ab_honeybee)) { thisrow$ab_honeybee <- 0 }
        if (is.na(thisrow$ab_bombus)) { thisrow$ab_bombus <- 0 }
        if (is.na(thisrow$ab_wildbees)) { thisrow$ab_wildbees <- 0 }			
        if (is.na(thisrow$ab_syrphids)) { thisrow$ab_syrphids <- 0 }	
        thisrow$ab_humbleflies <- NA
        thisrow$ab_other_flies <- NA
        thisrow$ab_beetles <- NA
        thisrow$ab_lepidoptera <- NA
        if (is.na(thisrow$ab_nonbee_hymenoptera)) { thisrow$ab_nonbee_hymenoptera <- 0 }
        thisrow$ab_others <- NA					
      }
      field_data2c = rbind(thisrow, field_data2c)
    }
    # if this is an visitation_rate study, we need to do this for the visit rate columns (visit_)	
  } else if (thisrow$abund_or_visits_edited=="use_visits") {
    if (!is.na(thisrow$visitation_rate)) {
      if (thisrow$taxa_recorded_edited=="bees") {
        if (is.na(thisrow$visit_honeybee)) { thisrow$visit_honeybee <- 0 }
        if (is.na(thisrow$visit_bombus)) { thisrow$visit_bombus <- 0 }
        if (is.na(thisrow$visit_wildbees)) { thisrow$visit_wildbees <- 0 }
        thisrow$visit_syrphids <- NA
        thisrow$visit_humbleflies <- NA
        thisrow$visit_other_flies <- NA
        thisrow$visit_beetles <- NA
        thisrow$visit_lepidoptera <- NA
        thisrow$visit_nonbee_hymenoptera <- NA
        thisrow$visit_others <- NA
      } else if (thisrow$taxa_recorded_edited=="all visitors") {
        if (is.na(thisrow$visit_honeybee)) { thisrow$visit_honeybee <- 0 }
        if (is.na(thisrow$visit_bombus)) { thisrow$visit_bombus <- 0 }
        if (is.na(thisrow$visit_wildbees)) { thisrow$visit_wildbees <- 0 }			
        if (is.na(thisrow$visit_syrphids)) { thisrow$visit_syrphids <- 0 }	
        if (is.na(thisrow$visit_humbleflies)) { thisrow$visit_humbleflies <- 0 }	
        if (is.na(thisrow$visit_other_flies)) { thisrow$visit_other_flies <- 0 }	
        if (is.na(thisrow$visit_beetles)) { thisrow$visit_beetles <- 0 }	
        if (is.na(thisrow$visit_lepidoptera)) { thisrow$visit_lepidoptera <- 0 }	
        if (is.na(thisrow$visit_nonbee_hymenoptera)) { thisrow$visit_nonbee_hymenoptera <- 0 }
        if (is.na(thisrow$visit_others)) { thisrow$visit_others <- 0 }	
      } else if (thisrow$taxa_recorded_edited=="bees and leps") { 
        if (is.na(thisrow$visit_honeybee)) { thisrow$visit_honeybee <- 0 }
        if (is.na(thisrow$visit_bombus)) { thisrow$visit_bombus <- 0 }
        if (is.na(thisrow$visit_wildbees)) { thisrow$visit_wildbees <- 0 }			
        thisrow$visit_syrphids <- NA
        thisrow$visit_humbleflies <- NA
        thisrow$visit_other_flies <- NA
        thisrow$visit_beetles <- NA
        if (is.na(thisrow$visit_lepidoptera)) { thisrow$visit_lepidoptera <- 0 }	
        thisrow$visit_nonbee_hymenoptera <- NA
        thisrow$visit_others <- NA			
      } else if (thisrow$taxa_recorded_edited=="bees and syrphids") { 
        if (is.na(thisrow$visit_honeybee)) { thisrow$visit_honeybee <- 0 }
        if (is.na(thisrow$visit_bombus)) { thisrow$visit_bombus <- 0 }
        if (is.na(thisrow$visit_wildbees)) { thisrow$visit_wildbees <- 0 }			
        if (is.na(thisrow$visit_syrphids)) { thisrow$visit_syrphids <- 0 }	
        thisrow$visit_humbleflies <- NA
        thisrow$visit_other_flies <- NA
        thisrow$visit_beetles <- NA
        thisrow$visit_lepidoptera <- NA
        thisrow$visit_nonbee_hymenoptera <- NA
        thisrow$visit_others <- NA			
      } else if (thisrow$taxa_recorded_edited=="bees, syrphids and wasps") { 
        if (is.na(thisrow$visit_honeybee)) { thisrow$visit_honeybee <- 0 }
        if (is.na(thisrow$visit_bombus)) { thisrow$visit_bombus <- 0 }
        if (is.na(thisrow$visit_wildbees)) { thisrow$visit_wildbees <- 0 }			
        if (is.na(thisrow$visit_syrphids)) { thisrow$visit_syrphids <- 0 }	
        thisrow$visit_humbleflies <- NA
        thisrow$visit_other_flies <- NA
        thisrow$visit_beetles <- NA
        thisrow$visit_lepidoptera <- NA
        if (is.na(thisrow$visit_nonbee_hymenoptera)) { thisrow$visit_nonbee_hymenoptera <- 0 }
        thisrow$visit_others <- NA					
      }
      field_data2c = rbind(thisrow, field_data2c)
    }
  }
}

# note that about 100 rows were removed due to NAs in either visitation rate 
# or abundance at this step, which is as expected
# most of the missing data are NAs scattered across many of the studies.  
# Of special note is "Thijs_Fijen_Allium_porrum_France" which has more NAs 
# than usual in the visit data, but it still looks better than the abundance 
# data which appears to have an issue with duplicated values.

# checking which studies it was:
unique(subset(field_data2b, site_id %in% setdiff(field_data2b$site_id, field_data2c$site_id))$study_id2)

nrow(field_data2b)
nrow(field_data2c)	# 2854


# data processing: Wild insect and all insect sums-----

# in the analysis, we are summing up the various recorded groups into "wild insects" 
# and maybe some other groups
# group sums that we might want:
# honeybee
# wild_insects
# wild_bees
# non_bees
# all_bees
# all_insects

# functions to do the sums for abundance data (note na.rm=T, which will ignore NAs
# when making the sums):
wild_insects_sums = function(x) {sum(x$ab_wildbees, x$ab_bombus, x$ab_syrphids, x$ab_humbleflies, x$ab_other_flies, x$ab_beetles, x$ab_lepidoptera, x$ab_nonbee_hymenoptera, x$ab_others, na.rm=T)}
wild_bees_sums = function(x) {sum(x$ab_wildbees, x$ab_bombus, na.rm=T)}
non_bee_sums = function(x) {sum(x$ab_syrphids, x$ab_humbleflies, x$ab_other_flies, x$ab_beetles, x$ab_lepidoptera, x$ab_nonbee_hymenoptera, x$ab_others, na.rm=T)}
all_bees_sums = function(x) {sum(x$ab_honeybee, x$ab_wildbees, x$ab_bombus, na.rm=T)}
all_insects_sums = function(x) {sum(x$ab_honeybee, x$ab_wildbees, x$ab_bombus, x$ab_syrphids, x$ab_humbleflies, x$ab_other_flies, x$ab_beetles, x$ab_lepidoptera, x$ab_nonbee_hymenoptera, x$ab_others, na.rm=T)}

# functions to do the sums for visit rate data:
wild_insects_sums_v = function(x) {sum(x$visit_wildbees, x$visit_bombus, x$visit_syrphids, x$visit_humbleflies, x$visit_other_flies, x$visit_beetles, x$visit_lepidoptera, x$visit_nonbee_hymenoptera, x$visit_others, na.rm=T)}
wild_bees_sums_v = function(x) {sum(x$visit_wildbees, x$visit_bombus, na.rm=T)}
non_bee_sums_v = function(x) {sum(x$visit_syrphids, x$visit_humbleflies, x$visit_other_flies, x$visit_beetles, x$visit_lepidoptera, x$visit_nonbee_hymenoptera, x$visit_others, na.rm=T)}
all_bees_sums_v = function(x) {sum(x$visit_honeybee, x$visit_wildbees, x$visit_bombus, na.rm=T)}
all_insects_sums_v = function(x) {sum(x$visit_honeybee, x$visit_wildbees, x$visit_bombus, x$visit_syrphids, x$visit_humbleflies, x$visit_other_flies, x$visit_beetles, x$visit_lepidoptera, x$visit_nonbee_hymenoptera, x$visit_others, na.rm=T)}

# apply functions to sum abundance and visit rates by row
for (i in 1:nrow(field_data2c)) {
  
  field_data2c$ab_wild_insects[i] = wild_insects_sums(field_data2c[i,])
  field_data2c$ab_wild_bees[i] = wild_bees_sums(field_data2c[i,])
  field_data2c$ab_non_bees[i] = non_bee_sums(field_data2c[i,])	
  field_data2c$ab_all_bees[i] = all_bees_sums(field_data2c[i,])
  field_data2c$ab_all_insects[i] = all_insects_sums(field_data2c[i,])
  
  field_data2c$visit_wild_insects[i] = wild_insects_sums_v(field_data2c[i,])
  field_data2c$visit_wild_bees[i] = wild_bees_sums_v(field_data2c[i,])
  field_data2c$visit_non_bees[i] = non_bee_sums_v(field_data2c[i,])
  field_data2c$visit_all_bees[i] = all_bees_sums_v(field_data2c[i,])
  field_data2c$visit_all_insects[i] = all_insects_sums_v(field_data2c[i,])
}

# make a new column for each insect group that pulls either visits or abundance 
# depending on what was chosen above

# define the list of studies
studies = unique(field_data2c$study_id2)

# ***JR code note: again, looping through subsets of the data, processing them, 
# and rebuilding the dataset using rbind

# create new empty dataframe
field_data2d = data.frame()

for (i in 1:length(studies)) {
  
  fd_sub = subset(field_data2c, study_id2==studies[i])
  
  # visit rate studies
  if (subset(field_data_summary_edited, study_id2==studies[i])$abund_or_visits_edited == "use_visits") {
    
    fd_sub$honeybee = fd_sub$visit_honeybee
    fd_sub$wild_bees = fd_sub$visit_wild_bees
    fd_sub$all_bees = fd_sub$visit_all_bees
    fd_sub$non_bees = fd_sub$visit_non_bees
    fd_sub$wild_insects = fd_sub$visit_wild_insects
    fd_sub$all_insects = fd_sub$visit_all_insects
    fd_sub$count_type = "visits"
    
    # abundance studies	
  } else {
    fd_sub$honeybee = fd_sub$ab_honeybee
    fd_sub$wild_bees = fd_sub$ab_wild_bees
    fd_sub$all_bees = fd_sub$ab_all_bees
    fd_sub$non_bees = fd_sub$ab_non_bees
    fd_sub$wild_insects = fd_sub$ab_wild_insects
    fd_sub$all_insects = fd_sub$ab_all_insects
    fd_sub$count_type = "abundance"
  }
  
  field_data2d = rbind(fd_sub, field_data2d)
  
}

nrow(field_data2c)
nrow(field_data2d)

# data processing: yield data-----

# fix yield data for ICP datasets to match Reilly et al 2020 - it seems at some point the yield vs yield2 were reversed for these:
# we can drop this part if we make the edits to CropPol

##IB note: it would be great to fix that in CropPol. Maybe is already fixed. Alfonso?


# sweet cherry should be fruit weight per branch
field_data2d[which(field_data2d$study_id2=="Theresa_PittsSinger_Prunus_dulcis_USA"),]$yield <- field_data2d[which(field_data2d$study_id2=="Theresa_PittsSinger_Prunus_dulcis_USA"),]$yield2
field_data2d[which(field_data2d$study_id2=="Shelby_Fleischer_Cucurbita_pepo_USA"),]$yield <- field_data2d[which(field_data2d$study_id2=="Shelby_Fleischer_Cucurbita_pepo_USA"),]$yield2
field_data2d[which(field_data2d$study_id2=="Neal_Williams_Citrullus_lanatus_USA"),]$yield <- field_data2d[which(field_data2d$study_id2=="Neal_Williams_Citrullus_lanatus_USA"),]$yield2
field_data2d[which(field_data2d$study_id2=="Jamie_Ellis_Citrullus_lanatus_USA"),]$yield <- field_data2d[which(field_data2d$study_id2=="Jamie_Ellis_Citrullus_lanatus_USA"),]$yield2
field_data2d[which(field_data2d$study_id2=="Nikki_Rothwell_Prunus_cerasus_USA"),]$yield <- field_data2d[which(field_data2d$study_id2=="Nikki_Rothwell_Prunus_cerasus_USA"),]$yield2
field_data2d[which(field_data2d$study_id2=="David_Biddinger_Prunus_cerasus_USA"),]$yield <- field_data2d[which(field_data2d$study_id2=="David_Biddinger_Prunus_cerasus_USA"),]$yield2

field_data2d[which(field_data2d$study_id2=="Theresa_PittsSinger_Prunus_dulcis_USA"),]$yield_units <- field_data2d[which(field_data2d$study_id2=="Theresa_PittsSinger_Prunus_dulcis_USA"),]$yield2_units
field_data2d[which(field_data2d$study_id2=="Shelby_Fleischer_Cucurbita_pepo_USA"),]$yield_units <- field_data2d[which(field_data2d$study_id2=="Shelby_Fleischer_Cucurbita_pepo_USA"),]$yield2_units
field_data2d[which(field_data2d$study_id2=="Neal_Williams_Citrullus_lanatus_USA"),]$yield_units <- field_data2d[which(field_data2d$study_id2=="Neal_Williams_Citrullus_lanatus_USA"),]$yield2_units
field_data2d[which(field_data2d$study_id2=="Jamie_Ellis_Citrullus_lanatus_USA"),]$yield_units <- field_data2d[which(field_data2d$study_id2=="Jamie_Ellis_Citrullus_lanatus_USA"),]$yield2_units
field_data2d[which(field_data2d$study_id2=="Nikki_Rothwell_Prunus_cerasus_USA"),]$yield_units <- field_data2d[which(field_data2d$study_id2=="Nikki_Rothwell_Prunus_cerasus_USA"),]$yield2_units
field_data2d[which(field_data2d$study_id2=="David_Biddinger_Prunus_cerasus_USA"),]$yield_units <- field_data2d[which(field_data2d$study_id2=="David_Biddinger_Prunus_cerasus_USA"),]$yield2_units

# berry wt per branch (the preferred yield var) column apparently didn't get uploaded, so need to fix this manually for sweet cherry in WA  
field_data2d[which(field_data2d$study_id2=="Robert_Gillespie_Prunus_avium_USA"),]$yield <- 
  c(
    587.506,
    710.985,
    571.16222,
    NA,
    140.12,
    173.98889,
    133.535,
    288.7975,
    195.05714,
    NA,
    NA,
    NA,
    NA,
    237.9425,
    196.06125,
    109.57,
    95.42,
    142.686,
    98.75286,
    94.55143,
    122.701,
    132.87143,
    223.94667,
    NA,
    NA,
    NA,
    216.61333,
    127.70333,
    161.42143,
    54.16333,
    114.76667,
    274.435,
    162.767,
    148.112,
    240.29
  )
field_data2d[which(field_data2d$study_id2=="Robert_Gillespie_Prunus_avium_USA"),]$yield_units <- "fruit weight per branch"

# remove rows with no yield data
field_data3 = subset(field_data2d, !is.na(field_data2d$yield))

# check number of rows
nrow(field_data2d)
nrow(field_data3)

# check number of studies
length(unique(field_data2d$study_id2))
length(unique(field_data3$study_id2))	#97


# data processing: transform to z-scores----

# note: for data that is all zeros, (e.g. when no non-bees were recorded for a crop), 
# it is not possible to calculate z-scores because there is no variance (divide by zero).
# This happens in several crops:
#ICP_Prunus_dulcis_USA-CA (wild bees)
#ICP_Vaccinium_corymbosum_USA-OR (non-bees)
#Mark_Otieno_Cajanus_cajan_Kenya (non-bees)
#Virginie_Boreux_Coffea_canephora_India (non-bees)

# note on NAs for richness z-scores: scale() removes NAs when calculating the means.  
# I checked that it was not introducing NAs for every row in a crop when there are only 
# a couple NAs--it was not.  This was only an issue for Natacha_Chacoff_Citrus_paradisi_Argentina, 
# which had NAs for richness in 4/12 fields.  Most or all of the other NAs in
# richness were already removed for missing other data. 

# define new list of studies that remain
studies = unique(field_data3$study_id2)

# create empty dataframe	
field_data4 = data.frame()

# loop, process, and rebuild
for (i in 1:length(studies)) {
  
  fd_sub = subset(field_data3, study_id2==studies[i])
  
  fd_sub$yield_z = scale(fd_sub$yield)[,1]
  fd_sub$honeybee_z = scale(fd_sub$honeybee)[,1]
  fd_sub$wild_bees_z = scale(fd_sub$wild_bees)[,1]
  fd_sub$all_bees_z = scale(fd_sub$all_bees)[,1]
  fd_sub$non_bees_z = scale(fd_sub$non_bees)[,1]
  fd_sub$wild_insects_z = scale(fd_sub$wild_insects)[,1]
  fd_sub$all_insects_z = scale(fd_sub$all_insects)[,1]
  fd_sub$richness_z = scale(fd_sub$observed_pollinator_richness)[,1]
  
  # also calculate SDs for comparing sd of HB to WI
  fd_sub$honeybee_sd = sd(fd_sub$honeybee)
  fd_sub$wild_insects_sd = sd(fd_sub$wild_insects)
  fd_sub$wild_bees_sd = sd(fd_sub$wild_bees)
  fd_sub$non_bees_sd = sd(fd_sub$non_bees)
  
  field_data4 = rbind(fd_sub, field_data4)
  
}

# remove rows with NA in the yield z-score
field_data5 = subset(field_data4, !is.na(field_data4$yield_z))

# check number of rows

nrow(field_data3)
nrow(field_data4)
nrow(field_data5) 

# check number of studies
length(unique(field_data4$study_id2))
length(unique(field_data5$study_id2)) # Annona squamosa and muricata dropped because 
# only 1 field each and can't calculate z scores (creates NAs) for yield (or insects 
# either, actually)
# 94

####################
# replace NAs in HB due to divide by zero in z-score with manually entered 0's.  
# I am currently thinking that this is better than just dropping them from model runs.
# this occurred in 2 studies:
# [1] "Finbarr_G_Horgan_Abelmoschus_esculentus_Philippines"		# appears HB recorded, but were 0
# [2] "Breno_M_Freitas_Malpighia_emarginata_Brazil" 			# appears HB recorded, but were 0

field_data5a = data.frame()
for (i in 1:nrow(field_data5)) {
  thisrow = field_data5[i,]
  if (is.na(thisrow$honeybee_z)) { thisrow$honeybee_z <- 0 }
  field_data5a = rbind(thisrow, field_data5a)
}
nrow(field_data5a)


# data processing: sample size-----

# drop studies with less than 3 fields, or more as suggested by Garibaldi

studies = names(subset(table(field_data5a$study_id2), table(field_data5a$study_id2)>=3)) 	
studies	#93

names(subset(table(field_data5a$study_id2), table(field_data5a$study_id2)>=8))
# or 76 with 8 studies

field_data6a = subset(field_data5a, study_id2 %in% studies)
nrow(field_data6a)

# which studies were deleted due to sample size?
setdiff(field_data5a$study_id2, field_data6a$study_id2)
# [1] "Bryony_Willcox_Macadamia_integrifolia_Australia"			# only 2 fields


# data processing: data subsets for main analyses----

# subset of studies that we can use for the main analysis (HB vs WI)
field_data8a = subset(field_data6a, !is.na(field_data6a$honeybee_z) & 
                        !is.na(field_data6a$wild_insects_z))
length(unique(field_data8a$study_id2))	# 93

nrow(field_data6a)
nrow(field_data8a) # no rows were removed at this step because NAs already dealt with

####################
# a subset of studies that provide richness data
field_data9a = subset(field_data8a, !is.na(richness_z))
length(unique(field_data9a$study_id2))	# 63

nrow(field_data8a)
nrow(field_data9a)
