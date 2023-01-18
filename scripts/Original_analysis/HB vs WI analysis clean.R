
# Note: this code file contains the main analyses for the "Wild insects and honey bees are equally important to crop yields in a global analysis" manuscript
# The goal is to use the CropPol database, which contains observations of insects visits and crop yields across many crops, to compare the importance of honey bee and wild insect visits to yield, and also compare the importance of total visits to richness
# this version was last edited on 12-5-2022 for the current version of the manuscript

# 2 csv files are required:
# 1) CropPol_field_level_data downloaded 10-25-2022.csv
# 2) field_data_summary_edited_v3.csv

# The general outline is:
#1) data processing of the CropPol database to get it ready for analysis (this version based on a downloaded csv of the database)
#2) plotting the proportion of visits by each group (figure 1) (starting on line 630)
#3) running mixed models, performing model selection (table 1,2,3), and plotting mixed model results (fig 2), starting on line 650
#4) sensitivity of results by subsampling and making funnel plots (fig 3), starting on line 1008

####################
# load necessary packages 

#IMPORTANT: I am commenting out tidyverse and MuMIn packages
# to avoid loading them as dependencies when using Renv. 
# You need to uncomment them to load them manually.
library(lme4)
#library(tidyverse)
#library(MuMIn)			# for AICc and model averaging
#library(optimx)		# for alternative convergence algorithms in lme4 (if needed)

####################
# load csv of data downloaded from CropPol project on github.  If we get the data directly from an internet link, we must consider whether any changes would break the processing steps below.

#field_data_old = read.csv(file="Reilly_models/CropPol_field_level_data downloaded 10-25-2022.csv", header=T) #file not uploaded to github!
#str(field_data)
# The version used in the published analysis should be: https://raw.githubusercontent.com/ibartomeus/OBservData/bf1deaababecb8c3e6ac2c674ae57f55c38cd69f/Final_Data/CropPol_field_level_data.csv
field_data <- read.csv("https://raw.githubusercontent.com/ibartomeus/OBservData/bf1deaababecb8c3e6ac2c674ae57f55c38cd69f/Final_Data/CropPol_field_level_data.csv")

#dim(field_data_old)
#dim(field_data)
#all.equal(field_data_old, field_data) #Yes, equal versions.

##################################################
# data processing: study IDs
##################################################

# in order to maximize the sample size within each study, I keep all the years of each study together.  CropPol currently gives each study_year a separate study_id, so I make a "study_id2" column that lumps years within study, and then fix various inconsistencies
field_data$study_id2 = ""
for (i in 1:nrow(field_data)) {
	thisrow = field_data[i,]
	# get rid of the year using strsplit() to separate the differnt parts of studyID and then subsetting away the year part
	field_data$study_id2[i] = paste(subset(unlist(strsplit(thisrow$study_id, split='_')), unlist(strsplit(thisrow$study_id, split='_'))!=thisrow$sampling_year), collapse="_")
}

# fix naming issues that didn't process correctly above
field_data[which(field_data$study_id2=="Alexandra_Maria_Klein_Coffea_arabica_Indonesia_2000_2001"),]$study_id2 <- "Alexandra_Maria_Klein_Coffea_arabica_Indonesia"
field_data[which(field_data$study_id2=="Alexandra_Maria_Klein_Coffea_canephora_Indonesia_2000_2001"),]$study_id2 <- "Alexandra_Maria_Klein_Coffea_canephora_Indonesia"
field_data[which(field_data$study_id2=="Alice_Classen_Coffea_arabica_Tanzania_2011_2012"),]$study_id2 <- "Alice_Classen_Coffea_arabica_Tanzania"
field_data[which(field_data$study_id2=="Bryony_Willcox_Mangifera_indica_Australia_2"),]$study_id2 <- "Bryony_Willcox_Mangifera_indica_Australia"
field_data[which(field_data$study_id2=="Davi_L_Ramos_Phaseolus_vulgaris L_Brazil_2015_2016"),]$study_id2 <- "Davi_L_Ramos_Phaseolus_vulgaris L_Brazil"
field_data[which(field_data$study_id2=="Margaret_Mayfield_Actinidia_deliciosa_New_Zealand_NA"),]$study_id2 <- "Margaret_Mayfield_Actinidia_deliciosa_New_Zealand"

# there are some "studies" where the different years come from different papers.
# keep these separate if the methods differed by year, otherwise lump them

# Dara_Stanley_Brassica_napus_Ireland
# 2009 and 2010 data come from different papers using differnt methods, so I will kept separate
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
# 2004-2007 are net data, 2008-12 are the visit data.  Note: these won't be used in the main analysis anyway since there is no associated yield data.
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

##################################################
# data processing: drop data that we can't use
##################################################

# drop pan trap data
field_data2a = subset(field_data, !study_id2 %in% c(
	"Breno_M_Freitas_Gossypium_hirsutum_Brazil",
	"Georg_Andersson_Fragaria_ananassa_Sweden",
	"Heather_Lee_Grab_Fragaria_ananassa_USA",
	"Rachel_Mallinger_Malus_domestica_USA",
	"Yi_Zou_Brassica_napus_China"
))

# drop enclosure data
field_data2a = subset(field_data2a, !study_id2 %in% c(	
	"Willis_Chan_Raine_Cucurbita_pepo_Canada"
))	

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
))

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

##################################################
# data processing: visit rates and abundance
##################################################

# In the CropPol database, insect visitors are recorded either as "abundance" or "visitation_rate".  Abundance is typically = net data, and visitation_rate = timed observations.  
# For this analysis, we will use visitation rate if it is available and the data looks reasonable, and otherwise use the abundance data.

# I wrote a loop to do this (not included in this file), but then I had to change a few studies manually after checking the data, so it is currently being done with a manually created csv.
# In particular:
# I use abundance instead of visits for the Garratt studies because it looks like the visit data might not include all pollinators, whereas the abundance data does
# I use visits instead of abundance for Taylor_Ricketts_Coffea_arabica_Costa_Rica because the visit data has only 4 fewer rows and that is probably worth it
# I use visits instead of abundance for  Thijs_Fijen_Allium_porrum_France, Thijs_Fijen_Allium_porrum_Italy, and Virginie_Boreux_Coffea_canephora_India because of potential data issues in the visit data

# read in a csv that lists whether to use the visitation_rate or abundance column for each study
field_data_summary_edited = read.csv(file="scripts/field_data_summary_edited_v3.csv", header=T)
field_data_summary_edited

# After going back to the relevant papers, I also made some manual edits to the richness_restriction column in the CropPol database, which lists which insect groups were recorded for each study. These edits are in the "taxa_recorded_edited" column in this csv.  
# In particular:
# Saul_A_Cunningham_Brassica_napus_Australia
	# from the paper: "Pollinators were counted in three categories: feral honeybees, hover flies and native bees."  I changed the category from NA to "bees and syrphids"
# Blande_Viana_Passiflora_edulis_Brazil
	# this data comes from a thesis and was also in Ricketts et al 2008.  It is clearly bees only, so I changed the category from NA to "bees"
# Alejandro_Trillo_Fragaria_ananassa_Spain:
	# the methods in the paper indicate that they used different methods to sample bees and other insects. This is a problem for comparisons, so I changed the category from "all visitors considered" to "bees", so that only the bees (including honey bee) part will be used. 

# merge this csv with the field data
field_data2b = merge(field_data2a, field_data_summary_edited, by="study_id2")

# check to make sure we still have the correct number of rows
nrow(field_data2a)	# 2993
nrow(field_data2b)	# 2993

#setdiff(unique(field_data2a$study_id2), unique(field_data2b$study_id2))


##################################################
# data processing: clean up 0's and NAs
##################################################

# there is some inconsistency in the database about whether missing insect data is recorded as 0s or NAs
# This can be a serious issue, especially when using na.rm=T, as we do in some of the functions below.
# To be safe, the following code will go through and manually reset NAs and 0's depending on which taxa were actually recorded based on the "taxa_recorded_edited" column

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

# ***JR code note: what I do in the loop below is subsetting a single row from the data, doing some processing on it, then rebuilding the dataframe from the processed rows using rbind.  There are various other ways this can be accomplished, but I find this method more intuitive and easier to troubleshoot.

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

####################

# note that about 100 rows were removed due to NAs in either visitation rate or abundance at this step, which is as expected
# most of the missing data are NAs scattered across many of the studies.  Of special note is "Thijs_Fijen_Allium_porrum_France" which has more NAs than usual in the visit data, but it still looks better than the abundance data which appears to have an issue with duplicated values.

# checking which studies it was:
unique(subset(field_data2b, site_id %in% setdiff(field_data2b$site_id, field_data2c$site_id))$study_id2)

nrow(field_data2b)
nrow(field_data2c)	# 2854


###################################################
# data processing: Wild insect and all insect sums
###################################################

# in the analysis, we are summing up the various recorded groups into "wild insects" and maybe some other groups
# group sums that we might want:
# honeybee
# wild_insects
# wild_bees
# non_bees
# all_bees
# all_insects

# functions to do the sums for abundance data (note na.rm=T, which will ignore NAs when making the sums):
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

####################
# make a new column for each insect group that pulls either visits or abundance depending on what was chosen above

# define the list of studies
studies = unique(field_data2c$study_id2)

# ***JR code note: again, looping through subsets of the data, processing them, and rebuilding the dataset using rbind

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

##################################################
# data processing: yield data
##################################################

# fix yield data for ICP datasets to match Reilly et al 2020 - it seems at some point the yield vs yield2 were reversed for these:
# we can drop this part if we make the edits to CropPol

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


##################################################
# data processing: transform to z-scores
##################################################

# note: for data that is all zeros, (e.g. when no non-bees were recorded for a crop), it is not possible to calculate z-scores because there is no variance (divide by zero).  This happens in several crops:
#ICP_Prunus_dulcis_USA-CA (wild bees)
#ICP_Vaccinium_corymbosum_USA-OR (non-bees)
#Mark_Otieno_Cajanus_cajan_Kenya (non-bees)
#Virginie_Boreux_Coffea_canephora_India (non-bees)

# note on NAs for richness z-scores: scale() removes NAs when calculating the means.  I checked that it was not introducing NAs for every row in a crop when there are only a couple NAs--it was not.  This was only an issue for Natacha_Chacoff_Citrus_paradisi_Argentina, which had NAs for richness in 4/12 fields.  Most or all of the other NAs in richness were already removed for missing other data. 

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

####################
# remove rows with NA in the yield z-score
field_data5 = subset(field_data4, !is.na(field_data4$yield_z))

# check number of rows

nrow(field_data3)
nrow(field_data4)
nrow(field_data5) 

# check number of studies
length(unique(field_data4$study_id2))
length(unique(field_data5$study_id2)) # Annona squamosa and muricata dropped because only 1 field each and can't calculate z scores (creates NAs) for yield (or insects either, actually)
# 94

####################
# replace NAs in HB due to divide by zero in z-score with manually entered 0's.  I am currently thinking that this is better than just dropping them from model runs.
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


##################################################
# data processing: sample size
##################################################

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


##################################################
# data processing: data subsets for main analyses
##################################################

# subset of studies that we can use for the main analysis (HB vs WI)
field_data8a = subset(field_data6a, !is.na(field_data6a$honeybee_z) & !is.na(field_data6a$wild_insects_z))
length(unique(field_data8a$study_id2))	# 93

nrow(field_data6a)
nrow(field_data8a) # no rows were removed at this step because NAs already dealt with

####################
# a subset of studies that provide richness data
field_data9a = subset(field_data8a, !is.na(richness_z))
length(unique(field_data9a$study_id2))	# 63

nrow(field_data8a)
nrow(field_data9a)


#############################################################################
# analysis: proportion of visits by HB vs WI, or HB, wild bees, and non bees
#############################################################################

# fraction of visits by WI boxplot (currently manuscript figure 1)

# calculate fraction of "visits" by each group
field_data8a$fracWI = field_data8a$wild_insects/field_data8a$all_insects
field_data8a$fracHB = field_data8a$honeybee/field_data8a$all_insects

# reorder factor levels of study so they will plot in order of increasing median when we run a boxplot
field_data8a$study_id2 <- reorder(field_data8a$study_id2, field_data8a$fracWI, median, na.rm=T)

#png(filename = "wild insect visits by crop boxplots 11-16-2022.png", width = 9, height = 5, units = "in", bg = "white", res = 600)
par(mar=c(4,4,4,4))
boxplot(data=field_data8a, fracWI ~ study_id2, col=c("light gray"), xaxt = 'n', xlab="study", ylab="proportion of visits by wild insects")
axis(1,seq(5,90,5))
#dev.off()


##################################################
# analysis: mixed models for HB vs WI
##################################################

# for this we are using field_data8a, which has HB and WI, and where we have entered manual 0's for HB z-scores in 2 systems where they were =0.
# note that some of these studies are bees only, and many do not have richness data

# note: all of the following models will return a warning of "boundary (singular) fit: see ?isSingular".  I think this is due to using z-scores, which make the mean of all studies =0, so it should not be a problem for the models.
# note: REML=F versions for AIC values in table 1

mod1 = lmer(data=field_data8a, yield_z ~ (1|study_id2), na.action=na.fail, REML=F)
AICc(mod1)

mod2 = lmer(data=field_data8a, yield_z ~ wild_insects_z + (1 + wild_insects_z|study_id2), na.action=na.fail, REML=F)
AICc(mod2)

mod3 = lmer(data=field_data8a, yield_z ~ honeybee_z + (1 + honeybee_z|study_id2), na.action=na.fail, REML=F)
AICc(mod3)

mod5 = lmer(data=field_data8a, yield_z ~ wild_insects_z + honeybee_z + (1 + wild_insects_z|study_id2) + (1 + honeybee_z|study_id2), na.action=na.fail, REML=F)
AICc(mod5)

mod9 = lmer(data=field_data8a, yield_z ~ wild_insects_z*honeybee_z + (1 + wild_insects_z|study_id2) + (1 + honeybee_z|study_id2) + (1 + wild_insects_z*honeybee_z|study_id2), na.action=na.fail, REML=F)
AICc(mod9)

# model selection table:
model.sel(mod1,mod2,mod3,mod5,mod9)

# REML versions for coefficient estimates in table 1
update(mod1, REML=T)
update(mod2, REML=T)
update(mod3, REML=T)
update(mod5, REML=T)
update(mod9, REML=T)


####################
# Akaike weights and relative importance values for table 2.

# deltaAIC values
bestAIC = AICc(mod5)
d1 = AICc(mod1)-bestAIC
d2 = AICc(mod2)-bestAIC
d3 = AICc(mod3)-bestAIC
d5 = AICc(mod5)-bestAIC
d9 = AICc(mod9)-bestAIC

# likelihood values
L1 = exp(-.5*d1)
L2 = exp(-.5*d2)
L3 = exp(-.5*d3)
L5 = exp(-.5*d5)
L9 = exp(-.5*d9)

# Akaike weights - these should match the weights in the model selection table above
w1 = L1/(L1+L2+L3+L5+L9)
w2 = L2/(L1+L2+L3+L5+L9)
w3 = L3/(L1+L2+L3+L5+L9)
w5 = L5/(L1+L2+L3+L5+L9)
w9 = L9/(L1+L2+L3+L5+L9)

# relative importance values for each variable = sum of akaike weights of models containing that variable

# WI
w2+w5+w9		# 0.999853

# HB
w3+w5+w9		# 0.9999984

##################################################
# analysis: mixed models for HB vs WI vs Richness
##################################################

# for this we are using field_data9a, which is the subset of studies providing data for HB, WI, and Richness

# I didn't get any convergence warnings when I ran these with the current data, but if it does happen, try adding one of these to the lmer statement:
# control = lmerControl(optimizer="Nelder_Mead")
# control = lmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B'))

# note that I'm using 0+ in the random effect instead of 1+ since all studies are z-scored and thus constrained to have a mean of 0. My logic here is that there is no need to try to estimate separate intercepts for each study. We are still estimating the random slopes as usual.

modR_1 = lmer(data=field_data9a, yield_z ~  (1|study_id2), na.action=na.fail, REML=F)
AICc(modR_1)

modR_2 = lmer(data=field_data9a, yield_z ~ wild_insects_z + (0 + wild_insects_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_2)

modR_3 = lmer(data=field_data9a, yield_z ~ honeybee_z + (0 + honeybee_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_3)

modR_4 = lmer(data=field_data9a, yield_z ~ richness_z + (0 + richness_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_4)

modR_5 = lmer(data=field_data9a, yield_z ~ wild_insects_z + honeybee_z + (0 + wild_insects_z|study_id2) + (0 + honeybee_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_5)

modR_6 = lmer(data=field_data9a, yield_z ~ wild_insects_z + richness_z + (0 + wild_insects_z|study_id2) + (0 + richness_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_6)

modR_7 = lmer(data=field_data9a, yield_z ~ honeybee_z + richness_z + (0 + honeybee_z|study_id2) + (0 + richness_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_7)

modR_8 = lmer(data=field_data9a, yield_z ~ wild_insects_z + honeybee_z + richness_z + (0 + wild_insects_z|study_id2) + (0 + honeybee_z|study_id2) + (0 + richness_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_8)

modR_9 = lmer(data=field_data9a, yield_z ~ wild_insects_z*honeybee_z + (0 + wild_insects_z*honeybee_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_9)

modR_10 = lmer(data=field_data9a, yield_z ~ wild_insects_z*richness_z + (0 + wild_insects_z*richness_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_10)

modR_11 = lmer(data=field_data9a, yield_z ~ honeybee_z*richness_z + (0 + honeybee_z*richness_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_11)

modR_12 = lmer(data=field_data9a, yield_z ~ wild_insects_z*honeybee_z + richness_z + (0 + wild_insects_z*honeybee_z|study_id2) + (0 + richness_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_12)

modR_13 = lmer(data=field_data9a, yield_z ~ wild_insects_z*richness_z + honeybee_z + (0 + wild_insects_z*richness_z|study_id2) + (0 + honeybee_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_13)

modR_14 = lmer(data=field_data9a, yield_z ~ honeybee_z*richness_z + wild_insects_z + (0 + honeybee_z*richness_z|study_id2) + (0 + wild_insects_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_14)

modR_15 = lmer(data=field_data9a, yield_z ~ honeybee_z*richness_z + wild_insects_z*richness_z + (0 + honeybee_z*richness_z|study_id2) + (0 + wild_insects_z*richness_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_15)

modR_16 = lmer(data=field_data9a, yield_z ~ wild_insects_z*honeybee_z + honeybee_z*richness_z + (0 + wild_insects_z*honeybee_z|study_id2) + (0 + honeybee_z*richness_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_16)

modR_17 = lmer(data=field_data9a, yield_z ~ wild_insects_z*honeybee_z + wild_insects_z*richness_z + (0 + wild_insects_z*honeybee_z|study_id2) + (0 + wild_insects_z*richness_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_17)

modR_18 = lmer(data=field_data9a, yield_z ~ wild_insects_z*honeybee_z + wild_insects_z*richness_z + honeybee_z*richness_z + (0 + wild_insects_z*honeybee_z|study_id2) + (0 + wild_insects_z*richness_z|study_id2) +  (0 + honeybee_z*richness_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_18)


# all insects models
modR_8b = lmer(data=field_data9a, yield_z ~ all_insects_z + richness_z + (0 + all_insects_z|study_id2) + (0 + richness_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_8b)
summary(modR_8b)

modR_5b = lmer(data=field_data9a, yield_z ~ all_insects_z + (0 + all_insects_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_5b)
summary(modR_5b)

modR_15b = lmer(data=field_data9a, yield_z ~ all_insects_z + richness_z + all_insects_z*richness_z + (0 + all_insects_z|study_id2) + (0 + richness_z|study_id2) + (0 + all_insects_z*richness_z|study_id2), na.action=na.fail, REML=F)
AICc(modR_15b)
summary(modR_15b)

# model selection table:
model.sel(modR_1, modR_2, modR_3,modR_4,modR_5,modR_5b,modR_6,modR_7,modR_8,modR_8b,modR_9,modR_10,modR_11,modR_12,modR_13,modR_14,modR_15,modR_15b,modR_16,modR_17,modR_18)

# REML versions for coefficient estimates in table 3
fixef(update(modR_1, REML=T))
fixef(update(modR_2, REML=T))
fixef(update(modR_3, REML=T))
fixef(update(modR_4, REML=T))
fixef(update(modR_5, REML=T))
fixef(update(modR_5b, REML=T))
fixef(update(modR_6, REML=T))
fixef(update(modR_7, REML=T))
fixef(update(modR_8, REML=T))
fixef(update(modR_8b, REML=T))
fixef(update(modR_9, REML=T))
fixef(update(modR_10, REML=T))
fixef(update(modR_11, REML=T))
fixef(update(modR_12, REML=T))
fixef(update(modR_13, REML=T))
fixef(update(modR_14, REML=T))
fixef(update(modR_15, REML=T))
fixef(update(modR_15b, REML=T))
fixef(update(modR_16, REML=T))
fixef(update(modR_17, REML=T))
fixef(update(modR_18, REML=T))


####################
# Akaike weights and relative importance values for table 2.

bestAIC = AICc(modR_8b)
d1 = AICc(modR_1)-bestAIC
d2 = AICc(modR_2)-bestAIC
d3 = AICc(modR_3)-bestAIC
d4 = AICc(modR_4)-bestAIC
d5 = AICc(modR_5)-bestAIC
d6 = AICc(modR_6)-bestAIC
d7 = AICc(modR_7)-bestAIC
d8 = AICc(modR_8)-bestAIC
d9 = AICc(modR_9)-bestAIC
d10 = AICc(modR_10)-bestAIC
d11 = AICc(modR_11)-bestAIC
d12 = AICc(modR_12)-bestAIC
d13 = AICc(modR_13)-bestAIC
d14 = AICc(modR_14)-bestAIC
d15 = AICc(modR_15)-bestAIC
d16 = AICc(modR_16)-bestAIC
d17 = AICc(modR_17)-bestAIC
d18 = AICc(modR_18)-bestAIC

bestAIC = AICc(modR_8b)
d1b = AICc(modR_1)-bestAIC
d4b = AICc(modR_4)-bestAIC
d5b = AICc(modR_5b)-bestAIC
d8b = AICc(modR_8b)-bestAIC
d15b = AICc(modR_15b)-bestAIC

# likelihood values
L1 = exp(-.5*d1)
L2 = exp(-.5*d2)
L3 = exp(-.5*d3)
L4 = exp(-.5*d4)
L5 = exp(-.5*d5)
L6 = exp(-.5*d6)
L7 = exp(-.5*d7)
L8 = exp(-.5*d8)
L9 = exp(-.5*d9)
L10 = exp(-.5*d10)
L11 = exp(-.5*d11)
L12 = exp(-.5*d12)
L13 = exp(-.5*d13)
L14 = exp(-.5*d14)
L15 = exp(-.5*d15)
L16 = exp(-.5*d16)
L17 = exp(-.5*d17)
L18 = exp(-.5*d18)

L1b = exp(-.5*d1b)
L4b = exp(-.5*d4b)
L5b = exp(-.5*d5b)
L8b = exp(-.5*d8b)
L15b = exp(-.5*d15b)

# Akaike weights
w1 = L1/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w2 = L2/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w3 = L3/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w4 = L4/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w5 = L5/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w6 = L6/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w7 = L7/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w8 = L8/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w9 = L9/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w10 = L10/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w11 = L11/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w12 = L12/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w13 = L13/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w14 = L14/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w15 = L15/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w16 = L16/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w17 = L17/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)
w18 = L18/(L1+L2+L3+L4+L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15+L16+L17+L18)

w1b = L1b/(L1b+L4b+L5b+L8b+L15b)
w4b = L4b/(L1b+L4b+L5b+L8b+L15b)
w5b = L5b/(L1b+L4b+L5b+L8b+L15b)
w8b = L8b/(L1b+L4b+L5b+L8b+L15b)
w15b = L15b/(L1b+L4b+L5b+L8b+L15b)


# relative importance values for each variable = sum of akaike weights of models containing that variable

# model set 2a
# WI
w2+w5+w6+w8+w9+w10+w12+w13+w14+w15+w16+w17+w18	
# HB
w3+w5+w7+w8+w9+w11+w12+w13+w14+w15+w16+w17+w18	
# RICH
w4+w6+w7+w8+w10+w11+w12+w13+w14+w15+w16+w17+w18

# model set 2b
# AI
w5b + w8b + w15b
# RICH
w4b + w8b + w15b


####################
# model averaging 

# model set 2a
modav = model.avg(list(modR_1, modR_2, modR_3,modR_4,modR_5,modR_6,modR_7,modR_8,modR_9,modR_10,modR_11,modR_12,modR_13,modR_14,modR_15,modR_16,modR_17,modR_18))

# output--look at the "full average" coefficients which include shrinkage
summary(modav)

# model set 2b
# the all insects models
modavb = model.avg(list(modR_1, modR_4, modR_5b,modR_8b,modR_15b))

# output--look at the "full average" coefficients which include shrinkage
summary(modavb)



################################################################################
# plots: random effect estimates by study, similar to Garibaldi et al 2013 fig 2b, but using model 5
################################################################################

mod5_REML = update(mod5, REML=T)
summary(mod5_REML)

a=fixef(mod5_REML)
a_wi = a[2]		#overall slope of wild insects 
a_hb = a[3]		#overall slope of honey bee 
b=ranef(mod5_REML)
qq <- attr(ranef(mod5_REML, condVar = TRUE)[[1]], "postVar") # Extract the variances of the random effects
qq_wi = qq[[1]]
qq_hb = qq[[2]]
e_wi=(sqrt(qq_wi))	
e_wi=e_wi[2,2,]		
e_hb=(sqrt(qq_hb))	
e_hb=e_hb[2,2,]		
b_wi = b[[1]][2]	# slope of wild insects by crop?  lmer calls this the conditional mean
b_hb = b[[1]][4]	# slope of hb by crop?  lmer calls this the conditional mean
mean_wi=(b_wi+a_wi)			# mean ...  add by crop part to overall slope
ci_lower_wi=(b_wi+a_wi)-(e_wi*2)	# lower CI
ci_upper_wi=(b_wi+a_wi)+(e_wi*2)	# upper CI
mean_hb=(b_hb+a_hb)			# mean ...  add by crop part to overall slope
ci_lower_hb=(b_hb+a_hb)-(e_hb*2)	# lower CI
ci_upper_hb=(b_hb+a_hb)+(e_hb*2)	# upper CI

###
modelout = cbind(mean_wi, ci_lower_wi, ci_upper_wi, mean_hb, ci_lower_hb, ci_upper_hb)
colnames(modelout)[1] = "mean_wi"
colnames(modelout)[2] = "ci_lower_wi"
colnames(modelout)[3] = "ci_upper_wi"
colnames(modelout)[4] = "mean_hb"
colnames(modelout)[5] = "ci_lower_hb"
colnames(modelout)[6] = "ci_upper_hb"

modelout$crop = rownames(modelout)
modelout

modelout_mod5 = modelout


# plot figure 2
#png(filename = "slopes by crop MODEL 5 10-25-2022.png", width = 10, height = 5, units = "in", bg = "white", res = 600)

par(mfrow=c(1,2))

modelout_mod5 = arrange(modelout_mod5, mean_wi)
plot(1:nrow(modelout_mod5), modelout_mod5$mean_wi, pch=19, col="#00bb5c", ylim=c(-.75,.75), ylab="wild insect slope (mean and 95% CI)", cex=1.5, xaxt = "n", xlab="crop study", main="Wild Insects")
axis(1, at=seq(0,nrow(modelout_mod5),10), labels=seq(0,nrow(modelout_mod5),10))
arrows(1:nrow(modelout_mod5), modelout_mod5$ci_lower_wi, 1:nrow(modelout_mod5), modelout_mod5$ci_upper_wi, length=0.05, angle=90, code=3, col="#00bb5c", lwd=1.5)
abline(h=0.060, lty="dashed", lwd=2)
abline(h=0, lwd=1, col="black", lty="dotted")

modelout_mod5 = arrange(modelout_mod5, mean_hb)
plot(1:nrow(modelout_mod5), modelout_mod5$mean_hb, pch=19, col="#c96d00", ylim=c(-.75,.75), ylab="honey bee slope (mean and 95% CI)", cex=1.5, xaxt = "n", xlab="crop study", main="Honey Bee")
axis(1, at=seq(0,nrow(modelout_mod5),10), labels=seq(0,nrow(modelout_mod5),10))
arrows(1:nrow(modelout_mod5), modelout_mod5$ci_lower_hb, 1:nrow(modelout_mod5), modelout_mod5$ci_upper_hb, length=0.05, angle=90, code=3, col="#c96d00", lwd=1.5)
abline(h=0.070, lty="dashed", lwd=2)
abline(h=0, lwd=1, col="black", lty="dotted")

#dev.off()


################################################################################
# sensitivity plots: subsampling and making funnel plots
################################################################################

# this section creates funnel plots showing how the overall mean slope of wild insects and HB converges as the number of studies increases

studylist = unique(field_data8a$study_id2)

# create empty dataframe
finalout = data.frame()

# number of replicates for subsampling
# warning, looks really nice with 1000 reps, but this will take quite a while
for (rep in 1:100) {

	# loop over subsamples of 2 up to all 89 studies
	for (ncrops in 2:93) {

		# draw the subsampled studies
		thissample = sample(studylist, ncrops)
		field_data_this = subset(field_data8a, study_id2 %in% thissample)

		# using model 5 here, which has both WI and HB
		mod5 = lmer(data=field_data_this, yield_z ~ wild_insects_z + honeybee_z + (1 + wild_insects_z|study_id2) + (1 + honeybee_z|study_id2), na.action=na.fail, REML=T)

		########################
		# extract random effects and calculate CIs
		a=fixef(mod5)
		a_wi = a[2]		#overall slope of wild insects 
		a_hb = a[3]		#overall slope of honey bee 
		b=ranef(mod5)
		qq <- attr(ranef(mod5, condVar = TRUE)[[1]], "postVar") # Extract the variances of the random effects
		qq_wi = qq[[1]]
		qq_hb = qq[[2]]
		e_wi=(sqrt(qq_wi))	
		e_wi=e_wi[2,2,]		
		e_hb=(sqrt(qq_hb))	
		e_hb=e_hb[2,2,]		
		b_wi = b[[1]][2]	# slope of wild insects by crop
		b_hb = b[[1]][4]	# slope of hb by crop
		mean_wi=(b_wi+a_wi)			# mean ...  add by crop part to overall slope
		ci_lower_wi=(b_wi+a_wi)-(e_wi*2)	# lower CI
		ci_upper_wi=(b_wi+a_wi)+(e_wi*2)	# upper CI
		mean_hb=(b_hb+a_hb)			# mean ...  add by crop part to overall slope
		ci_lower_hb=(b_hb+a_hb)-(e_hb*2)	# lower CI
		ci_upper_hb=(b_hb+a_hb)+(e_hb*2)	# upper CI

		modelout = cbind(mean_wi, ci_lower_wi, ci_upper_wi, mean_hb, ci_lower_hb, ci_upper_hb)
		colnames(modelout)[1] = "mean_wi"
		colnames(modelout)[2] = "ci_lower_wi"
		colnames(modelout)[3] = "ci_upper_wi"
		colnames(modelout)[4] = "mean_hb"
		colnames(modelout)[5] = "ci_lower_hb"
		colnames(modelout)[6] = "ci_upper_hb"

		modelout$crop = rownames(modelout)
		modelout

		modelout$range_wi = modelout$ci_upper_wi - modelout$ci_lower_wi
		modelout$range_hb = modelout$ci_upper_hb - modelout$ci_lower_hb

		out1 = data.frame(model="mod5", rep = rep, ncrops = ncrops, wi_slope = a_wi, hb_slope=a_hb, error_width_wi = mean(modelout$range_wi), error_width_hb = mean(modelout$range_hb))
		finalout = rbind(out1, finalout)

	}
}

#finalout


################

# this section creates funnel plots showing how the overall mean slope of ALL INSECTS and RICHNESS converges as the number of studies increases

studylist = unique(field_data9a$study_id2)

# create empty dataframe
finalout2 = data.frame()

# number of replicates for subsampling
# warning, looks really nice with 1000 reps, but this will take quite a while
for (rep in 1:100) {

	# loop over subsamples of 2 up to all 89 studies
	for (ncrops in 2:63) {

		# draw the subsampled studies
		thissample = sample(studylist, ncrops)
		field_data_this = subset(field_data9a, study_id2 %in% thissample)

		# using model 5 here, which has both WI and HB
		modR_8b = lmer(data=field_data_this, yield_z ~ all_insects_z + richness_z + (1 + all_insects_z|study_id2) + (1 + richness_z|study_id2), na.action=na.fail, REML=T)
		

		########################
		# extract random effects and calculate CIs
		a=fixef(modR_8b)
		a_ai = a[2]		#overall slope of aild insects 
		a_rich = a[3]		#overall slope of honey bee 
		b=ranef(modR_8b)
		qq <- attr(ranef(modR_8b, condVar = TRUE)[[1]], "postVar") # Extract the variances of the random effects
		qq_ai = qq[[1]]
		qq_rich = qq[[2]]
		e_ai=(sqrt(qq_ai))	
		e_ai=e_ai[2,2,]		
		e_rich=(sqrt(qq_rich))	
		e_rich=e_rich[2,2,]		
		b_ai = b[[1]][2]	# slope of aild insects by crop
		b_rich = b[[1]][4]	# slope of rich by crop
		mean_ai=(b_ai+a_ai)			# mean ...  add by crop part to overall slope
		ci_lower_ai=(b_ai+a_ai)-(e_ai*2)	# lower CI
		ci_upper_ai=(b_ai+a_ai)+(e_ai*2)	# upper CI
		mean_rich=(b_rich+a_rich)			# mean ...  add by crop part to overall slope
		ci_lower_rich=(b_rich+a_rich)-(e_rich*2)	# lower CI
		ci_upper_rich=(b_rich+a_rich)+(e_rich*2)	# upper CI

		modelout = cbind(mean_ai, ci_lower_ai, ci_upper_ai, mean_rich, ci_lower_rich, ci_upper_rich)
		colnames(modelout)[1] = "mean_ai"
		colnames(modelout)[2] = "ci_lower_ai"
		colnames(modelout)[3] = "ci_upper_ai"
		colnames(modelout)[4] = "mean_rich"
		colnames(modelout)[5] = "ci_lower_rich"
		colnames(modelout)[6] = "ci_upper_rich"

		modelout$crop = rownames(modelout)
		modelout

		modelout$range_ai = modelout$ci_upper_ai - modelout$ci_lower_ai
		modelout$range_rich = modelout$ci_upper_rich - modelout$ci_lower_rich

		out1 = data.frame(model="modR_8b", rep = rep, ncrops = ncrops, ai_slope = a_ai, rich_slope=a_rich, error_width_ai = mean(modelout$range_ai), error_width_rich = mean(modelout$range_rich))
		finalout2 = rbind(out1, finalout2)

	}
}

#finalout2


################
# plotting for figure 3

#png(filename = "funnel plots 11-2-2022.png", width = 7, height = 7, units = "in", bg = "white", res = 600)

par(mfrow=c(2,2), mar = c(6, 4, 0.1, 1), oma=c(0,1,1,1))

########################
# summarize results by number of studies
summary_by_ncrops <- finalout %>%
		group_by(ncrops) %>%
        summarize(
				slope_mean_wi = mean(wi_slope),
                slope_CIlower_wi = quantile(wi_slope, .025),
                slope_CIupper_wi = quantile(wi_slope, .975),	
				width_mean_wi = mean(error_width_wi),
	            width_CIlower_wi = quantile(error_width_wi, .025),
                width_CIupper_wi = quantile(error_width_wi, .975),
				
				slope_mean_hb = mean(hb_slope),
                slope_CIlower_hb = quantile(hb_slope, .025),
                slope_CIupper_hb = quantile(hb_slope, .975),	
				width_mean_hb = mean(error_width_hb),
	            width_CIlower_hb = quantile(error_width_hb, .025),
                width_CIupper_hb = quantile(error_width_hb, .975),
                ) 
data.frame(summary_by_ncrops)

########################
# make funnel plots

plot(finalout$ncrops, finalout$wi_slope, ylab="wild insects effect size (slope)", xlab="number of studies", col="#00bb5c", ylim=c(-.4, .6))
abline(h=0, lwd=2, col="gray", lty="dotted")
points(summary_by_ncrops$ncrops, summary_by_ncrops$slope_CIlower_wi, type="l", lwd=2, lty="dotted")
points(summary_by_ncrops$ncrops, summary_by_ncrops$slope_CIupper_wi, type="l", lwd=2, lty="dotted")

# G 2013
points(30, 0.322819, type="p", pch=19, cex=1.25, col="black")	# WI by itself model F
points(20, 0.288398, type="p", pch=19, cex=1.25, col="black")	# WI and HB model P (best model)

# rader 2016
points(19, 0.187, type="p", pch=19, cex=1.25, col="black")		# other bees
points(19, 0.12, type="p", pch=19, cex=1.25, col="black")		# non bees

abline(h=0.050, lty="dashed", lwd=2, col="black") 				# our estimate from this model 5


########

plot(finalout$ncrops, finalout$hb_slope, ylab="honey bee effect size (slope)", xlab="number of studies", col="#c96d00", ylim=c(-.4, .6))
abline(h=0, lwd=2, col="gray", lty="dotted")
points(summary_by_ncrops$ncrops, summary_by_ncrops$slope_CIlower_hb, type="l", lwd=2, lty="dotted")
points(summary_by_ncrops$ncrops, summary_by_ncrops$slope_CIupper_hb, type="l", lwd=2, lty="dotted")

# G 2013
points(22, 0.15127, type="p", pch=19, cex=1.25, col="black")		# HB by itself (model M)
points(20, 0.149216, type="p", pch=19, cex=1.25, col="black")	# WI and HB model P (best model)

# rader 2016
points(19, -.019, type="p", pch=19, cex=1.25, col="black")		# HB

abline(h=0.078, lty="dashed", lwd=2, col="black") 				# estimate from this model 5


########################
# summarize results by number of studies
summary_by_ncrops <- finalout2 %>%
		group_by(ncrops) %>%
        summarize(
				slope_mean_ai = mean(ai_slope),
                slope_CIlower_ai = quantile(ai_slope, .025),
                slope_CIupper_ai = quantile(ai_slope, .975),	
				width_mean_ai = mean(error_width_ai),
	            width_CIlower_ai = quantile(error_width_ai, .025),
                width_CIupper_ai = quantile(error_width_ai, .975),
				
				slope_mean_rich = mean(rich_slope),
                slope_CIlower_rich = quantile(rich_slope, .025),
                slope_CIupper_rich = quantile(rich_slope, .975),	
				width_mean_rich = mean(error_width_rich),
	            width_CIlower_rich = quantile(error_width_rich, .025),
                width_CIupper_rich = quantile(error_width_rich, .975),
                ) 
data.frame(summary_by_ncrops)


########################
# make funnel plots

plot(finalout2$ncrops, finalout2$ai_slope, ylab="all insects effect size (slope)", xlab="number of studies", col="red", ylim=c(-.4, .6))
abline(h=0, lwd=2, col="gray", lty="dotted")
points(summary_by_ncrops$ncrops, summary_by_ncrops$slope_CIlower_ai, type="l", lwd=2, lty="dotted")
points(summary_by_ncrops$ncrops, summary_by_ncrops$slope_CIupper_ai, type="l", lwd=2, lty="dotted")

# G 2015
points(33, 0.3, type="p", pch=19, cex=1.25, col="black")		# visits effect in "best model" table S3

# Dainese 2019
points(42, 0.08, type="p", pch=19, cex=1.25, col="black")	# total abundance --> pollination

abline(h=0.104, lty="dashed", lwd=2, col="black") 			# our estimate from this model R8b

########

plot(finalout2$ncrops, finalout2$rich_slope, ylab="richness effect size (slope)", xlab="number of studies", col="#00a2ca", ylim=c(-.4, .6))
abline(h=0, lwd=2, col="gray", lty="dotted")
points(summary_by_ncrops$ncrops, summary_by_ncrops$slope_CIlower_rich, type="l", lwd=2, lty="dotted")
points(summary_by_ncrops$ncrops, summary_by_ncrops$slope_CIupper_rich, type="l", lwd=2, lty="dotted")

# G 2013
points(18, 0.17, type="p", pch=19, cex=1.25, col="black")	# richness by itself (model B)

# G 2015
points(33, 0.07, type="p", pch=19, cex=1.25, col="black")	# richness effect in "best model" table S3

# Dainese 2019
points(42, 0.106, type="p", pch=19, cex=1.25, col="black")	# richness --> pollination

abline(h=0.043, lty="dashed", lwd=2, col="black") 			# our estimate from this model R8b


#dev.off()


















