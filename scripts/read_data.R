# This script reads data from CropPol and cleans it for analysis.
# clean data is stored in /data

# read data from CropPol----

field_data <- read.csv("https://raw.githubusercontent.com/ibartomeus/OBservData/master/Final_Data/CropPol_field_level_data.csv")

# data processing: study IDs----

# in order to maximize the sample size within each study, we keep all the years of each study together.  
# CropPol currently gives each study_year a separate study_id, so I make a "study_id2" column that lumps years within study, and then fix various inconsistencies
## IB note: I wonder if keeping years separated makes also sense

unique(field_data$study_id2)	# 144 studies as of 1/27/2023

# data processing: drop data that we can't use----

unique(field_data2a$use_visits_or_abundance)

# drop pan trap data
field_data2a = subset(field_data, !use_visits_or_abundance %in% c(
  "drop (pan trap)"
)) 

# drop enclosure data
field_data2a = subset(field_data2a, !use_visits_or_abundance %in% c(
  "drop (enclosure)"
))	

# drop for missing HB data (HB data not recorded)	
field_data2a = subset(field_data2a, !use_visits_or_abundance %in% c(	
  "drop (HB not recorded)"
)) 

# drop studies with unclear methods for now, unless we learn that they were ok--see below
field_data2a = subset(field_data2a, !use_visits_or_abundance %in% c(
  "drop (methods)"
))

# notes on why these 3 studies were dropped:
# Breno_M_Freitas_Bixa_orellana_Brazil
# abundance method listed as "other" in the database, so it might be pan trap data, etc, and it also might be bees only. The paper listed as Dainese et al 2019, but I can't find it in the data provided
# Johan_Ekroos_Vicia_faba_Sweden
# taxa listed as "only bees (bumblebees)"  does this mean that HB were not sampled, or were there just not any?  Were other wild bees than bumble bees sampled?  There is no paper associated with this data
# Rebecca_Steward_Fragaria_ananassa_Sweden
# there is only data for bees in the wild bees column.  The paper says HB were sampled, but most bees were bumblebees.  But there are 0 HB in the data--where are they?  

unique(field_data2a$study_id2)	# 125 as of 1/27/2023

# data processing: visit rates and abundance----

# In the CropPol database, insect visitors are recorded either as "abundance" or "visitation_rate".  
# Abundance is typically = net data, and visitation_rate = timed observations.  
# For this analysis, we will use visitation rate if it is available and the data looks reasonable, and otherwise use the abundance data.

# This decision is recorded in the "use_visits_or_abundance" column rather than being made automatically to allow some subjectivity.
# In particular:
# I use abundance instead of visits for the Garratt studies because it looks like the visit data might not include all pollinators, whereas the abundance data does
# I use visits instead of abundance for Taylor_Ricketts_Coffea_arabica_Costa_Rica because the visit data has only 4 fewer rows and that is probably worth it
# I use visits instead of abundance for  Thijs_Fijen_Allium_porrum_France, Thijs_Fijen_Allium_porrum_Italy, and Virginie_Boreux_Coffea_canephora_India because of potential data issues in the visit data
unique(field_data2a$use_visits_or_abundance)

# check to make sure we still have the correct number of rows
nrow(field_data2a) # 3009 as of 1/27/2023

# data processing: clean up 0's and NAs----

# there is some inconsistency in the database about whether missing insect data is recorded as 0s or NAs
# This can be a serious issue, especially when using na.rm=T, as we do in some of the functions below.
# To be safe, the following code will go through and manually reset NAs and 0's depending on which taxa 
# were actually recorded based on the "taxa_recorded" column

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

# ***JR code note: what I do in the loop below is subsetting a single row from the data, doing some processing on it, then rebuilding the dataframe from the processed rows using rbind. 
# There are various other ways this can be accomplished, but I find this method more intuitive and easier to troubleshoot.

# set up a new empty dataframe
field_data2c = data.frame()

# loop through each row, replacing NAs and 0's depending on which taxa were listed as being recorded
for (i in 1:nrow(field_data2a)) {
  thisrow = field_data2a[i,]
  # if this is an abundance study, we need to do this for the abundance columns (ab_)
  if (thisrow$use_visits_or_abundance=="use abundance") {
    if (!is.na(thisrow$abundance)) {
      if (thisrow$taxa_recorded=="bees") {
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
      } else if (thisrow$taxa_recorded=="all visitors") {
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
      } else if (thisrow$taxa_recorded=="bees and leps") { 
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
      } else if (thisrow$taxa_recorded=="bees and syrphids") { 
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
      } else if (thisrow$taxa_recorded=="bees, syrphids and wasps") { 
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
  } else if (thisrow$use_visits_or_abundance=="use visits") {
    if (!is.na(thisrow$visitation_rate)) {
      if (thisrow$taxa_recorded=="bees") {
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
      } else if (thisrow$taxa_recorded=="all visitors") {
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
      } else if (thisrow$taxa_recorded=="bees and leps") { 
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
      } else if (thisrow$taxa_recorded=="bees and syrphids") { 
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
      } else if (thisrow$taxa_recorded=="bees, syrphids and wasps") { 
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

# Note that about 100 rows were removed due to NAs in either visitation rate or abundance at this step, which is as expected
# Most of the missing data are NAs scattered across many of the studies.  
# Of special note is "Thijs_Fijen_Allium_porrum_France" which has more NAs than usual in the visit data, but it still looks better than the abundance data which appears to have an issue with duplicated values.

# list which studies had NAs dropped:
unique(subset(field_data2a, site_id %in% setdiff(field_data2a$site_id, field_data2c$site_id))$study_id2)

nrow(field_data2a)	# 3009 as of 1/27/2023
nrow(field_data2c)	# 2854 as of 1/27/2023


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

# make a new column for each insect group that pulls either visits or abundance depending on what is listed in the "use_visits_or_abundance" column
# ***JR code note: again, looping through subsets of the data, processing them, and rebuilding the dataset using rbind

# Create a field level data summary
field_data_summary = unique( field_data[,c("study_id2","taxa_recorded","use_visits_or_abundance")] )
field_data_summary 

# define the list of studies
studies = unique(field_data2c$study_id2)

# create new empty dataframe
field_data2d = data.frame()

for (i in 1:length(studies)) {
  
  fd_sub = subset(field_data2c, study_id2==studies[i])
  
  # visit rate studies
  if (subset(field_data_summary, study_id2==studies[i])$use_visits_or_abundance == "use visits") {
    
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
nrow(field_data2d) # rows should be the same

# data processing: yield data-----

# remove rows with no yield data
field_data3 = subset(field_data2d, !is.na(field_data2d$yield))

# check number of rows
nrow(field_data2d) # 2892 as of 1/27/2023
nrow(field_data3) # 2326 as of 1/27/2023

# check number of studies
length(unique(field_data2d$study_id2))
length(unique(field_data3$study_id2))	# 97 as of 1/27/2023


# data processing: transform to z-scores----

# note on NAs for richness z-scores: scale() removes NAs when calculating the means.  
# I checked that it was not introducing NAs for every row in a crop when there are only a couple NAs--it was not.  
# This was only an issue for Natacha_Chacoff_Citrus_paradisi_Argentina, which had NAs for richness in 4/12 fields.  Most or all of the other NAs in richness were already removed for missing other data. 

# define new list of studies that remain
studies = unique(field_data3$study_id2)

# create empty dataframe	
field_data4 = data.frame()

# loop, process, and rebuild
for (i in 1:length(studies)) {
  
  fd_sub = subset(field_data3, study_id2==studies[i])
  
  fd_sub$yield_z = scale(fd_sub$yield)[,1]
  fd_sub$honeybee_z = scale(fd_sub$honeybee)[,1]
  fd_sub$wild_insects_z = scale(fd_sub$wild_insects)[,1]
  fd_sub$all_insects_z = scale(fd_sub$all_insects)[,1]
  fd_sub$richness_z = scale(fd_sub$observed_pollinator_richness)[,1]
  
  # also calculate SDs for comparing sd of HB to WI
  fd_sub$honeybee_sd = sd(fd_sub$honeybee)
  fd_sub$wild_insects_sd = sd(fd_sub$wild_insects)
  
  field_data4 = rbind(fd_sub, field_data4)
  
}

nrow(field_data3)
nrow(field_data4) # rows should be the same

# remove rows with NA in the yield z-score
field_data5 = subset(field_data4, !is.na(field_data4$yield_z))

nrow(field_data4)
nrow(field_data5) # 2323 as of 1/27/2023

# check number of studies
length(unique(field_data4$study_id2)) # 97 as of 1/27/2023
length(unique(field_data5$study_id2)) # 94 as of 1/27/2023

# These studies were dropped because there was only 1 field each and we can't calculate z scores (creates NAs) for yield (or insects either, actually)
setdiff(unique(field_data4$study_id2), unique(field_data5$study_id2))

# "Timothy_Weekers_Malus_domestica_UK"
# "Breno_M_Freitas_Annona_squamosa_Brazil"
# "Breno_M_Freitas_Annona_muricata_Brazil"

# Replace NAs in HB due to divide by zero in z-score with manually entered 0's.  
# I currently think this is better than just dropping them from model runs.
# this occurred in 2 studies:
# [1] "Finbarr_G_Horgan_Abelmoschus_esculentus_Philippines"		# appears HB recorded, but were 0
# [2] "Breno_M_Freitas_Malpighia_emarginata_Brazil" 			# appears HB recorded, but were 0

field_data5a = data.frame()
for (i in 1:nrow(field_data5)) {
  thisrow = field_data5[i,]
  if (is.na(thisrow$honeybee_z)) { thisrow$honeybee_z <- 0 }
  field_data5a = rbind(thisrow, field_data5a)
}

nrow(field_data5)
nrow(field_data5a) # rows should be the same


# data processing: sample size-----

# drop studies with less than 3 fields (or more as suggested by LG)
studies = names(subset(table(field_data5a$study_id2), table(field_data5a$study_id2)>=3)) 	
studies	# 93 as of 1/27/2023

# names(subset(table(field_data5a$study_id2), table(field_data5a$study_id2)>=8))
# 76 as of 1/27/2023

field_data6a = subset(field_data5a, study_id2 %in% studies)
nrow(field_data5a)
nrow(field_data6a) # 2321 as of 1/27/2023

# which studies were deleted due to sample size?
setdiff(field_data5a$study_id2, field_data6a$study_id2)
# [1] "Bryony_Willcox_Macadamia_integrifolia_Australia"			# only 2 fields


# data processing: data subsets for main analyses----

# subset of studies that we can use for the main analysis (HB vs WI)
field_data8a = subset(field_data6a, !is.na(field_data6a$honeybee_z) & !is.na(field_data6a$wild_insects_z))
length(unique(field_data8a$study_id2)) # 93 as of 1/27/2023

nrow(field_data6a)
nrow(field_data8a) # no rows should be removed at this step because NAs already dealt with

# a subset of studies that provide richness data
field_data9a = subset(field_data8a, !is.na(richness_z))
length(unique(field_data9a$study_id2))	# 63 as of 1/27/2023 (many studies did not provide richness data)

nrow(field_data8a)
nrow(field_data9a) # 1129 as of 1/27/2023 

#Save data to make it accesible to other scripts down the chain.
save(field_data8a, field_data9a, file = "scripts/temp_data.RData")
