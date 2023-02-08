#This script keeps track of previous estimates in a .csv for future plotting

#build original csv----
# (run ONLY once, hence commented)

# estimates_data <- data.frame(
	# timepoint =            c("2013","2013","2013","2013","2015","2016","2016","2019"),
	# estimate_HB =          c(0.149, 0.151, NA,    NA,    NA,    -0.019,NA,    NA),
	# estimate_WI =          c(0.288, NA,    0.323, NA,    NA,    0.187, NA,    NA),
	# estimate_all_insects = c(NA,    NA,    NA,    NA,    0.3,   NA,    NA,    0.08),
	# estimate_richness =    c(NA,    NA,    NA,    0.17,  0.07,  NA,    NA,    0.106),
	# n =                    c(20,    22,    30,    18,    33,    19,    19,    42))
# estimates_data
# write.table(estimates_data, "scripts/estimates_data.csv",
#          row.names = FALSE, sep = ",")

#Add new estimates----

old_table <- read.csv("scripts/estimates_data.csv")
source("scripts/read_data.R") #Once ran here, data will be made available to other scripts.
library(lme4)

# HB and WI
mod5 = lmer(data=field_data8a, yield_z ~ wild_insects_z + honeybee_z + (1 + wild_insects_z|study_id2) + (1 + honeybee_z|study_id2), na.action=na.fail, REML=T)
#summary(mod5)

# all insects and richness
modR_8b = lmer(data=field_data9a, yield_z ~ all_insects_z + richness_z + (0 + all_insects_z|study_id2) + (0 + richness_z|study_id2), na.action=na.fail, REML=T)
#summary(modR_8b)


estimates_data1 <- data.frame(timepoint = c(Sys.Date()),
                             estimate_HB = c(fixef(mod5)["honeybee_z"]),
                             estimate_WI = c(fixef(mod5)["wild_insects_z"]),
                             n = c(length(unique(field_data8a$study_id2))))

if(round(old_table$estimate_WI[nrow(old_table)-1],4) != round(estimates_data1$estimate_WI[1],4) | 
   round(old_table$estimate_HB[nrow(old_table)-1],4) != round(estimates_data1$estimate_HB[1],4)){
  write.table(estimates_data1, "scripts/estimates_data.csv", 
              append = TRUE, sep = ",", row.names = FALSE,
              col.names = FALSE)
}  

estimates_data2 <- data.frame(timepoint = c(Sys.Date()),
                             estimate_all_insects = c(fixef(modR_8b)["all_insects_z"]),
                             estimate_richness = c(fixef(modR_8b)["richness_z"]),
                             n = c(length(unique(field_data9a$study_id2))))
							 
if(round(old_table$estimate_all_insects[nrow(old_table)],4) != round(estimates_data2$estimate_all_insects[1],4) | 
   round(old_table$estimate_richness[nrow(old_table)],4) != round(estimates_data2$estimate_richness[1],4)){
  write.table(estimates_data2, "scripts/estimates_data.csv", 
              append = TRUE, sep = ",", row.names = FALSE,
              col.names = FALSE)
}  


