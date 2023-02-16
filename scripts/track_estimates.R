#This script keeps track of previous estimates in a .csv for future plotting

#build original csv----
# (run ONLY once, hence commented)

# estimates_data <- data.frame(
	# timepoint =            c("2013-03-29", "2015-08-02","2016-01-05","2019-10-16", "2022-12-10"),
	# estimate_HB =          c(0.128,  NA,    -0.019,NA, 0.085),
	# estimate_WI =          c(0.302,  NA,    0.187, NA, 0.078),
	# estimate_richness =    c(0.201,  0.07,  NA,    0.106, 0.036),
	# estimate_HBxWI =       c(-0.030, NA,    NA,    NA, -0.015),	
	# estimate_HBxrichness = c(NA,     NA,    NA,    NA, 0.006),	
	# estimate_WIxrichness = c(NA,     NA,    NA,    NA, -0.033),		
	# n =                    c(17,     33,    19,    42, 93))
# estimates_data
# write.table(estimates_data, "scripts/estimates_data.csv", row.names = FALSE, sep = ",")

#Add new estimates----

old_table <- read.csv("scripts/estimates_data.csv")
source("scripts/read_data.R") #Once ran here, data will be made available to other scripts.
library(lme4)

# Full model (including all 2-way interactions)
mod_R18 = lmer(data=field_data9a, yield_z ~ wild_insects_z + honeybee_z + richness_z + wild_insects_z:honeybee_z + wild_insects_z:richness_z + honeybee_z:richness_z 
+ (1 + wild_insects_z|study_id2)  + (1 + honeybee_z|study_id2) + (1 + richness_z|study_id2) + (1 + wild_insects_z:honeybee_z|study_id2) + (1 + wild_insects_z:richness_z|study_id2) + (1 + honeybee_z:richness_z|study_id2), na.action=na.fail, REML=T)
summary(mod_R18)

#summary(mod_R18)


estimates_data <- data.frame(timepoint = c(Sys.Date()),
                             estimate_HB = c(fixef(mod_R18)["honeybee_z"]),
                             estimate_WI = c(fixef(mod_R18)["wild_insects_z"]),
							 estimate_richness = c(fixef(mod_R18)["richness_z"]),
							 estimate_HBxWI = c(fixef(mod_R18)["wild_insects_z:honeybee_z"]),
							 estimate_HBxrichness = c(fixef(mod_R18)["honeybee_z:richness_z"]),
							 estimate_WIxrichness = c(fixef(mod_R18)["wild_insects_z:richness_z"]),
                             n = c(length(unique(field_data8a$study_id2))))

if (
	is.na(old_table$estimate_WI[nrow(old_table)]) | 
	is.na(old_table$estimate_HB[nrow(old_table)]) | 
	is.na(old_table$estimate_richness[nrow(old_table)]) | 
	is.na(old_table$estimate_HBxWI[nrow(old_table)]) | 
	is.na(old_table$estimate_HBxrichness[nrow(old_table)]) | 
	is.na(old_table$estimate_WIxrichness[nrow(old_table)])
) {
  write.table(estimates_data, "scripts/estimates_data.csv", 
              append = TRUE, sep = ",", row.names = FALSE,
              col.names = FALSE)
} else if (
	round(old_table$estimate_WI[nrow(old_table)],4) != round(estimates_data$estimate_WI[1],4) | 
	round(old_table$estimate_HB[nrow(old_table)],4) != round(estimates_data$estimate_HB[1],4) |
	round(old_table$estimate_richness[nrow(old_table)],4) != round(estimates_data$estimate_richness[1],4) |
	round(old_table$estimate_HBxWI[nrow(old_table)],4) != round(estimates_data$estimate_HBxWI[1],4) |
	round(old_table$estimate_HBxrichness[nrow(old_table)],4) != round(estimates_data$estimate_HBxrichness[1],4) |
	round(old_table$estimate_WIxrichness[nrow(old_table)],4) != round(estimates_data$estimate_WIxrichness[1],4) 
) {
  write.table(estimates_data, "scripts/estimates_data.csv", 
              append = TRUE, sep = ",", row.names = FALSE,
              col.names = FALSE)
}  

