#This script keeps track of previous estimates in a .csv for future plotting

#build original csv----
# (run ONLY once, hence commented)

# estimates_data <- data.frame(timepoint = c("2013", "2016", "2019"),
#                              estimate_HB = c(0.1, 0.1, 0.1),
#                              estimate_WB = c(0.1, 0.1, 0.1),
#                              n = c(20,30,40))
# write.table(estimates_data, "scripts/estimates_data.csv",
#          row.names = FALSE, sep = ",")

#IMPORTANT: Dates and estimates are placeholders now. FETCH THE RIGHT DATA

#Add new estimates----

old_table <- read.csv("scripts/estimates_data.csv")
source("scripts/read_data.R") #If ran here, I think the best approach 
#is saving the model outputs somewhere so results qmd can use it.

library(lme4)
mod9 = lmer(data=field_data8a, 
            yield_z ~ wild_insects_z*honeybee_z + 
              (1 + wild_insects_z|study_id2) + 
              (1 + honeybee_z|study_id2) + 
              (1 + wild_insects_z*honeybee_z|study_id2), 
            na.action=na.fail, REML=T)
#summary(mod9)

estimates_data <- data.frame(timepoint = c(Sys.Date()),
                             estimate_HB = c(fixef(mod9)["honeybee_z"]),
                             estimate_WB = c(fixef(mod9)["wild_insects_z"]),
                             n = c(length(unique(field_data8a$study_id2))))

if(round(old_table$estimate_WB[nrow(old_table)],4) != round(estimates_data$estimate_WB[1],4) | 
   round(old_table$estimate_HB[nrow(old_table)],4) != round(estimates_data$estimate_HB[1],4)){
  write.table(estimates_data, "scripts/estimates_data.csv", 
              append = TRUE, sep = ",", row.names = FALSE,
              col.names = FALSE)
}  


