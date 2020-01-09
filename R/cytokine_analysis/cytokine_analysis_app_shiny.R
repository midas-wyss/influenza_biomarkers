#########------------ Cytokine analysis app ------------######### 

library(tidyverse) # for organizing data structures
library(readr) # for data parsing / loading in data files
library(readxl) # for reading Excel files
library(hash)
library(shiny)
library(ggpubr)# runExample("01_hello")
options(warn = -1) 
# devtools::install_github("kassambara/ggpubr")
# install.packages("shiny")
cytokine_raw_tibb <- readxl::read_excel('../../../data/cytokine_data/Human-10plex.xlsx',
                                      sheet="Conc in Range",
                                      skip=7
                                      )

list_of_well_types <- cytokine_raw_tibb$'...1'
start_tibb_index <- which(list_of_well_types %in% 'Type')[2]

cytokine_concs_tibb <- cytokine_raw_tibb[start_tibb_index+4:(dim(cytokine_raw_tibb)[1]-30), ]
cytokine_concs_tibb <- cytokine_concs_tibb %>% rename(Type = ...1,
                                                  Well = ...2,
                                                  Description = ...3)

well_names <- cytokine_concs_tibb$Well 

cytokine_concs_tibb <- cytokine_concs_tibb[, !(names(cytokine_concs_tibb) %in% "Well")]
cytokine_concs_df <- data.frame(cytokine_concs_tibb)
rownames(cytokine_concs_df) <- well_names
plate_map_conditions <- c("NHB control N3",
                          "NHB 0.1 MOI N3",
                          "NHB 0.05 MOI N2",
                          "NHB 0.05 MOI N3",
                          "NHB 0.01 MOI N2",
                          "NHB 0.01 MOI N3",
                          "COPD control N3",
                          "COPD 0.1 MOI N3",
                          "COPD 0.05 MOI N2",
                          "COPD 0.05 MOI N3",
                          "COPD 0.01 MOI N2",
                          "COPD 0.01 MOI N3") 

plate_map_wells <- list(c('A4', 'A5', 'A6'),
                     c('B4', 'B5', 'B6'),
                     c('C4', 'C5', 'C6'),
                     c('D4', 'D5', 'D6'),
                     c('E4', 'E5', 'E6'),
                     c('F4', 'F5', 'F6'),
                     c('A7', 'A8', 'A9'),
                     c('B7', 'B8', 'B9'),
                     c('C7', 'C8', 'C9'),
                     c('D7', 'D8', 'D9'),
                     c('E7', 'E8', 'E9'),
                     c('F7', 'F8', 'F9'))
plate_map_dict <- hash(keys=plate_map_conditions,
                       values=plate_map_wells)


for(condition in ls(plate_map_dict)){
  condition_all_info <- unlist(strsplit(condition, " "))
  well_1 <- plate_map_dict[[condition]][1]
  well_2 <- plate_map_dict[[condition]][2]
  well_3 <- plate_map_dict[[condition]][3]
  
  if(!grepl("control", condition))
    {
    cytokine_concs_df[well_1, "Condition"] <- condition
    cytokine_concs_df[well_1, "MOI"] <- condition_all_info[2]
    cytokine_concs_df[well_1, "transwell"] <- tail(condition_all_info, n=1)
    cytokine_concs_df[well_1, "Condition"] <- condition_all_info[1]
    cytokine_concs_df[well_1, "Condition_MOI"] <- paste(condition_all_info[1],condition_all_info[2], sep="_")
    
    cytokine_concs_df[well_2, "Condition"] <- condition
    cytokine_concs_df[well_2, "MOI"] <- condition_all_info[2]
    cytokine_concs_df[well_2, "transwell"] <- tail(condition_all_info, n=1)
    cytokine_concs_df[well_2, "Condition"] <- condition_all_info[1]
    cytokine_concs_df[well_2, "Condition_MOI"] <- paste(condition_all_info[1],condition_all_info[2], sep="_")
    
    cytokine_concs_df[well_3, "Condition"] <- condition
    cytokine_concs_df[well_3, "MOI"] <- condition_all_info[2]
    cytokine_concs_df[well_3, "transwell"] <- tail(condition_all_info, n=1)
    cytokine_concs_df[well_3, "Condition"] <- condition_all_info[1]
    cytokine_concs_df[well_3, "Condition_MOI"] <- paste(condition_all_info[1],condition_all_info[2], sep="_")
  }
  else{
    cytokine_concs_df[well_1, "Condition"] <- condition
    cytokine_concs_df[well_1, "MOI"] <- "control"
    cytokine_concs_df[well_1, "transwell"] <- tail(condition_all_info, n=1)
    cytokine_concs_df[well_1, "Condition"] <- condition_all_info[1]
    cytokine_concs_df[well_1, "Condition_MOI"] <- paste(condition_all_info[1],condition_all_info[2], sep="_")
    
    cytokine_concs_df[well_2, "Condition"] <- condition
    cytokine_concs_df[well_2, "MOI"] <- "control"
    cytokine_concs_df[well_2, "transwell"] <- tail(condition_all_info, n=1)
    cytokine_concs_df[well_2, "Condition"] <- condition_all_info[1]
    cytokine_concs_df[well_2, "Condition_MOI"] <- paste(condition_all_info[1],condition_all_info[2], sep="_")
    
    cytokine_concs_df[well_3, "Condition"] <- condition
    cytokine_concs_df[well_3, "MOI"] <- "control"
    cytokine_concs_df[well_3, "transwell"] <- tail(condition_all_info, n=1)
    cytokine_concs_df[well_3, "Condition"] <- condition_all_info[1]
    cytokine_concs_df[well_3, "Condition_MOI"] <- paste(condition_all_info[1],condition_all_info[2], sep="_")
  }
  
}

rename_cytokines <- function(x){
  new_cytokine_name <- unlist(strsplit(x, "\\.\\."))[1]
  return(new_cytokine_name)
}
cytokine_concs_df[is.na(cytokine_concs_df)] <- "standard"
cytokine_concs_df[ cytokine_concs_df == "OOR >" ] <- NA
cytokine_concs_df[ cytokine_concs_df == "OOR <" ] <- NA
cytokine_names_raw <- colnames(cytokine_concs_df)[3:(length(colnames(cytokine_concs_df))-4)]
cytokines <- unlist(lapply(cytokine_names_raw, rename_cytokines))
cytokines <- sapply(list(cytokines),function(x) gsub("\\.", "-", as.character(x)))

cytokines_rename_cols <- unlist(lapply(cytokine_names_raw, rename_cytokines))
cytokines_rename_cols <- sapply(list(cytokines_rename_cols),function(x) gsub("\\.", "_", as.character(x)))
colnames(cytokine_concs_df)[3:(length(colnames(cytokine_concs_df))-4)] <- cytokines_rename_cols
plate_map_list_final <- substr(plate_map_conditions,1,nchar(plate_map_conditions)-3)
potential_conds_to_compare_pm <- sort(unique(plate_map_list_final), decreasing = TRUE)
potential_conds_to_compare <- append("None",potential_conds_to_compare_pm)

cytokine_column_conversion <- hash(keys = cytokines, values = cytokines_rename_cols)

######## Transwell specified here ######## 
transwell = "N3"
cytokine = "IL-6"
cytokine_concs_remove_standards_df<-cytokine_concs_df[!(cytokine_concs_df["Condition"]=="standard"),]
cytokine_concs_specify_transwell_df<-cytokine_concs_remove_standards_df[(cytokine_concs_remove_standards_df["transwell"]==transwell),]
num_columns <- dim(cytokine_concs_specify_transwell_df)[2]
cytokine_concs_specify_transwell_df[3:(num_columns-4)] <- lapply(cytokine_concs_specify_transwell_df[3:(num_columns-4)], as.numeric)
cytokine_concs_specify_transwell_df$MOI <- factor(cytokine_concs_specify_transwell_df$MOI,
                       levels = c('control','0.01', '0.05', '0.1'),ordered = TRUE)

# geom_dotplot(binaxis='y', stackdir='center',
# position=position_dodge(1), dotsize = 0.5)
cytokine_concs_specify_transwell_df$Condition <- factor(cytokine_concs_specify_transwell_df$Condition,
                                                  levels = c('NHB','COPD'),ordered = TRUE)
my_comp <- list(c("COPD_control", "COPD_0.01"))
cytokine_conc_vals <- cytokine_concs_specify_transwell_df[cytokine_column_conversion[[cytokine]]]
cytokine_conc_vals <- cytokine_conc_vals[!is.na(cytokine_conc_vals)]
min_val <- min(cytokine_conc_vals) #+ 0.1*min(cytokine_conc_vals)
max_val <- max(cytokine_conc_vals) #+ 0.1*max(cytokine_conc_vals) 
ggplot(cytokine_concs_specify_transwell_df, aes_string(x = "MOI", y=cytokine_column_conversion[[cytokine]],fill="Condition")) + geom_boxplot()+
  stat_summary(fun.y=mean, geom="point",aes_string(group="Condition"), position=position_dodge(0.75), 
               color="black", size=4)+
  theme_classic() + ylim(min_val- 0.1*min_val, max_val + 0.1*max_val)+
  labs(fill="", x="MOI", y=paste(cytokine," (pg/uL)")) + scale_fill_manual(values=c("#B25116", "#FB84D1"))+
  stat_compare_means(aes(group = Condition),
    label = "p.format",
    method = "t.test",
    paired=FALSE, method.args = list(var.equal = FALSE))  +geom_point(shape= 21,position=position_dodge(0.75)) +  theme(legend.text=element_text(size=20),legend.key.width=unit(2,"cm"), legend.key.height=unit(2,"cm"))



