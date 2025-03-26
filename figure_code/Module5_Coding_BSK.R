# Placenta database cleaning
# Data retrieved 11-18-24
# By: Bailey Kane

# ----------------
# import libraries and data
# ----------------

# import libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# read in .csv
rawdata <- read.csv("Placenta Dbase11-18-24 - Studies.csv")

# ----------------
# create function to clean TRUE/FALSE/NAs
# ----------------

# logical.cleaner() function will clean all:
#   YES;TRUE into 1
#   NO;NONE;FALSE into 0
#   NA;N/A;-;blank;UNKNOWN;UNK;NULL into NA
logical.cleaner <- function(x) {
  x %>%
  gsub(".*\\b(yes|true|yyes)\\b.*", 1, ., ignore.case = TRUE) %>% 
  gsub(".*\\b(no|not|none|false|nno|noo|nl)\\b.*", 0, ., ignore.case = TRUE) %>% 
  gsub(".*\\b(na|n/a|^-$|^$|unknown|unk|null)\\b.*", NA, ., ignore.case = TRUE)
}

numeric.cleaner <- function(x) { #need this one because the data wasn't read in as numeric in some places (sample size)
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  return (x)
}

boolean.to.nat.lang <- function(x) { # for easy reading of Y/N/NA answers
  x %>%
  gsub(0, "No", .) %>%
  gsub(1, "Yes", .)
}
# ----------------
# clean data
# ----------------

# replace all periods in colnames with underscores, remove trailing underscores
names(rawdata) <- gsub("\\.+","_", names(rawdata))
names(rawdata) <- gsub("_+$","", names(rawdata))

# remove all rows that have not yet been annotated or don't have a data type
data <- rawdata[rawdata$GEO_Series_ID_GSE!="",]
data <- data[data$Data_type_from_CURE_list!="",]

# clean additional data types column
data$Additional_data_types_included_in_the_entry_list_if_any <- logical.cleaner(data$Additional_data_types_included_in_the_entry_list_if_any)
# indicate multiple data types if additional entry in 'additional data types'
data$Data_type_inc_multiple <- ifelse(data$Additional_data_types_included_in_the_entry_list_if_any %in% c(0,NA),
                                      data$Data_type_from_CURE_list,
                                      "Multiple")

# clean 'Organism' column with custom pipeline
data$Organism <- lapply(data$Organism, function(x){
  x %>%
  gsub(".*(;|,).*", "Multiple", .) %>%
  gsub(".*\\b(homo sapiens|human)\\b.*", "Homo sapiens", ., ignore.case = "TRUE") %>%
  gsub(".*\\b(mus musculus|mice|mouse)\\b.*", "Mus musculus", ., ignore.case = "TRUE") %>%
  gsub(".*\\b(Rattus norvegicus|rat|rattus)\\b.*", "Rattus norvegicus", ., ignore.case = "TRUE") %>%
  gsub(".*\\b(Bos taurus|cow)\\b.*", "Bos taurus", ., ignore.case = "TRUE") %>%
  gsub(".*\\b(Equus caballus|horse)\\b.*", "Equus caballus", ., ignore.case = "TRUE") %>%
  gsub(".*\\b(Sus scrofa|boar)\\b.*", "Sus scrofa", ., ignore.case = "TRUE") %>%
  logical.cleaner(.)
  
})
# flatten because this becomes a list for some reason....
data$Organism <- unlist(data$Organism)

# identify common organisms (n=>5)
Common_Organisms <- data %>% count(data$Organism)
Common_Organisms <- Common_Organisms[Common_Organisms$n>4,]
data$Common_Organisms <- ifelse(data$Organism %in% Common_Organisms$'data$Organism', 
                                data$Organism, 
                                "Other")

# clean y/n columns [43:64]
data[43:64] <- lapply(data[43:64], function(x) {
  logical.cleaner(x)
})

# clean numeric columns
data$Sample_size_placenta <- numeric.cleaner(data$Sample_size_placenta)
data$Total_GEO_sample_size <- numeric.cleaner(data$Total_GEO_sample_size)
data$Sample_size_decidua <- numeric.cleaner(data$Sample_size_decidua)

# pull SuperSerieses into dataframe 'SuperSeries'
SuperSeries <- data[data$SuperSeries_check_if_yes==TRUE,]
data <- data[data$SuperSeries_check_if_yes==FALSE,]

# ----------------
# ggplots
# ----------------

# pie chart of data type frequencies
datatable <- as.data.frame(table(data$Data_type_inc_multiple), useNA = "ifany")
ggplot(datatable, aes(x="", y=Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Frequency of Data Type",
       subtitle = "of non-SuperSeries GEO entries",
       fill = "Data Type") +
  geom_text(aes(label = Freq), position = position_stack(vjust=0.5), size = 2)
ggsave("Frequency of Data Type.png",
       width = 6,
       height = 4,
       bg = "white")

# same thing but with SuperSeries
superdatatable <- as.data.frame(table(SuperSeries$Data_type_inc_multiple))
superdatatable$Var1 <- as.character(superdatatable$Var1)
superdatatable[nrow(superdatatable)+1,] = c("Proteomics", as.numeric(0)) # zero proteomics
superdatatable$Var1 <- as.factor(superdatatable$Var1)
superdatatable$Freq <- as.numeric(superdatatable$Freq)

ggplot(superdatatable, aes(x="", y=Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Frequency of Data Type",
       subtitle = "of SuperSeries GEO entries only",
       fill = "Data Type") +
  geom_text(aes(label = Freq), position = position_stack(vjust=0.5), size = 2)
ggsave("Frequency of Data Type Superseries.png",
       width = 6,
       height = 4,
       bg = "white")

# histogram of GEO Sample Size frequency (of those <100), colored by Data Type
ggplot(data[data$Total_GEO_sample_size<100,], (aes(x=Total_GEO_sample_size, fill=Data_type_inc_multiple))) +
  geom_histogram(bins=25) +
  theme_minimal() +
  labs(title = "Sample Size Frequency by Data Type",
       subtitle = "of non-SuperSeries GEO entries with n < 100",
       x = "GEO Sample Size (n)",
       y = "Frequency",
       fill = "Data Type")
ggsave("Sample Size Frequency by Data Type.png",
       width = 6,
       height = 4,
       bg = "white")

# pie chart of organism frequencies
organismtable <- as.data.frame(table(data$Common_Organisms))
ggplot(organismtable, aes(x="", y=Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Organism Frequency",
       subtitle = "of non-SuperSeries GEO entries",
       fill = "Organism") +
  geom_text(aes(label = Freq), position = position_stack(vjust=0.5), size = 2) +
  scale_fill_brewer(palette = "RdYlBu")
ggsave("Organism Frequency.png",
       width = 6,
       height = 4,
       bg = "white")

# bar plot of data types by organism name
ggplot(data[is.na(data$Data_type_inc_multiple)==0,], aes(x=Data_type_inc_multiple, fill = Common_Organisms)) +
  geom_bar() + 
  theme_minimal() +
  labs(title = "Organism Prevalence by Data Type",
       subtitle = "of non-SuperSeries GEO entries",
       x = "Data Type",
       y = "Frequency",
       fill = "Organism") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.key.size = unit(0.5, 'cm')) +
  scale_fill_brewer(palette = "RdYlBu")
ggsave("Organism Prevalence by Data Type.png",
       width = 8,
       height = 4,
       bg = "white")

# of human data with decidual samples, how many studies report maternal features?
# subset data into organism = homo sapiens, sample size (decidua) > 0
humandecidualdata <- data[data$Sample_size_decidua>0, ] %>% .[.$Organism=="Homo sapiens",]

# maternal feature 1: Maternal Age
maternal_feature_1 <- as.data.frame(table(humandecidualdata$Maternal_age_at_sample_collection_provided_yes_no, useNA = "ifany"))
maternal_feature_1$Var1 <- boolean.to.nat.lang(maternal_feature_1$Var1)
ggplot(maternal_feature_1, aes(x = "", y = Freq, fill = Var1)) +
      geom_bar(width = 1, stat = "identity") +
      coord_polar(theta = "y") +
      theme_void() +
      labs(title = "Maternal Age Reported",
           subtitle = "In Homo sapiens datasets containing decidual samples",
           fill = "") +
      geom_text(aes(label = Freq), position = position_stack(vjust=0.5))

# maternal feature 2: Fetal complications reported
maternal_feature_2 <- as.data.frame(table(humandecidualdata$Fetal_complications_listed_yes_no, useNA = "ifany"))
maternal_feature_2$Var1 <- boolean.to.nat.lang(maternal_feature_2$Var1)
ggplot(maternal_feature_2, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Fetal Complications Reported",
       subtitle = "In Homo sapiens datasets containing decidual samples",
       fill = "") +
  geom_text(aes(label = Freq), position = position_stack(vjust=0.5))

# maternal feature 3: Pregnancy complications reported
maternal_feature_3 <- as.data.frame(table(humandecidualdata$Samples_from_pregnancy_complications_collected, useNA = "ifany"))
maternal_feature_3$Var1 <- boolean.to.nat.lang(maternal_feature_3$Var1)
ggplot(maternal_feature_3, aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Pregancy Complications Reported",
       subtitle = "In Homo sapiens datasets containing decidual samples",
       fill = "") +
  geom_text(aes(label = Freq), position = position_stack(vjust=0.5))

# a greater proportion of superseries consist of multiple data types than non-superseries
ggplot(as.data.frame(table(SuperSeries$Data_type_inc_multiple), useNA = "ifany"), aes(x="", y=Freq, fill = Var1)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  theme_void() +
  labs(title = "Frequency of Data Type",
       subtitle = "of SuperSeries GEO entries only",
       fill = "Data Type") +
  geom_text(aes(label = Freq), position = position_stack(vjust=0.5), size = 2)
ggsave("Frequency of Data Type (SuperSeries).png",
       width = 6,
       height = 4,
       bg = "white")

# what features are reported for all datasets?
# subset dataframe into features only
datafeatures <- data[c(44:60,62)]
names(datafeatures) <- c("Birthweight",
                         "Gestational Age at Delivery",
                         "Gestational Age at Sample Collection",
                         "Sex of Offspring",
                         "Parity",
                         "Gravidity",
                         "Number of Offspring per Pregnancy",
                         "Self-Reported Race/Ethnicity of Mother",
                         "Genetic Ancestry",
                         "Maternal Height",
                         "Maternal Pre-Pregnancy Weight",
                         "Paternal Height",
                         "Paternal Weight",
                         "Maternal Age at Sample Collection",
                         "Paternal Age",
                         "Samples from Pregnancy Complications",
                         "Mode of Delivery",
                         "Samples from Fetal Complications")
# clean age columns
datafeatures$`Gestational Age at Delivery` <- datafeatures$`Gestational Age at Delivery` %>%
  gsub(".*day.?.*", 1, ., ignore.case = TRUE) %>%
  gsub(".*week.?.*", 1, ., ignore.case = TRUE)
datafeatures$`Gestational Age at Sample Collection` <- datafeatures$`Gestational Age at Sample Collection` %>%
  gsub(".*(day.?|E).*", 1, ., ignore.case = TRUE) %>%
  gsub(".*week.?.*", 1, ., ignore.case = TRUE)
datafeatures$`Maternal Age at Sample Collection` <- datafeatures$`Maternal Age at Sample Collection` %>%
  gsub(".*(day.?|E).*", 1, ., ignore.case = TRUE) %>%
  gsub(".*week.?.*", 1, ., ignore.case = TRUE)

# clean columns where people accidentally described data instead of reporting features
datafeatures$`Samples from Fetal Complications` <- datafeatures$`Samples from Fetal Complications` %>%
  gsub(".*[A-Za-z].*", 1, ., ignore.case = TRUE)
datafeatures$`Mode of Delivery` <- datafeatures$`Mode of Delivery` %>%
  gsub(".*[A-Za-z].*", 1, ., ignore.case = TRUE)
datafeatures$`Samples from Pregnancy Complications` <- datafeatures$`Samples from Pregnancy Complications` %>%
  gsub(".*[A-Za-z].*", 1, ., ignore.case = TRUE)
datafeatures$`Maternal Height` <- datafeatures$`Maternal Height` %>%
  gsub(".*[A-Za-z].*", 0, ., ignore.case = TRUE)
datafeatures$`Maternal Pre-Pregnancy Weight` <- datafeatures$`Maternal Pre-Pregnancy Weight` %>%
  gsub(".*[A-Za-z].*", 0, ., ignore.case = TRUE)
datafeatures$`Genetic Ancestry` <- datafeatures$`Genetic Ancestry` %>%
  gsub(".*[A-Za-z].*", 1, ., ignore.case = TRUE)

# make long-format data so I can represent it with multiple bars
response_counts <- datafeatures %>%
  pivot_longer(cols = everything(), names_to = "Feature", values_to = "Response") %>%
  count(Feature, Response)
response_counts$Response <- lapply(response_counts$Response, function(x){
  boolean.to.nat.lang(x)
})
response_counts$Response <- unlist(response_counts$Response)

# create plot (finally)
ggplot(response_counts, aes(x=Feature, y = n, fill = Response)) + 
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Metadata Provided in Dataset?",
       subtitle = "in non-SuperSeries GEO entries",
       x = "",
       y = "Frequency",
       fill = "Data Provided") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.key.size = unit(0.5, 'cm'))
ggsave("Metadata Provided.png",
       width = 8,
       height = 5,
       bg = "white")
