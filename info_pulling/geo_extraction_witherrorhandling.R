#load packages
library(GEOquery)
library(readxl)
library(writexl)
library(dplyr)
#------------------------------------------------------------------------

#load .csv into a dataframe
gse_df <- read.csv("gse.csv", header = FALSE, stringsAsFactors = FALSE)
geo_ids <- as.character(gse_df[[1]])
print(geo_ids)

#limit processing to a subset for testing
geo_ids <- geo_ids[1:10]
geo_ids


#------------------------------------------------------------------------

#function to convert numeric GEO ID to standard GSE format
convert_to_gse <- function(numeric_ids) {
  sapply(numeric_ids, function(numeric_id) {
    numeric_id <- as.character(numeric_id)
    
    #if ID starts with "2000", remove "2000" prefix
    if (grepl("^2000", numeric_id)) {
      return(paste0("GSE", sub("^2000", "", numeric_id)))
    }
    
    #if ID starts with "200", remove "200" prefix
    if (grepl("^200", numeric_id)) {
      return(paste0("GSE", sub("^200", "", numeric_id)))
    }
    
    #otherwise, return as "GSE" + numeric_id
    return(paste0("GSE", numeric_id))
  }, USE.NAMES = FALSE)  
}

geo_ids <- convert_to_gse(geo_ids)
geo_ids
#------------------------------------------------------------------------


#function to extract GEO information with error handling
get_gse_metadata <- function(geo_ids) {
  failed_ids <- c()
  gse_dataframes <- list()
  
  for (geo_id in geo_ids) {
    result <- tryCatch({
      gse <- getGEO(geo_id, GSEMatrix = FALSE)
      meta_data <- Meta(gse)
      
      data.frame(
        GEO_ID = geo_id,
        contact_country = ifelse(!is.null(meta_data$contact_country), meta_data$contact_country, NA),
        contact_email = ifelse(!is.null(meta_data$contact_email), meta_data$contact_email, NA),
        contact_institute = ifelse(!is.null(meta_data$contact_institute), meta_data$contact_institute, NA),
        contact_name = ifelse(!is.null(meta_data$contact_name), meta_data$contact_name, NA),
        geo_accession = ifelse(!is.null(meta_data$geo_accession), meta_data$geo_accession, NA),
        overall_design = ifelse(!is.null(meta_data$overall_design), meta_data$overall_design, NA),
        platform_id = ifelse(!is.null(meta_data$platform_id), meta_data$platform_id, NA),
        pubmed_id = ifelse(!is.null(meta_data$pubmed_id), meta_data$pubmed_id, NA),
        relation = ifelse(!is.null(meta_data$relation), paste(meta_data$relation, collapse = ", "), NA),
        sample_number = ifelse(!is.null(meta_data$sample_id), paste(meta_data$sample_id, collapse = ", "), NA),
        sample_taxid = ifelse(!is.null(meta_data$sample_taxid), meta_data$sample_taxid, NA),
        submission_date = ifelse(!is.null(meta_data$submission_date), meta_data$submission_date, NA),
        last_update_date = ifelse(!is.null(meta_data$last_update_date), meta_data$last_update_date, NA),
        main_topic = ifelse(!is.null(meta_data$summary), meta_data$summary, NA),
        supplementary_file = ifelse(!is.null(meta_data$supplementary_file), paste(meta_data$supplementary_file, collapse = ", "), NA),
        title = ifelse(!is.null(meta_data$title), meta_data$title, NA),
        type = ifelse(!is.null(meta_data$type), meta_data$type, NA),
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      failed_ids <<- c(failed_ids, geo_id)
      return(NULL)
    })
    
    if (!is.null(result)) {
      gse_dataframes <- append(gse_dataframes, list(result))
    }
  }
  
  #combine all individual GSE data frames into one
  combined_df <- do.call(rbind, gse_dataframes)
  return(list(data = combined_df, failed = failed_ids))
}
#------------------------------------------------------------------------


#run it
metadata_result <- get_gse_metadata(geo_ids)
metadata_df <- metadata_result$data
failed_ids <- metadata_result$failed

#replace ",," and "," with a space in the contact_name column
metadata_df$contact_name <- gsub(",+", " ", metadata_df$contact_name)

#extract SubSeries, BioProject, and SRA IDs
metadata_df$SubSeries <- ifelse(grepl("SubSeries of: GSE[0-9]+", metadata_df$relation),
                                sub(".*SubSeries of: (GSE[0-9]+).*", "\\1", metadata_df$relation), 
                                NA)
metadata_df$BioProject <- ifelse(grepl("BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA[0-9]+", metadata_df$relation),
                                 sub(".*BioProject: https://www.ncbi.nlm.nih.gov/bioproject/(PRJNA[0-9]+).*", "\\1", metadata_df$relation), 
                                 NA)
metadata_df$SRA <- ifelse(grepl("SRA: https://www.ncbi.nlm.nih.gov/sra\\?term=SRP[0-9]+", metadata_df$relation),
                          sub(".*SRA: https://www.ncbi.nlm.nih.gov/sra\\?term=(SRP[0-9]+).*", "\\1", metadata_df$relation), 
                          NA)

metadata_df$sample_number <- sapply(strsplit(metadata_df$sample_number, ", "), length)
#------------------------------------------------------------------------

metadata_df$sample_taxid

taxonomy <- c(
  "10090" = "Mus musculus",
  "9315" = "Notamacropus eugenii",
  "9606" = "Homo sapiens",
  "9823" = "Sus scrofa",
  "10116" = "Rattus norvegicus",
  "9913" = "Bos taurus",
  "9361" = "Dasypus novemcinctus",  
  "10092" = "Mus musculus domesticus",
  "9544" = "Macaca mulatta",  
  "9915" = "Bos indicus",
  "9545"="Macaca nemestrina",
  "9986"="Oryctolagus cuniculus"
)
library(dplyr)
metadata_df$organism <- taxonomy[as.character(metadata_df$sample_taxid)]
metadata_df$sample_taxid
metadata_df$organism
#------------------------------------------------------------------------


#convert dates
metadata_df$last_update_date <- format(as.Date(metadata_df$last_update_date, format = "%b %d %Y"), "%m/%d/%Y")
metadata_df$submission_date <- format(as.Date(metadata_df$submission_date, format = "%b %d %Y"), "%m/%d/%Y")

#save output
write.csv(metadata_df, "gse_metadata.csv", row.names = FALSE)
write.csv(data.frame(Failed_GSE_IDs = failed_ids), "failed_gse_ids.csv", row.names = FALSE)

print("Processing complete. Failed GSE IDs saved to failed_gse_ids.csv")
