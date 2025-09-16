#################
# Load packages #
#################
library(GEOquery)
library(readxl)
library(writexl)
library(dplyr)
library(stringr)
library(purrr)
library(lubridate)
library(progress)  

#####################
# Utility functions #
#####################

#Null-coalescing operator
`%||%` = function(x, y) if (!is.null(x)) x else y

# Convert numeric IDs to GSE format
convert_to_gse = function(numeric_ids) {
  sapply(numeric_ids, function(numeric_id) {
    numeric_id = as.character(numeric_id)
    if (grepl("^2000", numeric_id)) {
      return(paste0("GSE", sub("^2000", "", numeric_id)))
    }
    if (grepl("^200", numeric_id)) {
      return(paste0("GSE", sub("^200", "", numeric_id)))
    }
    return(paste0("GSE", numeric_id))
  }, USE.NAMES = FALSE)  
}

#Extract GEO metadata with error handling
get_single_gse = function(geo_id) {
  tryCatch({
    gse = getGEO(geo_id, GSEMatrix = FALSE)
    meta_data = Meta(gse)
    
    tibble(
      GEO_ID = geo_id,
      contact_country = meta_data$contact_country %||% NA,
      contact_email = meta_data$contact_email %||% NA,
      contact_institute = meta_data$contact_institute %||% NA,
      contact_name = gsub(",+", " ", meta_data$contact_name %||% NA),
      geo_accession = meta_data$geo_accession %||% NA,
      overall_design = meta_data$overall_design %||% NA,
      platform_id = meta_data$platform_id %||% NA,
      pubmed_id = meta_data$pubmed_id %||% NA,
      relation = if (!is.null(meta_data$relation)) paste(meta_data$relation, collapse = ", ") else NA,
      sample_number = if (!is.null(meta_data$sample_id)) length(meta_data$sample_id) else NA,
      sample_taxid = meta_data$sample_taxid %||% NA,
      submission_date = meta_data$submission_date %||% NA,
      last_update_date = meta_data$last_update_date %||% NA,
      main_topic = meta_data$summary %||% NA,
      supplementary_file = if (!is.null(meta_data$supplementary_file)) paste(meta_data$supplementary_file, collapse = ", ") else NA,
      title = meta_data$title %||% NA,
      type = meta_data$type %||% NA
    )
  }, error = function(e) {
    tibble(GEO_ID = geo_id, error = as.character(e))
  })
}

#################
# Main workflow #
#################

#Load GEO IDs from CSV
gse_df = read.csv("gse.csv", header = FALSE, stringsAsFactors = FALSE)
geo_ids = as.character(gse_df[[1]])

#Limit for testing (change as needed)
geo_ids = geo_ids[1:100]

#Convert to GSE format
geo_ids = convert_to_gse(geo_ids)

#################################################
# Fetch metadata SEQUENTIALLY with progress bar #
#################################################
pb = progress_bar$new(
  format = "  Fetching [:bar] :percent ETA: :eta",
  total = length(geo_ids), clear = FALSE, width=60
)

results = vector("list", length(geo_ids))

for (i in seq_along(geo_ids)) {
  results[[i]] = get_single_gse(geo_ids[i])
  pb$tick()
}

#Separate successes and failures
metadata_df = bind_rows(keep(results, ~ !"error" %in% names(.x)))
failed_ids = map_chr(keep(results, ~ "error" %in% names(.x)), "GEO_ID")

###################
# Post-processing #
###################

#Extract SubSeries, BioProject, and SRA IDs with regex
metadata_df = metadata_df %>%
  mutate(
    SubSeries   = str_extract(relation, "GSE\\d+"),
    BioProject  = str_extract(relation, "PRJNA\\d+"),
    SRA         = str_extract(relation, "SRP\\d+")
  )

#Map taxonomy IDs
taxonomy = c(
  "10090" = "Mus musculus",
  "9315"  = "Notamacropus eugenii",
  "9606"  = "Homo sapiens",
  "9823"  = "Sus scrofa",
  "10116" = "Rattus norvegicus",
  "9913"  = "Bos taurus",
  "9361"  = "Dasypus novemcinctus",  
  "10092" = "Mus musculus domesticus",
  "9544"  = "Macaca mulatta",  
  "9915"  = "Bos indicus",
  "9545"  = "Macaca nemestrina",
  "9986"  = "Oryctolagus cuniculus"
)

metadata_df = metadata_df %>%
  mutate(organism = taxonomy[as.character(sample_taxid)])

#######################
# Fix dates for Excel #
#######################

metadata_df = metadata_df %>%
  mutate(
    # Convert to Date class
    last_update_date = as.Date(parse_date_time(last_update_date, orders = c("b d Y"))),
    submission_date  = as.Date(parse_date_time(submission_date, orders = c("b d Y"))),
    # Format as mm/dd/yyyy for Excel 
    last_update_date = format(last_update_date, "%m/%d/%Y"),
    submission_date  = format(submission_date, "%m/%d/%Y")
  )

################
#Save results #
################

write_xlsx(list(
  Metadata = metadata_df,
  Failed   = tibble(Failed_GSE_IDs = failed_ids)
), "gse_metadata.xlsx")

print("Processing complete. Results saved to gse_metadata.xlsx")


