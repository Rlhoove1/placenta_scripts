setwd("C:/Users/White/OneDrive/Desktop")

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

# Extract GSE-level metadata, plus GSM-level library/seq info
get_single_gse_full = function(geo_id) {
  tryCatch({
    gse = getGEO(geo_id, GSEMatrix = FALSE)
    meta_data = Meta(gse)
    
    # Extract GSM-level info
    gsm_list = GSMList(gse)
    
    # Define placenta-related keywords
    placenta_keywords <- c("placenta", "chorionic", "decidua", "amnion", "trophoblast", "chorion", "fetal membrane", "umbilical")
    placenta_pattern <- paste(placenta_keywords, collapse = "|")
    
    # Count placenta-related samples
    placenta_samples_count = sum(
      map_lgl(gsm_list, ~ str_detect(tolower(.x@header$title %||% ""), placenta_pattern) |
                str_detect(tolower(.x@header$source_name_ch1 %||% ""), placenta_pattern))
    )
    
    # Extract per-GSM library info and instrument models
    library_strategy = gsm_list %>%
      map_chr(~ .x@header$library_strategy %||% NA) %>% unique() %>% paste(collapse = ", ")
    library_selection = gsm_list %>%
      map_chr(~ .x@header$library_selection %||% NA) %>% unique() %>% paste(collapse = ", ")
    instrument_model = gsm_list %>%
      map_chr(~ .x@header$instrument_model %||% NA) %>% unique() %>% paste(collapse = ", ")
    data_processing = gsm_list %>%
      map_chr(~ .x@header$data_processing %||% NA) %>% unique() %>% paste(collapse = ", ")
    extraction_protocol = gsm_list %>%
      map_chr(~ .x@header$extract_protocol_ch1 %||% .x@header$extract_protocol %||% NA) %>% unique() %>% paste(collapse = ", ")
    
    tibble(
      GEO_ID = geo_id,
      contact_country = meta_data$contact_country %||% NA,
      contact_email = meta_data$contact_email %||% NA,
      contact_institute = meta_data$contact_institute %||% NA,
      contact_name = gsub(",+", " ", meta_data$contact_name %||% NA),
      overall_design = meta_data$overall_design %||% NA,
      platform_id = meta_data$platform_id %||% NA,
      pubmed_id = meta_data$pubmed_id %||% NA,
      doi = meta_data$relation %>% str_extract("10\\.\\d{4,9}/[-._;()/:A-Z0-9]+") %||% NA,
      relation = if (!is.null(meta_data$relation)) paste(meta_data$relation, collapse = ", ") else NA,
      sample_number = if (!is.null(meta_data$sample_id)) length(meta_data$sample_id) else NA,
      sample_taxid = meta_data$sample_taxid %||% NA,
      submission_date = meta_data$submission_date %||% NA,
      last_update_date = meta_data$last_update_date %||% NA,
      main_topic = meta_data$summary %||% NA,
      supplementary_file = if (!is.null(meta_data$supplementary_file)) paste(meta_data$supplementary_file, collapse = ", ") else NA,
      title = meta_data$title %||% NA,
      experiment_type = meta_data$type %||% NA,
      placenta_samples = placenta_samples_count,   # NEW: number of placenta-related samples
      extraction_protocol = extraction_protocol,
      library_strategy = library_strategy,
      library_selection = library_selection,
      instrument_model = instrument_model,
      data_processing = data_processing
    )
    
  }, error = function(e) {
    tibble(GEO_ID = geo_id, error = as.character(e))
  })
}


#################
# Main workflow #
#################

# Load GEO IDs
gse_df = read.csv("gse.csv", header = FALSE, stringsAsFactors = FALSE)
geo_ids = as.character(gse_df[[1]])

# Limit for testing
geo_ids = geo_ids[1:10]

# Convert to GSE format
geo_ids = convert_to_gse(geo_ids)

# Progress bar
pb = progress_bar$new(
  format = "  Fetching [:bar] :percent ETA: :eta",
  total = length(geo_ids), clear = FALSE, width = 60
)

results = vector("list", length(geo_ids))

for (i in seq_along(geo_ids)) {
  results[[i]] = get_single_gse_full(geo_ids[i])
  pb$tick()
}

# Separate successes and failures
metadata_df = bind_rows(keep(results, ~ !"error" %in% names(.x)))
failed_ids = map_chr(keep(results, ~ "error" %in% names(.x)), "GEO_ID")

# Map taxonomy IDs
taxonomy = c(
  "10090" = "Mus musculus",
  "9315"  = "Notamacropus eugenii",
  "9606"  = "Homo sapiens",
  "9823"  = "Sus scrofa",
  "10116" = "Rattus norvegicus",
  "9913"  = "Bos taurus",
  "9361"  = "Dasypus novemcinctus",  
  "10092" = "Mus musculus domesticus",
  "9544" = "Macaca mulatta",  
  "9915"  = "Bos indicus",
  "9545"  = "Macaca nemestrina",
  "9986"  = "Oryctolagus cuniculus"
)

# Map organism first
metadata_df = metadata_df %>%
  mutate(organism = taxonomy[as.character(sample_taxid)]) %>%
  select(-sample_taxid)  

# Add superseries TRUE/FALSE column
metadata_df <- metadata_df %>%
  mutate(
    part_of_superseries = str_detect(relation, "GSE\\d+") & 
      !str_detect(relation, paste0("^", GEO_ID, "$"))
  )

# Fix dates for Excel
metadata_df = metadata_df %>%
  mutate(
    last_update_date = as.Date(parse_date_time(last_update_date, orders = c("b d Y"))),
    submission_date  = as.Date(parse_date_time(submission_date, orders = c("b d Y"))),
    last_update_date = format(last_update_date, "%m/%d/%Y"),
    submission_date  = format(submission_date, "%m/%d/%Y")
  )

# Save results
write_xlsx(list(
  Metadata = metadata_df,
  Failed   = tibble(Failed_GSE_IDs = failed_ids)
), "gse_metadata_full.xlsx")

print("Processing complete. Results saved to gse_metadata_full.xlsx")

