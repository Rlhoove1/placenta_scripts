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
library(progressr)
library(tibble)
library(furrr)
library(data.table)

#####################
# Utility functions #
#####################

`%||%` = function(x, y) if (!is.null(x)) x else y

convert_to_gse = function(ids) {
  paste0("GSE", sub("^2000?|^200?", "", as.character(ids)))
}

#--------------------------
# Fetch single GSE
#--------------------------
get_single_gse_full = function(geo_id) {
  tryCatch({
    gse = getGEO(geo_id, GSEMatrix = FALSE, AnnotGPL = FALSE)
    meta_data = Meta(gse)
    gsm_list = GSMList(gse)
    
    placenta_keywords = c("placenta","chorionic","decidua","trophoblast","chorion")
    placenta_pattern = paste(placenta_keywords, collapse = "|")
    
    placenta_samples_count = sum(
      map_lgl(gsm_list, ~ str_detect(tolower(.x@header$title %||% ""), placenta_pattern) |
                str_detect(tolower(.x@header$source_name_ch1 %||% ""), placenta_pattern))
    )
    
    # Deduplicated GSM-level columns
    library_strategy = gsm_list %>% map_chr(~ .x@header$library_strategy %||% NA) %>% unique() %>% paste(collapse=", ")
    library_selection = gsm_list %>% map_chr(~ .x@header$library_selection %||% NA) %>% unique() %>% paste(collapse=", ")
    instrument_model = gsm_list %>% map_chr(~ .x@header$instrument_model %||% .x@header$instrument_model_ch1 %||% NA) %>% unique() %>% paste(collapse=", ")
    library_source = gsm_list %>% map_chr(~ .x@header$library_source %||% NA) %>% unique() %>% paste(collapse=", ")
    extracted_molecule = gsm_list %>% map_chr(~ .x@header$extracted_molecule_ch1 %||% NA) %>% unique() %>% paste(collapse=", ")
    
    # GSE-level metadata
    relation = if (!is.null(meta_data$relation)) paste(meta_data$relation, collapse = ", ") else NA
    supplementary_file = if (!is.null(meta_data$supplementary_file)) paste(meta_data$supplementary_file, collapse = ", ") else NA
    platform_id = if (!is.null(meta_data$platform_id)) paste(meta_data$platform_id, collapse = ", ") else NA
    sample_taxid = if (!is.null(meta_data$sample_taxid)) paste(meta_data$sample_taxid, collapse = ", ") else NA
    experiment_type = if (!is.null(meta_data$type)) paste(meta_data$type, collapse = ", ") else NA
    
    tibble(
      GEO_ID = geo_id,
      contact_country = meta_data$contact_country %||% NA,
      contact_email = meta_data$contact_email %||% NA,
      contact_institute = meta_data$contact_institute %||% NA,
      contact_name = gsub(",+", " ", meta_data$contact_name %||% NA),
      overall_design = meta_data$overall_design %||% NA,
      platform_id = platform_id,
      pubmed_id = meta_data$pubmed_id %||% NA,
      doi = {
        tmp = str_extract(relation %||% "", "10\\.\\d{4,9}/[-._;()/:A-Z0-9]+")
        if (length(tmp) == 0 || is.na(tmp) || tmp == "") NA_character_ else tmp
      },
      relation = relation,
      sample_number = if (!is.null(meta_data$sample_id)) length(meta_data$sample_id) else NA,
      sample_taxid = sample_taxid,
      submission_date = meta_data$submission_date %||% NA,
      last_update_date = meta_data$last_update_date %||% NA,
      main_topic = meta_data$summary %||% NA,
      supplementary_file = supplementary_file,
      title = meta_data$title %||% NA,
      experiment_type = experiment_type,
      placenta_samples = placenta_samples_count,
      library_strategy = library_strategy,
      library_selection = library_selection,
      library_source = library_source,
      instrument_model = instrument_model,
      extracted_molecule = extracted_molecule
    )
    
  }, error = function(e) {
    tibble(GEO_ID = geo_id, error = as.character(e))
  })
}

#--------------------------
# Cached wrapper with retries
#--------------------------
cache_dir = "geo_cache"
dir.create(cache_dir, showWarnings = FALSE)

get_single_gse_full_cached = function(geo_id, tries = 3) {
  cache_file = file.path(cache_dir, paste0(geo_id, ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  
  attempt = 1
  out = NULL
  while(attempt <= tries) {
    out = tryCatch(get_single_gse_full(geo_id), error = function(e) NULL)
    if (!is.null(out)) break
    Sys.sleep(3)
    attempt = attempt + 1
  }
  if (is.null(out)) out = tibble(GEO_ID = geo_id, error = "Failed after retries")
  saveRDS(out, cache_file)
  return(out)
}

#--------------------------
# Main workflow
#--------------------------
gse_df = read.csv("gse.csv", header = FALSE, stringsAsFactors = FALSE)
geo_ids = as.character(gse_df[[1]])
geo_ids = geo_ids[1:25] # adjust as needed
geo_ids = convert_to_gse(geo_ids)

plan(multisession, workers = 4)  # adjust cores

handlers(global = TRUE)
with_progress({
  p = progressor(along = geo_ids)
  results = future_map(geo_ids, function(id) {
    out = get_single_gse_full_cached(id)
    p()
    out
  }, .options = furrr_options(seed = TRUE))
})

#--------------------------
# Combine results
#--------------------------
metadata_df = rbindlist(keep(results, ~ !"error" %in% names(.x)), fill=TRUE)
failed_ids  = map_chr(keep(results, ~ "error" %in% names(.x)), "GEO_ID")

#--------------------------
# Superseries handling 
#--------------------------
metadata_df = metadata_df %>%
  mutate(
    Superseries = map_chr(
      relation,
      ~ {
        if (is.null(.x) || is.na(.x) || .x == "") {
          NA_character_
        } else {
          hits = str_extract_all(.x, "GSE\\d+")[[1]]
          if (length(hits) == 0) NA_character_ else paste(unique(hits), collapse = ", ")
        }
      }
    )
  )

#--------------------------
# Map taxonomy
#--------------------------
taxonomy = c(
  "10090" = "Mus musculus","9315"="Notamacropus eugenii","9606"="Homo sapiens",
  "9823"="Sus scrofa","10116"="Rattus norvegicus","9913"="Bos taurus",
  "9361"="Dasypus novemcinctus","10092"="Mus musculus domesticus","9544"="Macaca mulatta",
  "9915"="Bos indicus","9545"="Macaca nemestrina","9986"="Oryctolagus cuniculus"
)
metadata_df$organism = taxonomy[as.character(metadata_df$sample_taxid %||% NA)]
# Remove sample_taxid column
metadata_df = metadata_df %>% select(-sample_taxid)
#--------------------------
# Format dates
#--------------------------
metadata_df = metadata_df %>%
  mutate(
    last_update_date = as.Date(parse_date_time(last_update_date, orders = c("b d Y"))),
    submission_date  = as.Date(parse_date_time(submission_date, orders = c("b d Y"))),
    last_update_date = format(last_update_date, "%m/%d/%Y"),
    submission_date  = format(submission_date, "%m/%d/%Y")
  )

#--------------------------
# Save results
#--------------------------
write_xlsx(list(
  Metadata = metadata_df,
  Failed   = tibble(Failed_GSE_IDs = failed_ids)
), "gse_metadata_full.xlsx")

print("Processing complete. Results saved to gse_metadata_full.xlsx")


