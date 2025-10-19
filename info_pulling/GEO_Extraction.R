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
library(rvest)

#---------------------------------
# Clear GEOquery cache helper, helpful for better testing
#---------------------------------
unlink("geo_cache", recursive = TRUE)
dir.create("geo_cache", showWarnings = FALSE)

#---------------------------------
# Utility functions
#---------------------------------
`%||%` = function(x, y) if (!is.null(x)) x else y

convert_to_gse = function(ids) {
  paste0("GSE", sub("^2000?|^200?", "", as.character(ids)))
}

#---------------------------------
# GSM fields extractor (all in one pass)
#---------------------------------
gsm_cache_dir = file.path("geo_cache", "gsm_pages")
dir.create(gsm_cache_dir, showWarnings = FALSE, recursive = TRUE)

extract_gsm_fields = function(gsm_id) {
  cache_file = file.path(gsm_cache_dir, paste0(gsm_id, ".rds"))
  
  if (file.exists(cache_file)) {
    text = readRDS(cache_file)
  } else {
    page = tryCatch(
      read_html(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gsm_id)),
      error = function(e) return(NA)
    )
    if (is.na(page)) return(rep(NA_character_, 5))
    text = html_text2(page)
    saveRDS(text, cache_file)
    Sys.sleep(0.2)  # small delay
  }
  
  keywords = c("Library strategy", "Library selection", "Library source",
                "Instrument model", "Extracted molecule")
  
  sapply(keywords, function(k) {
    pattern = paste0("(?i)", k, "\\s*[:\\s]?\\s*([^;,\n]+)")
    match = str_match_all(text, pattern)[[1]]
    values = unique(str_trim(match[,2]))
    values = values[values != ""]
    if (length(values) == 0) NA_character_ else paste(values, collapse=", ")
  })
}

#---------------------------------
# Fetch single GSE
#---------------------------------
get_single_gse_full = function(geo_id) {
  cache_file = file.path("geo_cache", paste0(geo_id, ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  
  out = tryCatch({
    gse = getGEO(geo_id, GSEMatrix = FALSE, AnnotGPL = FALSE)
    meta_data = Meta(gse)
    gsm_list = GSMList(gse)
    gsm_ids = names(gsm_list)
    
    # placenta sample count
    placenta_keywords = c("placenta","chorionic","decidua","trophoblast","chorion")
    placenta_pattern = paste(placenta_keywords, collapse = "|")
    placenta_samples_count = sum(
      map_lgl(gsm_list, ~ str_detect(tolower(.x@header$title %||% ""), placenta_pattern) |
                str_detect(tolower(.x@header$source_name_ch1 %||% ""), placenta_pattern))
    )
    
    # Scrape GSMs for all five fields in parallel
    gsm_values = future_map(gsm_ids, extract_gsm_fields)
    
    library_strategy   = sapply(gsm_values, `[[`, 1) %>% na.omit() %>% unique() %>% paste(collapse=", ")
    library_selection  = sapply(gsm_values, `[[`, 2) %>% na.omit() %>% unique() %>% paste(collapse=", ")
    library_source     = sapply(gsm_values, `[[`, 3) %>% na.omit() %>% unique() %>% paste(collapse=", ")
    instrument_model   = sapply(gsm_values, `[[`, 4) %>% na.omit() %>% unique() %>% paste(collapse=", ")
    extracted_molecule = sapply(gsm_values, `[[`, 5) %>% na.omit() %>% unique() %>% paste(collapse=", ")
    
    # Build tibble with GSE-level metadata
    tibble(
      GEO_ID = geo_id,
      contact_country = meta_data$contact_country %||% NA,
      contact_email = meta_data$contact_email %||% NA,
      contact_institute = meta_data$contact_institute %||% NA,
      contact_name = gsub(",+", " ", meta_data$contact_name %||% NA),
      overall_design = meta_data$overall_design %||% NA,
      platform_id = if (!is.null(meta_data$platform_id)) paste(meta_data$platform_id, collapse=", ") else NA,
      pubmed_id = meta_data$pubmed_id %||% NA,
      relation = if (!is.null(meta_data$relation)) paste(meta_data$relation, collapse=", ") else NA,
      sample_number = if (!is.null(meta_data$sample_id)) length(meta_data$sample_id) else NA,
      sample_taxid = if (!is.null(meta_data$sample_taxid)) paste(meta_data$sample_taxid, collapse=", ") else NA,  # =- added back
      submission_date = meta_data$submission_date %||% NA,
      last_update_date = meta_data$last_update_date %||% NA,
      main_topic = meta_data$summary %||% NA,
      supplementary_file = if (!is.null(meta_data$supplementary_file)) paste(meta_data$supplementary_file, collapse=", ") else NA,
      title = meta_data$title %||% NA,
      experiment_type = if (!is.null(meta_data$type)) paste(meta_data$type, collapse=", ") else NA,
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
  
  saveRDS(out, cache_file)
  out
}

#---------------------------------
# Main workflow
#---------------------------------
gse_df = read.csv("gse.csv", header = FALSE, stringsAsFactors = FALSE)
geo_ids = convert_to_gse(as.character(gse_df[[1]]))
geo_ids = geo_ids[1:200]  # adjust for testing

plan(multisession, workers = 4)  # adjust cores
handlers(global = TRUE)

with_progress({
  p = progressor(along = geo_ids)
  results = future_map(geo_ids, function(id) {
    out = get_single_gse_full(id)
    p()
    out
  }, .options = furrr_options(seed = TRUE))
})

#---------------------------------
# Combine results
#---------------------------------
metadata_df = rbindlist(keep(results, ~ !"error" %in% names(.x)), fill = TRUE)
failed_ids = map_chr(keep(results, ~ "error" %in% names(.x)), "GEO_ID")

#---------------------------------
# Superseries handling
#---------------------------------
metadata_df = metadata_df %>%
  mutate(
    Superseries = map_chr(
      relation,
      ~ {
        if (is.null(.x) || is.na(.x) || .x == "") return(NA_character_)
        hits = str_extract_all(.x, "GSE\\d+")[[1]]
        if (length(hits) == 0) NA_character_ else paste(unique(hits), collapse=", ")
      }
    )
  )

#---------------------------------
# Map taxonomy
#---------------------------------
taxonomy = c(
  "10090" = "Mus musculus","9315"="Notamacropus eugenii","9606"="Homo sapiens",
  "9823"="Sus scrofa","10116"="Rattus norvegicus","9913"="Bos taurus",
  "9361"="Dasypus novemcinctus","10092"="Mus musculus domesticus","9544"="Macaca mulatta",
  "9915"="Bos indicus","9545"="Macaca nemestrina","9986"="Oryctolagus cuniculus"
)
metadata_df$organism = taxonomy[as.character(metadata_df$sample_taxid %||% NA)]

#---------------------------------
# Remove sample_taxid only at the end
#---------------------------------
metadata_df = metadata_df %>% select(-sample_taxid)

#---------------------------------
# Format dates
#---------------------------------
metadata_df = metadata_df %>%
  mutate(
    last_update_date = format(as.Date(parse_date_time(last_update_date, orders = c("b d Y"))), "%m/%d/%Y"),
    submission_date  = format(as.Date(parse_date_time(submission_date, orders = c("b d Y"))), "%m/%d/%Y")
  )

#---------------------------------
# Save results
#---------------------------------
write_xlsx(list(
  Metadata = metadata_df,
  Failed   = tibble(Failed_GSE_IDs = failed_ids)
), "gse_metadata_full.xlsx")

print("Processing complete. Results saved to gse_metadata_full.xlsx")
