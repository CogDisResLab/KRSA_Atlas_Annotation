# Get Data from the site

library(tidyverse)
library(httr)
library(jsonlite)

template <-
  "https://kinase-library.phosphosite.org/api/scorer/score-site/{modified}/Homo%20sapiens%20(Human)/OCHOA/true/false"

data <- read_csv("data/input_sequence_data.csv") |>
  mutate(request = str_glue(template, modified = new_sequence))

process_request <- function(peptide_id, sequence, site, url) {
  output_file <- file.path("data",
                           "individual",
                           str_glue("{peptide_id}_{site}.csv"))

  col_spec <- cols(
    .default = col_character(),
    sitePosition = col_double(),
    score = col_double(),
    scoreQuantile = col_double(),
    scorePosition = col_double(),
    scoreDistributionSize = col_double(),
    scoreRank = col_double(),
    quantileRank = col_double(),
    id = col_double(),
    weakCatalyticActivity = col_logical()
  )

  if (!file.exists(output_file)) {
    message(str_glue("Downloading data for {peptide_id} for site {site}"))
    Sys.sleep(1)
    res <- GET(url)

    data <- content(res, as = "text", encoding = "utf8") |>
      fromJSON() |>
      as_tibble() |>
      unnest(scores) |>
      unnest(motif) |>
      mutate(peptide_id = peptide_id,
             sequence = sequence) |>
      select(-visibility) |>
      write_csv(output_file)
  } else {
    message(str_glue("Data for {peptide_id} for site {site} already downloaded"))
    data <-
      read_csv(output_file,
               col_types = col_spec,
               show_col_types = FALSE)
  }
}

all_data <- data |>
  pmap_dfr( ~ process_request(..1, ..2, ..4, ..5)) |>
  write_csv("results/complete_kinase_specificity_raw.csv") |>
  select(
    peptide_id,
    sequence,
    siteLabel,
    kinase_name = name,
    gene_name = geneName,
    kinase_group = motifGroup,
    kinase_type = type,
    score,
    score_quantile = scoreQuantile,
    score_rank = scoreRank,
    quantile_rank = quantileRank,
    score_position = scorePosition,
    score_distribution_size = scoreDistributionSize,
    weak_activity = weakCatalyticActivity
  ) |>
  write_csv("results/complete_stk_specificity_map.csv")
