# Get Data from the site

suppressPackageStartupMessages({
    library(tidyverse)
    library(httr2)
})

data <- read_csv(file.path("data", "input_sequence_data.csv"))

process_request <- function(id, peptide_id, sequence, chip) {
    output_file <- file.path(
        "data",
        "individual",
        chip,
        peptide_id,
        str_glue("{id}.csv")
    )

    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

    if (file.exists(output_file)) {
        message(str_glue("Data for {peptide_id} already downloaded"))

        read_csv(output_file,
            show_col_types = FALSE
        )
    } else {
        message(str_glue("Downloading data for {id}"))
        user_agent <- str_glue("KRSA Annotation; https://github.com/CogDisResLab/stk_annotation_update; {gitr::gitr_current_sha()} ")
        Sys.sleep(1L)
        graphql_url <- "https://kinase-library.phosphosite.org/kinase-library-api/graphql"
        graphql_query <- "
            query ScoreSite($sequence: String!,
            $pp: Boolean,
            $non_canonical: Boolean) {
  scoreSite(sequence: $sequence, pp: $pp, nonCanonical: $non_canonical) {
    processedSequence
    kinases {
      family
      geneName
      name
      percentile
      percentileRank
      score
      scoreRank
      uniprotId
      displayName
    }
    __typename
  }
}"

        graphql_vars <- list(
            sequence = sequence,
            pp = TRUE,
            non_canonical = FALSE
        )
        req <- request(graphql_url) |>
            req_headers(
                "Content-Type" = "application/json", # nolint: nonportable_path_linter.
            ) |>
            req_body_json(list(
                query = graphql_query,
                variables = graphql_vars,
                operationName = "ScoreSite"
            )) |>
            req_user_agent(user_agent)
        res <- req_perform(req)
        resp <- resp_body_json(res, simplifyVector = TRUE)
        # Defensive: check for errors and expected structure
        if (!is.null(resp[["errors"]])) {
            stop("GraphQL error:", resp[["errors"]][[1L]][["message"]],
                call. = FALSE
            )
        }
        site <- resp[["data"]][["scoreSite"]]
        kinases <- site[["kinases"]]
        kinases |>
            mutate(
                id = id,
                peptide_id = peptide_id,
                sequence = sequence,
                processed_sequence = site[["processedSequence"]]
            ) |>
            write_csv(output_file)
    }
}

all_data <- data |>
    select(ID, PeptideID, prepared_sequence, chip_type) |>
    pmap(~ process_request(..1, ..2, ..3, ..4)) |>
    bind_rows() |>
    write_csv(file.path("results", "complete_kinase_specificity_map_raw.csv.gz")) |>
    select(
        ID = id,
        PeptideID = peptide_id,
        SourceSequence = sequence,
        ProcessedSequence = processed_sequence,
        Family = family,
        GeneName = geneName,
        DisplayName = displayName,
        Percentile = percentile,
        PercentileRank = percentileRank,
        Score = score,
        ScoreRank = scoreRank,
        UniprotID = uniprotId
    ) |>
    write_csv(file.path("results", "complete_kinase_specificity_map.csv.gz"))
