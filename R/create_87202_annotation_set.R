# Create an annotation for the 87012 version of the chip.

library(tidyverse)

data_files <- list.files(file.path("data", "individual", "STK", "ART87202"),
    recursive = TRUE,
    full.names = TRUE
)

# Load the simple annotation data into one frame
annotation_data <- data_files |>
    set_names(\(x) basename(dirname(x))) |>
    map(read_csv) |>
    map(\(df) filter(df, percentile > 80L)) |>
    bind_rows(.id = "Peptide")

# Diagnostics for each peptide by kinase and annotations
clean_annotation_data <- annotation_data |>
    select(Peptide, displayName, scoreRank) |>
    summarise(scoreRank = mean(scoreRank), .by = c(Peptide, displayName)) |>
    nest(.by = Peptide) |>
    mutate(
        n_kinase = map_dbl(data, nrow),
        min_rank = map_dbl(data, \(df) min(df[["scoreRank"]])),
        max_rank = map_dbl(data, \(df) max(df[["scoreRank"]]))
    ) |>
    select(-data)

# Clean the dataset by averaging multiple phosphosites
# and selecting the relevant columsn
clean_annotation_data <- annotation_data |>
    select(Peptide, displayName, scoreRank) |>
    summarise(scoreRank = mean(scoreRank), .by = c(Peptide, displayName)) |>
    group_by(Peptide) |>
    arrange(scoreRank, .by_group = TRUE) |>
    slice_min(order_by = scoreRank, prop = 0.25, with_ties = TRUE) |>
    ungroup()

annotation_mapping <- clean_annotation_data |>
    select(-scoreRank) |>
    nest(.by = Peptide) |>
    mutate(
        Kinase = map(data, \(df) {
            df |>
                pull(displayName) |>
                str_flatten(collapse = " ")
        })
    ) |>
    select(-data) |>
    unnest(Kinase) |>
    select(Substrates = Peptide, Kinases = Kinase)

annotation_mapping |>
    write_csv("results/ART87202-STK-Mapping.csv")

annotation_coverage <- clean_annotation_data |>
    select(Substrates = Peptide, Kin = displayName) |>
    unique()

annotation_coverage |>
    write_csv("results/ART87202-STK-Coverage.csv")
