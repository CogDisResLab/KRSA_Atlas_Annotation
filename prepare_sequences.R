# Comprehensive Peptide Sequence Preparation Script
# Processes peptide sequences with multiple modifications:
# - Adds * next to S, T, and Y residues
# - Generates variants with different phosphorylation sites
# - Handles ambiguous amino acids (B and Z)

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringi)
})

encode_phosphoprimer_location <- function(sequence) {
  cleaned <- str_remove_all(sequence, "\\W")
  str_split(cleaned, "") |>
    map_chr(
      \(split_string) {
        if_else(str_detect(split_string, "[a-z]"), "1", "0") |>
          str_c(collapse = "") |>
          strtoi(base = 2) |>
          as.character()
      }
    )
}

encode_phosphosite_location <- function(sequence) {
  index <- str_locate(sequence, "[STY]\\*")[[1]]
  residue <- str_sub(sequence, index, index)
  str_c(residue, index)
}

encode_aa_disambiguation <- function(sequence) {
  cleaned <- str_remove_all(sequence, "\\W")
  case_when(
    str_detect(cleaned, "SENT") ~ "B3xN",
    str_detect(cleaned, "SEDT") ~ "B3xD",
    .default = "O0xO"
  )
}

# Function to handle ambiguous amino acids (B and Z)
handle_ambiguous_aa <- function(seq) {
  # Check for B (asparagine or aspartic acid)
  if (str_detect(seq, "B")) {
    seq_n <- str_replace(seq, "B", "N")
    seq_d <- str_replace(seq, "B", "D")
    return(list(seq_n, seq_d))
  }

  # Check for Z (glutamine or glutamic acid)
  if (str_detect(seq, "Z")) {
    seq_e <- str_replace(seq, "Z", "E")
    seq_q <- str_replace(seq, "Z", "Q")
    return(list(seq_e, seq_q))
  }

  # If no ambiguous amino acids, return original sequence
  list(seq)
}

# Function to mark all phosphosites with a *
mark_phosphosites <- function(seq, chip) {
  if (chip == "STK") {
    s_positions <- str_locate_all(seq, "[ST]")[[1]][, 1]
  } else if (chip == "PTK") {
    s_positions <- str_locate_all(seq, "[Y]")[[1]][, 1]
  }

  sequences <- rep(NA, length(s_positions))

  for (idx in 1:length(sequences)) {
    replacement <- paste0(str_sub(seq, s_positions[idx], s_positions[idx]), "*")
    sequences[idx] <- stri_sub_replace(
      seq,
      from = s_positions[idx],
      to = s_positions[idx],
      replacement = replacement
    )
  }

  sequences
}

# Generate a permutation of primed (prephosphorylated) residues
mark_phosphoprimers <- function(seq) {
  non_query_phosphosites <- str_locate_all(seq, "[STY](?!\\*)")[[1]][, 1]

  primer_combinations <- map(
    0:length(non_query_phosphosites),
    ~ combn(non_query_phosphosites, .x, simplify = FALSE)
  ) |>
    unlist(recursive = FALSE)

  primer_variants <- map(
    primer_combinations,
    \(x) {
      modified_sequence <- seq
      for (loc in x) {
        extracted <- str_sub(modified_sequence, loc, loc)
        replacement <- str_to_lower(extracted)
        modified_sequence <- stri_sub_replace(modified_sequence,
          replacement = replacement, from = loc, to = loc
        )
      }
      modified_sequence
    }
  )

  primer_variants
}

# Main processing function
prepare_peptide_sequences <- function(input_data, chip) {
  # Prepare sequences with modifications
  processed_data <- input_data |>
    # Remove rows with NA sequences and those starting with 'p'
    filter(
      !is.na(Sequence),
      str_detect(Sequence, "\\(p[STY]\\)", negate = TRUE),
      str_detect(ID, fixed("Y253F"), negate = TRUE)
    ) |>
    # Expand sequences with ambiguous amino acids
    mutate(
      expanded_sequences = map(Sequence, handle_ambiguous_aa)
    ) |>
    unnest(expanded_sequences) |>
    unnest(expanded_sequences) |>
    # Generate phosphorylation variants
    mutate(
      phosphosite_variants = map(expanded_sequences, ~ mark_phosphosites(.x, chip))
    ) |>
    unnest(phosphosite_variants) |>
    mutate(phosphoprime_variants = map(phosphosite_variants, mark_phosphoprimers)) |>
    unnest(phosphoprime_variants) |>
    unnest(phosphoprime_variants) |>
    mutate(
      ambiguous_aa = encode_aa_disambiguation(phosphoprime_variants),
      phosphosite_location = encode_phosphosite_location(phosphoprime_variants),
      phosphoprime_location = str_c("0b", encode_phosphoprimer_location(phosphoprime_variants)),
      Old_ID = ID,
      ID = str_c(Old_ID, phosphosite_location, phosphoprime_location, ambiguous_aa, sep = "_")
    ) |>
    select(
      ID,
      source_sequence = Sequence,
      prepared_sequence = phosphoprime_variants,
      phosphosite = phosphosite_location,
      priming_status = phosphoprime_location,
      disambiguation = ambiguous_aa
    )

  # Return processed data for further inspection if needed
  processed_data |>
    unique() |>
    mutate(chip_type = chip)
}

# Load input data
ptk_peptides <- read_tsv(file.path("raw", "86412_Array_Layout.txt")) |>
  prepare_peptide_sequences(chip = "PTK")
stk_peptides <- read_tsv(file.path("raw", "87102_Array_Layout.txt")) |>
  prepare_peptide_sequences(chip = "STK")

# Process sequences
final_processed_sequences <- ptk_peptides |>
  bind_rows(stk_peptides) |>
  write_csv(file.path("data", "input_sequence_data.csv"))
