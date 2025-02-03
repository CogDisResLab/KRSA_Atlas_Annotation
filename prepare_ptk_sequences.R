# Extract, clean and prepare the peptide sequences

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringi)
})

layout <- read_tsv(file.path("raw", "86412_Array_Layout.txt"))

prepare_sequence <- function(seq) {
  modified <- str_locate_all(seq, "Y") |>
    as.data.frame() |>
    mutate(
      char = str_sub(seq, start, end),
      replacement = str_c(char, "*"),
      new_sequence = stri_sub_replace(
        seq,
        from = start,
        to = end,
        replacement = replacement
      ),
      site = str_c(char, start)
    )

  out <- modified

  out
}


final_list <- layout |>
  select(ID, Sequence) |>
  filter(!is.na(Sequence), str_detect(ID, "^p", negate = TRUE)) |>
  mutate(modified = map(Sequence, ~ prepare_sequence(.x))) |>
  unnest(modified) |>
  select(-start, -end, -char, -replacement) |>
  write_csv(file.path("data", "ptk_input_sequence_data.csv"))
