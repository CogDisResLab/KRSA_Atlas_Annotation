# Extract, clean and prepare the peptide sequences

library(tidyverse)
library(stringi)

layout <- read_tsv("raw/87102_Array_Layout.txt")

prepare_sequence <- function(seq) {

  modified <- str_locate_all(seq, "S|T") |>
    as.data.frame() |>
    mutate(
      char = str_sub(seq, start, end),
      replacement = str_c(char, "*"),
      new_sequence = stri_sub_replace(
        seq,
        from = start,
        to = end,
        replacement = replacement
      )
    )

  out <- modified |> pull(new_sequence) |> str_c(collapse = ", ")
}


final_list <- layout |>
  select(ID, Sequence) |>
  filter(!is.na(Sequence), str_detect(ID, '^p', negate = TRUE)) |>
  mutate(modified = map(Sequence, ~ prepare_sequence(.x))) |>
  separate_rows(modified, sep = ", ") |>
  write_csv("results/input_sequences")
