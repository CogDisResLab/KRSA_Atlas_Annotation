# Process the results and generate specifictiy metrics
suppressPackageStartupMessages({
    library(tidyverse)
    library(UpSetR)
})

mapkinases <- c(
    "JNK1", "JNK2", "JNK3",
    "ERK1", "ERK2", "ERK5", "ERK7",
    "P38A", "P38B", "P38D", "P38G"
    # "TAK1", "DYRK1A", "DYRK1B", "DYRK2",
    # "DYRK3", "DYRK4"
)

mapkinase_data <- read_csv(file.path("results", "complete_kinase_specificity_map.csv.gz")) |>
    filter(str_detect(ID, fixed("0b0")), DisplayName %in% mapkinases)

clean_data <- mapkinase_data |>
    select(ID, PeptideID, GeneName, DisplayName, Percentile) |>
    filter(Percentile >= 80L) |>
    mutate(Group = case_when(
        str_detect(DisplayName, "ERK") ~ "ERK",
        str_detect(DisplayName, "JNK") ~ "JNK",
        str_detect(DisplayName, "P38") ~ "P38",
        str_detect(DisplayName, "TAK") ~ "TAK",
        str_detect(DisplayName, "DYRK1") ~ "DYRK1",
        str_detect(DisplayName, "DYRK2") ~ "DYRK2",
        str_detect(DisplayName, "DYRK3") ~ "DYRK3",
        str_detect(DisplayName, "DYRK4") ~ "DYRK4",
    )) |>
    select(Group, DisplayName, PeptideID) |>
    unique() |>
    write_csv(file.path("results", "MAPK-mapped-kinases.csv"))

clean_lists <- clean_data |>
    group_by(Group, DisplayName) |>
    group_split()

common_erk12 <- clean_lists |>
    keep(~ .x[["Group"]][[1L]] == "ERK") |>
    set_names(c("ERK1", "ERK2", "ERK5", "ERK7")) |>
    keep_at(at = c(rep(TRUE, 2L), rep(FALSE, 2L))) |>
    map(~ pull(.x, PeptideID)) |>
    reduce(intersect)

common_erk_all <- clean_lists |>
    keep(~ .x[["Group"]][[1L]] == "ERK") |>
    set_names(c("ERK1", "ERK2", "ERK5", "ERK7")) |>
    map(~ pull(.x, PeptideID)) |>
    reduce(intersect)

common_jnk <- clean_lists |>
    keep(~ .x[["Group"]][[1L]] == "JNK") |>
    set_names(c("JNK1", "JNK2", "JNK3")) |>
    map(~ pull(.x, PeptideID)) |>
    reduce(intersect)

common_p38 <- clean_lists |>
    keep(~ .x[["Group"]][[1L]] == "P38") |>
    set_names(c("P38A", "P38B", "P38D", "P38D")) |>
    map(~ pull(.x, PeptideID)) |>
    reduce(intersect)


clean_groups <- clean_data |>
    select(DisplayName, PeptideID) |>
    filter(!DisplayName %in% c("ERK5", "ERK7")) |>
    nest(data = PeptideID) |>
    mutate(Peptide = map(data, ~ .x[["PeptideID"]])) |>
    select(-data) |>
    deframe()

sets <- upset(fromList(clean_groups),
    nsets = 11L,
    order.by = "freq", nintersects = 100L,
    empty.intersections = "on"
)

clean_groups <- clean_data |>
    select(Group, PeptideID) |>
    unique() |>
    nest(data = PeptideID) |>
    mutate(Peptide = map(data, ~ .x[["PeptideID"]])) |>
    select(-data) |>
    deframe()

sets <- upset(fromList(clean_groups),
    nsets = 11L,
    order.by = "freq", nintersects = 100L,
    empty.intersections = "on"
)
