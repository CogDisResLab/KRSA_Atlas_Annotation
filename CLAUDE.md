# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This project provides scripts and workflows for updating kinase specificity annotations, focusing on the STK and PTK kinome. It leverages the KRSA (Kinome Random Sampling Analyzer) platform and the latest kinome data resources.

### Key Resources
- **KRSA**: Kinome Random Sampling Analyzer - https://github.com/CogDisResLab/KRSA
- **reKRSA**: reproducible KRSA - https://github.com/CogDisResLab/reKRSA
- **STK Kinome Atlas Paper**: A Kinome Atlas for Serine/Threonine Kinases (PubMed: 35012345)
- **PTK Kinome Atlas Paper**: A Kinome Atlas for Protein Tyrosine Kinases (PubMed: 36098765)
- **GraphQL API**: https://kinase-library.phosphosite.org/kinase-library-api/graphql

## Architecture

The project has a simple R-based pipeline with the following structure:

```
├── raw/                           # Raw layout files from KRSA (Layout.txt)
├── raw_backup/                    # Backup of raw files after UTF-8 conversion
├── data/
│   ├── input_sequence_data.csv.bz2    # Prepared peptide sequences
│   └── individual/                 # Downloaded kinase specificity data (per peptide/chip)
├── results/
│   ├── complete_kinase_specificity_map.csv.gz       # Full kinase map
│   ├── complete_kinase_specificity_map_raw.csv.gz   # Raw kinase map
│   └── MAPK-mapped-kinases.csv                 # MAPK-focused subset
├── scripts/
│   └── utf8_converter.sh            # Batch UTF-8 conversion script
└── R/
    ├── prepare_sequences.R          # Prepare peptide sequences with modifications
    ├── download_specificity_data.R  # Download data from GraphQL API
    └── create_specificity_metrics.R # Process results and generate metrics
```

## Data Pipeline Workflow

The annotation pipeline follows these steps:

1. **Layout Processing**: Read Layout.txt files from `raw/` directory
2. **UTF-8 Conversion**: Convert legacy files to UTF-8 using `scripts/utf8_converter.sh`
3. **Sequence Preparation**: Run `prepare_sequences.R` to:
   - Add * markers to S, T, Y residues
   - Generate phosphorylation variants
   - Handle ambiguous amino acids (B/Z)
   - Create phosphoprimer/phosphosite encodings
   - Output: `data/input_sequence_data.csv.bz2`
4. **Data Download**: Run `download_specificity_data.R` to:
   - Query GraphQL API for each prepared sequence
   - Download kinase specificity scores
   - Output: `results/complete_kinase_specificity_map.csv.gz`
5. **Metrics Generation**: Run `create_specificity_metrics.R` to:
   - Filter kinases (MAPK, JNK, P38, etc.)
   - Create intersection analysis
   - Generate UpSet plots

## Common Development Tasks

### Running the Full Pipeline

```bash
# Step 1: UTF-8 conversion (if needed)
bash scripts/utf8_converter.sh

# Step 2: Prepare sequences
Rscript prepare_sequences.R

# Step 3: Download data from API
Rscript download_specificity_data.R

# Step 4: Generate metrics
Rscript create_specificity_metrics.R
```

### Running Individual R Scripts

```bash
# Prepare sequences only
Rscript prepare_sequences.R

# Download data only (uses cached data when files already exist)
Rscript download_specificity_data.R

# Generate metrics
Rscript create_specificity_metrics.R
```

### Running a Single Test Script

To run a specific analysis script:
```bash
Rscript [script_name].R
```

### Key R Packages Used

- tidyverse (includes: dplyr, tidyr, purrr, tibble, readr, stringr)
- httr2 (for HTTP requests to GraphQL API)
- glue (for string interpolation)
- stringi (for string operations)
- UpSetR (for intersection visualization)

## Configuration

### .Rprofile

The project uses renv for package management. The `.Rprofile` automatically sources `renv/activate.R`.

### .gitignore

- `.Rproj.user/`, `.Rhistory`, `.RData`, `.Ruserdata`
- `*.log`
- `x-*.R` (downloaded data from API)
- `.vscode/`

**Note**: Data files are excluded:
- `data/*.bz2` and `results/*.gz` are tracked (compression makes them reasonable to version)
- `raw/` files are excluded (they're regenerated)
- `raw_backup/` is excluded

### renv

The project uses renv for reproducible R environments. Run `renv::init()` or `renv::restore()` as needed.

## File Formats

### Input: Layout.txt (from raw directory)
```
PeptideID<TAB>Sequence<TAB>GeneName<TAB>...
```

### Output: input_sequence_data.csv.bz2 columns
- `ID` - Composite identifier (Old_ID + phosphosite + priming + disambiguation)
- `PeptideID` - Original peptide identifier
- `source_sequence` - Original peptide sequence
- `prepared_sequence` - Modified sequence with phosphorylation markers
- `phosphosite` - Phosphorylation site location (e.g., "Y4")
- `priming_status` - Priming location encoding (e.g., "0b0", "0b1")
- `disambiguation` - Amino acid disambiguation code (e.g., "O0xO")
- `ChipArticleID` - KRSA article identifier

### GraphQL API Response
The GraphQL query returns kinase data including:
- `processedSequence` - Processed sequence
- `kinases` array with: family, geneName, name, percentile, percentileRank, score, scoreRank, uniprotId, displayName

## Kinase Families

The project focuses on these kinase families:
- **ERK**: ERK1, ERK2, ERK5, ERK7
- **JNK**: JNK1, JNK2, JNK3
- **P38**: P38A, P38B, P38D, P38G
- **TAK**: TAK1
- **DYRK**: DYRK1A, DYRK1B, DYRK2, DYRK3, DYRK4

## Stale Data Handling

When the GraphQL API has been deprecated (as indicated in recent commits), data must be sourced from alternative repositories. See commits:
- `526eea2` - Cleanup of old stale data
- `ee2f412` - Include new prepared sequences from alternative sources

## Chip Types

- **STK**: Serine/Threonine kinase chip (targets S/T residues)
- **PTK**: Protein Tyrosine kinase chip (targets Y residues)

## gstack

Use the **/browse** skill (from gstack) for all web browsing tasks.

### Never Use
- Never use `mcp__claude-in-chrome__*` tools — always use `/browse` instead.

### Available Skills

Use these skills for planning, design, testing, and development tasks:

- **/office-hours** — Brainstorming, idea exploration, startup thinking
- **/plan-ceo-review** — CEO/founder-mode strategic review, scope expansion
- **/plan-eng-review** — Engineering review, architecture, data flow, edge cases
- **/plan-design-review** — Design critique, visual audits, layout review
- **/design-consultation** — Design system creation, brand guidelines
- **/review** — Pre-landing PR review, code review
- **/ship** — Ship workflow: merge, CI, deploy, verify
- **/browse** — Headless browser for QA testing and site dogfooding
- **/qa** — Systematic QA testing, test and fix bugs
- **/qa-only** — QA report mode: test without fixing
- **/design-review** — Visual QA, design polish, fix visual issues
- **/setup-browser-cookies** — Import cookies for authenticated testing
- **/retro** — Engineering retrospective, commit history analysis
- **/investigate** — Systematic debugging with root cause analysis
- **/document-release** — Update documentation post-ship
- **/codex** — Code review, challenge mode, second opinion
- **/careful** — Safety mode: destructive command warnings
- **/freeze** — Restrict edits to specific directory
- **/guard** — Full safety: warnings + directory-scoped edits
- **/unfreeze** — Clear edit restrictions
- **/gstack-upgrade** — Upgrade gstack to latest version

### Example Usage

```bash
# To view a website:
/browse

# To test a feature on a live site:
/browse

# For authentication before testing:
/setup-browser-cookies
```
