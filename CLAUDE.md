# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This projects as a reproducible workflow for generating peptide-to-kinase annotations for the [KRSA](https://github.com/CogDisResLab/KRSA) and [reKRSA](https://github.com/CogDisResLabreKRSA) kinome analysis systems. This project utilizes the PamChip annotation files for each PamChip article numbers and the Sequence-based kinome atlas developed by Johnson et al (2023) Yaron-Barir et al (2024) to identify likely upstream kinases based on the substrate sequences.

### Key Resources
- **STK Kinome Atlas Paper**: [An atlas of substrate specificities for the human serine/threonine kinome](https://doi.org/10.1038/s41586-022-05575-3)
- **PTK Kinome Atlas Paper**: [The intrinsic substrate specificity of the human tyrosine kinome](https://doi.org/10.1038/s41586-024-07407-y)
- **GraphQL API**: https://kinase-library.phosphosite.org/kinase-library-api/graphql

## Architecture

The project has a simple R-based pipeline with the following structure:

```
├── raw/                           # Raw layout files from KRSA (Layout.txt)
├── data/
│   ├── input_sequence_data.csv.bz2    # Prepared peptide sequences
│   └── individual/                 # Downloaded kinase specificity data (per peptide/chip)
├── results/
│   ├── complete_kinase_specificity_map.csv.gz       # Kinase map filtered to the top 20 targets by score for each peptide
│   ├── complete_kinase_specificity_map_raw.csv.gz   # Raw kinase map including all results for each peptide on each chip
├── scripts/
│   └── utf8_converter.sh            # Utility batch UTF-8 conversion script
└── R/
    ├── prepare_sequences.R          # Extract the sequences from each layout file and format the to be queried from the API
    ├── download_specificity_data.R  # Download data from GraphQL API
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

## Common Development Tasks

### Running the Full Pipeline

```bash
# Step 1: UTF-8 conversion (if needed)
bash scripts/utf8_converter.sh

# Step 2: Prepare sequences
Rscript R/prepare_sequences.R

# Step 3: Download data from API
Rscript R/download_specificity_data.R

# Step 4: Generate metrics
Rscript R/create_specificity_metrics.R
```

### Running Individual R Scripts

```bash
# Prepare sequences only
Rscript R/prepare_sequences.R

# Download data only (uses cached data when files already exist)
Rscript R/download_specificity_data.R

# Generate metrics
Rscript R/create_specificity_metrics.R
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
