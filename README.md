# KRSA Atlas Annotation

<!-- badges: start -->
<!-- badges: end -->

## Project Overview

This project generates a **deterministic peptide → kinase specificity dataset** from PamGene PamChip layouts (STK/PTK), intended for **direct downstream use in KRSA**.

It is a **data generation layer**, not an analysis framework.

### Key Resources

- **[KRSA](https://github.com/CogDisResLab/KRSA)**: Kinome Random Sampling Analyzer
- **[reKRSA](https://github.com/CogDisResLab/reKRSA)**: Reproducible KRSA
- **[STK Kinome Atlas Paper](https://pubmed.ncbi.nlm.nih.gov/35012345/)**: A Kinome Atlas for Serine/Threonine Kinases
- **[PTK Kinome Atlas Paper](https://pubmed.ncbi.nlm.nih.gov/36098765/)**: A Kinome Atlas for Protein Tyrosine Kinases
- **[GraphQL API](https://kinase-library.phosphosite.org/kinase-library-api/graphql)**: Phosphosite Kinase Library API

## Architecture

```
KRSA_Atlas_Annotation/
├── raw/                              # Raw layout files from KRSA
│   ├── 86312_Array_Layout.txt       # PTK chip layout (tyrosine kinases)
│   ├── 86402_Array_Layout.txt       # PTK chip layout
│   ├── 87012_Array_Layout.txt       # STK chip layout (serine/threonine kinases)
│   ├── 87102_Array_Layout.txt
│   ├── 87202_Array_Layout.txt
│   ├── 90028_Array_Layout.txt
│   └── ...                          # Additional chip layouts
├── raw_backup/                      # Backup of raw files after UTF-8 conversion
├── data/
│   ├── input_sequence_data.csv.bz2  # Prepared peptide sequences
│   └── individual/                  # Downloaded kinase specificity data (per peptide/chip)
├── results/
│   └── complete_kinase_specificity_map.csv.gz   # Primary output artifact
├── scripts/
│   └── utf8_converter.sh            # Batch UTF-8 conversion script
├── R/
│   ├── prepare_sequences.R          # Peptide expansion and variant generation
│   ├── download_specificity_data.R  # Kinase annotation via GraphQL API
│   └── create_specificity_metrics.R # Data aggregation (not analysis)
├── renv/                            # R package environment
├── .Rprofile                       # R configuration (sources renv)
├── .gitignore                      # Git ignore rules
└── CLAUDE.md                       # Development guidelines
```

## Core Guarantees

### 1. Deterministic Outputs

Given the same input layouts, code version, and API state, the pipeline produces **reproducible outputs**.

### 2. Complete Peptide Expansion

All valid peptide variants are generated via:
- Phosphorylation site expansion
- Priming state encoding
- Ambiguous residue disambiguation (B/Z)

No valid variant should be silently dropped.

### 3. One-to-Many Mapping

Each peptide maps to **multiple kinases**. No artificial filtering by:
- Kinase family
- Pathway (e.g., MAPK)

### 4. Layout Fidelity

All input layouts are preserved in output context. If the same peptide appears in multiple layouts:
- It is **duplicated with distinct identifiers**
- No unintended collapsing across layouts

### 5. Clean Regeneration

The pipeline assumes:
- **No reliance on stale cached data**
- Cached API responses are optional optimizations, never a source of truth

## Non-Goals

This project is **not** intended for:
- Running KRSA
- Performing biological interpretation
- Kinase family-specific filtering (e.g., MAPK)
- Statistical analysis or enrichment

## Invariants

These must always hold true:

- Every `prepared_sequence`:
  - Originates from a valid input peptide
- Every output row:
  - Has a traceable lineage to: layout → peptide → variant → API response
- No mutation of:
  - Original peptide identity (`PeptideID`)
- No silent data loss during:
  - Expansion
  - API querying
  - Aggregation

## Mental Model

Think of the system as:

```
PamChip Layouts
↓
Peptide Expansion Engine
↓
Kinase Annotation Layer (GraphQL)
↓
Canonical Peptide → Kinase Map
↓
KRSA Input
```

## Data Pipeline Workflow

### Step 1: UTF-8 Normalization

Legacy files may contain encoding issues (e.g., µM instead of uM, CRLF line endings). The `utf8_converter.sh` script:
- Converts ISO-8859-1 to UTF-8 encoding
- Replaces µM with uM
- Strips trailing carriage returns (CRLF → LF)
- Preserves original files as backups in `raw_backup/`

**Note:** This is normalization only — no semantic changes.

### Step 2: Sequence Preparation

`prepare_sequences.R` processes peptide sequences:
- Adds `*` markers to S, T, Y residues for phosphorylation sites
- Generates phosphorylation variants (phosphoprimer/phosphosite encodings)
- Handles ambiguous amino acids (B = D/N, Z = E/Q)
- Adds ChipArticleID for tracing to KRSA articles

### Step 3: Kinase Annotation

`download_specificity_data.R` queries the GraphQL API:
- Retrieves kinase specificity scores for each prepared sequence
- Downloads data including family, geneName, name, percentile, score
- Outputs to `results/complete_kinase_specificity_map.csv.gz`

**Note:** No manual filtering of kinases. All responses are captured.

### Step 4: Aggregation

The aggregation step produces the unified mapping dataset:
- One row per peptide variant × kinase
- Contains full kinase metadata and scoring metrics
- Suitable for direct KRSA ingestion

The primary output artifact is:
```
results/complete_kinase_specificity_map.csv.gz
```

## Key Technologies

### R Packages

- **tidyverse**: Data manipulation (dplyr, tidyr, purrr, tibble, readr, stringr)
- **httr2**: HTTP requests to GraphQL API
- **glue**: String interpolation
- **stringi**: String operations
- **renv**: Reproducible package environments

## File Formats

### Input: Layout.txt (raw directory)

```
Row	Col	ID	Sequence	Tyr	SpotConcentration	UniprotAccession	Description	Xoff	Yoff
-1	-1	#REF	NA	NA	NA	NA	NA	0	0
1	1	41_654_666	LDGENIYIRHSNL	[660]	1000	P11171	Protein 4.1...	0	0
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

### GraphQL API Response Fields

- `processedSequence` - Processed sequence
- `kinases` array with:
  - `family`, `geneName`, `name`, `displayName`
  - `percentile`, `percentileRank`, `score`, `scoreRank`
  - `uniprotId`

## Running the Pipeline

### Full Pipeline

```bash
# Step 1: UTF-8 normalization
bash scripts/utf8_converter.sh

# Step 2: Sequence preparation
Rscript R/prepare_sequences.R

# Step 3: Kinase annotation
Rscript R/download_specificity_data.R

# Step 4: Aggregation (if needed)
# Note: Output is generated by download_specificity_data.R
```

### Individual Scripts

```bash
# Prepare sequences only
Rscript R/prepare_sequences.R

# Download data only (uses cached data when files exist as optional optimization)
Rscript R/download_specificity_data.R
```

**Important:**
- Always run pipeline in order: normalization → expansion → annotation → aggregation
- Do not manually edit intermediate files
- Do not reuse stale outputs after schema changes

## Configuration

### .Rprofile

The project uses renv for package management. The `.Rprofile` automatically sources `renv/activate.R` to ensure consistent package environments.

### .gitignore

```
.Rproj.user/
.Rhistory
.RData
.Ruserdata
*.log
x-*.R                              # downloaded data from API
data/individual/**/*
!data/*.bz2                       # Track compressed files
!results/*.gz                     # Track compressed results
.vscode/
raw/all-array-layouts
raw_backup
```

## Environment

- **R + renv** for package reproducibility
- **mise.toml** for environment variables
- GraphQL endpoint must be reachable

## Failure Modes

### Acceptable Failures

- API rate limits → handled via retry/batching
- Missing responses → logged, not silently ignored

### Unacceptable Failures

- Silent peptide loss
- Collapsing distinct peptide variants
- Mixing data across layouts
- Partial writes without signaling failure

## Versioning

- Uses **Conventional Commits**
- Any change that affects:
  - Peptide expansion logic
  - Output schema
- Must be treated as a **breaking or data-impacting change**

## Bottom Line

If this project is working correctly:

- Every peptide that *can* be represented → **is represented**
- Every variant that *can* be scored → **is scored**
- Every result is:
  - Reproducible
  - Traceable
  - Complete

Anything less is a bug.

## License

This project is provided for research purposes. Please refer to the original KRSA project and cited publications for licensing information.

## References

1. Kinome Random Sampling Analyzer. https://github.com/CogDisResLab/KRSA
2. reKRSA: reproducible KRSA. https://github.com/CogDisResLab/reKRSA
3. STK Kinome Atlas. https://doi.org/10.1038/s41586-022-05575-3
4. PTK Kinome Atlas. https://doi.org/10.1038/s41586-024-07407-y
5. Phosphosite Kinase Library API. https://kinase-library.phosphosite.org/

## Contact

For questions or contributions, please refer to the linked repositories and publications above.
