# Installation Guide

**R 4.6.0 / Bioconductor 3.23**

Four options are provided. Choose based on your environment:

| Option | Best for | Reproducibility |
|--------|----------|-----------------|
| [A. renv](#a-renv-recommended) | R users with exact version lock | ★★★★★ |
| [B. Docker](#b-docker) | Any OS, fully isolated | ★★★★★ |
| [C. Singularity](#c-singularity-hpc) | HPC / cluster environments | ★★★★★ |
| [D. conda](#d-conda) | conda users (may lag Bioc releases) | ★★★☆☆ |
| [E. Manual R install](#e-manual-r-install) | Quick local setup | ★★★☆☆ |

---

## A. renv (Recommended)

`renv.lock` in this repo pins every package to the exact version used in the original analysis.

```bash
# 1. Open R in the project root
R

# 2. Install renv (if not already installed)
install.packages("renv")

# 3. Restore the exact environment from the lock file
renv::restore()

# 4. Run the pipeline
source("scripts/run_pipeline.R")
```

> If R asks to install packages to a personal library, answer **yes**.  
> Restoration may take 10–20 minutes on first run.

---

## B. Docker

No R installation required. All dependencies are bundled in the image.

### Prerequisites
- [Docker Desktop](https://www.docker.com/products/docker-desktop/) (Windows/macOS) or Docker Engine (Linux)

### Build and run

```bash
# Build the image (one-time, ~15 min)
docker build -t alopecia-fmt .

# Run the full pipeline with demo data
docker run --rm \
  -v "$(pwd)/sim_inputs:/project/sim_inputs" \
  -v "$(pwd)/outputs:/project/outputs" \
  alopecia-fmt

# Interactive shell (for development)
docker run --rm -it \
  -v "$(pwd)/sim_inputs:/project/sim_inputs" \
  -v "$(pwd)/outputs:/project/outputs" \
  alopecia-fmt bash
```

### Mount external mapping files (Steps 6–7)

Steps 6 and 7 require MetaCyc/KEGG files that are not bundled in the image.
Mount them at runtime:

```bash
docker run --rm \
  -v "$(pwd)/sim_inputs:/project/sim_inputs" \
  -v "$(pwd)/outputs:/project/outputs" \
  -v "/path/to/map_metacyc-pwy_name.txt.gz:/project/map_metacyc-pwy_name.txt.gz" \
  -v "/path/to/map_metacyc-pwy_lineage.tsv:/project/map_metacyc-pwy_lineage.tsv" \
  -v "/path/to/ko_ec_cache.rds:/project/ko_ec_cache.rds" \
  alopecia-fmt
```

---

## C. Singularity (HPC)

For cluster environments where Docker is not available.

### Build from Docker image

```bash
# Pull and convert Docker image to SIF
singularity pull alopecia-fmt.sif docker://bioconductor/bioconductor_docker:RELEASE_3_23

# Or build locally from Dockerfile (requires Docker installed first)
docker build -t alopecia-fmt .
singularity build alopecia-fmt.sif docker-daemon://alopecia-fmt:latest
```

### Run

```bash
# Run the pipeline
singularity exec \
  --bind $(pwd)/sim_inputs:/project/sim_inputs \
  --bind $(pwd)/outputs:/project/outputs \
  alopecia-fmt.sif \
  Rscript /project/scripts/run_pipeline.R

# Interactive shell
singularity shell \
  --bind $(pwd)/sim_inputs:/project/sim_inputs \
  --bind $(pwd)/outputs:/project/outputs \
  alopecia-fmt.sif
```

### Submit to SLURM

```bash
#!/bin/bash
#SBATCH --job-name=alopecia-fmt
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00

singularity exec \
  --bind $SLURM_SUBMIT_DIR/sim_inputs:/project/sim_inputs \
  --bind $SLURM_SUBMIT_DIR/outputs:/project/outputs \
  alopecia-fmt.sif \
  Rscript /project/scripts/run_pipeline.R
```

---

## D. conda

### Prerequisites
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge)

### Create environment

```bash
# Create from environment.yml (mamba is faster than conda)
mamba env create -f environment.yml

# Or with conda
conda env create -f environment.yml

# Activate
conda activate alopecia-fmt

# Run R
Rscript scripts/run_pipeline.R
```

> **Note:** conda-forge and bioconda may lag behind new Bioconductor releases.
> If `r-base=4.6.0` or a Bioconductor 3.23 package is unavailable,
> use Docker (Option B) instead.

### Update environment

```bash
conda env update -f environment.yml --prune
```

---

## E. Manual R install

Run the provided install script from within R:

```r
source("scripts/packages.R")
```

This installs all CRAN and Bioconductor packages automatically.

### Required versions

| Component | Version |
|-----------|---------|
| R | 4.6.0 |
| Bioconductor | 3.23 |
| ANCOMBC | 2.14.0 |
| lefser | 1.22.0 |
| ALDEx2 | 1.44.0 |

---

## External mapping files (Steps 6–7)

The functional enrichment steps require files not included in this repo:

| File | Source |
|------|--------|
| `map_metacyc-pwy_name.txt.gz` | [MetaCyc downloads](https://metacyc.org/downloads.shtml) → pathway names |
| `map_metacyc-pwy_lineage.tsv` | [MetaCyc downloads](https://metacyc.org/downloads.shtml) → pathway hierarchy |
| `ko_ec_cache.rds` | Generate via `KEGGREST`: `KEGGREST::keggLink("enzyme","ko")` |
| `ko_module_cache.rds` | Generate via `KEGGREST::keggLink("module","ko")` |
| `ko_pathway_cache.rds` | Generate via `KEGGREST::keggLink("pathway","ko")` |

Place these files in the project root before running Steps 6–7.
