# Alopecia FMT Microbiome Pipeline
# Base image: official Bioconductor Docker (R 4.6.0 + Bioconductor 3.23)
#
# Build:   docker build -t alopecia-fmt .
# Run:     docker run --rm -v $(pwd)/sim_inputs:/project/sim_inputs \
#                         -v $(pwd)/outputs:/project/outputs \
#                         alopecia-fmt
# Shell:   docker run --rm -it alopecia-fmt bash

FROM bioconductor/bioconductor_docker:RELEASE_3_23

LABEL description="Alopecia oral FMT — 16S + pathway microbiome analysis pipeline" \
      r.version="4.6.0" \
      bioconductor.version="3.23"

WORKDIR /project

# ── CRAN packages ─────────────────────────────────────────────────────────────
RUN R -e "install.packages(c( \
    'vegan', 'ggrepel', 'patchwork', 'lme4', 'lmerTest', \
    'emmeans', 'broom', 'RColorBrewer', 'scales', 'readxl', \
    'Cairo', 'rstatix' \
  ), repos='https://cloud.r-project.org', Ncpus=4)"

# ── Bioconductor packages ─────────────────────────────────────────────────────
RUN R -e "BiocManager::install(c( \
    'ANCOMBC', 'lefser', 'ALDEx2', 'Maaslin2', \
    'fgsea', 'KEGGREST', 'ssizeRNA' \
  ), ask=FALSE, Ncpus=4)"

# ── Copy project files ────────────────────────────────────────────────────────
COPY scripts/ scripts/
COPY sim_inputs/ sim_inputs/

# External mapping files (MetaCyc + KEGG) must be mounted at runtime if needed:
#   docker run ... -v /path/to/map_metacyc-pwy_name.txt.gz:/project/map_metacyc-pwy_name.txt.gz

CMD ["Rscript", "scripts/run_pipeline.R"]
