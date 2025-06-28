## Causal Inference for Latent Factor-Modeled Outcomes

This directory contains the code needed to reproduce all simulation studies and the data application from [the paper "Causal Inference for Latent Factor-Modeled Outcomes"](https://arxiv.org/abs/2506.20549). For the R package, see the [causalLFO](https://github.com/jennalandy/causalLFO) repository. 

We are unable to share the dataset or simulation parameters publicly, as they rely on private ICGC ARGO data, but we provide all processing code here that can be used by individuals with approval to access the private ICGC ARGO data. Publicly available mutational counts data for the 96-alphabet mutation classification can be accessed from ICGC ARGO (`WGS_PCAWG_2018_02_09/WGS_PCAWG.96.csv`, [instructions here](https://docs.icgc-argo.org/docs/data-access/icgc-25k-data#open-release-data---object-bucket-details)). If granted access to the private ICGC ARGO data, users can access the legacy PCAWG data through the SFTP server. The `germline_variations` subdirectory contains the results of germline mutation calling (`pcawg8.snps.indels.svs.phased.icgc.v2.controlled.vcf.gz`). The `donors_and_biospecimens` subdirectory contains age information as well as a mapping between tumor WGS ICGC specimen IDs and the normal tissue WGS aliquot IDs (`pcawg-data-releases.xlsx`).

## Directory Structure

- `visualization.R`: function and constants used in visualizations throughout this repository
- `methods.R`: functions to perform baseline and novel algorithms
- `application`:
    - Data
        - `application_data.qmd`: processes data for use in other application files. **This notebook must be run through before anything else in application is done.**
    -  Simulation parameters
        - `application_based_simulations_BA.qmd`: creates simulation parameters for simulations based on PCAWG breast-adenocarcinoma. **This notebook must be run through before anything in simulations is done.**
        - `application_based_simulations_BA_sub45.qmd`: chooses rank range for early onset breast-adenocarcinoma
        - figures from both are stored in `application_based_simulations_figures/` subdirectory
    - Application
        - `application.R`: R script to run algorithms on application data, utilized by slurm script `application.slurm`
        - `application.qmd`: Processes and analyzes results
        - `application_ref_heatmap.qmd`: cosine similarity between COSMIC reference signatures discovered in data application
        - `application_res.png`: application results figure used in paper
        - `application_sim.png`: cosine similairty of each algorithm's estimated signature to reference signature, used in paper appendix
- `simulations`:
    - Setup
        - `simulate.R`: functions to simulate datasets
        - `datasets.qmd`: simulate 100 datasets for each of $\pi = 0.2, 0.8$. **This notebook must be run through before anthing else in simulations is done.**
        - `define_realizations.R`: fixing 20 realizations of treatment assingment vectors for $\pi = 0.2, 0.8$ for use in indirect effects analyses, stores realizations in the `processed/` subdirectory. **This script must be run before anything in simulations/indirect effects is done.**
    - Methods
        - `indirect_effects.R`: functions to compute indirect effects / measures of interference
        - `direct_effects/`: directory containing slurm and R scripts to run direct effect simulation studies
        - `indirect_effects/`: directory containing slurm and R scripts to run indirect effect simulation studies
    - Results
        - `direct_effects.qmd`: create plots for direct effect results, places figures in `all_figures/` subdirectory
        - `indirect_effects.qmd`: create plots for indirect effects results, places figures in `all_figures/` subdirectory
        - `paper_figures.qmd`: creates final plots for use in the paper, places figures in `paper_figures/` subdirectory
