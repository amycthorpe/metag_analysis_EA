# Metagenome - Taxonomy and Assembly 

## About
Repository containing [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow for 
- running [Kraken2](https://ccb.jhu.edu/software/kraken/) taxonomy on preprocessed reads
- assembly with [Megahit](https://github.com/voutcn/megahit)

# Setup

## Databases
Downloaded dbs on `2023-09-13` as indicated below 
- from the [Kraken2 database page](https://benlangmead.github.io/aws-indexes/k2)
    - `/hdd0/susbus/databases/kraken2/pluspfp`: [pluspfp](https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20230605.tar.g)
    - `/hdd0/susbus/databases/kraken2/standard_20230913`: [standard](https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230605.tar.gz)

## BRACKEN
- Prior to running, install [bracken](https://github.com/jenniferlu717/Bracken)
- After installation, set *execute* permissions for the `bracken` & `combine_bracken_outputs.py` script
```bash
chmod +x PATH_to_INSTALL_DIR/Bracken/bracken
chmod +x PATH_to_INSTALL_DIR/Bracken/analyses_scripts/combine_bracken_outputs.py
```

## Conda

[Conda user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

```bash
# install miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod u+x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh # follow the instructions
```

Getting the repository including sub-modules
```bash
git clone git@github.com:UKCEH-MolecularEcology/metag_analyses.git
```

Create the main `snakemake` environment

```bash
# create venv
conda env create -f envs/snakemake.yaml -n "snakemake"
conda activate snakemake
```

## How to run

The workflow can be launched using one of the option as follows

```bash
snake_cores=48    # adjust as needed
snake_jobs=4   # adjust as needed
conda_prefix="/hdd0/susbus/tools/conda_env"

# dry-run
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --conda-prefix ${conda_prefix} --cores ${snake_cores} --jobs ${snake_jobs} --conda-frontend conda --rerun-trigger mtime -rpn  

# full run
snakemake -s workflow/Snakefile --configfile config/config.yaml --use-conda --conda-prefix ${conda_prefix} --cores ${snake_cores} --jobs ${snake_jobs} --conda-frontend conda --rerun-trigger mtime -rp
```


## Configs

All config files are stored in the folder `config/`

**Important Note(s)**: 
- Edit the paths to 
    - `data_dir`: `/ssd0/susbus/socd/data/preprocessed`
    - `results_dir`: `/ssd0/susbus/socd/results`
    - `env_dir`: `/ssd0/susbus/socd/kraken2/envs`
    - ***`databases`***: `/hdd0/susbus/databases/kraken2/pluspfp`

- Provide a `socd_sample_list` in the *config* folder
