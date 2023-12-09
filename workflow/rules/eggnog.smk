"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-09-23]
Run: snakemake -s workflow/rules/eggnog.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run eggnog on proteins
"""


############################################
rule eggnog:
    input:
        expand(os.path.join(RESULTS_DIR, "eggnog/{sid}/{sid}.emapper.annotations"), sid=SAMPLES.index)
    output:
        touch("status/eggnog.done")


############################################
localrules: download_eggnogDB


############################################
# EggNOG database installation
rule download_eggnogDB:
    output:
        done=os.path.join(RESULTS_DIR, "eggnog/db_download.done")
    log:
        out=os.path.join(RESULTS_DIR, "logs/setup.eggnog_DB.log")
    params:
        path=config["eggnog"]["db"]
    conda:
        os.path.join(ENV_DIR, "eggnog.yaml")
    message:
        "Setup: EggNOG database"
    shell:
        "(date && mkdir -p {params.path} && "
        "download_eggnog_data.py -y --data_dir {params.path} && "
        "touch {output} && date) &> >(tee {log})"

# EGGNOG mapping to annotations
rule emapper:
    input:
        dummy=os.path.join(RESULTS_DIR, "eggnog/db_download.done"),
        fasta=os.path.join(RESULTS_DIR, "prodigal/{sid}/{sid}.faa")
    output:
        os.path.join(RESULTS_DIR, "eggnog/{sid}/{sid}_emapper.seed_orthologs")
    conda:
        os.path.join(ENV_DIR, "eggnog.yaml")
    threads:
        config["eggnog"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/eggnog/emapper_{sid}.log")
    params:
        itype=config["eggnog"]["itype"],
        genepred=config["eggnog"]["genepred"],
        db=config["eggnog"]["db"]
    message:
        "Running EggNog-mapper on {wildcards.sid}"
    shell:
        "(date && mkdir -p {output} && "
        "emapper.py -m diamond --data_dir {params.db} --itype {params.itype} --no_file_comments --cpu {threads} -i {input} -o {wildcards.sid} --output_dir $(dirname {output}) && "
        "date) &> >(tee {log})"

# Final annotations
rule emapper_final:
    input:
        rules.emapper.output[0]
    output:
        os.path.join(RESULTS_DIR, "eggnog/{sid}/{sid}.emapper.annotations")
    conda:
        os.path.join(ENV_DIR, "eggnog.yaml")
    threads:
        config["eggnog"]["final_threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/eggnog/{sid}.final_emapper.out.log")
    params:
        itype=config["eggnog"]["itype"],
        genepred=config["eggnog"]["genepred"],
        db=config["eggnog"]["db"]
    message:
        "Running EggNog-mapper annotations on {wildcards.sid}"
    shell:
        "(date && "
        "emapper.py --data_dir {params.db} --annotate_hits_table {input} --no_file_comments -o $(echo {output} | sed 's/.emapper.annotations//g' ) --cpu {threads} --dbmem && "
        "date) &> >(tee {log})"
