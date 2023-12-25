"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-12-25]
Run: snakemake -s workflow/rules/amr.smk --use-conda --cores 64 -rp
Latest modification:
Purpose: To run AMR on concatenated assembly
"""


############################################
rule amr:
    input:
        os.path.join(RESULTS_DIR, "amr/rgi.txt")
    output:
        touch("status/amr.done")


############################################
localrules: download_rgi_db, setup_rgi_db


############################################
# Download RGI data
rule download_rgi_db:
    output:
        archive=temp(os.path.join(DB_DIR, "rgi/card-data.tar.bz2")),
        json=os.path.join(DB_DIR, "rgi/card.json")
    log:
        os.path.join(RESULTS_DIR, "logs/setup.rgi.db.log")
    params:
        db_url=config["rgi"]["db_url"]
    message:
        "Setup: download RGI data"
    shell:
        "(date && "
        "wget -O {output.archive} {params.db_url} --no-check-certificate && "
        "tar -C $(dirname {output.archive}) -xvf {output.archive} && "
        "date) &> >(tee {log})"

# Setup RGI: load required DB
# NOTE: to make sure that the same DB is used for all targets
rule setup_rgi_db:
    input:
        os.path.join(DB_DIR, "rgi/card.json")
    output:
        "status/rgi_setup.done"
    log:
        os.path.join(RESULTS_DIR, "logs/rgi.setup.log")
    conda:
        os.path.join(ENV_DIR, "rgi.yaml")
    message:
        "Setup: load RGI DB"
    shell:
        "(date && "
        "rgi clean --local && "
        "rgi load --card_json {input} --local && "
        "rgi database --version --local && "
        "touch {output} && date) &> >(tee {log})"

# Run RGI: Assembly (DNA)
rule annotation_rgi:
    input:
        fna=os.path.join(RESULTS_DIR, "assembly/cat_assembly_filter.fasta"),
        db=os.path.join(DB_DIR, "rgi/card.json"),
        setup="status/rgi_setup.done" # NOTE: to make sure that the same DB is used for all targets
    output:
        txt=os.path.join(RESULTS_DIR, "amr/rgi.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/rgi.annotation.log")
    threads:
        config["rgi"]["threads"]
    params:
        alignment_tool="DIAMOND"
    conda:
        os.path.join(ENV_DIR, "rgi.yaml")
    message:
        "Annotation: RGI: mmseqs_rep_seq_DNA"
    shell:
        "(date && "
        "rgi database --version --local && "
        # NOTE: https://github.com/arpcard/rgi/issues/93: KeyError: 'snp' --> re-run
        "rgi main --input_sequence {input.fna} --output_file {output.txt} --local -a {params.alignment_tool} --clean --low_quality -n {threads} || "
        "rgi main --input_sequence {input.fna} --output_file {output.txt} --local -a {params.alignment_tool} --clean --low_quality -n {threads} && "
        "date) &> >(tee {log})"

