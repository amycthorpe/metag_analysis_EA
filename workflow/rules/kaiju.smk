"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-09-23]
Run: snakemake -s workflow/rules/kaiju.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To classify eukaryotes from contigs using KAIJU
"""


############################################
rule kaiju:
    input:
        os.path.join(RESULTS_DIR, "kaiju/kaiju_names.txt")
    output:
        touch("status/kaiju.done")


############################################
# localrules:


############################################
# Preparing the files for KAIJU classification
rule kaiju_classify:
    input:
        FASTA=",".join([os.path.join(RESULTS_DIR, f"assembly/{sid}/{sid}.fasta") for sid in SAMPLES]) # FASTA=os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta")
    output:
        out=os.path.join(RESULTS_DIR, "kaiju/kaiju_out.txt"),
        names=os.path.join(RESULTS_DIR, "kaiju/kaiju_names.txt")
    conda:
        os.path.join(ENV_DIR, "kaiju.yaml")
    threads:
        config["kaiju"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/kaiju.log")
    params:
        fmi=config["kaiju"]["fmi"],
        nodes=config["kaiju"]["nodes"],
        names=config["kaiju"]["names"],
    message:
        "Running KAIJU against the nr_protein database"
    shell:
        "(date && "
        "kaiju-multi -v -z {threads} -t {params.nodes} -f {params.fmi} -i {input} > {output.out} && "
        "kaiju-addTaxonNames -t {params.nodes}  -n {params.names} -i {output.out} -o {output.names} && "
        "date) &> >(tee {log})"
