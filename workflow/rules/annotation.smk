"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-09-23]
Run: snakemake -s workflow/rules/annotation.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run gene calling via PRODIGAL on contigs
"""


############################################
rule annotation:
    input:
        expand(os.path.join(RESULTS_DIR, "prodigal/{sid}/{sid}.faa"), sid=SAMPLES)
    output:
        touch("status/annotation.done")


############################################
# localrules:


############################################
# Annotating the genes in assemblies
rule prodigal:
    input:
        FASTA=os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta")
    output:
        FAA=os.path.join(RESULTS_DIR, "prodigal/{sid}/{sid}.faa")
    conda:
        os.path.join(ENV_DIR, "prodigal.yaml")
    threads:
        config['prodigal']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/prodigal.{sid}.log")
    message:
        "Running PRODIGAL on {wildcards.sid}"
    shell:
        "(date && prodigal -a {output.FAA} -p meta -i {input.FASTA} && date) &> >(tee {log})"

