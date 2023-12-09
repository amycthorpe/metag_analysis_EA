"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-09-23]
Run: snakemake -s workflow/rules/seed.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run SEED annotations based on SAMSA2
"""


############################################
rule seed:
    input:
        expand(os.path.join(RESULTS_DIR, "seed/{sid}/{sid}.seed.reduced_counts.tsv"), sid=SAMPLES.index)
    output:
        touch("status/seed.done")


############################################
localrules:


############################################
# SEED mapping to annotations
rule seed_diamond:
    input:
        fasta=os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta")
    output:
        daa=os.path.join(RESULTS_DIR, "seed/{sid}/{sid}_seed.daa"),
        tsv=os.path.join(RESULTS_DIR, "seed/{sid}/{sid}_seed.tsv")
    conda:
        os.path.join(ENV_DIR, "seed.yaml")
    threads:
        config["seed"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/seed/seed_diamond_{sid}.log")
    params:
        dmnd_db=config["seed"]["dmnd_db"],
        tmpdir=temp(os.path.join(RESULTS_DIR, "seed/tmp"))
    message:
        "Running SEED-mapper on {wildcards.sid}"
    shell:
        "(date && mkdir -p $(dirname {output}) && "
        "diamond blastx --db {params.dmnd_db} -q {input.fasta} -a {output.daa} -t {params.tmpdir} -k 1 --threads {threads} && "
        "diamond view --daa {output.daa} -o {output.tsv} -f tab --threads {threads} && "
        "date) &> >(tee {log})"

# SEED gene counts
rule gene_counts:
    input:
        rules.seed_diamond.output.tsv
    output:
        counts=os.path.join(RESULTS_DIR, "seed/{sid}/{sid}.seed.counts.tsv"),
        hits=os.path.join(RESULTS_DIR, "seed/{sid}/{sid}.seed.per_hit.tsv")
    conda:
        os.path.join(ENV_DIR, "seed_python.yaml")
    threads:
        config["seed"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/seed/{sid}.gene_counts.log")
    params:
        fa_db=config["seed"]["fa_db"],
        src=os.path.join(SRC_DIR, "DIAMOND_subsystems_analysis_counter.py")
    message:
        "Getting gene counts for SEED: {wildcards.sid}"
    shell:
        "(date && "
        "python2 {params.src} -I {input} -D {params.fa_db} -O {output.counts} -P {output.hits} && "
        "date) &> >(tee {log})"

rule reduce_counts:
    input:
        rules.gene_counts.output.counts
    output:
        reduced=os.path.join(RESULTS_DIR, "seed/{sid}/{sid}.seed.reduced_counts.tsv")
    conda:
        os.path.join(ENV_DIR, "seed_python.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/seed/{sid}.reduce_counts.log")
    params:
        src=os.path.join(SRC_DIR, "subsys_reducer.py")
    message:
        "Reduce identical subsystems annotations for: {wildcards.sid}"
    shell:
        "(date && "
        "python {params.src} -I {input} -O {output.reduced} && "
        "date) &> >(tee {log})"

    
