"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-08-30]
Run: snakemake -s workflow/rules/assembly.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run MEGAHIT assembly on reads
"""


############################################
rule assembly:
    input:
        expand(os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta"), sid=SAMPLES.index)
    output:
        touch("status/assembly.done")


############################################
# localrules: phyloseq_input_kraken2


############################################
# Assembling the reads
rule megahit:
    input:
        sr1=os.path.join(RESULTS_DIR, "preprocessed/trimmed/{sid}/{sid}_val_1.fq.gz"),
        sr2=os.path.join(RESULTS_DIR, "preprocessed/trimmed/{sid}/{sid}_val_2.fq.gz")
    output:
        os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta")
    conda:
        os.path.join(ENV_DIR, "megahit.yaml")
    threads:
        config['megahit']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/megahit.{sid}.log")
    message:
        "Running MEGAHIT on {wildcards.sid}"
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    shell:
        "(date && megahit -1 {input.sr1} -2 {input.sr2} --kmin-1pass -m 0.9 --k-list 27,37,47,57,67,77,87 --min-contig-len 1000 -t {threads} -o $(dirname {output})/tmp && "
        "cd $(dirname {output}) && "
        "rsync -avP tmp/ . && "
        "ln -sf final.contigs.fa $(basename {output}) && "
        "rm -rf tmp/ && "
        "date) &> >(tee {log})"

