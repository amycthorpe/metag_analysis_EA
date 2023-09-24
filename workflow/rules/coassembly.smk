"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-08-30]
Run: snakemake -s workflow/rules/coassembly.smk --use-conda --cores 72 -rp
Latest modification:
Purpose: To run MEGAHIT assembly on all reads after deduplication
"""


############################################
rule coassembly:
    input:
        os.path.join(RESULTS_DIR, "coassembly/assembly/COASSEMBLY.fasta")
    output:
        touch("status/coassembly.done")


############################################
# localrules: 


############################################
# Concatenating the reads
rule concat:
    input:
        sr1=expand(os.path.join(DATA_DIR, "{sid}_R1.fq.gz"), sid=SAMPLES),
        sr2=expand(os.path.join(DATA_DIR, "{sid}_R2.fq.gz"), sid=SAMPLES)
    output:
        or1=os.path.join(RESULTS_DIR, "coassembly/reads/merged_R1.fq.gz"),
        or2=os.path.join(RESULTS_DIR, "coassembly/reads/merged_R2.fq.gz")
    log:
        os.path.join(RESULTS_DIR, "logs/concat_reads.log")
    message:
        "Concatenating all reads for coassembly"
    shell:
        "(date && cat {input.sr1} > {output.or1} && "
        "cat {input.sr2} > {output.or2} && date) &> >(tee {log})"

rule coassembly_dedup:
    input:
        r1=rules.concat.output.or1,
        r2=rules.concat.output.or2
    output:
        odup1=os.path.join(RESULTS_DIR, "coassembly/dedup/merged_R1.fq.gz"),
        odup2=os.path.join(RESULTS_DIR, "coassembly/dedup/merged_R2.fq.gz")
    log:
        os.path.join(RESULTS_DIR, "logs/deduplicate_coassembly_reads.log")
    threads:
        config["clumpify"]["threads"]
    conda:
        os.path.join(ENV_DIR, "bbmap.yaml")
    message:
        "Removing duplicate reads for easier downstream assembly"
    shell:
        "(date && clumpify.sh in={input.r1} in2={input.r2} out={output.odup1} out2={output.odup2} dupedist={config[clumpify][dupedist]} dedupe=t optical=t threads={threads} groups={config[clumpify][groups]} -Xmx{config[clumpify][memory]} && date) &> >(tee {log})"

rule coassembly_megahit:
    input:
        sr1=rules.coassembly_dedup.output.odup1,
        sr2=rules.coassembly_dedup.output.odup2
    output:
        os.path.join(RESULTS_DIR, "coassembly/assembly/COASSEMBLY.fasta")
    log:
        os.path.join(RESULTS_DIR, "logs/coassembly_megahit.log")
    threads:
        config["megahit"]["coassembly_threads"]
    conda:
        os.path.join(ENV_DIR, "megahit.yaml")
    message:
        "Coassemblye all the reads: MEGAHIT"
    shell:
        "(date && megahit -1 {input.sr1} -2 {input.sr2} --kmin-1pass -m 0.9 --k-list 27,37,47,57,67,77,87 --min-contig-len 1000 -t {threads} -o $(dirname {output})/tmp && "
        "cd $(dirname {output}) && "
        "rsync -avP tmp/ . && "
        "ln -sf final.contigs.fa $(basename {output}) && "
        "rm -rf tmp/ && "
        "date) &> >(tee {log})"

