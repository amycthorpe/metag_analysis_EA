"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-11-20]
Run: snakemake -s workflow/rules/bin_prep.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To prepare assemblies for binning
"""


############################################
rule bin_prep:
    input:
        expand(os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}_filter.fasta.sa"), sid=SAMPLES.index)
    output:
        touch("status/bin_prep.done")


############################################
# localrules: phyloseq_input_kraken2


############################################
# Preparing the assemblies for binning 
rule filter_length:
    input:
        os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta")
    output:
        os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}_filter.fasta")
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    threads:
        config['filter_length']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/filter/length.{sid}.log")
    message:
        "Removing contigs less than 1.5Kb from {wildcards.sid}"
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    shell:
        "(date && seqkit seq -j {threads} -o {output} -m 1499 {input} && date) &> >(tee {log})"

rule bin_mapping_index:
    input:
        rules.filter_length.output
    output:
        os.path.join(RESULTS_DIR,"assembly/{sid}/{sid}_filter.fasta.sa")
    log:
        os.path.join(RESULTS_DIR, "logs/mapping/{sid}_mapping.bwa.index.log")
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    message:
        "Mapping: BWA index for assembly mapping for {wildcards.sid}"
    shell:
        "(date && bwa index {input} && date) &> >(tee {log})"    


rule filter_mapping:
    input:
        read1=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid_2}/{sid_2}_filtered.R1.fq"),
        read2=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid_2}/{sid_2}_filtered.R2.fq"),
        cont=rules.filter_length.output,
        idx=rules.bin_mapping_index.output
    output:
        temp(os.path.join(RESULTS_DIR,"bam/{sid}/{sid}_{sid_2}.bam"))
    threads:
        config["mapping"]["threads"]
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/mapping/{sid}_{sid_2}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/mapping/{sid}_{sid_2}.err.log")
    message:
        "Running bwa to produce sorted bams"
    shell:
        "(date && bwa mem -t {threads} {input.cont} {input.read1} {input.read2} | samtools sort -@{threads} -o {output} - && " 
        "samtools index {output} && date) 2> {log.err} > {log.out}"
