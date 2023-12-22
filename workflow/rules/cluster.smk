"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-11-20]
Run: snakemake -s workflow/rules/cluster.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To prepare a concatenated assembly for binning
"""


############################################
rule bin_prep:
    input:
        expand(os.path.join(RESULTS_DIR, "assembly/cat_assembly_filter.fasta.sa"), sid=SAMPLES.index),
        expand(os.path.join(RESULTS_DIR,"bam/{sid}/cat_assembly_{sid}.bam"), sid=SAMPLES.index)
    output:
        touch("status/cluster.done")


############################################
# localrules: phyloseq_input_kraken2


############################################
# Clustering the assemblies for binning 
rule ass_cat:
    input:
        expand(os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta"), sid=SAMPLES.index)
    output:
        os.path.join(RESULTS_DIR, "assembly/cat_assembly.fasta")
    log:
        os.path.join(RESULTS_DIR, "logs/assembly/concatenation.log")
    message:
        "Concatenating all assemblies"
    shell:
        "(date && cat {input} > {output} && date) &> >(tee {log})"

rule cat_mmseqs2:
    input:
        rules.ass_cat.output
    output:
        os.path.join(RESULTS_DIR, "mmseqs/cat_assembly_rep_seq.fasta")
    log:
        os.path.join(RESULTS_DIR, "logs/mmseqs_cluster.log")
    conda:
        os.path.join(ENV_DIR, "mmseqs2.yaml")
    threads:
        config["mmseqs2"]["threads"]
    params:
        tmp=os.path.join(RESULTS_DIR, "tmpdir"),
        cov=config["mmseqs2"]["cov"],
        cov_mode=config["mmseqs2"]["cov_mode"],
        min_id=config["mmseqs2"]["min_id"]
    shell:
        "(date && "
        "mmseqs easy-cluster --threads {threads} --force-reuse 0 --cov-mode {params.cov_mode} --min-seq-id {params.min_id} -c {params.cov} {input} $(echo $(basename {output} | sed 's,_rep_seq.fasta,,g')) {params.tmp} && "
        "date) &> >(tee {log})"

rule cat_filter_length:
    input:
        rules.cat_mmseqs2.output
    output:
        os.path.join(RESULTS_DIR, "assembly/cat_assembly_filter.fasta")
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    threads:
        config['filter_length']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/filter/length.log")
    message:
        "Removing contigs less than 1.5K"
    shell:
        "(date && seqkit seq -j {threads} -o {output} -m 1499 {input} && date) &> >(tee {log})"

rule cat_bin_mapping_index:
    input:
        rules.cat_filter_length.output
    output:
        os.path.join(RESULTS_DIR,"assembly/cat_assembly_filter.fasta.sa")
    log:
        os.path.join(RESULTS_DIR, "logs/mapping/mapping.bwa.index.log")
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    message:
        "Mapping: BWA index for assembly mapping"
    shell:
        "(date && bwa index {input} && date) &> >(tee {log})"    

rule filter_mapping:
    input:
        read1=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R1.fq"),
        read2=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R2.fq"),
        cont=rules.cat_filter_length.output,
        idx=rules.cat_bin_mapping_index.output
    output:
        os.path.join(RESULTS_DIR,"bam/{sid}/cat_assembly_{sid}.bam")
    threads:
        config["mapping"]["threads"]
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/mapping/cat_assembly_{sid}.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/mapping/cat_assembly_{sid}.err.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Running bwa to produce sorted bams"
    shell:
        "(date && bwa mem -t {threads} {input.cont} {input.read1} {input.read2} | samtools sort -@{threads} -o {output} - && " 
        "samtools index {output} && date) 2> {log.err} > {log.out}"
