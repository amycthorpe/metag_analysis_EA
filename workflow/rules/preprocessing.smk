"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-10-15]
Run: snakemake -s workflow/rules/preprocessing.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To trim and filter reads
"""


############################################
# params:
samples=list(SAMPLES.index)


############################################
rule preprocessing:
    input:
        expand(os.path.join(RESULTS_DIR, "preprocessed/trimmed/{sid}/{sid}_val_{rid}.fq.gz"), sid=SAMPLES.index, rid=["1", "2"]),
        expand(os.path.join(RESULTS_DIR, "preprocessed/fastqc/{sid}/{sid}_{rid}_fastqc.zip"), sid=SAMPLES.index, rid=["R1", "R2"]),
        os.path.join(RESULTS_DIR, "preprocessed/multiqc/fastqc/multiqc_report.html")
    output:
        touch("status/preprocessing.done")


############################################
localrules:


############################################
# Trimming raw fastq reads
rule trim_galore_pe:
    input:
        [lambda wildcards: SAMPLES.loc[wildcards.sid, "sR1"], lambda wildcards: SAMPLES.loc[wildcards.sid, "sR2"]],
    output:
        fasta_fwd=os.path.join(RESULTS_DIR, "preprocessed/trimmed/{sid}/{sid}_val_1.fq.gz"),
#        report_fwd=os.path.join(RESULTS_DIR, "preprocessed/trimmed/reports/{sid}_1.fq.gz_trimming_report.txt"),
        fasta_rev=os.path.join(RESULTS_DIR, "preprocessed/trimmed/{sid}/{sid}_val_2.fq.gz"),
#        report_rev=os.path.join(RESULTS_DIR, "preprocessed/trimmed/reports/{sid}_2.fq.gz_trimming_report.txt")
    threads:
        config["trim_galore"]["threads"]
    params:
        extra="--illumina -q 25",
        debug=lambda wildcards: print(f"Wildcards for rule trim_galore_pe: {wildcards}")
    log:
        os.path.join(RESULTS_DIR, "logs/trim_galore/{sid}.log")
    conda:
        os.path.join(ENV_DIR, "trim_galore.yaml")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Trimming paired end reads for {wildcards.sid}"
    shell:
        "(date && "
         "trim_galore -j {threads} {params.extra} --basename {wildcards.sid} -o $(dirname {output.fasta_fwd}) --paired {input} &&"
        "date) &> {log}"

#    wrapper:
#        "v2.13.0/bio/trim_galore/pe"

# Indexing fasta file to be filtered from
rule index:
    input:
        config["filter"]["fasta"]
    output:
        expand(os.path.join(RESULTS_DIR, "preprocessed/filter/GRCh38_latest_genomic.fna.gz.{ext}"), ext=BWA_IDX_EXT)
    threads:
        config["bwa"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/filter/filter.index.log")
    params:
        algorithm="bwtsw"
    message:
        "Indexing the filtering fasta file"
    wrapper:
        "v2.6.0/bio/bwa/index"

# Mapping raw reads to filter
rule map_to_mask:
    input:
        r1=rules.trim_galore_pe.output.fasta_fwd,
        r2=rules.trim_galore_pe.output.fasta_rev,
        index=os.path.join(RESULTS_DIR, "preprocessed/filter/GRCh38_latest_genomic.fna.gz.sa")
    output:
        bam=temp(os.path.join(RESULTS_DIR, "preprocessed/bam/{sid}_filtered.bam"))
    threads:
        config["bwa"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/bam/{sid}.log")
    conda:
        os.path.join(ENV_DIR, "bwa.yaml")
    message:
        "Mapping trimmed {wildcards.sid} reads to filter"
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    shell:
        "(date && "
        "bwa mem -t {threads} {input.index} {input.r1} {input.r2} | samtools view -b -f 12 -@{threads} - | samtools sort -@{threads} - > {output.bam} && "
        "date) &> {log}"

# Convert bam to paired-end fq
rule bam_to_fastq:
    input:
        rules.map_to_mask.output.bam
    output:
        r1=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R1.fq.gz"),
        r2=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R2.fq.gz")
    threads:
        config["bwa"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/bam_to_fastq/{sid}.log")
    conda:
        os.path.join(ENV_DIR, "bedtools.yaml")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Getting fastq files from BAM for {wildcards.sid}"
    shell:
        "(date && bamToFastq -i {input} -fq {output.r1} -fq2 {output.r2} && date) &> {log}"

# Running QC
rule fastqc:
    input:
        os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.{rid}.fq.gz")
    output:
        zip=os.path.join(RESULTS_DIR, "preprocessed/fastqc/{sid}/{sid}_{rid}_fastqc.zip"),
        html=os.path.join(RESULTS_DIR, "preprocessed/fastqc/{sid}/{sid}_{rid}_fastqc.html")
    log:
        os.path.join(RESULTS_DIR, "logs/fastqc/{sid}_{rid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index),
        rid="R1|R2"
    threads:
        config["fastqc"]["threads"]
    conda:
        os.path.join(ENV_DIR, "fastqc.yaml")
    message:
        "FastQC: {wildcards.sid}"
    shell:
        "fastqc -q -f fastq -t {threads} -o $(dirname {output.zip}) {input} &> {log}"

# Collating QC results
#rule multiqc_fastqc:
#    input:
##         lambda wildcards: expand(os.path.join(RESULTS_DIR, "preprocessed/fastqc/{{sid}}/{sid}_{rid}_fastqc.zip"), sid="|".join(SAMPLES.index), rid=["R1", "R2"])
#        expand(os.path.join(RESULTS_DIR, "preprocessed/fastqc/{{sid}}/{sid}_{rid}_fastqc.zip"), sid=list(SAMPLES.index), rid=["R1", "R2"])
#    output:
#        html=os.path.join(RESULTS_DIR, "preprocessed/multiqc/fastqc/multiqc_report.html"),
#        stat=os.path.join(RESULTS_DIR, "preprocessed/multiqc/fastqc/multiqc_data/multiqc_fastqc.txt"),
#    log:
#        os.path.join(RESULTS_DIR, "multiqc/fastqc/multiqc.log")
#    threads:
#        config["fastqc"]["threads"]
#    conda:
#        os.path.join(ENV_DIR, "multiqc.yaml")
#    message:
#        "MultiQC (FastQC)"
#    wildcard_constraints:
#        sid="|".join(SAMPLES.index)
#    shell:
#        "multiqc --interactive -p -f -m fastqc -o $(dirname {output.html}) $(dirname {input[0]}) &> {log}"

rule multiqc_fastqc:
    input:
        lambda wildcards: expand(os.path.join(RESULTS_DIR, "preprocessed/fastqc/{sid}/{sid}_{rid}_fastqc.zip"), sid=samples, rid=["R1", "R2"])
    output:
        html=os.path.join(RESULTS_DIR, "preprocessed/multiqc/fastqc/multiqc_report.html"),
        stat=os.path.join(RESULTS_DIR, "preprocessed/multiqc/fastqc/multiqc_data/multiqc_fastqc.txt"),
    log:
        os.path.join(RESULTS_DIR, "multiqc/fastqc/multiqc.log")
    threads:
        config["fastqc"]["threads"]
    conda:
        os.path.join(ENV_DIR, "multiqc.yaml")
    wildcard_constraints:
        sid="|".join(samples)
    message:
        "MultiQC (FastQC)"
    shell:
        "multiqc --interactive -p -f -m fastqc -o $(dirname {output.html}) $(dirname {input[0]}) &> {log}"

