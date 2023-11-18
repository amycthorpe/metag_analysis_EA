"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-10-09]
Run: snakemake -s workflow/rules/coverage.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To estimate contig and gene coverage
"""


############################################
rule coverage:
    input:
        expand(os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}_depth.txt"), sid=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}.gene.len"), sid=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}_gene_coverage.txt"), sid=SAMPLES)
    output:
        touch("status/coverage.done")


#############################################
# Params
BWA_IDX_EXT = ["amb", "ann", "bwt", "pac", "sa"]  


############################################
localrules: gene_depth, contig_length, contig_gene_link


############################################
# Mapping reads to contigs
rule mapping_index:
    input:
        FASTA=os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta")
    output:
        expand(os.path.join(RESULTS_DIR, "bam/{{sid}}/{{sid}}.{ext}"), ext=BWA_IDX_EXT)
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    threads:
        config['bwa']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/index.{sid}.log")
    params:
        idx_prefix=lambda wildcards, output: os.path.splitext(output[0])[0]
    message:
        "Indexing assembly from {wildcards.sid}"
    shell:
        "(date && bwa index {input} -p {params.idx_prefix} && date) &> {log}"

rule mapping:
    input:
        asm=os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta"),
        sr1=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R1.fq"),
        sr2=os.path.join(RESULTS_DIR, "preprocessed/reads/{sid}/{sid}_filtered.R2.fq"),
        idx=expand(os.path.join(RESULTS_DIR, "bam/{{sid}}/{{sid}}.{ext}"), ext=BWA_IDX_EXT)
    output:
        os.path.join(RESULTS_DIR, "bam/{sid}/{sid}.sorted.bam")
    conda:
        os.path.join(ENV_DIR, "mapping.yaml")
    threads:
         config['bwa']['map_threads']
    log:
        os.path.join(RESULTS_DIR, "logs/mapping.{sid}.log")
    params:
        idx_prefix=lambda wildcards, input: os.path.splitext(input.idx[0])[0],
        bam_prefix=lambda wildcards, output: os.path.splitext(output[0])[0],
        chunk_size=config["samtools"]["sort"]["chunk_size"]
    message:
        "Mapping short reads to {wildcards.sid} assembly"
    shell:
        "(date && "
        "bwa mem -t {threads} {params.idx_prefix} {input.sr1} {input.sr2} | "
        "samtools view -@ {threads} -SbT {input.asm} | "
        "samtools sort -@ {threads} -m {params.chunk_size} -T {params.bam_prefix} -o {output} && "
        "date) &> {log}"

rule summarise_depth:
    input:
        rules.mapping.output[0]
    output:
        depth=os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}_depth.txt"),
        paired=os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}_paired.txt")
    conda:
        os.path.join(ENV_DIR, "metabat2.yaml")
    threads:
        config['bwa']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/depth.{sid}.log")
    message:
        "Getting coverage for {wildcards.sid}"
    shell:
        "(date && "
        "jgi_summarize_bam_contig_depths --outputDepth {output.depth} --pairedContigs {output.paired} {input} && "
        "date) &> {log}"

rule contig_length:
    input:
        asm=os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta")
    output:
        length=os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}.length.txt")
    conda:
        os.path.join(ENV_DIR, "perl.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/length.{sid}.log")
    params:
        calc=os.path.join(SRC_DIR, "fastaNamesSizes.pl")
    message:
        "Getting the contig lengths for {wildcards.sid}"
    shell:
        "(date && perl {params.calc} {input.asm} > {output.length} && date ) &> {log}"

rule gene_depth:
    input:
        gff=os.path.join(RESULTS_DIR, "prodigal/{sid}/{sid}.gff"),
        BAM=os.path.join(RESULTS_DIR, "bam/{sid}/{sid}.sorted.bam"),
        length=rules.contig_length.output.length
    output:
        hist=os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}.gene_depth.hist"),
        average=os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}.gene_depth.avg"), 
        len=os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}.gene.len")
    conda:
        os.path.join(ENV_DIR, "perl.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/gene_depth.{sid}.log")
    params:
        tmp=os.path.join(RESULTS_DIR, "coverage/{sid}/tmp"),
        avg_calc=os.path.join(SRC_DIR, "calcAvgCoverage.pl")
    message:
        "Estimating gene length and depth for {wildcards.sid}"
    shell:
        """
        TMP_FILE=$(mktemp --tmpdir={params.tmp} -t "annotation_cov_XXXXXX.bed")
        sortBed -g {input.length} -i {input.gff} | awk '$4 < $5' | coverageBed -hist -sorted -g {input.length} -b {input.BAM} -a "stdin" | grep -v "^all" > $TMP_FILE

        paste <(cat $TMP_FILE | cut -f9 | cut -f1 -d \";\" | sed -e \"s/ID=//g\") \
         <(cut -f10,11,12,13 $TMP_FILE) > {output.hist}
    
        rm $TMP_FILE
        # Record average gene coverage
        {params.avg_calc} {output.hist} > {output.average}

        # Record gene length file
        cut -f 1,4 {output.hist} | uniq > {output.len}
        """   

rule contig_gene_link:
    input:
        faa=os.path.join(RESULTS_DIR, "prodigal/{sid}/{sid}.faa"),
        average=rules.gene_depth.output.average
    output:
        txt=os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}_gene_contig.txt"),
        gene_cov=os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}_gene_coverage.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/gene_contig_{sid}.log")
    message:
        "Linking genes to contig IDs from FAA: {wildcards.sid}"
    shell:
        "(date && grep '>' {input.faa} | cut -f1,9 -d ' ' | sed 's/>//g' | cut -f1 -d ';' | sed -e 's/ID=//g' > {output.txt} && "
        "join -1 2 -2 1 <(sort -k2,2 {output.txt}) <(sort -k1,1 {input.average}) > {output.gene_cov} && " 
        "date) &> {log}"
