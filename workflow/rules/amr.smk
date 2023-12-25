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
# localrules: 


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
        os.path.join(RESULTS_DIR, "logs/setup.rgi.setup.log")
    conda:
        os.path.join(ENV_DIR, "rgi.yaml")
    message:
        "Setup: load RGI DB"
    shell:
        "(date && rgi clean --local && "
        "rgi load --card_json {input} --local && rgi database --version --local && "
        "date) &> >(tee {log}) && touch {output}"







        os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}_modified.fasta")
    log:
        os.path.join(RESULTS_DIR, "logs/modify_fasta/{sid}.log")
    message:
        "Adding filename to fasta headers: {wildcards.sid}"
    shell:
        "(date && filename={wildcards.sid} && "
        """awk -v filename="$filename" '/^>/ {{ print $1"_"filename; next }} 1' {input} > {output} && """
        "date) &> >(tee {log})"

# Concatetnating the assemblies for binning 
rule ass_cat:
    input:
        expand(os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}_modified.fasta"), sid=SAMPLES.index)
    output:
        os.path.join(RESULTS_DIR, "assembly/cat_assembly.fasta")
    log:
        os.path.join(RESULTS_DIR, "logs/assembly/concatenation.log")
    message:
        "Concatenating all assemblies"
    shell:
        "(date && cat {input} > {output} && date) &> >(tee {log})"

# Clustering the assemblies for binning
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
        "mmseqs easy-cluster --threads {threads} --force-reuse 0 --cov-mode {params.cov_mode} --min-seq-id {params.min_id} -c {params.cov} {input} $(echo {output} | sed 's,_rep_seq.fasta,,g') {params.tmp} && "
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
        "Running bwa to produce sorted bams: {wildcards.sid}"
    shell:
        "(date && bwa mem -t {threads} {input.cont} {input.read1} {input.read2} | samtools sort -@{threads} -o {output} - && " 
        "samtools index {output} && date) 2> {log.err} > {log.out}"
