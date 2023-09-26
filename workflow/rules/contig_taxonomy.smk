"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-09-25]
Run: snakemake -s workflow/rules/contig_taxonomy.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run Kraken2+BRACKEN on metagenome assemblies, i.e. contigs
"""


############################################
rule contig_taxonomy:
    input:
        os.path.join(RESULTS_DIR, "mpa_report/contig/combined_output.tsv"),
        os.path.join(RESULTS_DIR, "bracken/contig/combined_bracken.txt")
    output:
        touch("status/contig_taxonomy.done")


############################################
localrules: phyloseq_input_kraken2


############################################
# Taxonomic classification using KRAKEN2
rule contig_kraken2:
    input:
        os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta"),
    output:
        report=os.path.join(RESULTS_DIR, "kraken2/contig/{sid}_kraken.report"),
        summary=os.path.join(RESULTS_DIR, "kraken2/contig/{sid}_kraken.out")
    conda:
        os.path.join(ENV_DIR, "kraken2.yaml")
    threads:
        config['kraken2']['threads']
    params:
        db=config['kraken2']['db'],
        confidence=config['kraken2']['contig_confidence']
    log:
        os.path.join(RESULTS_DIR, "logs/contig/kraken2.{sid}.log")
    message:
        "Running kraken2 on {wildcards.sid}"
    shell:
        "(date && kraken2 --threads {threads} --db {params.db} --confidence {params.confidence} --output {output.summary} --report {output.report} {input} && date) &> >(tee {log})"

# Running Struo2 database
use rule contig_kraken2 as contig_struo2_kraken2 with:
    output:
        report=os.path.join(RESULTS_DIR, "kraken2/contig/{struo2_{sid}_kraken.report"),
        summary=os.path.join(RESULTS_DIR, "kraken2/contig/struo2_{sid}_kraken.out")
    threads:
        config['struo2_kraken2']['threads']
    params:
        db=config['struo2_kraken2']['db'],
        confidence=config['kraken2']['confidence']
    log:
        os.path.join(RESULTS_DIR, "logs/contig/struo2_kraken2.{sid}.log")
    message:
        "Running {wildcards.sid} with the Struo2_Kraken2 db"

# Running KRAKEN2+BRACKEN as suggested by KRAKEN2 website
rule contig_bracken:
    input:
        report=rules.contig_kraken2.output.report,
        r1=os.path.join(DATA_DIR, "{sid}_R1.fq.gz"),  # os.path.join(DATA_DIR, "00.RawData/{sid}/{sid}_R1.fastq.gz"),
        r2=os.path.join(DATA_DIR, "{sid}_R2.fq.gz") # os.path.join(DATA_DIR, "00.RawData/{sid}/{sid}_R2.fastq.gz")
    output:
        bracken=os.path.join(RESULTS_DIR, "bracken/contig/{sid}.bracken"),
        report=os.path.join(RESULTS_DIR, "bracken/contig/{sid}_bracken.report")
    threads:
        config['kraken2']['threads']
    conda:
        os.path.join(ENV_DIR, "bracken.yaml")
    params:
        db=config['kraken2']['db'],
        read=config['kraken2']['read'],
        level=config['kraken2']['contig_level'],
        bracken=config['bracken']['bin']
    log:
        os.path.join(RESULTS_DIR, "logs/contig/bracken.{sid}.log")
    message:
        "Running kraken & bracken for {wildcards.sid}"
    shell:
        "(date && {params.bracken} -d {params.db} -i {input.report} -o {output.bracken} -w {output.report} -r {params.read} -l {params.level} && date)  &> >(tee {log})"

rule contig_remove_uncultured:
    input:
        bracken=os.path.join(RESULTS_DIR, "bracken/contig/{sid}.bracken")
    output:
        edited=os.path.join(RESULTS_DIR, "bracken/contig/{sid}_edited.bracken")
    log:
        os.path.join(RESULTS_DIR, "logs/contig/edited_bracken_{sid}")
    message:
        "Removing 'uncultured' taxa from bracken output from {wildcards.sid} due to combining issues"
    shell:
        "(date && grep -v 'uncultured' {input.bracken} | grep -v 'endosymbionts' | grep -v 'Incertae Sedis' > {output.edited} && date) &> >(tee {log})"

rule contig_combine_bracken:
    input:
        bracken=expand(os.path.join(RESULTS_DIR, "bracken/contig/{sid}_edited.bracken"), sid=SAMPLES)
    output:
        out=os.path.join(RESULTS_DIR, "bracken/contig/combined_bracken.txt")
    conda:
        os.path.join(ENV_DIR, "python2.yaml")
    params:
        combine=config['bracken']['combine']
    log:
        os.path.join(RESULTS_DIR, "logs/contig/bracken_combine.log")
    message:
        "Combining all the output from BRACKEN"
    shell:
        "(date && python {params.combine} --files {input.bracken} -o {output.out} && date)  &> >(tee {log})"


#########################
### MPA-style report ###
rule contig_mpa_report:
    input:
        report=os.path.join(RESULTS_DIR, "bracken/contig/{sid}_bracken.report")
    output:
        mpa=os.path.join(RESULTS_DIR, "mpa_report/contig/{sid}_mpa.tsv")
    conda:
        os.path.join(ENV_DIR, "bracken_new.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/contig/mpa_{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES)
    message:
        "Creating mpa-style report for {wildcards.sid}"
    shell:
        "(date && kreport2mpa.py -r {input.report} -o {output.mpa} && date)  &> >(tee {log})"

rule contig_combine_mpa:
    input:
        mpa=expand(os.path.join(RESULTS_DIR, "mpa_report/contig/{sid}_mpa.tsv"), sid=SAMPLES)
    output:
        combined=os.path.join(RESULTS_DIR, "mpa_report/contig/combined_output.tsv")
    conda:
        os.path.join(ENV_DIR, "krakentools.yaml")
    params:
        combine=os.path.join(SRC_DIR, "combine_mpa_modified.py")
    log:
        os.path.join(RESULTS_DIR, "logs/contig/mpa_combine.log")
    message:
        "Creating a combined mpa-style report"
    shell:
        "(date && {params.combine} -i {input.mpa} -d $(dirname {output.combined}) && date)  &> >(tee {log})"


