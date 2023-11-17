"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-08-30]
Run: snakemake -s workflow/rules/taxonomy.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run Kraken2+BRACKEN on reads
"""


############################################
# params:
samples=list(SAMPLES.index)


############################################
rule taxonomy:
    input:
        expand(os.path.join(RESULTS_DIR, "kraken2/{sid}_kraken.report"), sid=SAMPLES.index), 
        expand(os.path.join(RESULTS_DIR, "bracken/{sid}.bracken"), sid=SAMPLES.index),
        expand(os.path.join(RESULTS_DIR, "mpa_report/{sid}_mpa.tsv"), sid=SAMPLES.index),
        os.path.join(RESULTS_DIR, "mpa_report/combined_output.tsv"),
        os.path.join(RESULTS_DIR, "bracken/combined_bracken.txt")
    output:
        touch("status/taxonomy.done")


############################################
localrules: phyloseq_input_kraken2


############################################
# Taxonomic classification using KRAKEN2
rule kraken2:
    input:
        r1=os.path.join(RESULTS_DIR, "preprocessed/trimmed/{sid}/{sid}_val_1.fq.gz"),
        r2=os.path.join(RESULTS_DIR, "preprocessed/trimmed/{sid}/{sid}_val_2.fq.gz")
#         [lambda wildcards: SAMPLES.loc[wildcards.sid, "sR1"], lambda wildcards: SAMPLES.loc[wildcards.sid, "sR2"]]        
#        os.path.join(DATA_DIR, "{sid}_R1.fq.gz"),  
#        os.path.join(DATA_DIR, "{sid}_R2.fq.gz") 
    output:
        report=os.path.join(RESULTS_DIR, "kraken2/{sid}_kraken.report"),
        summary=os.path.join(RESULTS_DIR, "kraken2/{sid}_kraken.out")
    conda:
        os.path.join(ENV_DIR, "kraken2.yaml")
    threads:
        config['kraken2']['threads']
    params:
        db=config['kraken2']['db'],
        confidence=config['kraken2']['confidence']
    log:
        os.path.join(RESULTS_DIR, "logs/kraken2.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Running kraken2 on {wildcards.sid}"
    shell:
        "(date && kraken2 --threads {threads} --db {params.db} --confidence {params.confidence} --paired --output {output.summary} --report {output.report} {input} && date) &> >(tee {log})"
#        "(date && kraken2 --threads {threads} --db {params.db} --paired --output {output.summary} --report {output.report} {input} && date) &> >(tee {log})"

# Running Struo2 database
use rule kraken2 as struo2_kraken2 with:
    output:
        report=os.path.join(RESULTS_DIR, "kraken2/struo2_{sid}_kraken.report"),
        summary=os.path.join(RESULTS_DIR, "kraken2/struo2_{sid}_kraken.out")
    threads:
        config['struo2_kraken2']['threads']
    params:
        db=config['struo2_kraken2']['db'],
        confidence=config['kraken2']['confidence']
    log:
        os.path.join(RESULTS_DIR, "logs/struo2_kraken2.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Running {wildcards.sid} with the Struo2_Kraken2 db"

# Running KRAKEN2+BRACKEN as suggested by KRAKEN2 website
rule bracken:
    input:
        report=rules.kraken2.output.report,
        r1=os.path.join(RESULTS_DIR, "preprocessed/trimmed/{sid}/{sid}_val_1.fq.gz"),  # os.path.join(DATA_DIR, "00.RawData/{sid}/{sid}_R1.fastq.gz"),
        r2=os.path.join(RESULTS_DIR, "preprocessed/trimmed/{sid}/{sid}_val_2.fq.gz") # os.path.join(DATA_DIR, "00.RawData/{sid}/{sid}_R2.fastq.gz")
    output:
        bracken=os.path.join(RESULTS_DIR, "bracken/{sid}.bracken"),
        report=os.path.join(RESULTS_DIR, "bracken/{sid}_bracken.report")
    threads:
        config['kraken2']['threads']
    conda:
        os.path.join(ENV_DIR, "bracken.yaml")
    params:
        db=config['kraken2']['db'],
        read=config['kraken2']['read'],
        level=config['kraken2']['level'],
        bracken=config['bracken']['bin']
    log:
        os.path.join(RESULTS_DIR, "logs/bracken.{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Running kraken & bracken for {wildcards.sid}"
    shell:
        "(date && {params.bracken} -d {params.db} -i {input.report} -o {output.bracken} -w {output.report} -r {params.read} -l {params.level} && date)  &> >(tee {log})"

rule remove_uncultured:
    input:
        bracken=os.path.join(RESULTS_DIR, "bracken/{sid}.bracken")
    output:
        edited=os.path.join(RESULTS_DIR, "bracken/{sid}_edited.bracken")
    log:
        os.path.join(RESULTS_DIR, "logs/edited_bracken_{sid}")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Removing 'uncultured' taxa from bracken output from {wildcards.sid} due to combining issues"
    shell:
        "(date && grep -v 'uncultured' {input.bracken} | grep -v 'endosymbionts' | grep -v 'Incertae Sedis' > {output.edited} && date) &> >(tee {log})"

rule combine_bracken:
    input:
        bracken=expand(os.path.join(RESULTS_DIR, "bracken/{sid}_edited.bracken"), sid=samples)
    output:
        out=os.path.join(RESULTS_DIR, "bracken/combined_bracken.txt")
    conda:
        os.path.join(ENV_DIR, "python2.yaml")
    params:
        combine=config['bracken']['combine']
    log:
        os.path.join(RESULTS_DIR, "logs/bracken_combine.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Combining all the output from BRACKEN"
    shell:
        "(date && python {params.combine} --files {input.bracken} -o {output.out} && date)  &> >(tee {log})"


#########################
### MPA-style report ###
rule mpa_report:
    input:
        report=os.path.join(RESULTS_DIR, "bracken/{sid}_bracken.report")
    output:
        mpa=os.path.join(RESULTS_DIR, "mpa_report/{sid}_mpa.tsv")
    conda:
        os.path.join(ENV_DIR, "bracken_new.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/mpa_{sid}.log")
    wildcard_constraints:
        sid="|".join(SAMPLES.index)
    message:
        "Creating mpa-style report for {wildcards.sid}"
    shell:
        "(date && kreport2mpa.py -r {input.report} -o {output.mpa} && date)  &> >(tee {log})"

rule combine_mpa:
    input:
        mpa=expand(os.path.join(RESULTS_DIR, "mpa_report/{sid}_mpa.tsv"), sid=samples)
    output:
        combined=os.path.join(RESULTS_DIR, "mpa_report/combined_output.tsv")
    conda:
        os.path.join(ENV_DIR, "krakentools.yaml")
    params:
        combine=os.path.join(SRC_DIR, "combine_mpa_modified.py")
    log:
        os.path.join(RESULTS_DIR, "logs/mpa_combine.log")
    message:
        "Creating a combined mpa-style report"
    shell:
        "(date && {params.combine} -i {input.mpa} -d $(dirname {output.combined}) && date)  &> >(tee {log})"



##########################################
### LEGACY CODE - to be edited/removed ###
##########################################
### Phyloseq conversion ###
# Convert kraken output to phyloseq input
rule phyloseq_input_kraken2_sample:
    input:
        rules.kraken2.output.report
    output:
        biom=temp(os.path.join(RESULTS_DIR, "kraken2/{sid}.report.biom")),
        tsv=os.path.join(RESULTS_DIR, "kraken2/{sid}.report.tsv")
    conda:
        os.path.join(ENV_DIR, "biom.yaml")
    message:
        "Phyloseq input for Kraken2: {wildcards.sid}"
    shell:
        # kraken-report -> biom -> TSV
        "kraken-biom {input} -o {output.biom} --fmt hdf5 && "
        "biom convert -i {output.biom} -o {output.tsv} --to-tsv --header-key taxonomy && "
        # add unclassified
        "echo -e \"\\n0\\t$(grep -P 'U\\t0\\t\s*unclassified$' {input} | cut -f3)\\tk__; p__; c__; o__; f__; g__; s__\" >> {output.tsv}"

# Convert Struo2 kraken output to phyloseq input
rule phyloseq_input_struo2_kraken2:
    input:
        rules.struo2_kraken2.output.report
    output:
        biom=temp(os.path.join(RESULTS_DIR, "kraken2/{sid}.struo2_report.biom")),
        tsv=os.path.join(RESULTS_DIR, "kraken2/{sid}.struo2_report.tsv")
    wildcard_constraints:
        sid="|".join(SAMPLES),
    conda:
        os.path.join(ENV_DIR, "biom.yaml")
    message:
        "Phyloseq input for Kraken2: {wildcards.sid}"
    shell:
        # kraken-report -> biom -> TSV
        "kraken-biom {input} -o {output.biom} --fmt hdf5 && "
        "biom convert -i {output.biom} -o {output.tsv} --to-tsv --header-key taxonomy && "
        # add unclassified
        "echo -e \"\\n0\\t$(grep -P 'U\\t0\\t\s*unclassified$' {input} | cut -f3)\\tk__; p__; c__; o__; f__; g__; s__\" >> {output.tsv}"        


############
# Merging all samples together into one table
rule phyloseq_input_kraken2:
    input:
        expand(os.path.join(RESULTS_DIR, "kraken2/{sid}.report.tsv"), sid=SAMPLES)
    output:
        abund=os.path.join(RESULTS_DIR, "kraken2/phyloseq.kraken_abundance.tsv"),
        tax=os.path.join(RESULTS_DIR, "kraken2/phyloseq.kraken_taxonomy.tsv")
    # wildcard_constraints:
    #     db="|".join(config["kraken2"]["db"].keys())
    message:
        "Phyloseq input for Kraken2"
    wildcard_constraints:
        sid="|".join(SAMPLES)
    run:
        import os
        import re
        import csv
        import pandas
        from collections import OrderedDict
        
        # Required functions
        from scripts.utils import read_kraken2_report2biom2tsv, lineage2dict, RANKS
#        RANKS = OrderedDict({'k': 'kingdom', 'p': 'phylum', 'c': 'class', 'o': 'order', 'f': 'family', 'g': 'genus', 's': 'species'})

        dfs   = None
        tax_d = dict()
        for ifile_path in input:
            sid = os.path.basename(ifile_path).split(".")[0]
            df  = read_kraken2_report2biom2tsv(ifile_path)
            df  = df.rename(columns={"count": sid})
            df.set_index("taxid", drop=True, inplace=True, verify_integrity=True)
            if dfs is None:
                dfs = df[[sid]].copy()
                tax_d.update(df["lineage"].to_dict())
            else:
                dfs = dfs.merge(
                    right=df[[sid]],
                    how="outer",
                    left_index=True, right_index=True
                )
                tax_d.update(df["lineage"].to_dict())
        # abundance
        df_a = dfs.copy()
        # taxonomy
        df_t = pandas.DataFrame(index=tax_d.keys(), columns=RANKS.values())
        for t_id, t_lin in tax_d.items():
            for t_r, t_n in lineage2dict(t_lin).items():
                df_t.loc[t_id, RANKS[t_r]] = t_n
        # save
        df_a.to_csv(output.abund, sep="\t", header=True, index=True, index_label="taxonID", na_rep=0)
        df_t.to_csv(output.tax, sep="\t", header=True, index=True, index_label="taxonID")


# Phyloseq input from Struo2 report files
rule phyloseq_input_struo2_kraken2_sample:
    input:
        expand(os.path.join(RESULTS_DIR, "kraken2/{sid}.struo2_report.tsv"), sid=SAMPLES)
    output:
        abund=os.path.join(RESULTS_DIR, "kraken2/phyloseq.struo2_kraken_abundance.tsv"),
        tax=os.path.join(RESULTS_DIR, "kraken2/phyloseq.struo2_kraken_taxonomy.tsv")
    message:
        "Phyloseq input for Struo2_Kraken2"
    run:
        import os
        import re
        import csv
        import pandas
        from collections import OrderedDict
        
        # Required functions
        from scripts.utils import read_kraken2_report2biom2tsv, lineage2dict, RANKS

        dfs   = None
        tax_d = dict()
        for ifile_path in input:
            sid = os.path.basename(ifile_path).split(".")[0]
            df  = read_kraken2_report2biom2tsv(ifile_path)
            df  = df.rename(columns={"count": sid})
            df.set_index("taxid", drop=True, inplace=True, verify_integrity=True)
            if dfs is None:
                dfs = df[[sid]].copy()
                tax_d.update(df["lineage"].to_dict())
            else:
                dfs = dfs.merge(
                    right=df[[sid]],
                    how="outer",
                    left_index=True, right_index=True
                )
                tax_d.update(df["lineage"].to_dict())
        # abundance
        df_a = dfs.copy()
        # taxonomy
        df_t = pandas.DataFrame(index=tax_d.keys(), columns=RANKS.values())
        for t_id, t_lin in tax_d.items():
            for t_r, t_n in lineage2dict(t_lin).items():
                df_t.loc[t_id, RANKS[t_r]] = t_n
        # save
        df_a.to_csv(output.abund, sep="\t", header=True, index=True, index_label="taxonID", na_rep=0)
        df_t.to_csv(output.tax, sep="\t", header=True, index=True, index_label="taxonID")        
