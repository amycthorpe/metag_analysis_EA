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
        expand(os.path.join(RESULTS_DIR, "seed/{sid}/{sid}.seed.reduced_counts.tsv"), sid=SAMPLES.index),
        os.path.join(RESULTS_DIR, "seed/abundance_tab.csv")
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

# Concatenating counts for all samples
rule concatenate_abundance:
    output:
        os.path.join(RESULTS_DIR, "seed/allsamples_anno_concatenated.tsv")
    params:
        counts=os.path.join(RESULTS_DIR, "seed/{sid}/{sid}.seed.reduced_counts.tsv")
    log:
        os.path.join(RESULTS_DIR, "logs/seed/concatenating_allsamples.log")
    message:
        "Concatenating reduced SEED counts"
    run:
        # Get the input and output paths
        input_path = os.path.dirname(output[0])
        output_path = output[0]

        with open(output_path, "w") as out:
            for root, dirs, filenames in os.walk(input_path):
                for filename in filenames:

                    # get files with ".seed.reduced_counts.tsv" ending
                    if filename.endswith('.seed.reduced_counts.tsv'):

                        # Get full file path of current file.
                        filename_path = os.path.join(root, filename)

                        # Get sample name from filename.
                        sample_name = filename.split(".seed.reduced_counts.tsv")[0]

                        # Open file and read it line by line.
                        with open(filename_path, "r") as samsa_reduced_abund:

                            # Iterate over each line.
                            for line in samsa_reduced_abund:

                                # Split line into list of fields.
                                fields = line.strip("\n").split("\t")

                                # Extract fields.
                                prop = fields[0]
                                no_reads = fields[1]
                                gene = fields[2]

                                # Check if gene name is blank.
                                if gene == "":

                                    # If gene name is blank, replace it with "unclassified_" and fourth field.
                                    gene = "unclassified_" + fields[3]

                                # Write sample name, gene name, and number of reads to output file.
                                out.write(sample_name + "\t" + gene + "\t" + no_reads + "\n")

# Converting concatenated abundance to matrix
rule matrify:
    input:
        rules.concatenate_abundance.output[0]
    output:
        os.path.join(RESULTS_DIR, "seed/abundance_tab.csv")
    message:
        "Converting concatenated SEED abundances to a matrix"
    run:
        import pandas as pd
        
        # Define the column names
        column_names = ['sample_name', 'gene', 'no_read']
        
        # Read the TSV file into a DataFrame
        all_anno_concatenated = pd.read_csv(input[0], sep='\t', quoting=3, names=column_names, header=None)

        # Assuming that 'sample_name', 'gene', and 'no_reads, i.e. number of reads' columns exist in the DataFrame
        abundance_tab = all_anno_concatenated.pivot(index='sample_name', columns='gene', values='no_reads')

        # Write the DataFrame to a CSV file
        abundance_tab.to_csv(output[0], index=True)

