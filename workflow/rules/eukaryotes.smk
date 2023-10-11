"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-09-23]
Run: snakemake -s workflow/rules/eukaryotes.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To classify eukaryotes from contigs
"""


############################################
rule eukaryotes:
    input:
#        os.path.join(RESULTS_DIR, "eukulele_output"),
        os.path.join(RESULTS_DIR, "euk_tax_cov/merged_eukaryote_abundances.txt"),
        os.path.join(RESULTS_DIR, "euk_tax_cov/merged_eukulele_all_abundances.txt")
    output:
        touch("status/eukaryotes.done")


############################################
localrules: merge_cov, euk_abundances, merge_cov_all, euk_abundances_all


############################################
# Preparing the files for eukaryote classification
rule file_prep:
    input:
        os.path.join(RESULTS_DIR, "prodigal/{sid}/{sid}.faa")
    output:
        os.path.join(RESULTS_DIR, "eukulele_input/{sid}.faa")
    log:
        os.path.join(RESULTS_DIR, "logs/eukulele/{sid}.log")
    message:
        "Symlinking {wildcards.sid} for eukulele"
    shell:
        "(date && ln -vs {input} {output} && date) &> >(tee {log})"


rule eukulele:
    input:
        os.path.join(RESULTS_DIR, "eukulele_input")
    output:
        protected(directory(os.path.join(RESULTS_DIR, "eukulele_output")))
    conda:
        os.path.join(ENV_DIR, "eukulele.yaml")
    threads:
        config["eukulele"]["threads"]
    log:
        os.path.join(RESULTS_DIR, "logs/eukulele.log")
    params:
        db=config["eukulele"]["db"],
        scratch=config["eukulele"]["scratch"],
        aligner=config["eukulele"]["aligner"],
        ref_dir=config["eukulele"]["ref_dir"]
    message:
        "Running EUKULELE"
    shell:
        "(date && EUKulele all -m mets --sample_dir {input} --out_dir {output[0]} --database {params.db} --scratch {params.scratch} --CPUs {threads} --alignment_choice {params.aligner} --reference_dir {params.ref_dir} && "
        "date) &> >(tee {log})"

rule merge_cov:
    input:
        cov=os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}_gene_coverage.txt")
    output:
        taxa=os.path.join(RESULTS_DIR, "euk_tax_cov/{sid}_all_taxa.txt"),
        euk=os.path.join(RESULTS_DIR, "euk_tax_cov/{sid}_eukaryotes.txt")
    params:
        tax=os.path.join(RESULTS_DIR, "eukulele_output/taxonomy_estimation/{sid}-estimated-taxonomy.out")
    log:
        os.path.join(RESULTS_DIR, "logs/cov_merge.{sid}.log")
    message:
        "Merging Eukaryotic taxonomy with coverage for {wildcards.sid}"
    run:
        tax=pd.read_csv(params.tax, header=0, sep="\t", index_col=0)
        cov=pd.read_csv(input.cov, header=None, sep=" ")
        cov.columns=['gene_id', 'transcript_name', 'coverage']

        # merging taxonomy with coverage
        merged=filt_euks.merge(cov, how="left", on="transcript_name")
        filt_merged=merged[['full_classification','classification', 'coverage']]

        # Grouping same taxonomy and getting sum of coverage
        final=filt_merged.groupby(['full_classification','classification'], as_index=False)['coverage'].sum()
        
        # writing to file
        final.to_csv(output.taxa, index=None, sep="\t")
        
        # Eukaryotes
        # keeping only rows that contain 'Eukary' in the 'full_classification' column
        euks=tax[tax['full_classification'].str.contains("Eukary", na=False)].drop(['counts'], axis=1)
        # keeping only rows that have at least 50% 'max_pid'
        filt_euks=euks.query("max_pid >=50")
        # merging taxonomy with coverage
        merged_euks=filt_euks.merge(cov, how="left", on="transcript_name")
        filt_merged_euk=merged_euks[['full_classification','classification', 'coverage']]

        # Grouping same taxonomy and getting sum of coverage
        final_euks=filt_merged_euk.groupby(['full_classification','classification'], as_index=False)['coverage'].sum()
        
        # writing to file
        final_euks.to_csv(output.euk, index=None, sep="\t")


rule euk_abundances:
    input:
        expand(os.path.join(RESULTS_DIR, "euk_tax_cov/{sid}_eukaryotes.txt"), sid=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "euk_tax_cov/merged_eukaryote_abundances.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/euk_abundances.log")
    message:
        "Merging all eukaryote abundances for all samples"
    run:
        # Collecting all files in folder 
        directory=os.path.dirname(input[0])
        os.chdir(directory)
	
        # verify the path using getcwd() 
        cwd = os.getcwd() 
  
        # print the current directory 
        print("Current working directory is:", cwd) 

        mylist=[f for f in glob.glob("*.txt")]
        mylist

        # making individual dataframes for each file
        dataframes= [ pd.read_csv( f, header=0, sep="\t", usecols=["full_classification", "coverage"]) for f in mylist ] # add arguments as necessary to the read_csv method

        # Merging all files based on common column
        merged=reduce(lambda left,right: pd.merge(left,right,on='full_classification', how='outer'), dataframes)

        # Giving appropriate column names
        names=['full_classification']+mylist
        new_cols=list(map(lambda x: x.replace('_eukaryotes.txt',''),names))
        merged.columns=new_cols

        # checking if any values are "NA"
        merged.isnull().values.any()
        # if "NA" run the following
        merged.fillna('', inplace=True)

        # Removing rows with all zeroes (0 or 0.0)
        merged.set_index('full_classification', inplace=True)  # first to make first column as rownames
        edited=merged.loc[~(merged==0).all(axis=1)]

        # Writing file without zeroes
        edited.to_csv(output[0], sep='\t', index=True, header=True)

rule merge_cov_all:
    input:
        cov=os.path.join(RESULTS_DIR, "coverage/{sid}/{sid}_gene_coverage.txt")
    output:
        os.path.join(RESULTS_DIR, "euk_tax_cov/{sid}_eukulele_ALL.txt")
    params:
        tax=os.path.join(RESULTS_DIR, "eukulele_output/taxonomy_estimation/{sid}-estimated-taxonomy.out")
    log:
        os.path.join(RESULTS_DIR, "logs/cov_merge_all.{sid}.log")
    message:
        "Merging all EUKulele taxonomy with coverage for {wildcards.sid}"
    run:
        tax=pd.read_csv(input.tax, header=0, sep="\t", index_col=0)
        cov=pd.read_csv(input.cov, header=None, sep="\t")
        cov.columns=['transcript_name', 'coverage']

        # keeping only rows that contain 'Eukary' in the 'full_classification' column
        euks=tax.drop(['counts'], axis=1)

        # merging taxonomy with coverage
        merged=euks.merge(cov, how="left", on="transcript_name")
        filt_merged=merged[['full_classification','classification', 'coverage']]

        # Grouping same taxonomy and getting sum of coverage
        final=filt_merged.groupby(['full_classification','classification'], as_index=False)['coverage'].sum()

        # writing to file
        final.to_csv(output[0], index=None, sep="\t")

rule euk_abundances_all:
    input:
        expand(os.path.join(RESULTS_DIR, "euk_tax_cov/{sid}_eukulele_ALL.txt"), sid=SAMPLES)
    output:
        os.path.join(RESULTS_DIR, "euk_tax_cov/merged_eukulele_all_abundances.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/euk_all_abundances.log")
    message:
        "Merging all taxa abundances"
    run:
        # Collecting all files in folder
        directory=os.path.dirname(input[0])
        os.chdir(directory)

        # verify the path using getcwd()
        cwd = os.getcwd()

        # print the current directory
        print("Current working directory is:", cwd)

        mylist=[f for f in glob.glob("*ALL.txt")]
        mylist

        # making individual dataframes for each file
        dataframes= [ pd.read_csv( f, sep="\t", usecols=['full_classification', 'coverage']) for f in mylist ] # add arguments as necessary to the read_csv method

        # Merging all files based on common column
        merged=reduce(lambda left,right: pd.merge(left,right,on='full_classification', how='outer'), dataframes)

        # Giving appropriate column names
        names=['full_classification']+mylist
        new_cols=list(map(lambda x: x.replace('_eukulele_all.txt',''),names))
        merged.columns=new_cols

        # checking if any values are "NA"
        merged.isnull().values.any()
        # if "NA" run the following
        merged.fillna('', inplace=True)

        # Removing rows with all zeroes (0 or 0.0)
        merged.set_index('full_classification', inplace=True)  # first to make first column as rownames
        edited=merged.loc[~(merged==0).all(axis=1)]

        # Writing file without zeroes
        edited.to_csv(output[0], sep='\t', index=True, header=True)
