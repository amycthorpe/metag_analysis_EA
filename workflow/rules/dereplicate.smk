"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-11-20]
Run: snakemake -s workflow/rules/dereplicate.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: Post-processing bins
"""


############################################
rule dereplicate:
    input:
        os.path.join(RESULTS_DIR, "bins/finalbins")
    output:
        touch("status/dereplicate.done")


############################################
# localrules: 


############################################
# Dereplicating the bins
rule prepare_dasTool:
    input:
        metabat=rules.metabat2.output,
        concoct=rules.concoct.output,
        metabinner=rules.metabinner.output
    output:
        metabatout=os.path.join(RESULTS_DIR, "bins/dastool/metabat_das.tsv"),
        concoctout=os.path.join(RESULTS_DIR, "bins/dastool/concoct_das.tsv"),
        metabinnerout=os.path.join(RESULTS_DIR, "bins/dastool/metabinner_das.tsv")
    params:
        value="bin"
    conda: 
        os.path.join(ENV_DIR, "dastool.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/prepare_dasTool.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/prepare_dasTool.err.log")
    message:
        "Preparing files for DasTool"
    shell:
        """
        (date && Fasta_to_Contig2Bin.sh -e fa -i {input.metabat} > {output.metabatout} && 
        perl -pe 's/metabat./{params.value}_metabat_/g' {output.metabatout} > t.txt && mv t.txt {output.metabatout} && 
        perl -pe 's/\t/\t{params.value}_metabinner_/g' {input.metabinner} > {output.metabinnerout} && 
        perl -pe 's/,/\t{params.value}_concoct_/g' {input.concoct} | tail -n +2 > {output.concoctout} && date) 2> {log.err} > {log.out}
        """

rule dasTool:
    input:
        metabat=rules.prepare_dasTool.output.metabatout,
        concoct=rules.prepare_dasTool.output.concoctout,
        metabinner=rules.prepare_dasTool.output.metabinnerout,
        cont=rules.cat_filter_length.output
    output:
        os.path.join(RESULTS_DIR, "bins/dastool/das_DASTool_summary.tsv")
    conda:
        os.path.join(ENV_DIR, "dastool.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/dasTool.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/dasTool.err.log")
    threads:
        config["dasTool"]["threads"]
    message:
        "Running DasTool"
    shell:
        """
        (date && 
        DAS_Tool -i {input.concoct},{input.metabat},{input.metabinner} -l concoct,metabat,metabinner -c {input.cont} -o $(dirname {output})/das --write_bins --search_engine diamond --threads {threads} && 
        date) 2> {log.err} > {log.out}        
        """

rule checkm:
    input:
        rules.dasTool.output
    output:
        os.path.join(RESULTS_DIR, "bins/dastool_checkm.tsv")
    conda:
        os.path.join(ENV_DIR, "checkm.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/dasTool_checkm.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/dasTool_checkm.err.log")
    threads:
        config["checkm"]["threads"]
    message:
        "Running CheckM on each DasTool output"
    shell:
        """
        (date && checkm lineage_wf -t {threads} -f {output} --tab_table -x fa $(dirname {input})/*bins $(dirname {output})/dastool_checkm && date) 2> {log.err} > {log.out}
        """

rule drep_prepare:
    input:
        check=rules.checkm.output
    output:
        final=os.path.join(RESULTS_DIR, "bins/checkmbeforedrep.tsv"),
        temp=temp(os.path.join(RESULTS_DIR, "bins/checkm_temp.tsv"))
    shell:
        """
        awk '{{print $1".fa,"$13","$14}}' {input.check} | perl -pe "s/Bin.*\n//g" > {output.temp} && 
        (echo "genome,completeness,contamination" && cat {output.temp}) > {output.final}
        """

def symlink_relative(files, input_dir, output_dir):
    """create symlink with and adjust for relative path"""

    input_dir_rel=os.path.relpath(input_dir, output_dir)

    for f in files:
        os.symlink(os.path.join(input_dir_rel, f), os.path.join(output_dir, f))
        
rule get_all_bins:
    input:
        bins=os.path.join(RESULTS_DIR, "bins/dastool/das_DASTool_summary.tsv")
    output:
        directory(os.path.join(RESULTS_DIR,"bins/binsbeforedrep")),
    run:
        os.mkdir(output[0])
        #folder=os.path.dirname(input.bins)
        for files in input.bins:
            bin_folder=os.path.dirname(files) + '/das_DASTool_bins'
            fasta_files=[f for f in os.listdir(bin_folder) if f.endswith(".fa")]
            symlink_relative(fasta_files, bin_folder, output[0])

rule drep:
    input:
        check=rules.drep_prepare.output.final,
        bins=os.path.join(RESULTS_DIR,"bins/binsbeforedrep")
    output:
        temp=directory(os.path.join(RESULTS_DIR, "bins/drep/dereplicated_genomes")),
        final=directory(os.path.join(RESULTS_DIR, "bins/finalbins"))
    conda:
        os.path.join(ENV_DIR, "drep.yaml")
    log:
        out=os.path.join(RESULTS_DIR, "logs/drep/drep.out.log"),
        err=os.path.join(RESULTS_DIR, "logs/drep/drep.err.log")
    params:
        comp=config["drep"]["comp"],
        cont=config["drep"]["cont"]
    threads:
        config["drep"]["threads"]
    message:
        "Running dRep on all the bins"
    shell:
        "(date && "
        "dRep dereplicate $(dirname {output.temp}) -p {threads} -comp {params.comp} -con {params.cont} --genomeInfo {input.check} -g {input.bins}/*fa && "
        "cp -r {output.temp} {output.final} && date) 2> {log.err} > {log.out}"
