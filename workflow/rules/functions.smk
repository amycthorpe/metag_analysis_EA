"""
Author: Susheel Bhanu BUSI
Affiliation: Molecular Ecology group, UKCEH
Date: [2023-09-23]
Run: snakemake -s workflow/rules/functions.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: To run MagicLamp and Lithogenies on contigs
"""


############################################
rule functions:
    input:
        expand(os.path.join(RESULTS_DIR, "magiclamp/{sid}/lithogenie_output"), sid=SAMPLES)
    output:
        touch("status/functions.done")


############################################
localrules: install_magiclamp


############################################
# MagicLamp Initial Setup 
rule install_magiclamp:
    output:
        done=os.path.join(RESULTS_DIR, "magiclamp/magiclamp.installed")
    log:
        out=os.path.join(RESULTS_DIR, "logs/setup.magiclamp.log")
    params:
        script=os.path.join(SRC_DIR, "install_magiclamp.sh"), 
        path=os.path.join(SUBMODULES, "magiclamp")
    conda:
        os.path.join(ENV_DIR, "magiclamp.yaml")
    message:
        "Setup: install Magiclamp"
    shell:
        "script=$(realpath {params.script}) && cd {params.path} && ${{script}} && touch {output}"

rule magiclamp:
    input:
        contigs=os.path.join(RESULTS_DIR, "assembly/{sid}/{sid}.fasta"),
        installed=os.path.join(RESULTS_DIR, "magiclamp/magiclamp.installed")
    output:
        directory(os.path.join(RESULTS_DIR, "magiclamp/{sid}/lithogenie_output"))
    params:
        path=os.path.join(SUBMODULES, "MagicLamp")
    conda:
        os.path.join(ENV_DIR, "magiclamp.yaml")
    threads:
        config['magiclamp']['threads']
    log:
        os.path.join(RESULTS_DIR, "logs/magiclamp/{sid}.log")
    message:
        "Running MagicLamps from Magiclamp on {wildcards.sid}"
    shell:
        "(date && export PATH=$PATH:{params.path} && "
        "MagicLamp.py LithoGenie -bin_dir $(dirname {input.contigs}) -bin_ext fasta -out {output} -t {threads} --norm && date) &> >(tee {log})" 

