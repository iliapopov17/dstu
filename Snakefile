import os

# Set default path if not provided
DATA_PATH = config.get('data_path', 'genbank_sequences')

rule all:
    input:
        "demo_data/tree_ufb.treefile"

rule concatenate_fasta:
    input:
        lambda wildcards: expand(os.path.join(DATA_PATH, "{{sample}}.fasta"), sample=glob_wildcards(os.path.join(DATA_PATH, "{sample}.fasta")).sample)
    output:
        "demo_data/all_seqs.fa"
    shell:
        "cat {input} > {output}"

rule mafft_alignment:
    input:
        "demo_data/all_seqs.fa"
    output:
        "demo_data/all_seqs_mafft.fa"
    shell:
        "mafft --auto {input} > {output}"

rule iqtree_model_finding:
    input:
        "demo_data/all_seqs_mafft.fa"
    output:
        tree="demo_data/tree_MF2.treefile",
        logfile="demo_data/tree_MF2.iqtree"
    shell:
        "iqtree2 -m MFP -s {input} --prefix demo_data/tree_MF2 -T AUTO"

rule extract_best_model:
    input:
        "demo_data/tree_MF2.iqtree"
    output:
        "demo_data/best_model.txt"
    shell:
        "grep 'Best-fit model according to BIC:' {input} | cut -d ':' -f 2 | xargs > {output}"

rule iqtree_specific_model:
    input:
        model="demo_data/best_model.txt",
        seqs="demo_data/all_seqs_mafft.fa"
    output:
        "demo_data/tree_ufb.treefile"
    shell:
        "iqtree2 -s {input.seqs} -m $(cat {input.model}) -pre demo_data/tree_ufb -bb 1000 -nt AUTO"
