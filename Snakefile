rule all:
    input:
        auspice_tree = "auspice/hpylori_tree.json",
        auspice_meta = "auspice/hpylori_meta.json"

# Config variables to be used by rules
# Parameters are defined within their own rules

rule config:
    params:
        seq = "data/Hp_J99_032020.vcf.gz",
        ref = "data/AE001439.fa",
        meta = "data/Hp_J99_metadata_032020.tsv",
        exclude = "config/dropped_strains.txt",
        generef = "config/Hp_J99.gff3",
        genes = "config/genes.txt",
        clades = "config/clades.tsv",
        colors = "config/color.tsv",
        config = "config/config.json",
        

config = rules.config.params #so we can use config.x rather than rules.config.params.x
#end of config definition

rule filter:
    input:
        seq = config.seq,
        meta = config.meta,
        exclude = config.exclude
    output:
        "results/filtered.vcf.gz"
    shell:
        """
        augur filter --sequences {input.seq} \
            --metadata {input.meta} \
            --exclude {input.exclude} \
            --output {output}
        """

rule tree:
    input:
        aln = rules.filter.output,
        ref = config.ref,
    output:
        "results/tree_raw.nwk"
    params:
        method = 'iqtree'
    shell:
        """
        augur tree --alignment {input.aln} \
            --vcf-reference {input.ref} \
            --method {params.method} \
            --output {output}
        """

rule refine:
    input:
        tree = rules.tree.output,
        aln = rules.filter.output,
        metadata = config.meta,
        ref = config.ref
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json",
    params:
        root = 'best',
        coal = 'opt'
    shell:
        """
        augur refine --tree {input.tree} \
            --alignment {input.aln} \
            --vcf-reference {input.ref} \
            --metadata {input.metadata} \
            --timetree \
            --root {params.root} \
            --coalescent {params.coal} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data}
        """

rule ancestral:
    input:
        tree = rules.refine.output.tree,
        alignment = rules.filter.output,
        ref = config.ref
    output:
        nt_data = "results/nt_muts.json",
        vcf_out = "results/nt_muts.vcf"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral --tree {input.tree} \
            --alignment {input.alignment} \
            --vcf-reference {input.ref} \
            --inference {params.inference} \
            --output {output.nt_data} \
            --output-vcf {output.vcf_out}
        """

rule translate:
    input:
        tree = rules.refine.output.tree,
        ref = config.ref,
        gene_ref = config.generef,
        vcf = rules.ancestral.output.vcf_out,
        genes = config.genes
    output:
        aa_data = "results/aa_muts.json",
        vcf_out = "results/translations.vcf",
        vcf_ref = "results/translations_reference.fasta"
    shell:
        """
        augur translate --tree {input.tree} \
            --vcf-reference {input.ref} \
            --ancestral-sequences {input.vcf} \
            --genes {input.genes} \
            --reference-sequence {input.gene_ref} \
            --output {output.aa_data} \
            --alignment-output {output.vcf_out} \
            --vcf-reference-output {output.vcf_ref}
        """

rule clades:
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.aa_data,
        nuc_muts = rules.ancestral.output.nt_data,
        clades = config.clades
    output:
        clade_data = "results/clades.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output {output.clade_data}
        """

rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = config.meta,
        branch_lengths = rules.refine.output.node_data,
        nt_muts = rules.ancestral.output.nt_data,
        aa_muts = rules.translate.output.aa_data,
        color_defs = config.colors,
        config = config.config,
        clades = rules.clades.output.clade_data
    output:
        tree = rules.all.input.auspice_tree,
        meta = rules.all.input.auspice_meta
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.aa_muts} {input.nt_muts} {input.clades} \
            --auspice-config {input.config} \
            --colors {input.color_defs} \
            --output-tree {output.tree} \
            --output-meta {output.meta}
        """
