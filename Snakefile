import os

configfile: "config.yaml"

with open(config["accessions_file"]) as f:
    SAMPLES = [line.strip() for line in f if line.strip()]

with open(config["strain_id_map"]) as f:
    next(f) # header
    STRAIN_ID = dict(line.strip().split("\t") for line in f if line.strip())

SRA_DIR =        config["sra_dir"]
FASTQ_DIR =      config["fastq_dir"]
TRIMMED_DIR =    config["trimmed_dir"]
ASSEMBLY_DIR =   config["assembly_dir"]
ANNOTATION_DIR = config["annotation_dir"]
MODEL_DIR =      config["model_dir"] 
LOGS_DIR =       config["logs_dir"]    
TMP_DIR =        config["tmp_dir"]
BIOLOG_DATA =    config["biolog_data"]

ALL_SAMPLES = list(SAMPLES)

def _pick(key):
	chosen = config.get(key)
	if chosen is None:
		return list(ALL_SAMPLES)
	unknown = [s for s in chosen if s not in ALL_SAMPLES]
	if unknown:
		raise ValueError(f"{key}: not in accessions_file: {unknown}")
	return list(chosen)

MODEL_SAMPLES = _pick("model_samples") 
SCAN_SAMPLES  = _pick("scan_samples")
BOTH_SAMPLES  = [s for s in MODEL_SAMPLES if s in SCAN_SAMPLES]

rule all:
    input:
#        expand(f"{FASTQ_DIR}/{{sample}}_1.fastq.gz", sample=SAMPLES),
#        expand(f"{FASTQ_DIR}/{{sample}}_2.fastq.gz", sample=SAMPLES),
#
#        expand(f"{TRIMMED_DIR}/{{sample}}_1.trim.fastq.gz", sample=SAMPLES),
#        expand(f"{TRIMMED_DIR}/{{sample}}_2.trim.fastq.gz", sample=SAMPLES),
#
        expand(f"{ASSEMBLY_DIR}/{{sample}}/contigs.fasta", sample=ALL_SAMPLES),
        expand(f"{ASSEMBLY_DIR}/{{sample}}/busco", sample=ALL_SAMPLES),
        expand(f"{ANNOTATION_DIR}/{{sample}}/bakta/{{sample}}.faa", sample=ALL_SAMPLES),
        expand(f"{MODEL_DIR}/{{sample}}/{{sample}}.draft.xml", sample=MODEL_SAMPLES),
        expand(f"{MODEL_DIR}/{{sample}}/{{sample}}.gapfilled.xml", sample=MODEL_SAMPLES),
        expand(f"{MODEL_DIR}/{{sample}}/{{sample}}.fba_media_test.tsv", sample=MODEL_SAMPLES),
        expand(f"{MODEL_DIR}/{{sample}}/{{sample}}.fba_biolog_compare.tsv", sample=MODEL_SAMPLES),
        expand(f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.operon_gene_check.tsv", sample=MODEL_SAMPLES),
        expand(f"{MODEL_DIR}/{{sample}}/{{sample}}.summary_report.pdf", sample=MODEL_SAMPLES),
        expand(f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.metabolic_gene_scan.tsv", sample=SCAN_SAMPLES),
        expand(f"{MODEL_DIR}/{{sample}}/{{sample}}.knockout_single.tsv", sample=BOTH_SAMPLES),
        expand(f"{MODEL_DIR}/{{sample}}/{{sample}}.knockout_combined.tsv", sample=BOTH_SAMPLES),
        expand(f"{MODEL_DIR}/{{sample}}/{{sample}}.knockout_crosscheck.tsv", sample=BOTH_SAMPLES),
        expand(f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.glc_locus_check.tsv", sample=SCAN_SAMPLES),
        [f"{MODEL_DIR}/cross_strain/candidate_table.tsv"] if BOTH_SAMPLES else []


# convert sra to fastq
#rule sra_to_fastq:
#    input:
#        sra=f"{SRA_DIR}/{{sample}}/{{sample}}.sra",
#    output:
#        r1=f"{FASTQ_DIR}/{{sample}}_1.fastq.gz",
#        r2=f"{FASTQ_DIR}/{{sample}}_2.fastq.gz",
#    conda:
#        "envs/sra_tools.yaml"
#    threads: 8
#    resources:
#        mem_mb=8000,
#        runtime="2h",
#    params:
#        tmp=TMP_DIR,
#        outdir=FASTQ_DIR,
#    log:
#        f"{LOGS_DIR}/{{sample}}/{{sample}}.fasterq.log",
#    shell:
#        """
#        mkdir -p {params.tmp} {params.outdir} {LOGS_DIR}
#
#        fasterq-dump {input.sra} \
#            --split-files \
#            --threads {threads} \
#            --temp {params.tmp} \
#            --outdir {params.outdir} \
#            > {log} 2>&1
#
#        pigz -p {threads} \
#            {params.outdir}/{wildcards.sample}_1.fastq \
#            {params.outdir}/{wildcards.sample}_2.fastq
#        """


# read trimming/QC
#rule fastp_trim:
#    input:
#        r1=f"{FASTQ_DIR}/{{sample}}_1.fastq.gz",
#        r2=f"{FASTQ_DIR}/{{sample}}_2.fastq.gz",
#    output:
#        r1=f"{TRIMMED_DIR}/{{sample}}_1.trim.fastq.gz",
#        r2=f"{TRIMMED_DIR}/{{sample}}_2.trim.fastq.gz",
#        json=f"{TRIMMED_DIR}/{{sample}}.fastp.json",
#        html=f"{TRIMMED_DIR}/{{sample}}.fastp.html",
#    conda:
#        "envs/fastp.yaml"
#    threads: 4
#    resources:
#        mem_mb=4000,
#        runtime="1h",
#    log:
#        f"{LOGS_DIR}/{{sample}}/{{sample}}.fastp.log",
#    shell:
#        """
#        mkdir -p {TRIMMED_DIR}
#
#        fastp \
#            -i {input.r1} -I {input.r2} \
#            -o {output.r1} -O {output.r2} \
#            --json {output.json} \
#            --html {output.html} \
#            --thread {threads} \
#            > {log} 2>&1
#        """

# genome assembly
#rule spades_assembly:
#    input:
#        r1=f"{TRIMMED_DIR}/{{sample}}_1.trim.fastq.gz",
#        r2=f"{TRIMMED_DIR}/{{sample}}_2.trim.fastq.gz",
#    output:
#        contigs=f"{ASSEMBLY_DIR}/{{sample}}/contigs.fasta",
#    conda:
#        "envs/spades.yaml"
#    threads: 16
#    resources:
#        mem_mb=128000,
#        runtime="6h",
#        partition="fat",
#    params:
#        outdir=f"{ASSEMBLY_DIR}/{{sample}}",
#        mem_gb=lambda wc, resources: resources.mem_mb // 1000,
#    log:
#        f"{LOGS_DIR}/{{sample}}/{{sample}}.spades.log",
#    shell:
#        """
#        spades.py \
#            --isolate \
#            -1 {input.r1} -2 {input.r2} \
#            -o {params.outdir} \
#            --threads {threads} \
#            --memory {params.mem_gb} \
#            > {log} 2>&1
#        """

# assembly completeness & QC
# must do before on head node - downloads to default busco_downloads
# conda activate busco
# busco --download enterobacterales_odb10
# conda deactivate
rule busco_qc:
    input:
        assembly=f"{ASSEMBLY_DIR}/{{sample}}/contigs.fasta",
    output:
        outdir=directory(f"{ASSEMBLY_DIR}/{{sample}}/busco"),
    conda:
        "envs/busco.yaml"
    threads: 8
    resources:
        mem_mb=8000,
        runtime="2h",
        partition="medium",
    params:
        lineage="enterobacterales_odb10",
        outpath=f"{ASSEMBLY_DIR}/{{sample}}",
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.busco.log",
    shell:
        """
        busco \
            -i {input.assembly} \
            -o busco \
            --out_path {params.outpath} \
            -m genome \
            -l {params.lineage} \
            -c {threads} \
            -f \
            > {log} 2>&1
        """

# gene annotation
rule bakta_annotate:
    input:
        assembly=f"{ASSEMBLY_DIR}/{{sample}}/contigs.fasta",
    output:
        faa=f"{ANNOTATION_DIR}/{{sample}}/bakta/{{sample}}.faa",
        txt=f"{ANNOTATION_DIR}/{{sample}}/bakta/{{sample}}.txt",
        gff3=f"{ANNOTATION_DIR}/{{sample}}/bakta/{{sample}}.gff3",
    conda:
        "envs/bakta.yaml"
    threads: 8
    resources:
        mem_mb=8000,
        runtime="2h",
        partition="medium",
    params:
        outdir=f"{ANNOTATION_DIR}/{{sample}}/bakta",
        db=config["baktadb"],
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.bakta.log",
    shell:
        """
        bakta \
            --db {params.db} \
            --output {params.outdir} \
            --prefix {wildcards.sample} \
            --genus Escherichia \
            --species coli \
            --strain {wildcards.sample} \
            --threads {threads} \
            --force \
            {input.assembly} \
            > {log} 2>&1
        """

# CarveMe draft reconstruction with SCIP solver
rule carveme_draft:
    input:
        faa=f"{ANNOTATION_DIR}/{{sample}}/bakta/{{sample}}.faa",
    output:
        model=f"{MODEL_DIR}/{{sample}}/{{sample}}.draft.xml",
    conda:
        "envs/carveme.yaml"
    threads: 4
    resources:
        mem_mb=8000,
        runtime="4h",
        partition="medium",
    params:
        outdir=f"{MODEL_DIR}/{{sample}}",
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.carveme_draft.log",
    shell:
        """
        carve \
            {input.faa} \
            -o {output.model} \
            -u gramneg \
            --solver scip \
            > {log} 2>&1
        """

# CarveMe gapfill with M9,LB generic media
rule carveme_gapfill:
    input:
        model=f"{MODEL_DIR}/{{sample}}/{{sample}}.draft.xml",
    output:
        model=f"{MODEL_DIR}/{{sample}}/{{sample}}.gapfilled.xml",
    conda:
        "envs/carveme.yaml"
    threads: 2
    resources:
        mem_mb=8000,
        runtime="2h",
        partition="medium",
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.carveme_gapfill.log",
    shell:
        """
        gapfill \
            {input.model} \
            -m M9,LB \
            -o {output.model} \
            > {log} 2>&1
        """

# FBA validation for gapfilled model with M9+glyoxylate
# and M9+rhamnose matching for 13E0634
rule fba_media_test:
    input:
        model=f"{MODEL_DIR}/{{sample}}/{{sample}}.gapfilled.xml",
        script="scripts/fba_media_test.py",
    output:
        results=f"{MODEL_DIR}/{{sample}}/{{sample}}.fba_media_test.tsv",
    conda:
        "envs/cobra.yaml"
    threads: 1
    resources:
        mem_mb=4000,
        runtime="30m",
        partition="medium",
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.fba_media_test.log",
    shell:
        """
        python {input.script} \
            {input.model} \
            -o {output.results} \
            > {log} 2>&1
        """

# Compare FBA validation with paper's results 
# against actual Biolog data in suplementary XLSX
rule fba_biolog_compare:
    input:
        fba_results=f"{MODEL_DIR}/{{sample}}/{{sample}}.fba_media_test.tsv",
        biolog=BIOLOG_DATA,
        script="scripts/compare_fba_biolog.py",
    output:
        compare=f"{MODEL_DIR}/{{sample}}/{{sample}}.fba_biolog_compare.tsv",
    conda:
        "envs/compare.yaml"
    threads: 1
    resources:
        mem_mb=2000,
        runtime="15m",
        partition="medium",
    params:
        labor_id=lambda wc: STRAIN_ID[wc.sample],
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.fba_biolog_compare.log",
    shell:
       """
        python {input.script} \
            --fba {input.fba_results} \
            --biolog {input.biolog} \
            --labor-id {params.labor_id} \
            -o {output.compare} \
            > {log} 2>&1
        """

# which E. coli K-12 MG165 rha/glc operon genes does this
# strain's annotated protein have a good hit for?
# !!! Run extract_reference_genes.py BEFORE running this step !!!
rule operon_gene_check:
    input:
        query="files/rha_glc_genes.faa",
        target_faa=f"{ANNOTATION_DIR}/{{sample}}/bakta/{{sample}}.faa",
        script="scripts/interpret_operon_hits.py",
    output:
        results=f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.operon_gene_check.tsv",
    conda:
        "envs/diamond.yaml"
    threads: 4
    resources:
        mem_mb=4000,
        runtime="15m",
        partition="medium",
    params:
        db=f"{ANNOTATION_DIR}/{{sample}}/bakta/{{sample}}_diamond_db",
        raw_hits=f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.operon_gene_hits.raw.tsv",
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.operon_gene_hits.raw.tsv",
    shell:
        """
        diamond makedb --in {input.target_faa} -d {params.db} > {log} 2>&1

        diamond blastp \
            --query {input.query} \
            --db {params.db} \
            --out {params.raw_hits} \
            --outfmt 6 qseqid sseqid pident qcovhsp evalue bitscore \
            --max-target-seqs 1 \
            --threads {threads} \
            >> {log} 2>&1

        python {input.script} \
            --query-fasta {input.query} \
            --raw-hits {params.raw_hits} \
            -o {output.results} \
            >> {log} 2>&1
        """

# !!! Run build_iml1515_reference.py BEFORE starting this step !!!
# Reciprocal-best-hit map: each iML1515 reference protein -> its best hit in the
# strain's Bakta proteins, and whether that hit is reciprocal (reverse search
# against the whole K-12 proteome). Adapted from get_rbh.py.
rule ortholog_map:
    input:
        query="files/iML1515_metabolic_genes.faa",
        k12_proteome="files/GCF_000005845.2_ASM584v2_protein.faa",
        feature_table="files/GCF_000005845.2_ASM584v2_feature_table.txt",
        target_faa=f"{ANNOTATION_DIR}/{{sample}}/bakta/{{sample}}.faa",
        script="scripts/ortholog_map.py",
    output:
        results=f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.ortholog_map.tsv",
    conda:
        "envs/diamond.yaml"
    threads: 4
    resources:
        mem_mb=4000,
        runtime="15m",
        partition="medium",
    params:
        db_dir=f"{ANNOTATION_DIR}/{{sample}}/ortholog_map_dbs",
        fwd=f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.ortholog_map.fwd.raw.tsv",
        rev=f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.ortholog_map.rev.raw.tsv",
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.ortholog_map.log",
    shell:
        """
        python {input.script} \
            --query {input.query} \
            --reference-proteome {input.k12_proteome} \
            --target {input.target_faa} \
            --feature-table {input.feature_table} \
            --out {output.results} \
            --db_dir {params.db_dir} \
            --fwd_hits {params.fwd} \
            --rev_hits {params.rev} \
            --threads {threads} \
            > {log} 2>&1
        """
# !!! fix_ref_isoforms.py was run after this step for the fix and this step
#     was re-run again. !!!
# retired script. patch implemented in build_iml1515_reference.py

# Genome-wide scan: intact / disrupted / absent per iML1515 metabolic gene,
# from the ortholog map (no DIAMOND call here any more).
rule metabolic_gene_scan:
    input:
        query="files/iML1515_metabolic_genes.faa",
        gene_table="files/iML1515_gene_table.tsv",
        ortholog_map=f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.ortholog_map.tsv",
        script="scripts/scan_metabolic_genes.py",
    output:
        results=f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.metabolic_gene_scan.tsv",
    conda:
        "envs/diamond.yaml"
    threads: 1
    resources:
        mem_mb=2000,
        runtime="10m",
        partition="medium",
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.metabolic_gene_scan.log",
    shell:
        """
        python {input.script} \
            --query-fasta {input.query} \
            --gene-table {input.gene_table} \
            --ortholog-map {input.ortholog_map} \
            -o {output.results} \
            > {log} 2>&1
        """

# Direct search of the reference operon proteins against the assembly (a tblastn
# equivalent), to tell a gene the annotation missed, a frameshifted pseudogene and
# a true absence apart. Place after operon_gene_check. Uses the existing diamond
# env.
rule glc_locus_check:
    input:
        query="files/rha_glc_genes.faa",
        assembly=f"{ASSEMBLY_DIR}/{{sample}}/contigs.fasta",
        gff3=f"{ANNOTATION_DIR}/{{sample}}/bakta/{{sample}}.gff3",
        script="scripts/locus_check.py",
    output:
        loci=f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.glc_locus_check.tsv",
        neighbourhood=f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.glc_neighbourhood.tsv",
    conda:
        "envs/diamond.yaml"
    threads: 4
    resources:
        mem_mb=4000,
        runtime="15m",
        partition="medium",
    params:
        db_dir=f"{ANNOTATION_DIR}/{{sample}}/locus_check",
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.glc_locus_check.log",
    shell:
        """
        python {input.script} \
            --query {input.query} \
            --assembly {input.assembly} \
            --gff3 {input.gff3} \
            --out-loci {output.loci} \
            --out-neighbourhood {output.neighbourhood} \
            --window-genes glcA,glcB,glcC,glcD,glcE,glcF,glcG \
            --db_dir {params.db_dir} \
            --threads {threads} \
            > {log} 2>&1
        """

# Add glyoxylate transport + exchange to the gap-filled model. CarveMe's
# universe has no EX_glx_e at all (confirmed from the universe build
# log, not gene-driven). Uses model.metabolites'
# real BiGG IDs (glx_c, h_e, h_c) directly, unlike an earlier
# standalone attempt at this that guessed wrong/differently-spelled
# IDs (glyox[c] vs the real glx_c) and ended up building an entire
# pathway disconnected from the model. GLXCL already exists natively
# in every model checked so far — this curation step only plugs the
# external-transport gap; whether that's enough to reach biomass is
# an empirical question the FBA rerun below actually answers, not
# something hand-built pathway reactions should presume upfront.
rule curate_model:
    input:
        model=f"{MODEL_DIR}/{{sample}}/{{sample}}.gapfilled.xml",
        script="scripts/curate_model.py",
    output:
        model=f"{MODEL_DIR}/{{sample}}/{{sample}}.curated.xml",
    conda:
        "envs/cobra.yaml"
    threads: 1
    resources:
        mem_mb=2000,
        runtime="15m",
        partition="medium",
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.curate_model.log",
    shell:
        """
        python {input.script} {input.model} -o {output.model} > {log} 2>&1
        """


# Same script as the original fba_media_test, just pointed at the
# curated model instead of gapfilled.xml — kept as a separate output
# (not overwriting fba_media_test.tsv) so before/after stays comparable.
rule fba_media_test_curated:
    input:
        model=f"{MODEL_DIR}/{{sample}}/{{sample}}.curated.xml",
        script="scripts/fba_media_test.py",
    output:
        results=f"{MODEL_DIR}/{{sample}}/{{sample}}.fba_media_test_curated.tsv",
    conda:
        "envs/cobra.yaml"
    threads: 1
    resources:
        mem_mb=4000,
        runtime="30m",
        partition="medium",
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.fba_media_test_curated.log",
    shell:
        """
        python {input.script} {input.model} -o {output.results} > {log} 2>&1
        """

rule fba_biolog_compare_curated:
    input:
        fba_results=f"{MODEL_DIR}/{{sample}}/{{sample}}.fba_media_test_curated.tsv",
        biolog=BIOLOG_DATA,
        script="scripts/compare_fba_biolog.py",
    output:
        compare=f"{MODEL_DIR}/{{sample}}/{{sample}}.fba_biolog_compare_curated.tsv",
    conda:
        "envs/compare.yaml"
    threads: 1
    resources:
        mem_mb=2000,
        runtime="15m",
        partition="medium",
    params:
        labor_id=lambda wc: STRAIN_ID[wc.sample],
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.fba_biolog_compare_curated.log",
    shell:
        """
        python {input.script} \
            --fba {input.fba_results} \
            --biolog {input.biolog} \
            --labor-id {params.labor_id} \
            -o {output.compare} \
            > {log} 2>&1
        """

# Knock out the metabolic genes the genome scan called disrupted or absent,
# via the iML1515 GPR rules (reaction-level), in the curated model, and
# rerun FBA on M9 + glucose / glyoxylate / rhamnose. One row per flagged
# gene on its own, plus one combined all-flagged knockout (table + model).
rule knockout_flagged_genes:
    input:
        model=f"{MODEL_DIR}/{{sample}}/{{sample}}.curated.xml",
        scan=f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.metabolic_gene_scan.tsv",
        iml1515="files/iML1515.json.gz",
        script="scripts/knockout_flagged_genes.py",
        media_script="scripts/fba_media_test.py",  # imported by the knockout script
    output:
        single=f"{MODEL_DIR}/{{sample}}/{{sample}}.knockout_single.tsv",
        combined=f"{MODEL_DIR}/{{sample}}/{{sample}}.knockout_combined.tsv",
        model=f"{MODEL_DIR}/{{sample}}/{{sample}}.knockout_combined.xml",
    conda:
        "envs/cobra.yaml"
    threads: 1
    resources:
        mem_mb=4000,
        runtime="30m",
        partition="medium",
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.knockout_flagged_genes.log",
    shell:
        """
        python {input.script} \
            --iml1515-json {input.iml1515} \
            --scan {input.scan} \
            --model {input.model} \
            --out-single {output.single} \
            --out-combined {output.combined} \
            --out-model {output.model} \
            > {log} 2>&1
        """

# Cross-check of knockout_flagged_genes: knocks out the strain's own locus
# tags with COBRApy (flagged genes whose best hit is a clean, unshared gene
# in the strain model) and compares the reactions removed with the
# iML1515-GPR knockout_single table. Genes that cannot be checked this way
# stay in the table with the reason.
rule knockout_crosscheck:
    input:
        model=f"{MODEL_DIR}/{{sample}}/{{sample}}.curated.xml",
        scan=f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.metabolic_gene_scan.tsv",
        single=f"{MODEL_DIR}/{{sample}}/{{sample}}.knockout_single.tsv",
        script="scripts/knockout_crosscheck.py",
        media_script="scripts/fba_media_test.py",              # imported
        knockout_script="scripts/knockout_flagged_genes.py",   # imported
    output:
        results=f"{MODEL_DIR}/{{sample}}/{{sample}}.knockout_crosscheck.tsv",
    conda:
        "envs/cobra.yaml"
    threads: 1
    resources:
        mem_mb=4000,
        runtime="30m",
        partition="medium",
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.knockout_crosscheck.log",
    shell:
        """
        python {input.script} \
            --scan {input.scan} \
            --model {input.model} \
            --knockout-single {input.single} \
            -o {output.results} \
            > {log} 2>&1
        """

# It runs once for all strains in BOTH_SAMPLES (the accession file), so with
# pruns_SRR.txt it merges all three strains into one table. Column headers use
# the LaborIDs from the strain map (13E0634_call, ...).
rule candidate_table:
    input:
        scan=expand(f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.metabolic_gene_scan.tsv", sample=BOTH_SAMPLES),
        single=expand(f"{MODEL_DIR}/{{sample}}/{{sample}}.knockout_single.tsv", sample=BOTH_SAMPLES),
        cross=expand(f"{MODEL_DIR}/{{sample}}/{{sample}}.knockout_crosscheck.tsv", sample=BOTH_SAMPLES),
        locus=expand(f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.glc_locus_check.tsv", sample=BOTH_SAMPLES),
        script="scripts/candidate_table.py",
    output:
        table=f"{MODEL_DIR}/cross_strain/candidate_table.tsv",
    params:
        strains=lambda wc, input: " ".join(
            f"--strain {STRAIN_ID[s]} {sc} {sg} {cr}"
            for s, sc, sg, cr in zip(BOTH_SAMPLES, input.scan, input.single, input.cross)
        ),
        locus=lambda wc, input: " ".join(
            f"--locus-check {STRAIN_ID[s]} {l}" for s, l in zip(BOTH_SAMPLES, input.locus)
        ),
        threshold=0.01,
    conda:
        "envs/cobra.yaml"
    threads: 1
    resources:
        mem_mb=1000,
        runtime="5m",
        partition="medium",
    log:
        f"{LOGS_DIR}/cross_strain/candidate_table.log",
    shell:
        """
        mkdir -p $(dirname {output.table}) $(dirname {log})
        python {input.script} {params.strains} {params.locus} \
            --growth-threshold {params.threshold} \
            --out {output.table} \
            > {log} 2>&1
        """

# Final per-strain deliverable: pulls BUSCO, Bakta, both CarveMe models,
# fba_media_test, fba_biolog_compare, operon_gene_check, AND the
# curated-model re-test into one PDF with plots.
rule summary_report:
    input:
        busco_dir=f"{ASSEMBLY_DIR}/{{sample}}/busco",
        bakta_txt=f"{ANNOTATION_DIR}/{{sample}}/bakta/{{sample}}.txt",
        draft_model=f"{MODEL_DIR}/{{sample}}/{{sample}}.draft.xml",
        gapfilled_model=f"{MODEL_DIR}/{{sample}}/{{sample}}.gapfilled.xml",
        fba_results=f"{MODEL_DIR}/{{sample}}/{{sample}}.fba_media_test.tsv",
        biolog_compare=f"{MODEL_DIR}/{{sample}}/{{sample}}.fba_biolog_compare.tsv",
        operon_check=f"{ANNOTATION_DIR}/{{sample}}/{{sample}}.operon_gene_check.tsv",
        biolog_xlsx=BIOLOG_DATA,
        curated_compare=f"{MODEL_DIR}/{{sample}}/{{sample}}.fba_biolog_compare_curated.tsv",
        script="scripts/generate_summary_report.py",
    output:
        pdf=f"{MODEL_DIR}/{{sample}}/{{sample}}.summary_report.pdf",
    conda:
        "envs/report.yaml"
    threads: 1
    resources:
        mem_mb=2000,
        runtime="15m",
        partition="medium",
    params:
        labor_id=lambda wc: STRAIN_ID[wc.sample],
    log:
        f"{LOGS_DIR}/{{sample}}/{{sample}}.summary_report.log",
    shell:
        """
        python {input.script} \
            --sample {wildcards.sample} \
            --labor-id {params.labor_id} \
            --busco-dir {input.busco_dir} \
            --bakta-txt {input.bakta_txt} \
            --draft-model {input.draft_model} \
            --gapfilled-model {input.gapfilled_model} \
            --fba-results {input.fba_results} \
            --biolog-compare {input.biolog_compare} \
            --operon-check {input.operon_check} \
            --biolog-xlsx {input.biolog_xlsx} \
            --curated-compare {input.curated_compare} \
            -o {output.pdf} \
            > {log} 2>&1
        """
