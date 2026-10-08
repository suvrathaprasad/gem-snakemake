# Project GEM

Genome-scale metabolic modeling of Shiga toxin-producing *Escherichia coli* (STEC)
strains, with model predictions validated against real Biolog phenotype microarray data.

This is a **pilot pipeline**. It takes a bacterial genome through assembly QC,
annotation, automated model reconstruction (CarveMe), flux balance analysis (COBRApy),
and a direct comparison of predicted growth against measured phenotypes, then writes a
one-page-per-topic PDF summary for each strain.

- Source data: BioProject [PRJNA559322](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA559322)
- Phenotype data and strain classification: Barth et al. 2020, *Toxins* 12(6):414
  ([doi:10.3390/toxins12060414](https://doi.org/10.3390/toxins12060414))
- Orchestration: Snakemake on Slurm, one small conda environment per rule

## Pipeline

```
SRA reads -> fastp trim -> SPAdes assembly -> BUSCO (QC) -> Bakta annotation
   -> CarveMe draft model -> CarveMe gap-fill (generic M9 + LB)
   -> FBA growth test (M9 + single carbon source)
   -> compare to Biolog -> operon gene check (diamond)
   -> model curation (glyoxylate transport) -> FBA + Biolog comparison again
   -> per-strain PDF summary report
```

| Stage | Snakemake rule | Tool | Main output |
|---|---|---|---|
| Reads to FASTQ | `sra_to_fastq` | SRA Toolkit, pigz | `fastq/{sample}_{1,2}.fastq.gz` |
| Trimming | `fastp_trim` | fastp | `trimmed/{sample}_{1,2}.trim.fastq.gz` |
| Assembly | `spades_assembly` | SPAdes (`--isolate`) | `assembly/{sample}/contigs.fasta` |
| Assembly QC | `busco_qc` | BUSCO (`enterobacterales_odb10`) | `assembly/{sample}/busco/` |
| Annotation | `bakta_annotate` (or `prokka_annotate`) | Bakta | `annotation/{sample}/bakta/{sample}.faa` |
| Draft model | `carveme_draft` | CarveMe (`-u gramneg`, SCIP) | `models/{sample}/{sample}.draft.xml` |
| Gap-fill | `carveme_gapfill` | CarveMe `gapfill` (M9, LB) | `models/{sample}/{sample}.gapfilled.xml` |
| FBA test | `fba_media_test` | COBRApy (`scripts/fba_media_test.py`) | `{sample}.fba_media_test.tsv` |
| FBA vs Biolog | `fba_biolog_compare` | `scripts/compare_fba_biolog.py` | `{sample}.fba_biolog_compare.tsv` |
| Gene check | `operon_gene_check` | DIAMOND blastp | `annotation/{sample}/{sample}.operon_gene_check.tsv` |
| Curation | `curate_model` | COBRApy (`scripts/curate_model.py`) | `{sample}.curated.xml` |
| Curated FBA | `fba_media_test_curated`, `fba_biolog_compare_curated` | as above | `*_curated.tsv` |
| Report | `summary_report` | matplotlib, reportlab | `models/{sample}/{sample}.summary_report.pdf` |

`spades_assembly` and the read-processing rules were superseded in this work by publicly
available assemblies (see "Input data" below), but are kept in the Snakefile.

## Repository layout

```
Snakefile                  workflow definition
config.yaml                pipeline paths and database locations
strain_id_map.tsv          SRR accession -> strain name (as used in the Biolog table)
environment-core.yml       driver environment (Snakemake + Slurm executor plugin)
envs/                      one small conda env per tool, built by Snakemake on demand
scripts/                   Python scripts called by the rules
workflow/profiles/slurm/   Snakemake profile (cluster-generic sbatch template)
bash_scripts/              original manual SRA download script (superseded)
other_environments/        early standalone conda envs (superseded by envs/)
```

Data and databases (`PRJNA559322/`, `reference/`, BUSCO and Bakta downloads, the Biolog
table) are created locally and are not part of the repository.

## Requirements

- A Slurm cluster, `mamba`/`conda`, and internet access on the login node
- Compute nodes without internet access are fine, as long as the one-time downloads
  below are done on the login node first

```bash
mamba env create -f environment-core.yml
mamba activate gem-core
```

All other tools are installed automatically by Snakemake into per-rule environments
(`use-conda: True`, cached under `$HOME/.cache/snakemake/conda`).

## Input data and one-time setup

Do these on the login node.

1. **Biolog phenotype table.** Table S1 from the supplementary material of
   Barth et al. 2020 (`Supplement_TableS1_Biolog_AUC.xlsx`). Set `biolog_data` in
   `config.yaml`.
2. **Assemblies.** Place one assembly per sample at
   `PRJNA559322/assembly/{SRR}/contigs.fasta`. Either use public assemblies for the
   BioProject, or run the read-based rules (`sra_to_fastq`, `fastp_trim`,
   `spades_assembly`) to build them.
3. **BUSCO lineage.** Pre-fetch so that compute nodes do not need internet:
   ```bash
   mamba env create -f envs/busco.yaml && mamba activate busco
   busco --download enterobacterales_odb10
   ```
4. **Bakta database.** The light database was used. At the time of writing,
   `bakta_db download` failed its MD5 check; downloading the file directly works:
   ```bash
   wget -c "https://zenodo.org/records/14916843/files/db-light.tar.xz?download=1" -O db-light.tar.xz
   md5sum db-light.tar.xz     # expect 4a6e059ded39e9c5537ef4137d2f5648
   tar -xJf db-light.tar.xz
   ```
   Set `baktadb` in `config.yaml` to the extracted folder.
5. **Reference genes for the operon check** (*E. coli* K-12 MG1655, RefSeq GCF_000005845.2):
   ```bash
   B=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2
   wget $B/GCF_000005845.2_ASM584v2_protein.faa.gz $B/GCF_000005845.2_ASM584v2_feature_table.txt.gz
   gunzip GCF_000005845.2_ASM584v2_*.gz
   mkdir -p reference
   python scripts/extract_reference_genes.py \
       --feature-table GCF_000005845.2_ASM584v2_feature_table.txt \
       --protein-faa   GCF_000005845.2_ASM584v2_protein.faa \
       -o reference/rha_glc_genes.faa
   ```
   This extracts 14 genes: *glcA-G* and *rhaA, B, D, M, R, S, T*. Check the output for
   `[warn]` lines.

## Configuration

- `config.yaml`: paths for each stage, `accessions_file` (one SRR per line),
  `baktadb`, `biolog_data`.
- `strain_id_map.tsv`: tab-separated `SRR` and `LaborID` (the strain name used in the
  Biolog table). Add a row for every strain you run.
- Slurm: set `partition=` in each rule's `resources:` block in the `Snakefile`, and the
  defaults in `workflow/profiles/slurm/config.yaml`. Add `--account=...` to the
  `sbatch` template there if your cluster needs one.

## Running

```bash
mamba activate gem-core
snakemake --profile workflow/profiles/slurm -n    # dry run
snakemake --profile workflow/profiles/slurm
```

Snakemake runs on the login node as a light driver and submits each job to Slurm.

**Note on the Snakefile state.** Stages that were already confirmed working are
commented out (rules and their `rule all` targets), so the active part currently covers
only the later stages. To run everything from scratch, uncomment the rules you need and
the matching targets in `rule all`. Any script a rule calls is declared as an `input:`,
so editing a script triggers a rerun automatically.

Logs are written to `PRJNA559322/logs/` (per-rule logs under `{sample}/`, Slurm
stdout/stderr under `{rule}/`).

## What it found (3 strains, 2 substrates)

Strains: 13E0634 (persistent colonizer), 12E0115 (unknown), 13E0659 (sporadic).
Substrates: glyoxylic acid and L-rhamnose, the two substrates Barth et al. identified
as best separating persistent from sporadic colonizers.

A Biolog result counts as positive when AUC minus the negative-control AUC is at least
500 (the paper's own rule).

| Strain | Glyoxylate (after curation) | Rhamnose |
|---|---|---|
| 13E0634 (persistent) | match | **mismatch**, explained by a truncated *rhaS* |
| 12E0115 (unknown) | match | match |
| 13E0659 (sporadic) | match | match |

- **Glyoxylate was a model gap, not biology.** The initial no-growth prediction was
  identical in all three strains, including those with an intact *glcA*. CarveMe's
  reaction universe (built from the BiGG database) contains no extracellular glyoxylate
  exchange reaction, so no genome could ever produce one. `scripts/curate_model.py` adds
  `glx_e`, a proton-symport transport reaction (`GLXt2r`) and `EX_glx_e`, applied to
  every strain and based on genomic evidence, not tuned to any strain's phenotype. A
  flux trace (`scripts/verify_curated_model.py`) shows growth routes through native
  reactions (`GLXCL`, then the tartronate-semialdehyde pathway), and `GLXt2r` is mass
  balanced.
- **Rhamnose is a real gene-level signal.** The model predicts growth in all three
  strains, but 13E0634 does not grow on rhamnose in the experiment. The operon check
  finds *rhaS* (the pathway's transcriptional activator) truncated in that strain only
  (100% identity, 53% coverage), with all structural genes intact. FBA has no
  representation of gene regulation, so it cannot see this. Barth et al. independently
  report the same *rhaS* frameshift in O26 persistent strains.

## Caveats

- Pilot scope: three strains and two substrates; the Biolog table covers 28 strains and
  380 substrates.
- CarveMe's `carve` and `gapfill` support CPLEX, Gurobi or SCIP, **not GLPK**. SCIP is
  used here; the COBRApy stages use GLPK.
- The single-carbon-source FBA medium leaves CO2 open for uptake. Real cells fix a small
  amount of CO2 (anaplerotic reactions), so neither fully open nor fully closed is exact;
  capping uptake at a small value would be more faithful. Not yet changed.
- The curation applies one uncapped transport bound to every strain, regardless of how
  divergent *glcA* is (63% identity in 13E0634 vs. about 99% in the others). Scaling the
  bound by sequence evidence is not implemented.
- The operon check classifies each reference gene as intact, disrupted or absent from
  identity and coverage (thresholds 70% and 80%) against the strain's Bakta proteins. It
  cannot distinguish a frameshift from a deletion or an assembly gap.

## References

- Barth SA, Weber M, Schaufler K, Berens C, Geue L, Menge C. Metabolic traits of bovine
  Shiga toxin-producing *Escherichia coli* (STEC) strains with different colonization
  properties. *Toxins* 2020;12(6):414. doi:10.3390/toxins12060414
- Machado D, et al. Fast automated reconstruction of genome-scale metabolic models for
  microbial species and communities. *Nucleic Acids Res* 2018 (gky537).
- Bakta: doi:10.1099/mgen.0.000685
- Tools: Snakemake, COBRApy, BUSCO, DIAMOND, SPAdes, fastp, SRA Toolkit, BiGG database.
- The per-rule conda environments and cluster-generic Slurm profile follow the pattern
  used in [snake_coli](https://gitlab.com/micweber/snake_coli).

## Author and license

Developed by Suvratha Jayaprasad, Postdoc at Friedrich-Loeffler-Institut (FLI) Jena.
