# Targeted Strain Assembly: MN-Col-091924-B vs. *S. elongatus* PCC 7942

## Background

This project analyzes 10 ONT metagenome samples processed through a Snakemake
pipeline (QC → trim → filter → Flye assembly → Racon polish → coverage →
MetaBAT2 binning). Polished assemblies live at `results/04_polish/{sample}.racon2.fasta`;
filtered reads live at `results/02_filter/{sample}.fastq.gz`.

One sample — MN-Col-091924-B (Columbia Lake, MN, ~168x coverage) — contains a
strain with a 16S rRNA sequence very closely related to *Synechococcus
elongatus* PCC 7942. This notebook performs a targeted assembly of that
strain's genome and compares it to the PCC 7942 lab reference to identify
indels and substitutions.

Prior to this analysis, the PCC 7942 16S sequence was already available at
`results/14_phylo_16S/cyano_refs_16S/`, and the GenBank record CP000100.1
(chromosome) had been downloaded to `results/ref/S_elongatus_PCC7942.gbk`.
The full reference genome FASTA (CP000100.1 + plasmid CP000101.1) did not yet
exist in the project and needed to be fetched.

## Approach (`notebooks/04_env7942_assembly.ipynb`)

1. **Fetch reference** — download the full PCC 7942 reference genome
   (CP000100.1 chromosome + CP000101.1 plasmid pANL) as FASTA via Biopython
   Entrez / NCBI efetch, cached to `results/ref/PCC7942_reference.fasta` so
   re-running skips the download.

2. **Map & extract** — map the MN-Col-091924-B filtered reads to the PCC 7942
   reference with minimap2 (`map-ont` preset), then extract only the mapped
   reads with samtools for reassembly.

3. **Targeted assembly** — assemble the extracted reads with Flye in
   `--nano-hq` mode (not metagenome mode, since this targets a single genome),
   output to `results/15_env7942/`.

4. **Polish** — two rounds of Racon polishing using the extracted reads.

5. **Variant calling** — align the polished assembly to the PCC 7942
   reference with minimap2 (`asm5`/`asm10` preset) and call variants with
   bcftools (`mpileup` → `call`) to produce a VCF of SNPs and indels, plus a
   summary table of variant counts by type.

6. **Synteny visualization** — a dot plot / whole-genome alignment figure
   (paftools or matplotlib) showing synteny between the environmental
   assembly and PCC 7942.

## Conventions

The notebook follows the style of the rest of the project: section headers
with `───` markers, tool-availability checks via `shutil.which()` at the top,
disk-cached results so cells are idempotent, and figures saved to
`results/figures/`. Shared utilities live in `notebooks/utils.py`
(`find_project_root`, `load_fasta`, `is_placeholder`, `parse_title`);
project configuration is in `config/config.yaml` (see `config/config.yaml.example`
for the template). Analysis was run in a `genomics_arm64` conda environment
with minimap2, samtools, flye, racon, bcftools, and the Biopython/pandas/matplotlib
stack.
