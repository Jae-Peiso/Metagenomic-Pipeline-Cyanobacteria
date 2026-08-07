# --- Metagenome Snakemake Pipeline (ONT long-read, reference-free) ---
import yaml
# After SAMPLES/FASTQS are loaded
import os
os.environ.setdefault("CONDA_SUBDIR", "osx-64")
configfile: "config/medaka_models.yml" # Can change main config file if needed

# ---- Upstream rerun guard ------------------------------------------------
# LOCK_UPSTREAM=1 (default): protect expensive upstream outputs from being overwritten
# LOCK_UPSTREAM=0: allow reruns/overwrites if you intentionally want to regenerate
LOCK_UPSTREAM = os.getenv("LOCK_UPSTREAM", "1").lower() not in ("0", "false", "no", "off")

def maybe_protected(path):
    """Wrap an output path as protected() when upstream outputs should not be overwritten."""
    return protected(path) if LOCK_UPSTREAM else path
# -------------------------------------------------------------------------

# ---- PATHS (single source of truth) ---------------------------------
DIR = {
    "flt":       "results/02_filter",
    "polish":    "results/04_polish",
    "map":       "results/05_map",
    "cov":       "results/06_cov",
    "bins":      "results/07_bins",
    "refine":    "results/08_bin_refine",
    "qc_bins":   "results/09_qc_bins",
    "tax":       "results/10_tax_profile",
}

def medaka_consensus(s):   return f"{DIR['polish']}/{s}.medaka/consensus.fasta"
def cov_mmi(s):            return f"{DIR['cov']}/{s}/contigs.mmi"
def cov_bam(s,t):          return f"{DIR['cov']}/{s}/maps/{t}.bam"
def cov_bai(s,t):          return f"{DIR['cov']}/{s}/maps/{t}.bam.bai"
def cov_depth(s):          return f"{DIR['cov']}/{s}/depth.tsv"
def bins_metabat2(s):      return f"{DIR['bins']}/{s}/metabat2"
def bins_maxbin2(s):       return f"{DIR['bins']}/{s}/maxbin2"
def refine_dastool(s):     return f"{DIR['refine']}/{s}/DASTool"
def refine_bins_dir(s):    return f"{DIR['refine']}/{s}/bins"
def bin_fa(s,b):           return f"{DIR['refine']}/{s}/bins/{b}.fa"
def bin_polished(s,b):     return f"{DIR['qc_bins']}/{s}/polished/{b}.medaka.fasta"
def tax_contigs_rep(s):    return f"{DIR['tax']}/{s}/contigs.kraken.report"
def tax_contigs_tsv(s):    return f"{DIR['tax']}/{s}/contigs.kraken.tsv"
def tax_bin_rep(s,b):      return f"{DIR['tax']}/{s}/bins/{b}.kraken.report"
def tax_bin_tsv(s,b):      return f"{DIR['tax']}/{s}/bins/{b}.kraken.tsv"

# Which contigs feed the meta branch:
# TEMP: start meta branch from Racon2
META_CONTIGS = lambda s: f"{DIR['polish']}/{s}.racon2.fasta"

# ---- SAMPLE Lists -------------------------------------------------
# Source of truth: config/samples.tsv (2 columns: sample_id <TAB> fastq_path)
print("ENV FocusSAMPLES =", os.getenv("FocusSAMPLES"))
print("ENV COV_READ_SAMPLES =", os.getenv("COV_READ_SAMPLES"))
SAMPLES, FASTQS = [], {}
with open("config/samples.tsv") as f:
    next(f)
    for line in f:
        sid, fq = line.rstrip("\n").split("\t")
        SAMPLES.append(sid) 
        FASTQS[sid] = fq

# Targets to process (assemblies/contiugs): optionally restrict to focusSAMPLES = "ID1, ID2, ..." 
_env_focus = os.getenv("FocusSAMPLES", "").strip() 
FOCUS = [s.strip() for s in _env_focus.split(",") if s.strip()] if _env_focus else []
if FOCUS: 
    unknown = sorted(set(FOCUS) - set(SAMPLES))
    if unknown:
        raise ValueError(f"Requested target(s) from 'FocusSAMPLE' not found in config/samples.tsv: {unknown}")
TARGET_SAMPLES = FOCUS if FOCUS else SAMPLES

# Read sets used for coverage estimation: defaults to *all* samples unless overridden
_env_cov = os.getenv("COV_READ_SAMPLES", "").strip()
COV_READ_SAMPLES = (
    [s.strip() for s in _env_cov.split(",") if s.strip()]
    if _env_cov else list(SAMPLES)
)

# Optional validation for COV_READ_SAMPLES
unknown_cov = sorted(set(COV_READ_SAMPLES) - set(SAMPLES))
if unknown_cov: 
    raise ValueError(f"Requested target(s) from 'COV_READ_SAMPLES' not found in config/samples.tsv: {unknown_cov}")

print("TARGET_SAMPLES:", TARGET_SAMPLES)
print("COV_READ_SAMPLES:", COV_READ_SAMPLES) 

# =====A Rule to Rule Them All=================================================================================================================

# Caution: this will create |TARGET_SAMPLES| × |COV_READ_SAMPLES| BAMs. If you focus to one 
# target but leave COV_READ_SAMPLES as all samples (defualt), you’ll get the “1×N” behavior.

rule all:
    input:
        # expand("results/04_polish/{s}.medaka/consensus.fasta", s=TARGET_SAMPLES),
        expand("results/06_cov/{s}/contigs.mmi", s=TARGET_SAMPLES),
        expand("results/06_cov/{s}/maps/{t}.bam", s=TARGET_SAMPLES, t=COV_READ_SAMPLES),
        expand("results/06_cov/{s}/depth.tsv", s=TARGET_SAMPLES),
        expand("results/07_bins/{s}/metabat2", s=TARGET_SAMPLES),
        # expand("results/07_bins/{s}/maxbin2",  s=TARGET_SAMPLES),      # OFF on macOS    
        # expand("results/08_bin_refine/{s}/DASTool", s=TARGET_SAMPLES), # OFF on macOS
        # expand("results/08_bin_refine/{s}/bins",    s=TARGET_SAMPLES), # OFF on macOS
        expand("results/10_tax_profile/{s}/contigs.kraken.report", s=TARGET_SAMPLES),
        expand("results/10_tax_profile/{s}/contigs.kraken.tsv",    s=TARGET_SAMPLES)
# =============================================================================================================================
# MAIN BRANCH (ASSEMBLY + POLISHING)
# ------- helpers -------

def medaka_model(sample):
    import yaml
    with open("config/medaka_models.yml") as yf:
        cfg = yaml.safe_load(yf)
    return cfg.get("models", {}).get(sample, cfg.get("default_model", "r1041_e82_400bps_sup"))

# =============================================================================================================================

# ------- ONLY these touch the original FASTQ paths -------
rule qc_nanoplot:
    input:
        fq = lambda wc: FASTQS[wc.s]
    output:
        rep = "results/00_qc/{s}/NanoPlot-report.html"
    threads: 4
    log:
        "logs/qc_nanoplot/{s}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/00_qc/{wildcards.s} logs/qc_nanoplot
        NanoPlot --fastq {input.fq} -o results/00_qc/{wildcards.s} --threads {threads} > {log} 2>&1
        """

rule qc_nanoplot_afterFilt:
    input:
        fq = "results/02_filter/{s}.flt.fq.gz"
    output:
        rep = "results/00_qc/{s}_afterFilt/NanoPlot-report.html"
    threads: 4
    log:
        "logs/qc_nanoplot_afterFilt/{s}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/00_qc/{wildcards.s}_afterFilt logs/qc_nanoplot_afterFilt
        NanoPlot --fastq {input.fq} -o results/00_qc/{wildcards.s}_afterFilt --threads {threads} > {log} 2>&1
        """

rule trim:
    input: lambda wc: FASTQS[wc.s]
    output: maybe_protected("results/01_trim/{s}.trim.fq.gz")
    threads: 8
    run:
        import shutil
        exe = shutil.which("porechop_abi") or shutil.which("porechop")
        if exe:
            shell(f"{exe} -i {input} -o {output} --threads {threads}")
        else:
            # no-op if no trimmer available (recent ONT often adapter-trimmed already)
            shell(f"ln -sf $(realpath {input}) {output} || cp {input} {output}")

rule filt:
    input: "results/01_trim/{s}.trim.fq.gz"
    output: maybe_protected("results/02_filter/{s}.flt.fq.gz")
    threads: 4
    shell: "filtlong --min_length 1000 --keep_percent 90 {input} | gzip > {output}"

# ------- From here on, consume explicit outputs (never FASTQS) -------
rule assemble_flye:
    input: "results/02_filter/{s}.flt.fq.gz"
    output: maybe_protected("results/03_assembly/{s}.flye/assembly.fasta")
    threads: 16
    params: 
        outdir = lambda wc: f"results/03_assembly/{wc.s}.flye"
    shell: 
         "flye --nano-raw {input} --out-dir {params.outdir} --meta --threads {threads} --min-overlap 2000"
    # Flye writes assembly.fasta to outdir; we declare that file as the rule output.


# Use unique suffixes to avoid wildcard collisions:
# --- FIRST MAPPING: make PAF for racon + BAM for indexing ---
rule map0_self_paf:
    input:
        asm="results/03_assembly/{s}.flye/assembly.fasta",
        reads="results/02_filter/{s}.flt.fq.gz"
    output:
        paf=maybe_protected("results/05_map/{s}.self.paf")
    threads: 8
    log:
        "logs/map0_self/{s}.paf.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/05_map logs/map0_self
        # default minimap2 output is PAF (perfect for racon)
        minimap2 -x map-ont -t {threads} {input.asm} {input.reads} > {output.paf} 2> {log}
        """

rule map0_self_bam:
    input:
        asm="results/03_assembly/{s}.flye/assembly.fasta",
        reads="results/02_filter/{s}.flt.fq.gz"
    output:
        bam=maybe_protected("results/05_map/{s}.self.bam"),
        bai=maybe_protected("results/05_map/{s}.self.bam.bai")
    threads: 8
    log:
        "logs/map0_self/{s}.bam.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/05_map logs/map0_self
        minimap2 -ax map-ont -t {threads} {input.asm} {input.reads} 2> {log} \
          | samtools view -b - \
          | samtools sort -@4 -o {output.bam}
        samtools index {output.bam} >> {log} 2>&1
        """
rule racon1:
    input:
        reads="results/02_filter/{s}.flt.fq.gz",
        ovl="results/05_map/{s}.self.paf",
        asm="results/03_assembly/{s}.flye/assembly.fasta"
    output:
        maybe_protected("results/04_polish/{s}.racon1.fasta")
    threads: 8
    log:
        "logs/racon1/{s}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/04_polish logs/racon1
        racon -m 8 -x -6 -g -8 -w 500 -t {threads} {input.reads} {input.ovl} {input.asm} > {output} 2> {log}
        """

# --- SECOND MAPPING against racon1 contigs: again PAF + BAM ---
rule map1_r1_paf:
    input:
        asm="results/04_polish/{s}.racon1.fasta",
        reads="results/02_filter/{s}.flt.fq.gz"
    output:
        paf=maybe_protected("results/05_map/{s}.r1.paf")
    threads: 8
    log:
        "logs/map1_r1/{s}.paf.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/05_map logs/map1_r1
        minimap2 -x map-ont -t {threads} {input.asm} {input.reads} > {output.paf} 2> {log}
        """

rule map1_r1_bam:
    input:
        asm="results/04_polish/{s}.racon1.fasta",
        reads="results/02_filter/{s}.flt.fq.gz"
    output:
        bam=maybe_protected("results/05_map/{s}.r1.bam"),
        bai=maybe_protected("results/05_map/{s}.r1.bam.bai")
    threads: 8
    log:
        "logs/map1_r1/{s}.bam.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/05_map logs/map1_r1
        minimap2 -ax map-ont -t {threads} {input.asm} {input.reads} 2> {log} \
          | samtools view -b - \
          | samtools sort -@4 -o {output.bam}
        samtools index {output.bam} >> {log} 2>&1
        """

rule racon2:
    input:
        reads="results/02_filter/{s}.flt.fq.gz",
        ovl="results/05_map/{s}.r1.paf",
        asm="results/04_polish/{s}.racon1.fasta"
    output:
        maybe_protected("results/04_polish/{s}.racon2.fasta")
    threads: 8
    log:
        "logs/racon2/{s}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/04_polish logs/racon2
        racon -m 8 -x -6 -g -8 -w 500 -t {threads} {input.reads} {input.ovl} {input.asm} > {output} 2> {log}
        """
from pathlib import Path

rule medaka_polish:
    input:
        reads = "results/02_filter/{s}.flt.fq.gz",
        asm   = "results/04_polish/{s}.racon2.fasta"
    output:
        fa  = medaka_consensus("{s}"),
        vcf=maybe_protected("results/04_polish/{s}.medaka/consensus.vcf.gz")
    params:
        outdir = lambda wc, output: str(Path(output.fa).parent),
        bacteria_flag = lambda wc: "--bacteria" if config.get("medaka", {}).get("bacteria", False) else ""
    threads: 8
    conda: "envs/medaka.yml"
    log: "logs/medaka/{s}.log"
    shell: r"""
        set -euo pipefail
        mkdir -p {params.outdir} "$(dirname {output.vcf})" "$(dirname {log})"
        medaka_consensus -i {input.reads} -d {input.asm} -o {params.outdir} -t {threads} {params.bacteria_flag} > {log} 2>&1
        bgzip -f {params.outdir}/consensus.vcf
        tabix -f {params.outdir}/consensus.vcf.gz
    """

# =============================================================================================================================
# META BRANCH (DIFF COVERAGE)
# =============================================================================================================================
# This bins each sample’s Medaka contigs using coverage across all 10 samples (differential coverage). 
# It then reconciles bins with DAS_Tool and produces a final non-redundant set of bins per sample.
# =============================================================================================================================

# META_CONTIGS = medaka_consensus  # <- if your medaka path differs, update here

# 1) Index each sample's contigs for fast mapping
rule idx_contigs:
    input:
        fa = META_CONTIGS("{s}") # <- triggers medaka_polish if needed
    output:
        mmi = "results/06_cov/{s}/contigs.mmi"
    threads: 2
    conda: "envs/minimap2.yml"
    shell: r"mkdir -p results/06_cov/{wildcards.s} && minimap2 -d {output.mmi} {input.fa}"

# 5) Map ALL samples' reads to EACH sample's contigs  (differential coverage)
#    Produces one BAM per (target assembly s, read-sample t).
# Purpose:
#   Create BAM files that quantify how much each sample's reads map to a target contig set.
#   For target contigs s and read sample t, the BAM encodes cov_{contig, t} for contigs in s.
# Why:
#   Binning benefits from "differential coverage": each contig gets a coverage vector across samples.
# Notes:
#   Off-diagonal mappings (t != s) are informative only if samples share genomes/lineages.
#   If samples are very different, consider restricting t to neighbors to reduce noise and compute.

rule map_reads_for_coverage:
    input:
        mmi = cov_mmi("{s}"),
        fq  = "results/02_filter/{t}.flt.fq.gz"
    output:
        bam = cov_bam("{s}","{t}"),
        bai = cov_bai("{s}","{t}")
    threads: 12
    conda: "envs/minimap2_samtools.yml"
    shell: r"""
        set -euo pipefail
        mkdir -p results/06_cov/{wildcards.s}/maps
        # allocate most threads to minimap2, rest to samtools sort
        mm2_t=$(( {threads} - 2 ))
        [ $mm2_t -lt 1 ] && mm2_t=1
        sam_t=$(( {threads} - mm2_t ))
        [ $sam_t -lt 1 ] && sam_t=1

        minimap2 -ax map-ont --secondary=no -t $mm2_t {input.mmi} {input.fq} \
        | samtools sort -@ $sam_t -m 1G -o {output.bam}
        samtools index {output.bam}
    """

# Convenience expansion for all s,t pairs
ALL_COV_BAMS = expand(
    "results/06_cov/{s}/maps/{t}.bam",
    s=TARGET_SAMPLES, t=TARGET_SAMPLES
)

# 6a) Summarize contig depths (MetaBAT2 helper) per target assembly s
# Purpose:
#   Convert one or more BAMs into a contig-by-sample coverage table.
# Output:
#   A table where each contig has a coverage vector across the chosen read samples.
# Why:
#   Binners cluster contigs using both composition and similarity of these coverage vectors
#   (often correlation/distance, i.e., "co-variation" across samples).
# Notes:
#   Use multiple BAMs (t over samples) for differential coverage; use only self BAM for debugging.

rule depth_table:
    input:
        fa   = META_CONTIGS("{s}"),
        bams = lambda wc: expand(cov_bam(wc.s,"{t}"), t = COV_READ_SAMPLES)
    output:
        tsv  = cov_depth("{s}")
    threads: 4
    # This methods call comes straight from the documentation for CoverM see https://wwood.github.io/CoverM/coverm-contig.html
    shell: r"""
    set -euo pipefail
    coverm contig --bam-files {input.bams} --methods metabat --threads {threads} --output-file {output.tsv}
    """

# 6b) Run multiple binners (MetaBAT2 + MaxBin2). Possibly add CONCOCT later if you like.
# Purpose:
#   Cluster contigs into draft genome bins using:
#     - sequence composition (k-mers), and
#     - differential coverage vectors from depth_table
# Why:
#   Contigs from the same genome tend to co-vary in abundance across samples.

rule bin_metabat2:
    input:
        fa  = META_CONTIGS("{s}"),
        dep = cov_depth("{s}")
    output:
        bins_dir = directory(bins_metabat2("{s}"))
    threads: 8
    conda: "envs/metabat2.yml"
    shell: r"""
        mkdir -p {output.bins_dir}
        metabat2 -i {input.fa} -a {input.dep} -o {output.bins_dir}/{wildcards.s} -t {threads}
    """
# bin_maxbin2 (optional)
# Purpose:
#   Alternative binner; can capture bins MetaBAT misses.
# Important:
#   MaxBin2 expects specific abundance inputs; ensure depth output is in the format it expects
#   (may require a conversion step from the coverm/metabat depth table).

rule bin_maxbin2:
    input:
        fa  = META_CONTIGS("{s}"),
        dep = cov_depth("{s}")
    output:
        bins_dir = directory(bins_maxbin2("{s}"))
    threads: 8
    conda: "envs/maxbin2.yml"
    shell: r"""
        mkdir -p {output.bins_dir}
        run_MaxBin.pl -thread {threads} -contig {input.fa} -abund {input.dep} -out {output.bins_dir}/{wildcards.s}
    """

# 7) Reconcile with DAS_Tool
# Purpose:
#   Refine/reconcile bins from multiple binners by scoring contig assignments and selecting a consensus.
# Why:
#   Different binners have different failure modes; DAS Tool often improves overall bin quality.

rule das_tool:
    input:
        fa   = META_CONTIGS("{s}"),
        mb2  = bins_metabat2("{s}"),
        mx2  = bins_maxbin2("{s}")
    output:
        out_dir = directory(refine_dastool("{s}"))
    threads: 8
    conda: "envs/dastool.yml"
    params:
        prefix = "{s}"
    shell: r"""
        mkdir -p {output.out_dir}
        DAS_Tool \
          -i {input.mb2},{input.mx2} \
          -l metabat2,maxbin2 \
          -c {input.fa} \
          -o {output.out_dir}/{params.prefix} \
          --threads {threads}
    """

# 8) Collect final reconciled bins (FASTA) from DAS_Tool
#    DAS_Tool usually writes *_DASTool_bins/ with multiple *.fa
rule collect_bins:
    input:
        dir = refine_dastool("{s}")
    output:
        bins_dir = directory(refine_bins_dir("{s}"))
    threads: 1
    conda: "envs/python.yml"
    shell: r"""
        mkdir -p {output.bins_dir}
        src=$(ls -d {input.dir}/*_DASTool_bins 2>/dev/null | head -n1)
        if [ -z "$src" ]; then echo "No DAS_Tool bins found for {wildcards.s}" >&2; exit 1; fi
        cp -r "$src"/* {output.bins_dir}/
    """

# 9) Per-bin Medaka polish using the original sample's reads
rule medaka_polish_bin:
    input:
        binfa = bin_fa("{s}","{b}"),
        reads = f"{DIR['flt']}/{{s}}.flt.fq.gz"
    output:
        fa    = bin_polished("{s}","{b}")
    params:
        outdir = f"{DIR['qc_bins']}/{{s}}/polished/{{b}}.medaka", 
        model  = os.environ.get("MEDAKA_MODEL", "r1041_e82_400bps_sup_v4.2.0")
    threads: 8
    conda: "envs/medaka.yml"
    log: "logs/medaka_bins/{s}_{b}.log"
    shell: r"""
        mkdir -p {params.outdir} logs/medaka_bins
        medaka_consensus -i {input.reads} -d {input.binfa} -o {params.outdir} -m {params.model} -t {threads} > {log} 2>&1
        cp {params.outdir}/consensus.fasta {output.fa}
    """
# =============================================================================================================================
# 10) Taxonomy on contigs (per sample; fast & always available)
rule kraken_contigs:
    input:
        fa = META_CONTIGS("{s}")
    output:
        rep = tax_contigs_rep("{s}"),
        tsv = tax_contigs_tsv("{s}")
    threads: 2
    conda: "envs/kraken2.yml"
    log: "logs/kraken_contigs/{s}.log"
    shell: r"""
        set -euo pipefail
        mkdir -p $(dirname {output.rep}) $(dirname {log})
        kraken2 --db db/kraken2-microbial-fatfree --threads {threads} \
            --report {output.rep} --output {output.tsv} \
            --use-names {input.fa} > {log} 2>&1
    """

# (Optional) Taxonomy on each reconciled bin FASTA
# Only include this in `rule all` once bins exist / are enumerated.
rule kraken_bins:
    input:
        fa = bin_fa("{s}","{b}")
    output:
        rep = tax_bin_rep("{s}","{b}"),
        tsv = tax_bin_tsv("{s}","{b}")
    threads: 2
    conda: "envs/kraken2.yml"
    log: "logs/kraken_bins/{s}_{b}.log"
    shell: r"""
        set -euo pipefail
        mkdir -p $(dirname {output.rep}) $(dirname {log})
        kraken2 --db db/kraken2-microbial-fatfree --threads {threads} \
            --report {output.rep} --output {output.tsv} \
            --use-names {input.fa} > {log} 2>&1
    """
# =============================================================================================================================