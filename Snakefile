# --- Metagenome Snakemake Pipeline (ONT long-read, reference-free) ---
# Active pipeline: QC → Trim → Filter → Assemble → Racon×2 → Coverage → MetaBAT2 → Kraken2
# See Snakefile_master.smk for full rule library (Medaka, MaxBin2, DAS_Tool, bin QC, etc.)
import os
os.environ.setdefault("CONDA_SUBDIR", "osx-64")

# ---- Upstream rerun guard ------------------------------------------------
# LOCK_UPSTREAM=1 (default): protect expensive upstream outputs from being overwritten
# LOCK_UPSTREAM=0: allow reruns/overwrites if you intentionally want to regenerate
LOCK_UPSTREAM = os.getenv("LOCK_UPSTREAM", "1").lower() not in ("0", "false", "no", "off")

def maybe_protected(path):
    return protected(path) if LOCK_UPSTREAM else path

# ---- PATHS ---------------------------------------------------------------
DIR = {
    "flt":    "results/02_filter",
    "polish": "results/04_polish",
    "map":    "results/05_map",
    "cov":    "results/06_cov",
    "bins":   "results/07_bins",
    "tax":    "results/10_tax_profile",
}

def cov_mmi(s):       return f"{DIR['cov']}/{s}/contigs.mmi"
def cov_bam(s, t):    return f"{DIR['cov']}/{s}/maps/{t}.bam"
def cov_bai(s, t):    return f"{DIR['cov']}/{s}/maps/{t}.bam.bai"
def cov_depth(s):     return f"{DIR['cov']}/{s}/depth.tsv"
def bins_metabat2(s): return f"{DIR['bins']}/{s}/metabat2"

# Contigs entering the meta branch (Racon2 output)
META_CONTIGS = lambda s: f"{DIR['polish']}/{s}.racon2.fasta"

# ---- SAMPLE LISTS --------------------------------------------------------
# Source of truth: config/samples.tsv (2 columns: sample_id <TAB> fastq_path)
print("ENV FocusSAMPLES =",    os.getenv("FocusSAMPLES"))
print("ENV COV_READ_SAMPLES =", os.getenv("COV_READ_SAMPLES"))

SAMPLES, FASTQS = [], {}
with open("config/samples.tsv") as f:
    next(f)
    for line in f:
        sid, fq = line.rstrip("\n").split("\t")
        SAMPLES.append(sid)
        FASTQS[sid] = fq

# Target assemblies: optionally restrict via FocusSAMPLES="ID1,ID2,..."
_env_focus = os.getenv("FocusSAMPLES", "").strip()
FOCUS = [s.strip() for s in _env_focus.split(",") if s.strip()] if _env_focus else []
if FOCUS:
    unknown = sorted(set(FOCUS) - set(SAMPLES))
    if unknown:
        raise ValueError(f"FocusSAMPLES not in samples.tsv: {unknown}")
TARGET_SAMPLES = FOCUS if FOCUS else SAMPLES

# Read sets for coverage estimation: defaults to all samples unless overridden
_env_cov = os.getenv("COV_READ_SAMPLES", "").strip()
COV_READ_SAMPLES = (
    [s.strip() for s in _env_cov.split(",") if s.strip()]
    if _env_cov else list(SAMPLES)
)
unknown_cov = sorted(set(COV_READ_SAMPLES) - set(SAMPLES))
if unknown_cov:
    raise ValueError(f"COV_READ_SAMPLES not in samples.tsv: {unknown_cov}")

print("TARGET_SAMPLES:  ", TARGET_SAMPLES)
print("COV_READ_SAMPLES:", COV_READ_SAMPLES)

# ==========================================================================
# Rule all
# NOTE: creates |TARGET_SAMPLES| × |COV_READ_SAMPLES| BAMs for coverage.
#       If COV_READ_SAMPLES is left as all samples but FocusSAMPLES is one,
#       you get 1×N BAMs. Restrict COV_READ_SAMPLES for related samples only.
# ==========================================================================
rule all:
    input:
        expand("results/06_cov/{s}/contigs.mmi",       s=TARGET_SAMPLES),
        expand("results/06_cov/{s}/maps/{t}.bam",      s=TARGET_SAMPLES, t=COV_READ_SAMPLES),
        expand("results/06_cov/{s}/depth.tsv",         s=TARGET_SAMPLES),
        expand("results/07_bins/{s}/metabat2",         s=TARGET_SAMPLES),
        expand("results/10_tax_profile/{s}/contigs.kraken.report", s=TARGET_SAMPLES),
        expand("results/10_tax_profile/{s}/contigs.kraken.tsv",    s=TARGET_SAMPLES),
        expand("results/10_tax_profile/{s}/all_16S.fna",           s=TARGET_SAMPLES),

# ==========================================================================
# QC
# ==========================================================================
rule qc_nanoplot:
    input:  fq  = lambda wc: FASTQS[wc.s]
    output: rep = "results/00_qc/{s}/NanoPlot-report.html"
    threads: 4
    log: "logs/qc_nanoplot/{s}.log"
    shell: r"""
        set -euo pipefail
        mkdir -p results/00_qc/{wildcards.s} logs/qc_nanoplot
        NanoPlot --fastq {input.fq} -o results/00_qc/{wildcards.s} --threads {threads} > {log} 2>&1
    """

rule qc_nanoplot_afterFilt:
    input:  fq  = "results/02_filter/{s}.flt.fq.gz"
    output: rep = "results/00_qc/{s}_afterFilt/NanoPlot-report.html"
    threads: 4
    log: "logs/qc_nanoplot_afterFilt/{s}.log"
    shell: r"""
        set -euo pipefail
        mkdir -p results/00_qc/{wildcards.s}_afterFilt logs/qc_nanoplot_afterFilt
        NanoPlot --fastq {input.fq} -o results/00_qc/{wildcards.s}_afterFilt --threads {threads} > {log} 2>&1
    """

# ==========================================================================
# TRIM + FILTER
# ==========================================================================
rule trim:
    input:  lambda wc: FASTQS[wc.s]
    output: maybe_protected("results/01_trim/{s}.trim.fq.gz")
    threads: 8
    run:
        import shutil
        exe = shutil.which("porechop_abi") or shutil.which("porechop")
        if exe:
            shell(f"{exe} -i {input} -o {output} --threads {threads}")
        else:
            shell(f"ln -sf $(realpath {input}) {output} || cp {input} {output}")

rule filt:
    input:  "results/01_trim/{s}.trim.fq.gz"
    output: maybe_protected("results/02_filter/{s}.flt.fq.gz")
    threads: 4
    shell: "filtlong --min_length 1000 --keep_percent 90 {input} | gzip > {output}"

# ==========================================================================
# ASSEMBLY
# ==========================================================================
rule assemble_flye:
    input:  "results/02_filter/{s}.flt.fq.gz"
    output: maybe_protected("results/03_assembly/{s}.flye/assembly.fasta")
    threads: 16
    params: outdir = lambda wc: f"results/03_assembly/{wc.s}.flye"
    shell:  "flye --nano-raw {input} --out-dir {params.outdir} --meta --threads {threads} --min-overlap 2000"

# ==========================================================================
# POLISH (Racon ×2)
# ==========================================================================
rule map0_self_paf:
    input:
        asm   = "results/03_assembly/{s}.flye/assembly.fasta",
        reads = "results/02_filter/{s}.flt.fq.gz"
    output:
        paf = maybe_protected("results/05_map/{s}.self.paf")
    threads: 8
    log: "logs/map0_self/{s}.paf.log"
    shell: r"""
        set -euo pipefail
        mkdir -p results/05_map logs/map0_self
        minimap2 -x map-ont -t {threads} {input.asm} {input.reads} > {output.paf} 2> {log}
    """

rule racon1:
    input:
        reads = "results/02_filter/{s}.flt.fq.gz",
        ovl   = "results/05_map/{s}.self.paf",
        asm   = "results/03_assembly/{s}.flye/assembly.fasta"
    output: maybe_protected("results/04_polish/{s}.racon1.fasta")
    threads: 8
    log: "logs/racon1/{s}.log"
    shell: r"""
        set -euo pipefail
        mkdir -p results/04_polish logs/racon1
        racon -m 8 -x -6 -g -8 -w 500 -t {threads} {input.reads} {input.ovl} {input.asm} > {output} 2> {log}
    """

rule map1_r1_paf:
    input:
        asm   = "results/04_polish/{s}.racon1.fasta",
        reads = "results/02_filter/{s}.flt.fq.gz"
    output:
        paf = maybe_protected("results/05_map/{s}.r1.paf")
    threads: 8
    log: "logs/map1_r1/{s}.paf.log"
    shell: r"""
        set -euo pipefail
        mkdir -p results/05_map logs/map1_r1
        minimap2 -x map-ont -t {threads} {input.asm} {input.reads} > {output.paf} 2> {log}
    """

rule racon2:
    input:
        reads = "results/02_filter/{s}.flt.fq.gz",
        ovl   = "results/05_map/{s}.r1.paf",
        asm   = "results/04_polish/{s}.racon1.fasta"
    output: maybe_protected("results/04_polish/{s}.racon2.fasta")
    threads: 8
    log: "logs/racon2/{s}.log"
    shell: r"""
        set -euo pipefail
        mkdir -p results/04_polish logs/racon2
        racon -m 8 -x -6 -g -8 -w 500 -t {threads} {input.reads} {input.ovl} {input.asm} > {output} 2> {log}
    """

# ==========================================================================
# COVERAGE (minimap2 index + per-sample BAMs + depth table)
# ==========================================================================
rule idx_contigs:
    input:  fa  = META_CONTIGS("{s}")
    output: mmi = "results/06_cov/{s}/contigs.mmi"
    threads: 2
    conda: "envs/minimap2.yml"
    shell: r"mkdir -p results/06_cov/{wildcards.s} && minimap2 -d {output.mmi} {input.fa}"

rule map_reads_for_coverage:
    input:
        mmi = cov_mmi("{s}"),
        fq  = "results/02_filter/{t}.flt.fq.gz"
    output:
        bam = cov_bam("{s}", "{t}"),
        bai = cov_bai("{s}", "{t}")
    threads: 12
    conda: "envs/minimap2_samtools.yml"
    shell: r"""
        set -euo pipefail
        mkdir -p results/06_cov/{wildcards.s}/maps
        mm2_t=$(( {threads} - 2 )); [ $mm2_t -lt 1 ] && mm2_t=1
        sam_t=$(( {threads} - $mm2_t ));  [ $sam_t -lt 1 ] && sam_t=1
        minimap2 -ax map-ont --secondary=no -t $mm2_t {input.mmi} {input.fq} \
          | samtools sort -@ $sam_t -m 1G -o {output.bam}
        samtools index {output.bam}
    """

rule depth_table:
    input:
        fa   = META_CONTIGS("{s}"),
        bams = lambda wc: expand(cov_bam(wc.s, "{t}"), t=COV_READ_SAMPLES)
    output:
        tsv = cov_depth("{s}")
    threads: 4
    shell: r"""
        set -euo pipefail
        COVERM=$(command -v coverm 2>/dev/null || echo "/Users/jaemac/miniforge3/envs/meta-ont-macos-arm64/bin/coverm")
        $COVERM contig --bam-files {input.bams} --methods metabat --threads {threads} --output-file {output.tsv}
    """

# ==========================================================================
# BINNING
# ==========================================================================
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

# ==========================================================================
# TAXONOMY (Kraken2 on contigs)
# ==========================================================================
rule kraken_contigs:
    input:  fa  = META_CONTIGS("{s}")
    output:
        rep = "results/10_tax_profile/{s}/contigs.kraken.report",
        tsv = "results/10_tax_profile/{s}/contigs.kraken.tsv"
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

# ==========================================================================
# 16S rRNA EXTRACTION (barrnap on full polished assembly — no Kraken2 pre-filter)
# Extracts ALL 16S sequences for downstream BLAST classification and community
# profiling in notebooks/02_all_samples_community.ipynb
# ==========================================================================
rule extract_16S:
    input:  fa  = META_CONTIGS("{s}")
    output:
        gff = "results/10_tax_profile/{s}/barrnap_full.gff",
        fna = "results/10_tax_profile/{s}/all_16S.fna"
    threads: 4
    log: "logs/extract_16S/{s}.log"
    shell: r"""
        set -euo pipefail
        mkdir -p $(dirname {output.gff}) logs/extract_16S
        barrnap --kingdom bac --threads {threads} --outseq {output.fna} \
            {input.fa} > {output.gff} 2> {log}
        # Keep only 16S sequences in the FASTA (barrnap --outseq includes all rRNA types)
        python3 -c "
import sys
keep = False
with open('{output.fna}') as f:
    for line in f:
        if line.startswith('>'):
            keep = '16S_rRNA' in line
        if keep:
            sys.stdout.write(line)
" > {output.fna}.tmp && mv {output.fna}.tmp {output.fna}
    """
