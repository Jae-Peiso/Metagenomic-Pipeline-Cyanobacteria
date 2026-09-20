# Metagenome (ONT) Pipeline — Project Scaffold

Long-read (ONT) metagenomics of freshwater cyanobacterial blooms — assembly, differential-coverage binning, and 16S profiling built around separating true cyanobacteria from algal chloroplasts.

## ⚠️ Data-sharing restriction

`notebooks/sample_map.ipynb` and the exported figures in `notebooks/Maps/` plot all field sites, including a small number collected without a valid collection permit. Those sites are marked as a separate category on the maps (not grouped with the three field expeditions) and are recorded as `No Permit -Destroyed` in `Cyano_FieldSamplesLog_All.xlsx`. **Sequencing reads and any other data tied to those samples must not be published or redistributed.** They are left visible here, clearly flagged, for internal tracking only.

## ⚠️ Methodological Flags

### FLAG 1 — Kraken2 community profiles are unreliable for these samples
Kraken2 k-mer classification on metagenome-assembled contigs produced **incorrect community profiles** for KY-Mam-100624-C and should not be used as primary taxonomy for any sample in this dataset.

**Root cause:** Low-coverage, unenriched metagenomes lack sufficient k-mer depth, and novel cyanobacteria lineages are underrepresented in the `kraken2-microbial-fatfree` database.

**When Kraken2 will work better:** After cyano enrichment — higher coverage across target genomes will make k-mer classification reliable.

### FLAG 2 — barrnap in `01_community_analysis.ipynb` was run on Kraken2-filtered contigs only
Barrnap was applied only to contigs pre-classified as cyanobacteria by Kraken2, not the full assembly. This compounds FLAG 1 — misclassified cyano contigs were invisible to barrnap.

**Fix:** In `02_all_samples_community.ipynb`, barrnap runs on the full polished assembly for all 10 samples. All 16S sequences are extracted and classified via BLAST with no Kraken2 pre-filtering.

### FLAG 3 — `01_community_analysis.ipynb` is a B+C pilot only
Cross-sample community analysis across all 10 samples lives in `02_all_samples_community.ipynb`.

### FLAG 4 — a provisional label leaked into results as an identification (RESOLVED)
Separating free-living cyanobacteria from algal chloroplasts was identified as the central
risk for this dataset **at the outset**. It is why `03a_contig_id_gene_context.ipynb` §3
exists as a flanking-gene test, and why `03b_16S_phylogeny.ipynb` was built around the
question "cyanobacteria vs. chloroplast origin" — that notebook's plan even predicted the
answer before the tree was run.

Contigs 225, 227, 427, 531, 1504 and 3677 from KY-Mam-100624-C and MN-Col-091924-B were
tagged `'Novel cyanobacterium'` as a **working marker**, meaning "16S falls in SILVA's
chloroplast bin with no named host — identity unresolved, keep these samples in view."

**The defect was that the provisional tag got consumed downstream as a positive
identification.** `is_cyano()` counted it as cyanobacterial, the tree viz relabelled SILVA's
`Chloroplast` node to `Novel cyanobacteria`, and the prose in `02` §10 and
`00_driver_2.0.ipynb` §13 hardened it into a claimed finding: a novel Oscillatoriales genus
shared between a Kentucky cave and a Minnesota lake. The purpose-built tests then answered
the question the tag was holding open — **they are chloroplast 16S from green algae.**

The classification logic now resolves identity rather than deferring it, and the counts and
conclusions in `02_all_samples_community.ipynb` §10 are updated.

**Where the tag hardened:** `00_driver_2.0.ipynb` §13 ran the 16S BLAST against NCBI nt
**with a Cyanobacteria filter**. In NCBI taxonomy a chloroplast sequence is filed under its
eukaryotic host, not under Cyanobacteria, so that filter makes a chloroplast hit impossible
to return — every plastid query is forced onto its nearest free-living cyanobacterial
relative, and the resulting 94–99% identity reads as a novel lineage. Because that BLAST
could not express "this is a plastid", it converted an open question into a false positive.
The SILVA labelling bug below then propagated the same conclusion through notebooks 02 and 03.

**Root cause:** SILVA files plastid SSU sequences in two places, and both are organelles:

| SILVA lineage | What it is |
|---|---|
| `Bacteria;Cyanobacteriota;Cyanobacteriia;Chloroplast;…` | SILVA's chloroplast bin |
| `Eukaryota;Archaeplastida;Chloroplastida;…` | the same organelles, on the host's branch |

The first sits under Cyanobacteriota because plastids **descend from** cyanobacteria by
primary endosymbiosis — it encodes evolutionary origin, not free-living status. When the
deepest rank is a placeholder (`uncultured`, `metagenome`, `Incertae Sedis`) it means the
**host alga is unidentified in SILVA**, not that the sequence is a novel cyanobacterium.
`utils.parse_title()` labelled that case `'Novel cyanobacterium'`, and every downstream
consumer inherited it. The 93.6–97.5% identities that read as excitingly divergent are just
the normal distance from plastid 16S to cyanobacterial references.

**Evidence (three independent lines):**
- **Phylogeny** — `results/14_phylo_16S/all_16S.treefile`, rebuilt with a proper reference set
  (IQ-TREE, 1000 UFBoot; **30 taxa: 6 query + 7 sample contigs, 10 cyanobacterial and 7
  chloroplast references**). The original tree had 21 taxa and only 2 plastid references,
  neither of them a taxon our data actually hit. All six query contigs now fall **inside** the
  Chlorophyta plastid radiation — and not as a clade of their own: ***Desmodesmus abundans***
  chloroplast sits **among** them, sister to contig_531 at **94.2/98**. Supporting splits:
  (contig_225, contig_227) 100/100; (contig_1504, contig_3677) 94.2/98; that group with
  contig_531 + *Desmodesmus* at 94.6/87; the whole assembly with *Coelastrella saipanensis*
  plastid at **100/100**, then *Bracteacoccus* 84.1/87, (*Chlamydomonas*, *Mychonastes*)
  100/100, *Chlorella vulgaris* 99.6/100, and *Cyanophora paradoxa* 80.7/92.
  Meanwhile every sample contig pairs with its own cyanobacterial reference — contig_3535 with
  *S. elongatus* PCC 7942 at 100/100, contig_134 with *Nodosilinea nodulosa* at 99.6/100,
  contig_35/112 with *Limnothrix redekei* at 100/100, contig_117 with *Spirulina major* at
  98.7/99. *Loriellopsis cavernicola* is in the tree and no query contig is near it.
- **Genomic context** — `03a_contig_id_gene_context.ipynb` §3 flanking-gene BLASTp: `Bacterial-like hits: 0
  | Plastid/eukaryotic hits: 3` for both samples. The genes adjacent to the 16S are psbC (PSII
  44 kDa, 92–96% id to *Chlamydomonas* / *Chlorella* chloroplast), atpF/atpH (ATP synthase CF0,
  92–99% id to *Desmodesmus* / *Bracteacoccus*), and plastid LAGLIDADG intron ORFs. Independent
  of any 16S database bias.
- **Genome architecture** — contig_1686 (78 kb) and contig_227 (72.7 kb) each carry the 16S in
  **two opposite-orientation copies** tens of kb apart: the chloroplast inverted repeat. Read
  against plastid genomes (103–204 kb in these lineages) rather than 5–8 Mb cyanobacterial
  chromosomes, these contigs are ~65–76% complete **organelle genomes**, not 1–2% fragments.

**Fix applied:**
- `utils.py` — added `is_plastid_lineage()`, `is_cyano_lineage()` and `silva_lineage()`.
  `parse_title()` now returns `'Chloroplast (unidentified host)'` where it previously
  returned `'Novel cyanobacterium'`. Classification keys off the **lineage**, never off the
  parsed label string.
- `02_all_samples_community.ipynb` — classification cell, `cyano_hits` filter, tree
  `ORDER_REMAP` (which relabelled `Chloroplast` → `Novel cyanobacteria`), and the bubble-plot
  `get_category()` all corrected. Findings summary rewritten with SILVA-derived counts.
- `03a_contig_id_gene_context.ipynb` — `is_novel_cyano()` renamed `is_unnamed_plastid_bin()`. The
  contig selection is unchanged, since it defines exactly the set the evidence covers.
- `03b_16S_phylogeny.ipynb` — reference fetch rewritten. The old query
  `'<organism>[orgn] 16S ribosomal RNA[title] chloroplast[filter]'` returned NOT FOUND for
  11 of 13 plastids and 7 of 16 cyanobacteria, because plastid rrn16 is a **gene inside a
  chloroplast genome record**, not a standalone titled entry — and when a genome record did
  come back, the length filter rejected it for being ~150 kb. `fetch_16S()` is now two-tier:
  standalone record first, then extract the annotated 16S rRNA feature from a complete genome
  record. Now **17/17 cyanobacteria and 15/15 plastids**. Also fixed: barrnap emits several
  partial 16S fragments per contig and all were written with identical FASTA headers
  (contig_3677 ×4) — only the longest is kept now; and the 800 bp alignment floor silently
  dropped contig_427 (775 bp) and contig_3677 (787 bp), so two of six query contigs never
  reached the tree. Floor lowered to 700 bp.
- `00_driver_2.0.ipynb` — §13 rewritten with the corrected assignments; §9 and §14 carry
  correction banners (see below).

**Downstream corrections in `00_driver_2.0.ipynb`:**
- **§13** — the origin. Corrected assignments from unfiltered SILVA: contig_1652 →
  ***Chlorella pyrenoidosa*** chloroplast (98.7%, 1473 bp), contig_1686 → ***Mychonastes
  jurisii*** chloroplast (99.7%, 1484 bp), contig_225/227 → chloroplast, unidentified host
  (93.6%). The "uncultured cyanobacterium Dpcom212" clones that contig_1686 matched at 100%
  are themselves plastid sequences deposited under a cyanobacterial name — in SILVA they sit
  in the chloroplast bin too. The *Loriellopsis* and "first genome for the Dpcom clade" claims
  do not hold.
- **§9 — KaiABC absence.** The absence is real; the interpretation is not. The three
  "well-covered cyanobacterial contigs" searched are plastid genomes, and **plastids do not
  carry kaiABC** — it was lost during endosymbiotic genome reduction. Zero Kai hits in a
  chloroplast genome is the expected result, not cave-adapted gene loss under relaxed
  selection. That hypothesis is untested by these data. The KaiC hits in Rhizobiaceae are
  unaffected and still stand.
- **§14 — completeness and FastANI.** Nothing here was measured wrong. The cyanobacterial
  genome size used as the completeness denominator was a deliberate fixed benchmark, applied
  uniformly so recovery was comparable across samples while the organisms' identity was still
  the open question. Now that identity is settled, the same lengths carry a second reading:
  against green algal plastid genomes (103–204 kb) these are ~65–76% complete organelle
  genomes rather than ~2% fragments. The ~74–76% ANI against three cyanobacterial genomes is
  likewise exactly what a chloroplast gives, so it does not indicate a novel cyanobacterial
  genus, and the "at least 3 novel lineages" summary does not hold. Same for
  `03a_contig_id_gene_context.ipynb` cell 8, which uses the same benchmark and now carries both readings.

**Effect on counts:** KY-Mam-100624-C 5 → **0** cyanobacterial 16S (18/39 plastid);
MN-Col-091924-B 10 → **2** (10/21 plastid). All other samples unchanged.

**Not affected:** MN-Col-091924-B's *Synechococcus* (contig_3535) is genuine — zero branch
length to PCC 7942 — so `04_env7942_assembly.ipynb` stands. The *Limnothrix*, *Cyanobium*,
*Cyanodictyon*, *Nodosilinea* and Phormidiaceae assignments are unaffected.

**Still open:**
- **`contig_3345` (KY-Mam-100624-C2) is a plastid that reached the tree as a confirmed-cyano
  anchor.** Its SILVA top hit is the chloroplast bin at 97.0% and it sat at zero branch length
  from contig_531 inside the *Desmodesmus* clade. It is now excluded in `03b` cell 4
  (`MISCLASSIFIED_ANCHORS`), but the root cause remains: it came from the stale
  `results/10_tax_profile/all_samples_cyano_16S_outgroup.fna`, built with the old
  classification. Regenerate that file from the corrected `02_all_samples_community.ipynb`
  and re-check the remaining sample contigs the same way, then drop the exclusion.
- Four green algal plastids could not be retrieved (*Micractinium conductrix*,
  *Auxenochlorella pyrenoidosa*, *Chlorella pyrenoidosa*, *Monoraphidium neglectum*) — their
  GenBank records carry no annotated `rrn16` feature, so neither fetch tier can reach them.
  *Chlorella vulgaris* / *variabilis* cover the same class. Not blocking.
- The corrected notebooks have **not been re-executed end-to-end**; cells whose source changed
  had their stale outputs cleared. `results/14_phylo_16S/` *has* been regenerated (references,
  alignment, trimmed alignment, tree). The figures in `results/figures/` still reflect the old
  classification and need a rerun of `02_all_samples_community.ipynb`.

**Generalisable lesson:** barrnap does not distinguish bacterial from organellar 16S. On any
photosynthetic community, filter plastids explicitly and on lineage, and include chloroplast
references in any cyanobacterial phylogeny.

---

## Quick start
1) Copy `config/samples.tsv.example` to `config/samples.tsv` and edit it with your sample IDs and FASTQ paths.
2) Copy `config/config.yaml.example` to `config/config.yaml` and fill in your NCBI Entrez email/API key (see comments in the file).
3) Edit `config/medaka_models.yml` to match your ONT basecalling model(s).
4) Create the conda env:
   conda env create -f envs/base.yml
   conda activate meta-ont
5) Dry-run Snakemake:
   snakemake -n
6) Run:
   snakemake --cores 16

Outputs will appear under `results/` (QC → trim → filter → assembly → polish).

## What is and isn't in this repo

`results/` is **not tracked**. Every stage in it is reproducible from the Snakemake workflow
plus the notebooks, so a clone rebuilds it rather than carrying ~120 MB of Kraken2 tables,
reference genomes and figures. Expect a fresh clone to have no `results/` directory at all
until you run the pipeline.

Two things this means in practice:

- **Figures and tables are build artifacts.** Run `02_all_samples_community.ipynb` and
  `03a` / `03b` to regenerate them. They are not shipped.
- **Manual, non-regenerable inputs live in `data/manual/`, which IS tracked.** Anything no
  code can fetch belongs there, not in `results/`. Currently that is
  `data/manual/WB7RKP8C014-Alignment.txt`, an NCBI web-BLAST download read by
  `03a_contig_id_gene_context.ipynb` §3b. If you add another hand-downloaded file, put it
  there or it will be lost on the next clone.

Also untracked: `db/` (SILVA, Kraken2 — publicly re-downloadable), `config/config.yaml` and
`config/samples.tsv` (machine-specific paths and an API key; copy the `.example` templates).

## Apple Silicon (macOS arm64) notes
- Use `envs/macos_arm64_min.yml` for local runs through polishing (NanoPlot → Porechop → Filtlong → Flye → Racon). 
- For **Medaka** and most **binning/QC** tools (MetaBAT2/MaxBin2/CheckM2/BUSCO/GTDB-Tk), prefer Linux (HPC/VM/Docker). 
  - Option A (recommended): run the whole pipeline on a Linux server using the provided `envs/linux_x86_full.yml`. 
  - Option B: run steps 0–4 locally on macOS, then copy `consensus/contigs` to Linux for binning/QC.
  - Option C: use Docker (e.g., colima) with biocontainers or `ontresearch/medaka` images for the missing tools.

The `notebooks/00_driver_2.0.ipynb` can orchestrate either local runs or remote Linux runs via SSH.
