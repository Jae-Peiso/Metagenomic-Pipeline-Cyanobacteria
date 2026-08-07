"""
utils.py — shared utilities for the cyanobacteria metagenomics notebooks.

Bootstrap pattern (use in every notebook setup cell):
    import sys
    from pathlib import Path
    _p = Path.cwd().resolve()
    for _parent in [_p] + list(_p.parents):
        if (_parent / 'Snakefile').exists():
            sys.path.insert(0, str(_parent / 'notebooks'))
            break
    from utils import find_project_root, load_config, ...
"""

import json
import subprocess
from pathlib import Path
import re
import yaml  # pip/conda: pyyaml


# ── Project root detection ────────────────────────────────────────────────────

def find_project_root(start=None):
    """Walk up from start (default: cwd) until a directory containing Snakefile is found.
    Returns the project root as a Path. Raises FileNotFoundError if not found.
    """
    p = (start or Path.cwd()).resolve()
    for parent in [p] + list(p.parents):
        if (parent / 'Snakefile').exists():
            return parent
    raise FileNotFoundError('No Snakefile found — is the working directory inside the project?')


# ── Config loading ────────────────────────────────────────────────────────────

def load_config(project_root=None):
    """Load config/config.yaml from the project root.

    Returns a dict with resolved Path values where applicable.
    Falls back gracefully if config.yaml is absent (returns empty dict).
    """
    root = project_root or find_project_root()
    cfg_path = root / 'config' / 'config.yaml'
    if not cfg_path.exists():
        return {}
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f) or {}

    # Resolve silva_db: use config value if set, otherwise derive from project root
    if cfg.get('silva_db'):
        cfg['silva_db'] = str(Path(cfg['silva_db']).expanduser())
    else:
        cfg['silva_db'] = str(root / 'db' / 'silva' / 'silva_SSU_NR99')

    return cfg


# ── FASTA I/O ─────────────────────────────────────────────────────────────────

def load_fasta(path):
    """Parse a FASTA file and return a dict of {header: sequence}.
    Header is the full description line minus the leading '>'.
    """
    seqs = {}
    header, seq = None, []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if header:
                    seqs[header] = ''.join(seq)
                header, seq = line[1:].strip(), []
            else:
                seq.append(line.strip())
        if header:
            seqs[header] = ''.join(seq)
    return seqs


# ── SILVA taxonomy helpers ────────────────────────────────────────────────────

PLACEHOLDERS = {
    'metagenome', 'uncultured', 'unclassified', 'unknown',
    'environmental sample', 'environmental samples',
    'uncultured bacterium', 'uncultured organism',
    'Incertae Sedis', 'incertae sedis',
}


def is_placeholder(name):
    """True if a SILVA taxon name is generic and carries no useful identity information.
    Catches exact matches in PLACEHOLDERS plus any name starting with
    'uncultured' or 'unclassified'.
    """
    n = name.lower().strip()
    return (n in {p.lower() for p in PLACEHOLDERS}
            or n.startswith('uncultured')
            or n.startswith('unclassified'))


def silva_lineage(title):
    """Split a SILVA BLAST hit title into its list of taxon names.

    SILVA titles look like '<accession> Bacteria;Cyanobacteriota;...;Chloroplast;...'.
    Returns [] for non-lineage titles (no ';').
    """
    if ';' not in title:
        return []
    lineage_str = title.split(' ', 1)[-1] if ' ' in title else title
    return [t.strip() for t in lineage_str.split(';') if t.strip()]


def is_plastid_lineage(title):
    """True if a SILVA hit is an organellar (chloroplast/plastid) 16S sequence.

    SILVA files plastid SSU sequences in TWO places, and both are organelles:

      1. Bacteria;Cyanobacteriota;Cyanobacteriia;Chloroplast;...
         This is SILVA's chloroplast bin. It sits under Cyanobacteriota because
         plastids are DESCENDED from cyanobacteria by primary endosymbiosis —
         it encodes evolutionary origin, not that the sequence came from a
         free-living cell. A placeholder at the deepest rank ('uncultured',
         'Incertae Sedis', 'metagenome') means the HOST alga is unidentified in
         SILVA; it does NOT mean the sequence is a novel free-living cyanobacterium.

      2. Eukaryota;Archaeplastida;Chloroplastida;...
         The same organelles, filed on the eukaryotic host's branch.

    Treating case (1) as novel cyanobacteria is the error corrected here; see
    FLAG 4 in README.md and the tree in results/14_phylo_16S/all_16S.treefile,
    where every such sequence resolved as a green algal chloroplast.
    """
    taxa = silva_lineage(title)
    if not taxa:
        return False
    if taxa[0] == 'Eukaryota':
        return True
    return 'Chloroplast' in taxa or 'Plastid' in taxa


def parse_title(title):
    """Extract the most informative taxon label from a SILVA BLAST hit title.

    Rules applied in order:
    - Non-lineage titles (no ';'): parse genus/species from the accession label.
    - Plastid lineages with a named host  → 'Putative chloroplast'.
    - Plastid lineages with no named host → 'Chloroplast (unidentified host)'.
    - All others: walk lineage from deepest node upward, return first non-placeholder.

    Both plastid labels are organelles. Use is_plastid_lineage() rather than
    comparing against these strings when you need to filter plastids out.
    """
    if ';' not in title:
        m = re.search(r'\|\s*([A-Z][a-z]+(?:\s+[a-z]+)?)', title)
        if m:
            return m.group(1)
        clean = title.split('>')[-1].strip().split()
        return f'{clean[0]} {clean[1]}' if len(clean) >= 2 else title[:40]

    taxa = silva_lineage(title)

    if is_plastid_lineage(title):
        deep = taxa[-1].strip() if taxa else ''
        if deep and not is_placeholder(deep) and deep.lower() != 'incertae sedis':
            return 'Putative chloroplast'
        # Host alga not identified in SILVA — still a chloroplast, not a novel cyano.
        return 'Chloroplast (unidentified host)'

    for taxon in reversed(taxa):
        if not is_placeholder(taxon):
            return taxon
    return taxa[0] if taxa else title[:40]


CYANO_KEYWORDS = [
    'Cyanobacteria', 'Cyanobacteri', 'Cyanobacteriota', 'Nostoc', 'Anabaena',
    'Microcystis', 'Synechococcus', 'Oscillatoria', 'Loriellopsis',
    'Aphanizomenon', 'Planktothrix', 'Cylindrospermopsis', 'Lyngbya',
    'Phormidium', 'Gloeobacter', 'Trichodesmium', 'Prochlorococcus',
    'Cyanothece', 'Calothrix', 'Scytonema', 'Tolypothrix', 'Fischerella',
    'Hapalosiphon',
]


def is_cyano_lineage(title, keywords=None):
    """True if a SILVA hit is a genuine (non-organellar) cyanobacterium.

    Plastids are excluded first, so a hit is only counted as cyanobacterial when
    it carries a cyanobacterial name AND is not in either SILVA plastid bin.
    """
    if is_plastid_lineage(title):
        return False
    kws = keywords if keywords is not None else CYANO_KEYWORDS
    text = (parse_title(title) + ' ' + title).lower()
    return any(kw.lower() in text for kw in kws)


# ── BLAST / 16S taxonomy ──────────────────────────────────────────────────────

def blast_16S_local(sample, overwrite=False, tax_dir=None, silva_db=None):
    """BLASTn all barrnap-extracted 16S sequences for a sample against the local SILVA SSU NR99 DB.

    Results are cached as blast_16S_silva.json in the sample's tax_profile directory.
    Set overwrite=True to force a rerun even if the cache exists.

    Parameters
    ----------
    sample    : sample ID string (e.g. 'KY-Mam-100624-C')
    overwrite : if True, re-run BLAST even when cached JSON exists
    tax_dir   : path to results/10_tax_profile (default: derived from project root)
    silva_db  : SILVA SSU NR99 database prefix (default: from config.yaml or project root)

    Returns the path to the JSON file, or None if the input 16S FASTA is missing/empty.
    """
    if tax_dir is None:
        tax_dir = find_project_root() / 'results' / '10_tax_profile'
    if silva_db is None:
        cfg = load_config()
        silva_db = cfg['silva_db']
    tax_dir = Path(tax_dir)

    fna_path = tax_dir / sample / 'all_16S.fna'
    out_json = tax_dir / sample / 'blast_16S_silva.json'

    if out_json.exists() and not overwrite:
        print(f'{sample}: cached, skipping')
        return out_json
    if not fna_path.exists() or fna_path.stat().st_size == 0:
        print(f'{sample}: no 16S file, skipping')
        return None

    out_tsv = tax_dir / sample / 'blast_16S_silva.tsv'
    cmd = [
        'blastn', '-db', str(silva_db),
        '-query', str(fna_path),
        '-out', str(out_tsv),
        '-outfmt', '6 qseqid sseqid stitle pident length evalue',
        '-max_target_seqs', '5',
        '-num_threads', '4',
    ]
    print(f'{sample}: running local BLAST...')
    subprocess.run(cmd, check=True)

    results = {}
    if out_tsv.exists():
        for line in out_tsv.read_text().splitlines():
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            qid, sid, stitle, pident, length, evalue = parts
            results.setdefault(qid, []).append({
                'title':     stitle,
                'pct_id':    float(pident),
                'evalue':    float(evalue),
                'align_len': int(length),
            })

    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2)
    print(f'{sample}: done → {out_json}')
    return out_json


def load_blast_results(sample, tax_dir=None):
    """Load cached SILVA BLASTn results for a sample from blast_16S_silva.json.

    Parameters
    ----------
    sample  : sample ID string
    tax_dir : path to results/10_tax_profile (default: derived from project root)

    Returns an empty dict if no results file exists.
    """
    if tax_dir is None:
        tax_dir = find_project_root() / 'results' / '10_tax_profile'
    path = Path(tax_dir) / sample / 'blast_16S_silva.json'
    if not path.exists():
        return {}
    with open(path) as f:
        return json.load(f)


def top_hit_table(sample, tax_dir=None):
    """Build a per-sequence DataFrame of top SILVA hits for a sample.

    Columns: sample, query (barrnap ID), top_hit (parsed taxon label),
    pct_id, evalue, raw_title (original SILVA lineage string).
    Queries with no BLAST result get top_hit='No hit'.

    Parameters
    ----------
    sample  : sample ID string
    tax_dir : path to results/10_tax_profile (default: derived from project root)
    """
    import pandas as pd
    results = load_blast_results(sample, tax_dir=tax_dir)
    rows = []
    for query_id, hits in results.items():
        if hits:
            top = hits[0]
            rows.append({
                'sample':    sample,
                'query':     query_id,
                'top_hit':   parse_title(top['title']),
                'pct_id':    top.get('pct_id'),
                'evalue':    top.get('evalue'),
                'raw_title': top['title'],
            })
        else:
            rows.append({
                'sample':    sample,
                'query':     query_id,
                'top_hit':   'No hit',
                'pct_id':    None,
                'evalue':    None,
                'raw_title': '',
            })
    return pd.DataFrame(rows)
