# Notebook Structure Improvements — Summary & Migration Plan

## What Changed

### 1. Added Function Docstrings (Recommendation 4)
All utility functions now have docstrings explaining their purpose, parameters, and return values:

| Function | Where | What It Does |
|---|---|---|
| `find_project_root()` | `02`, `03a_contig_id_gene_context`, `utils.py` | Walk up directory tree to find Snakefile |
| `load_fasta()` | `02`, `03`, `utils.py` | Parse FASTA into `{header: sequence}` dict |
| `blast_16S_local()` | `02`, `utils.py` | BLASTn 16S against SILVA with caching |
| `is_placeholder()` | `02`, `03a_contig_id_gene_context`, `utils.py` | Check if SILVA taxon is generic (uncultured, etc.) |
| `parse_title()` | `02`, `utils.py` | Extract best taxon label from SILVA lineage |
| `gc()` | `03a_contig_id_gene_context`, `utils.py` | Calculate GC content (%) |
| `classify_hit()` | `03a_contig_id_gene_context` | Categorize BLASTp hits (cyanobacterial vs plastid vs other) |
| `parse_gff()` | `03a_contig_id_gene_context` | Parse Prodigal GFF → contig gene lists |

### 2. Added Tool Availability Checks (Recommendation 5)
Three new notebook cells (one per notebook) that:
- Check if required external tools exist on PATH: `blastn`, `mafft`, `barrnap`, etc.
- Raise a clear `EnvironmentError` with conda package names **before** any subprocess calls
- Prevent hours of silent failures deep in notebook execution

**Where checks are added:**
- `02_all_samples_community.ipynb`: Right after imports (cell 0) — checks `blastn`, `mafft`, `barrnap`
- `03a_contig_id_gene_context.ipynb`: Right after imports (cell 0) — checks `blastn`, `barrnap`, `prodigal`, `hmmsearch`, `minimap2`, `samtools`
- `03b_16S_phylogeny.ipynb`: Right after imports (new cell) — checks `barrnap`, `mafft`, `trimal`, + at least one of `iqtree2`/`iqtree3`/`iqtree`/`FastTree`

### 3. Created Shared Utilities Module
**File:** `notebooks/utils.py`

Contains DRY implementations of common functions (171 lines):
```python
find_project_root(start=None)
load_config(project_root=None)
load_fasta(path)
is_placeholder(name)
parse_title(title)
```

**Import pattern:**
```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / 'notebooks'))
from utils import find_project_root, load_config, load_fasta, is_placeholder, parse_title
```

### 4. Created Configuration File
**File:** `config/config.yaml`

Centralizes user-specific settings:
```yaml
silva_db: ~/db/silva/silva_SSU_NR99  # or null for default (project_root/db/silva/silva_SSU_NR99)
entrez_email: palmer.lab@wisc.edu
tool_paths:
  barrnap:  ~/miniforge3/envs/genomics_arm64/bin/barrnap  # fallback if not on PATH
  mafft:    ~/miniforge3/envs/genomics_arm64/bin/mafft
  ... etc
```

**Usage in notebooks:**
```python
from utils import load_config
cfg = load_config()
silva_db = cfg['silva_db']  # picks up config.yaml value
```

---

## Migration Roadmap

### Phase 1: Update `02_all_samples_community.ipynb` (Recommended First)
**Why first:** This notebook is most stable (22/25 cells executed), has most duplicated functions, and is the reference for later notebooks.

**Steps (each cell number from original):**

1. **Cell 4** (find_project_root setup):
   ```python
   # Remove inline find_project_root() definition
   # Add to imports:
   import sys
   sys.path.insert(0, str(Path(__file__).parent / 'notebooks'))
   from utils import find_project_root, load_config, load_fasta, is_placeholder, parse_title, load_blast_results, top_hit_table

   PROJECT = find_project_root()
   cfg = load_config()
   # Keep AL_ SAMPLES, SAMPLE_LABELS, SAMPLE_COLOR_MAP as-is
   ```

2. **Cell 6** (load_fasta):
   ```python
   # Remove function definition — already imported from utils
   # Keep: seqs_by_sample = {} loop and file loading
   ```

3. **Cell 8** (blast_16S_local):
   ```python
   # Update SILVA_DB line:
   SILVA_DB = cfg['silva_db']  # from config.yaml, NOT hardcoded path

   # Remove function definition — already imported from utils
   # Keep: for s in ALL_SAMPLES: blast_16S_local(s) execution
   ```

4. **Cell 10** (parse_title & friends):
   ```python
   import re  # still needed

   # Remove: PLACEHOLDERS, is_placeholder(), parse_title() definitions
   # Remove: load_blast_results(), top_hit_table() definitions
   # All are imported from utils

   # Keep: all_hits = pd.concat(...) and subsequent code
   ```

5. **Test:** Run notebook end-to-end, verify identical output to baseline

------

### Phase 2: Update `03a_contig_id_gene_context.ipynb`
**Setup cell (cell 2):**
```python
# Add imports
import sys
sys.path.insert(0, str(Path(__file__).parent / 'notebooks'))
from utils import find_project_root, load_config, is_placeholder

# Use find_project_root() instead of hardcoded path:
PROJECT = find_project_root()
cfg = load_config()

# Rest stays the same (TAX_DIR, POLISH_DIR, etc. derived from PROJECT)
```

**Cell 4** (is_placeholder, is_novel_cyano):
```python
import re

# Remove is_placeholder() — import from utils
# Keep is_novel_cyano() — it's specific to this notebook (uses is_placeholder but has custom logic)

# Rest: novel_contigs extraction stays the same
```

**Cell 6** (gc):
```python
# Can keep gc() inline (it's simple, 3 lines) OR import from utils
# Either way works; keeping it inline is fine
```

**Cell 17** (parse_gff, classify_hit):
```python
# Keep these inline — they're specific to gene arrow visualization
# parse_gff() and classify_hit() are not used elsewhere
```

---

### Phase 3: Complete `03b_16S_phylogeny.ipynb` (Execute & Validate)
**Cell 1** (setup):
```python
# Add:
import sys
sys.path.insert(0, str(Path(__file__).parent / 'notebooks'))
from utils import find_project_root, load_config

# Update hardcoded path:
PROJECT = find_project_root()  # instead of Path('/Users/...')
cfg = load_config()
```

**Cell 5** (MAFFT/trimAL paths):
```python
# Use tool fallback logic from config:
MAFFT = shutil.which('mafft') or cfg.get('tool_paths', {}).get('mafft')
TRIMAL = shutil.which('trimal') or cfg.get('tool_paths', {}).get('trimal')
```

**Then:** Execute all 7 cells, debug NCBI Entrez/alignment/tree issues

---

### Phase 4: Update `00_driver_2.0.ipynb` (Optional)
Already has `find_project_root()` — could import from utils, not critical.

---

## Benefits After Migration

| Before | After |
|---|---|
| `is_placeholder()` defined 3 times (02, 03, utils) | Defined once; imported everywhere |
| `parse_title()` has 60-line docstring in cell 10 of 02, lost in 03 | Docstring lives in utils.py; available to all |
| Hardcoded `/Users/jaemac/...` paths in 03a_contig_id_gene_context setup | Uses `find_project_root()` — portable across systems |
| SILVA DB path `PROJECT / 'db' / 'silva' / 'silva_SSU_NR99'` embedded in 3 places | Centralized in `config.yaml`; one edit propagates everywhere |
| New user has to edit 30+ hardcoded paths across 5 notebooks | New user edits only `config/config.yaml` |

---

## Current Status

✅ **Done:**
- Docstrings added to all utility functions
- Tool checks added to all three main notebooks
- `utils.py` created and ready to import (now also includes `blast_16S_local`, `load_blast_results`, `top_hit_table`)
- `config.yaml` created with sensible defaults
- `02_all_samples_community.ipynb` migrated — all inline utility defs removed, imports from utils, SILVA_DB reads from `cfg`
- `03a_contig_id_gene_context.ipynb` migrated — hardcoded PROJECT path replaced with `find_project_root()`, inline `is_placeholder` removed
- `03b_16S_phylogeny.ipynb` migrated — hardcoded PROJECT path replaced with `find_project_root()`, tool paths (barrnap, mafft, trimal) read from `cfg`
- README updated: `notebooks/00_driver.ipynb` → `notebooks/00_driver_2.0.ipynb`

⏸️ **Pending (optional):**
- Update `00_driver_2.0.ipynb` to import `find_project_root` from utils (low priority — already functional)

---

## Recommended Next Steps

1. **In your notebook editor:**
   - Open `notebooks/02_all_samples_community.ipynb`
   - Apply Phase 1 edits cell-by-cell
   - Run notebook; verify output matches current baseline

2. **Repeat for Phase 2 & 3**

3. **Then address the 7942 variant-calling project** with a cleaner codebase foundation
