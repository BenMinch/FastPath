# FastPath-KEGG

Low-memory, high-speed KEGG pathway profiler for metagenomic reads.

Maps short reads directly to KEGG Orthology (KO) entries via k-mer matching, then scores KEGG pathways by KO coverage and harmonic-mean abundance. The hot path runs entirely in a compiled C extension; typical throughput is **~1–2 M reads/second** on a single core.

---

## How it works

```
Environmental FASTQ reads
        │
        ▼
  [C extension: sliding 2-bit k-mer window]
        │  canonical k-mer lookup in open-address hash table
        ▼
  KO integer ID  →  RPK accumulation
        │
        ▼
  Normalise to CPM  →  PathwayScorer
        │  coverage threshold + harmonic mean
        ▼
  {sample}_ko_abundances.tsv
  {sample}_pathways.tsv
  {sample}_stats.json
```

The k-mer index maps each canonical k-mer (default k=31) to a KO integer ID. k-mers shared by ≥5 KOs are marked as uninformative and ignored during profiling. Index RAM scales linearly: ~12 bytes per slot × 2× `max_kmers`.

---

## Installation

### Requirements

- Python ≥ 3.9
- GCC or Clang
- zlib development headers (`apt install zlib1g-dev` / `brew install zlib`)

### Build

```bash
git clone https://github.com/yourorg/fastpath-kegg
cd fastpath-kegg
python setup.py build_ext --inplace
pip install -e .
```

---

## Workflow

### 1. Download KEGG data

Fetches pathway→KO mappings and gene nucleotide sequences via the KEGG REST API. Results are cached locally; re-runs are incremental.

```bash
fastpath download \
    --out-dir  /data/kegg_db \
    --pathways map00,map01 \        # map00=metabolism, map01=overview (default)
    --max-genes-per-ko 50           # genes fetched per KO (default 50)
```

**Options:**

| Flag | Default | Description |
|---|---|---|
| `--out-dir` | *(required)* | Where to store downloaded data |
| `--pathways` | `map00,map01` | Comma-separated KEGG map prefixes |
| `--max-genes-per-ko` | `50` | Sequences per KO; lower → smaller index |
| `--organism` | *(all prokaryotes)* | Restrict to one KEGG organism code, e.g. `eco` |
| `--force` | off | Re-download even if cached |

KEGG allows ≤3 requests/second for academic use; the downloader respects this automatically.

---

### 2. Build k-mer index

```bash
fastpath build \
    --kegg-dir  /data/kegg_db \
    --out-dir   /data/kegg_index \
    --max-kmers 40000000
```

**RAM guide for `--max-kmers`:**

| Value | Peak RAM | Use case |
|---|---|---|
| `10_000_000` | ~200 MB | Testing / small datasets |
| `40_000_000` | ~800 MB | Recommended for most runs |
| `100_000_000` | ~2 GB | High-depth / broad coverage |

Outputs written to `--out-dir`:
- `kmer_index.bin.gz` — compressed k-mer hash table
- `pathways.json` — pathway→KO integer ID map
- `meta.json` — k, gene lengths, ko2id mapping

---

### 3. Profile reads

```bash
fastpath run \
    --input   sample.fastq.gz \
    --index   /data/kegg_index \
    --out-dir results/
```

Multiple input files are merged before profiling:

```bash
fastpath run \
    --input lane1.fastq.gz lane2.fastq.gz \
    --index /data/kegg_index \
    --out-dir results/
```

**Key flags:**

| Flag | Default | Description |
|---|---|---|
| `--min-hits` | `2` | Min k-mer hits to a KO for a read to count |
| `--min-identity` | `0.80` | Min fraction of sampled k-mers matching best KO |
| `--pathway-coverage` | `0.80` | Min KO coverage for a pathway to be reported |

---

### 4. Inspect results

```bash
fastpath inspect results/sample_pathways.tsv --head 30
```

---

## Output format

### `{sample}_ko_abundances.tsv`
```
# KO	Abundance (CPM)
K00001	4521.883200
K00024	3187.441100
...
```

### `{sample}_pathways.tsv`
```
# Pathway	Abundance
ko00010	2843.117000
ko00020	1204.553000
...
```
Pathway abundance = harmonic mean CPM of covered KOs. Only pathways meeting `--pathway-coverage` are reported.

### `{sample}_stats.json`
```json
{
  "total_reads": 5000000,
  "mapped_reads": 3214887,
  "mapping_rate": 0.6430
}
```

---

## Python API

```python
from fastpath import KmerIndex, AbundanceProfiler, PathwayScorer
import json, array

# Load index
idx = KmerIndex.load("kegg_index/kmer_index.bin.gz")
with open("kegg_index/meta.json") as f:
    meta = json.load(f)

gene_lengths = array.array('I', meta["gene_lengths"])
ko2id        = meta["ko2id"]
id2ko        = {v: k for k, v in ko2id.items()}

# Profile
profiler = AbundanceProfiler(index=idx, gene_lengths=gene_lengths)
profiler.profile("sample.fastq.gz")

ko_abund = profiler.normalized_abundances()   # {ko_int_id: CPM}
print(profiler.stats)

# Score pathways
with open("kegg_index/pathways.json") as f:
    pathways = json.load(f)
scorer   = PathwayScorer(pathways, coverage_threshold=0.80)
pw_abund = scorer.score({int(k): v for k, v in ko_abund.items()})

for pw, ab in sorted(pw_abund.items(), key=lambda x: -x[1])[:10]:
    print(f"{pw}\t{ab:.2f}")
```

---

## Running tests

```bash
python -m pytest tests/ -v
```

---

## Differences from the ChocoPhlAn version

| | ChocoPhlAn version | KEGG version |
|---|---|---|
| Reference DB | ChocoPhlAn pangenomes | KEGG nucleotide CDS |
| Gene family unit | UniRef90 cluster | KEGG Orthology (KO) |
| Pathway format | Custom JSON (reaction lists) | KEGG pathway→KO membership |
| Build input | `.ffn.gz` files from ChocoPhlAn | Downloaded via KEGG REST API |
| Output label | UniRef90 ID | KO ID (e.g. `K00001`) |
| C extension | Unchanged | Unchanged |

The C extension (`_kmercore.c`) is identical in both versions — it has no database-specific logic.
