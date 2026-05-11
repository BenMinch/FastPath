# FastPath

**Low-memory, high-speed metagenomic metabolic pathway profiler**

FastPath profiles the abundance of microbial metabolic pathways from raw metagenomic reads.  
It is inspired by [HUMAnN3](https://huttenhower.sph.harvard.edu/humann/) but is designed from the ground up for **speed and low RAM** — making it practical on laptops, HPC nodes with tight memory limits, and large cohort studies.

---

## Design philosophy

| Concern | HUMAnN3 | FastPath |
|---|---|---|
| Nucleotide search | Bowtie2 (full alignment) | k-mer exact match (hash lookup) |
| Translated search | DIAMOND | *(roadmap)* |
| Memory | 10–40 GB (pangenome + index) | Tunable budget (default ≈ 400 MB) |
| Speed | Minutes–hours | Seconds–minutes |
| Streaming | No (loads reads into RAM) | Yes (O(1) per read) |
| Pathway scoring | MinPath ILP | Harmonic-mean coverage scoring |
| Output | HUMAnN TSV | HUMAnN-compatible TSV |

### Core algorithms

1. **k-mer index** (`k=31` canonical k-mers, blake2b hashing, open-address dict)  
   Each database gene is broken into canonical 31-mers; each k-mer maps to a list of gene-family IDs.

2. **Streaming read profiler**  
   Reads are processed one at a time from FASTQ/FASTA (gzip supported).  
   Per read: k-mer hit counts are tallied per gene family → filtered by `--min-hits` and `--min-identity` → fractionally distributed across tied families → accumulated as RPK.

3. **Normalisation**  
   Gene RPK values are sum-normalised to **copies per million (CPM)** — analogous to HUMAnN's community-level abundances.

4. **Pathway scoring**  
   A pathway is scored only if ≥ `--pathway-coverage` fraction of its reactions are present.  
   Pathway abundance = harmonic mean of per-reaction abundances (penalises missing steps more strongly than arithmetic mean).

---

## Installation

```bash
# From source (pure Python, zero required dependencies)
git clone https://github.com/yourorg/fastpath.git
cd fastpath
pip install -e .

# Verify
fastpath --version
```

Optional faster hashing:
```bash
pip install mmh3
```

---

## Quick start (self-contained demo)

```bash
fastpath demo --out-dir ./demo_out --num-reads 50000
```

This will:
1. Generate a 12-gene synthetic pangenome (3 metabolic families)
2. Simulate 50 000 reads (70% from the pangenome, 30% random background)
3. Build a k-mer index, profile the reads, and score 4 toy pathways
4. Write TSV outputs to `demo_out/demo_results/`

Expected output:
```
[FastPath demo] Generating demo database …
  Done: 12 genes, ~7,000 unique k-mers
[FastPath demo] Generating 50,000 synthetic reads …
[FastPath run]  k=31  index=demo_out/demo_db  input=[...]
  Loaded 6,980 k-mers
  Profiling demo_out/demo_reads.fastq …
  3 pathways detected
  Mapping rate: 64.3% (32,179 / 50,000 reads)
✓ Done in 8.2s
```

---

## Full workflow

### Step 1 — Install humann (for its databases)

FastPath uses the same ChocoPhlAn pangenome and MetaCyc pathway databases as HUMAnN3.
Install humann to get access to its downloader and bundled pathway files:

```bash
pip install humann
```

### Step 2 — Download ChocoPhlAn and extract `pangenome.fna` + `gene_families.tsv`

```bash
python scripts/prepare_pangenome.py \
    --db-dir /data/chocophlan \
    --out-dir /data/fastpath_db
```

This downloads the full ChocoPhlAn database (~15 GB, 30–60 min on a fast connection)
and extracts two files into `--out-dir`:

- `pangenome.fna` — all gene nucleotide sequences concatenated
- `gene_families.tsv` — `gene_name<TAB>UniRef90_family` for every gene

If you have already run `humann_databases --download chocophlan full <dir>`, skip the
download with `--no-download`:

```bash
python scripts/prepare_pangenome.py \
    --db-dir /data/chocophlan \
    --out-dir /data/fastpath_db \
    --no-download
```

To work with a subset of species, filter by filename substring:

```bash
python scripts/prepare_pangenome.py \
    --db-dir /data/chocophlan \
    --out-dir /data/fastpath_db \
    --no-download \
    --species "Bacteroides_dorei" "Prevotella_copri"
```

### Step 3 — Build `pathways.json` from MetaCyc

```bash
python scripts/build_pathways_json.py \
    --out-dir /data/fastpath_db
```

This reads two files bundled with humann (auto-detected after `pip install humann`):

- `metacyc_pathways_structured_filtered_v24` — 887 MetaCyc pathway definitions
- `metacyc_reactions_level4ec_only.uniref.bz2` — reaction → UniRef90 family mapping

And writes to `--out-dir`:

- `pathways.json` — 880 pathways with integer-keyed gene-family reaction slots
- `family_index.tsv` — `integer_id<TAB>UniRef90_family` lookup table

If you have the files in non-standard locations, supply them explicitly:

```bash
python scripts/build_pathways_json.py \
    --pathways  /path/to/metacyc_pathways_structured_filtered \
    --reactions /path/to/metacyc_reactions_level4ec_only.uniref.bz2 \
    --out-dir   /data/fastpath_db
```

### Step 4 — Build the FastPath k-mer index

```bash
fastpath build \
    --fasta     /data/fastpath_db/pangenome.fna \
    --gene-map  /data/fastpath_db/gene_families.tsv \
    --pathways  /data/fastpath_db/pathways.json \
    --out-dir   /data/fastpath_index \
    --max-kmers 200000000    # ~1.6 GB RAM for full ChocoPhlAn; use 50000000 (≈400 MB) for a lighter build
```

### Step 5 — Profile reads

```bash
fastpath run \
    --input  sample_R1.fastq.gz sample_R2.fastq.gz \
    --index  /data/fastpath_index \
    --out-dir results/ \
    --min-hits        2     \
    --min-identity    0.80  \
    --pathway-coverage 0.80
```

### Step 6 — Inspect results

```bash
fastpath inspect results/sample_pathways.tsv
fastpath inspect results/sample_genefamilies.tsv --head 30
```

---

## Output files

### `<sample>_genefamilies.tsv`
```
# Gene Family	Abundance (CPM)
UniRef90_A0A000	58423.12
UniRef90_B0B001	31204.87
...
```

### `<sample>_pathways.tsv`
```
# Pathway	Abundance
PWY-5484	142310.44
PWY-6737	98204.11
...
```

### `<sample>_stats.json`
```json
{
  "total_reads": 1000000,
  "mapped_reads": 621034,
  "mapping_rate": 0.621
}
```

---

## `pathways.json` format

FastPath uses a simple JSON structure where each pathway lists its reactions and the gene family IDs that can catalyse each reaction:

```json
{
  "PWY-5484": [
    [42, 107],
    [88],
    [12, 13, 14]
  ]
}
```

- Outer array → reactions (ordered steps)
- Inner array → gene-family IDs that can catalyse that reaction (OR logic: any one suffices)
- A reaction is "present" if any of its gene families has CPM > 0

---

## Python API

```python
from fastpath import KmerIndex, AbundanceProfiler, PathwayScorer

# Build index
idx = KmerIndex(k=31, max_kmers=50_000_000)
idx.add_sequence(gene_sequence, gene_id=42)
idx.save("index.bin.gz")

# Load and profile
idx = KmerIndex.load("index.bin.gz")
profiler = AbundanceProfiler(idx, gene_lengths={42: 900})
profiler.profile("sample.fastq.gz")

gene_cpm = profiler.normalized_abundances()

# Score pathways
scorer = PathwayScorer(pathways={"PWY-A": [[42], [7, 8]]})
pathway_abund = scorer.score(gene_cpm)
```

---

## Memory tuning

The dominant memory cost is the k-mer index. Use `--max-kmers` to set a hard budget:

| `--max-kmers` | Approx. RAM | Suitable for |
|---|---|---|
| 10 000 000 | ~80 MB | Small custom databases |
| 50 000 000 | ~400 MB | Standard use (default) |
| 200 000 000 | ~1.6 GB | Full ChocoPhlAn-scale |

When the budget is reached, additional k-mers are silently skipped. Sensitivity decreases gracefully; you can inspect how many k-mers were indexed with the build log.

---

## Development

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

---

## Roadmap

- [ ] Translated search (protein k-mers / minimisers) for unclassified reads
- [ ] Multi-threaded read processing (`--threads N`)
- [ ] Streaming merge of paired-end reads
- [ ] Integration with MetaPhlAn4 taxonomic profiles for per-species pathway stratification
- [ ] Rust extension for k-mer hashing hot path

---

## Citation

If you use FastPath in your research, please also cite HUMAnN3, whose design this tool is inspired by:

> Beghini et al. (2021). *Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3.* eLife 10:e65088. https://doi.org/10.7554/eLife.65088
