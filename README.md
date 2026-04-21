# bio_tools

A small, practical Python toolbox for common bioinformatics tasks used in this repo: FASTA wrangling, phylogenetic tree utilities (ETE4/Biopython), taxonomy helpers, and a few workflow-oriented helpers (e.g. collapsing gene trees to higher-level representatives, detecting redundant paralog clades).

This is primarily a **Python package** (importable modules). 

## Requirements

- Python `>= 3.13` (as declared in `pyproject.toml`)
- Key third-party deps used throughout the package: `ete4`, `biopython`, `pandas`, `pyyaml`, `requests`, `matplotlib`/`seaborn`, `pycirclize`

## Install

This repo is set up to work well with [`uv`](https://github.com/astral-sh/uv).

Create a project environment and install dependencies:

```bash
uv sync
```

Run tests inside the `uv` environment:

```bash
uv run pytest -q
```

If you specifically want an editable install (PEP 660):

```bash
uv venv
source .venv/bin/activate
uv pip install -e .
```

## What’s inside

The code lives under `src/bio_tools/`.

- `bio_tools.files` — FASTA helpers (filtering, merging, DataFrame ↔ FASTA conversions)
- `bio_tools.taxa` — NCBI taxonomy helpers (mapping scientific names ↔ taxids, fetching higher-rank taxids)
- `bio_tools.phylo` — tree reduction / collapsing utilities
- `bio_tools.homology` — homology-oriented helpers (e.g. redundant paralog clade detection)
- `bio_tools.viz` — tree/heatmap visualization helpers (ETEs SmartView, pycirclize)

## Usage

### Reduce a gene tree to family-level representatives

`bio_tools.phylo.reduce.reduce_tree_to_family_representatives` finds clades annotated with taxonomic ranks (family/subfamily/tribe/genus), chooses a representative sequence per clade (currently: **longest sequence**), writes:

- a mapping JSON: `representative -> all clade members`
- a filtered FASTA where redundant members are removed (while preserving “characterized bait” sequences)

Example (mirrors the test setup):

```python
from ete4 import PhyloTree
from bio_tools.phylo.reduce import reduce_tree_to_family_representatives
from pathlib import Path

tree_path = "tests/data/test_reduce_to_fam_level.nwk"
fasta_path = "tests/data/test_reduce_to_fam_level.fasta"

t = PhyloTree(tree_path, sp_naming_function=lambda name: name.split("__")[-1])
t.annotate_ncbi_taxa(taxid_attr="species")

out_dir = Path("./out")
out_dir.mkdir(parents=True, exist_ok=True)

mapping, out_fasta, out_json = reduce_tree_to_family_representatives(
	tree=t,
	fasta_path=fasta_path,
	output_dir=str(out_dir),
)

print(len(mapping), out_fasta, out_json)
```

Notes:

- `annotate_ncbi_taxa(...)` populates per-node taxonomy props such as `rank`, which the reducer relies on.
- The first time you use `ete4.NCBITaxa`, it may download/build a local taxonomy DB under `~/.etetoolkit/taxa.sqlite`.

### Detect redundant paralog clades (same-species clades)

`bio_tools.homology.paralog.detect_redundant_paralog_clades` scans a gene tree to find clades that contain leaves from a *single species* (using a caller-provided `species_extractor`).

```python
from ete4 import Tree
from bio_tools.homology.paralog import detect_redundant_paralog_clades

t = Tree("(gene1__9606:0.1,gene2__9606:0.1,gene3__10090:0.2);")
species_extractor = lambda name: name.split("__")[-1]

clades = detect_redundant_paralog_clades(
	tree=t,
	species_extractor=species_extractor,
	return_as_strings=True,
)

print(clades)
```

To turn clades into an explicit “keep one representative, drop the rest” mapping and filter a FASTA accordingly:

```python
from pathlib import Path
from bio_tools.homology.paralog import (
	map_representative_paralog_to_all_redundant_paralogs,
	reduce_seq_collection_to_non_redundancy,
)

fasta = Path("input.fasta")
grouped = [["seqA__9606", "seqB__9606"], ["seqC__10090", "seqD__10090"]]

out_dir = Path("./out")
out_dir.mkdir(parents=True, exist_ok=True)

mapping = map_representative_paralog_to_all_redundant_paralogs(
	input_fasta=fasta,
	grouped_redundant_paralogs_as_str=grouped,
	output_dir=out_dir,
)

filtered_fasta = reduce_seq_collection_to_non_redundancy(
	input_fasta=fasta,
	redundancy_mapping=mapping,
	output_dir=out_dir,
)

print(filtered_fasta)
```

### FASTA utilities

Common helpers live in `bio_tools.files.fasta`.

```python
from bio_tools.files.fasta import filter_fasta, merge_fastas

filter_fasta(
	input_fasta="in.fasta",
	output_fasta="out.fasta",
	accessions={"XP_12345", "my_seq_id"},
	mode="remove",  # or "keep"
)

merge_fastas(["a.fasta", "b.fasta"], "merged.fasta")
```

### Taxonomy helpers

The taxonomy module uses `ete4.NCBITaxa`.

```python
from bio_tools.taxa.taxonomy import map_scientific_notation_to_tax_id, get_taxonomic_ranks

taxids = map_scientific_notation_to_tax_id([
	"Arabidopsis thaliana",
	"Oryza sativa",
])

# Map species taxids to their family taxids
species_taxids = list(taxids.values())
family_taxids = get_taxonomic_ranks(species_taxids, rank="family")
print(family_taxids)
```

## Development

### Run tests

```bash
uv run pytest -q
```

### Build / package

This project is a standard setuptools build (see `pyproject.toml`).

```bash
uv run python -m pip install build
uv run python -m build
```



