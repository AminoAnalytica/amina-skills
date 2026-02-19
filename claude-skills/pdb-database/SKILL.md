---
name: pdb-database
description: Query and retrieve protein/nucleic acid structures from RCSB PDB. Supports text, sequence, and structure-based searches, coordinate downloads, and metadata retrieval for structural biology workflows.
---

# RCSB Protein Data Bank

The Protein Data Bank hosts over 200,000 experimentally determined macromolecular structures plus computed models from AlphaFold and ModelArchive. This skill provides programmatic access to search, download, and analyze structural data.

## Applicable Scenarios

| Task | Examples |
|------|----------|
| Structure Retrieval | Download coordinates for a known PDB ID |
| Similarity Search | Find structures similar by sequence or 3D fold |
| Metadata Access | Get resolution, method, organism, ligands |
| Dataset Building | Compile structures for ML training or analysis |
| Drug Discovery | Identify ligand-bound structures for a target |
| Quality Filtering | Select high-resolution, well-refined structures |

## Setup

```bash
uv pip install rcsb-api requests
```

The `rcsb-api` package provides:
- `rcsbapi.search` - Query construction and execution
- `rcsbapi.data` - Metadata retrieval via REST/GraphQL

## Quick Reference

### Search Queries

```python
from rcsbapi.search import TextQuery, AttributeQuery, SequenceQuery, StructSimilarityQuery
from rcsbapi.search.attrs import rcsb_entity_source_organism, rcsb_entry_info, exptl

# Text search
results = list(TextQuery("kinase inhibitor")())

# Filter by organism
human = AttributeQuery(
    attribute=rcsb_entity_source_organism.scientific_name,
    operator="exact_match",
    value="Homo sapiens"
)

# Filter by resolution
high_res = AttributeQuery(
    attribute=rcsb_entry_info.resolution_combined,
    operator="less",
    value=2.0
)

# Combine queries: & (AND), | (OR), ~ (NOT)
results = list(TextQuery("kinase") & human & high_res)()

# Sequence similarity (MMseqs2)
seq_query = SequenceQuery(
    value="MKTAYIAKQRQISFVK...",
    evalue_cutoff=1e-5,
    identity_cutoff=0.7
)

# Structure similarity (3D fold)
struct_query = StructSimilarityQuery(
    structure_search_type="entry",
    entry_id="4HHB"
)
```

### Data Retrieval

```python
from rcsbapi.data import fetch, Schema

# Entry metadata
data = fetch("4HHB", schema=Schema.ENTRY)
print(data["struct"]["title"])
print(data["rcsb_entry_info"]["resolution_combined"])

# Polymer entity (chain info + sequence)
entity = fetch("4HHB_1", schema=Schema.POLYMER_ENTITY)
sequence = entity["entity_poly"]["pdbx_seq_one_letter_code"]
```

### File Downloads

| Format | URL Pattern |
|--------|-------------|
| PDB | `https://files.rcsb.org/download/{ID}.pdb` |
| mmCIF | `https://files.rcsb.org/download/{ID}.cif` |
| Assembly | `https://files.rcsb.org/download/{ID}.pdb1` |
| FASTA | `https://www.rcsb.org/fasta/entry/{ID}` |

```python
import requests
from pathlib import Path

def download_structure(pdb_id: str, fmt: str = "cif", outdir: str = ".") -> Path:
    url = f"https://files.rcsb.org/download/{pdb_id}.{fmt}"
    resp = requests.get(url)
    resp.raise_for_status()
    outpath = Path(outdir) / f"{pdb_id}.{fmt}"
    outpath.write_text(resp.text)
    return outpath
```

## Common Workflows

### Find High-Quality Human Structures

```python
from rcsbapi.search import TextQuery, AttributeQuery
from rcsbapi.search.attrs import rcsb_entity_source_organism, rcsb_entry_info, exptl

query = (
    TextQuery("receptor") &
    AttributeQuery(
        attribute=rcsb_entity_source_organism.scientific_name,
        operator="exact_match",
        value="Homo sapiens"
    ) &
    AttributeQuery(
        attribute=rcsb_entry_info.resolution_combined,
        operator="less",
        value=2.5
    ) &
    AttributeQuery(
        attribute=exptl.method,
        operator="exact_match",
        value="X-RAY DIFFRACTION"
    )
)
results = list(query())
```

### Batch Metadata Retrieval

```python
import time
from rcsbapi.data import fetch, Schema

def fetch_batch(pdb_ids: list, delay: float = 0.3) -> dict:
    """Fetch metadata with rate limiting."""
    results = {}
    for pdb_id in pdb_ids:
        try:
            data = fetch(pdb_id, schema=Schema.ENTRY)
            results[pdb_id] = {
                "title": data["struct"]["title"],
                "resolution": data.get("rcsb_entry_info", {}).get("resolution_combined"),
                "method": data.get("exptl", [{}])[0].get("method"),
            }
        except Exception as e:
            results[pdb_id] = {"error": str(e)}
        time.sleep(delay)
    return results
```

### Find Drug-Bound Structures

```python
from rcsbapi.search import AttributeQuery
from rcsbapi.search.attrs import rcsb_nonpolymer_entity_instance_container_identifiers

# Find structures containing imatinib (ligand ID: STI)
query = AttributeQuery(
    attribute=rcsb_nonpolymer_entity_instance_container_identifiers.comp_id,
    operator="exact_match",
    value="STI"
)
drug_complexes = list(query())
```

### GraphQL for Complex Queries

```python
import requests

query = """
{
  entry(entry_id: "4HHB") {
    struct { title }
    rcsb_entry_info {
      resolution_combined
      deposited_atom_count
    }
    polymer_entities {
      rcsb_polymer_entity { pdbx_description }
      entity_poly { pdbx_seq_one_letter_code }
    }
  }
}
"""

response = requests.post(
    "https://data.rcsb.org/graphql",
    json={"query": query}
)
result = response.json()["data"]["entry"]
```

## Key Concepts

| Term | Definition |
|------|------------|
| **PDB ID** | 4-character alphanumeric code (e.g., "4HHB"). AlphaFold uses "AF_" prefix |
| **Entity** | Distinct molecular species. A homodimer has one entity appearing twice |
| **Resolution** | Quality metric in angstroms. Lower is better; <2.0 A is high quality |
| **Biological Assembly** | Functional oligomeric state (may differ from asymmetric unit) |
| **mmCIF** | Modern format replacing legacy PDB; required for large structures |

## Best Practices

| Practice | Rationale |
|----------|-----------|
| Use mmCIF format | PDB format has atom count limits |
| Filter by resolution | <2.5 A for most analyses; <2.0 A for detailed work |
| Check experimental method | X-ray vs cryo-EM vs NMR have different quality metrics |
| Rate limit requests | 2-3 req/s to avoid 429 errors |
| Cache downloads | Structures rarely change after release |
| Prefer GraphQL | Reduces requests for complex data needs |

## Troubleshooting

| Issue | Resolution |
|-------|------------|
| 404 on entry fetch | Entry may be obsoleted; check RCSB website for superseding ID |
| 429 Too Many Requests | Implement exponential backoff; reduce request rate |
| Empty search results | Check query syntax; use `query.to_dict()` to debug |
| Large structure fails | Use mmCIF format instead of PDB |
| Missing sequence data | Query polymer entity (`{ID}_1`) not entry |

## References

See `references/api-reference.md` for:
- Complete REST endpoint documentation
- All searchable attributes and operators
- Advanced query patterns
- Rate limiting strategies

## External Links

- RCSB PDB: https://www.rcsb.org
- API Overview: https://www.rcsb.org/docs/programmatic-access/web-apis-overview
- Python Package: https://rcsbapi.readthedocs.io/
- GraphQL Explorer: https://data.rcsb.org/graphql
