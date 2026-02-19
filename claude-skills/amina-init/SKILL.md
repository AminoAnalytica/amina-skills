---
name: amina-init
description: Foundation skill for autonomous protein engineering via the Amina CLI. Use at session start when working with AminoAnalytica's Amina CLI. Triggers on: protein design, structure prediction, binder design, docking, molecular dynamics, enzyme engineering, or any mention of Amina/amina-cli. Related tool-specific skills (e.g., amina-rfdiffusion, amina-proteinmpnn, amina-boltz-2, amina-diffdock) and workflow skills (e.g., workflow-binder-design, workflow-enzyme-engineering) provide deeper guidance.
---

# Amina CLI

Run autonomous protein engineering workflows. Docs: [CLI](https://app.aminoanalytica.com/docs/cli) | [Tools](https://app.aminoanalytica.com/docs/tools)

## Session Setup

Verify silently at session start. Only prompt user if something fails.

```bash
amina --version                 # Install: pip install amina-cli
amina auth status               # If unauthenticated: amina auth set-key "{key}"
```

If user lacks a Python environment, offer to create one with `python -m venv .venv`.

## Agent Behavior

**Autonomy first**: Use defaults and infer parameters from context. Only ask for truly blocking inputs (target PDB for binder design, ligand for docking).

**Adapt to expertise**: Jargon (scRMSD, pLDDT, contigs) → be terse. Plain English → explain steps.

**Parse requests for**: Target (PDB ID/sequence/file), Goal, Constraints (hotspots, chains, lengths), Scale (number of designs).

**Classify workflow**:
- "binder", "bind to" → De Novo Binder Design
- "design protein" → De Novo Monomer Design
- "motif", "scaffold", "contigs" → Motif Scaffolding
- "enzyme", "catalytic", "mutant" → Enzyme Engineering
- "dock", "ligand", "SMILES" → Small Molecule Docking
- "validate", "assess" → Structure Validation
- "interface", "PPI" → PPI Analysis
- "embedding", "cluster" → Sequence Space Exploration

**Research when useful**: For unfamiliar proteins, use `WebSearch` for PDB IDs, binding interfaces, known hotspots. Skip for self-contained requests.

## CLI Reference

```bash
# Discovery (always check before running)
amina tools list                # List all tools
amina run <tool> --help         # Show parameters (never guess)

# Execution
amina run <tool> [params] -o <dir>           # Sync
amina run <tool> --background --job-name <n> # Async (>30s jobs)

# Job management
amina jobs list --status running
amina jobs status <job-id>
amina jobs wait <job-id> && amina jobs download <job-id> -o ./
amina jobs logs <job-id>        # Debug failures
```

**Use `--background` for**: Structure prediction (ESMFold, Boltz), RFDiffusion, MD simulations, any job >30s.

## Execution Patterns

### Project Organization

```bash
# Create project directory
~/amina-projects/{workflow_name}_{YYYYMMDD_HHMMSS}/
```

Single tool: `./results/<tool>_<jobname>/`
Multi-stage: `./project/01_structure/ 02_design/ 03_validation/`
Iterative: `./project/round_1/ round_2/ round_3/`

### Parallelization

- **`--background`**: Long tools (e.g., protein design, simulations, structure prediction). Wait with `amina jobs wait`.
- **Bash `run_in_background: true`**: Many short parallel jobs. Collect with `TaskOutput(block=true)`.

Always parallelize independent tasks within a stage.

### Tool Chaining

| From | Output | To | Input |
|------|--------|-----|-------|
| RFdiffusion | `*.pdb` | ProteinMPNN | `--pdb` |
| ProteinMPNN | `*.fasta` | ESMFold | `--sequence` (parse FASTA) |
| ESMFold | `*.pdb` | Simple RMSD | `--mobile` |
| P2Rank | `*_predictions.csv` | Vina | `--center-x/y/z` (parse CSV) |
| PDB Cleaner | `*_cleaned.pdb` | Any PDB tool | `--pdb` |

After `amina run`, use `Glob` to find outputs. Never guess metrics—wait for actual results.

## Results & Iteration

Present results naturally for the task. Read output files and interpret for user. Suggest 1-2 follow-ups conversationally.

**Key principle**: Never rerun expensive stages. Check what exists before re-executing.

## Error Handling

```bash
amina auth status               # Auth issues
amina jobs logs <job-id>        # Job failures
amina tools <tool>              # Format/parameter issues
```
