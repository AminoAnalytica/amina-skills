---
name: amina-init
description: Foundation skill for autonomous protein engineering via the Amina CLI. Use at session start when working with AminoAnalytica's Amina CLI or when you need help using the Amina CLI. Triggers on: protein design, structure prediction, binder design, docking, molecular dynamics, enzyme engineering, or any mention of Amina/amina-cli. Related tool-specific skills (e.g., amina-rfdiffusion, amina-proteinmpnn, amina-boltz-2, amina-diffdock) and workflow skills (e.g., workflow-binder-design, workflow-enzyme-engineering) provide deeper guidance.
---

# Amina CLI

Run autonomous protein engineering workflows.
Docs for you to read to get started: [CLI](https://app.aminoanalytica.com/docs/cli) | [Tools](https://app.aminoanalytica.com/docs/tools)

## Session Setup

**CRITICAL: Work in the user's current directory.** The working directory where Claude Code was invoked is the project root. Never create files in `~/` or other arbitrary locations.

### 1. Understand Project Context First

Before creating any files or directories, examine what already exists:

```bash
pwd                             # Confirm working directory
ls -la                          # See existing structure
```

Look for: existing project organization, input files (PDBs, FASTAs), config files, previous Amina outputs. Adapt to the existing structure rather than imposing a new one.

### 2. Dependency Management

**All packages must be installed in the same environment as amina-cli (Python 3.11+).**

```bash
# Verify environment
python --version                          # Must be 3.11+
which python                              # Confirm correct interpreter
amina --version                           # Confirm amina accessible
```

**Install dependencies using the same method as amina-cli:**

```bash
# If amina was installed with pip:
pip install biopython matplotlib

# If amina was installed with uv:
uv pip install biopython matplotlib
```

If Python is <3.11 or amina is not found, the user needs to activate the correct environment or install amina-cli first:

```bash
# Example setup (if needed)
uv venv --python 3.11                     # or: python3.11 -m venv .venv
source .venv/bin/activate
pip install amina-cli                     # or: uv pip install amina-cli
```

### 3. Verify Amina Authentication (silently)

Only prompt user if something fails:

```bash
amina auth status               # If unauthenticated: amina auth set-key "{key}"
```

### 4. Output Session Summary

Detect the active environment:

```bash
echo $VIRTUAL_ENV                         # Set by venv/uv venv when activated
which python                              # Shows interpreter path
```

- **venv/uv venv**: `$VIRTUAL_ENV` is set (e.g., `/path/to/project/.venv`)
- **uv run** (ephemeral): `$VIRTUAL_ENV` empty but python points to `.venv` or uv cache
- **System Python**: `$VIRTUAL_ENV` empty, python points to `/usr/bin/python` or similar

After completing setup checks, output a brief summary:

```
**Amina Session Ready**
- **Directory**: {cwd}
- **Environment**: {$VIRTUAL_ENV path, "uv-managed", or "system Python"}
- **Python**: {version} | **Amina**: {version} ✓
- **Auth**: {authenticated/needs setup}
- **Project files**: {e.g., "3 PDB files, 1 FASTA" or "empty directory"}
```

Keep it concise. Only mention issues if action is needed.

## Agent Behavior

**Autonomy first**: Use defaults and infer parameters from context. Only ask for truly blocking inputs (target PDB for binder design, ligand for docking).

**Adapt to expertise**: Jargon (scRMSD, pLDDT, contigs) → be terse. Plain English → explain steps.

**Clarity in communication**: Output information in a way that LLMs can easily parse and understand. For user facing information, use markdown formatting and emojis sparingly.

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

**Always use relative paths from the current working directory.** Never use absolute paths like `~/` or `/Users/...`.

```bash
# All outputs go in current directory or subdirectories
./                              # Current working directory is project root
./downloads/                    # Downloaded files (e.g., PDBs, FASTAs)
./results/<tool>_<jobname>/     # Single tool runs
./01_structure/ 02_design/      # Multi-stage workflows
./round_1/ round_2/             # Iterative refinement
```

If the directory already has structure, follow it. If empty, create minimal organization as needed.

### Parallelization

Use `--background` for ALL parallel jobs, both long and short:

1. **Submit**: Launch all jobs with `--background` in a single message (multiple Bash tool calls — each returns a job ID instantly)
2. **Wait**: `amina jobs wait {id1} {id2} {id3} ... --poll-interval 10` (blocks until all complete)
3. **Download**: `amina jobs download {id} -o {dir}/` for each job

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
