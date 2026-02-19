---
name: amina-init
description: Initialize Claude Code for Amina CLI. Use this skill at the start of every conversation when working with AminoAnalytica's Amina CLI in Claude Code. Provides guidance on setup and getting started with key commands.
---

# Amina CLI

Foundation skill for running protein engineering tools via the Amina command-line interface.

For latest documentation, see:
- **CLI documentation**: [Amina CLI Docs](https://app.aminoanalytica.com/docs/cli)
- **Tool catalog**: [Amina Tools](https://app.aminoanalytica.com/docs/tools)

## Skill Types

This skill is part of a network of protein engineering skills:

| Type | Purpose | Examples |
|------|---------|----------|
| **Amina Tool Skills** | Best practices for specific Amina CLI tools | rfdiffusion, proteinmpnn, boltz, diffdock |
| **General Skills** | External tools and databases for protein work | PDB database, UniProt, PyMOL, BioPython |
| **Workflow Skills** | End-to-end pipelines combining multiple tools | binder-design, enzyme-engineering |

This foundation skill covers Amina CLI usage; tool-specific and workflow skills provide deeper guidance.

## Environment Setup

**At the start of every session, verify the setup:**

### 1. Check Amina CLI Installation

```bash
amina --version
```

If missing, install with `pip install amina-cli` (try `pip3` or `python3 -m pip` as fallbacks).

### 2. Check Authentication

```bash
amina auth status
```

If not authenticated, ask for their API key (sign up free at https://aminoanalytica.com), then set with:
```bash
amina auth set-key "{key}"
```

### 3. Ask About Python Environment

Since the user installed `amina-cli` via pip, they likely already have a Python environment. Ask:

- "Do you have a Python virtual environment set up for this project? If so, what's the activation command?"

If no environment exists, offer to create one:

```bash
# Create virtual environment
python -m venv .venv

# Activate (macOS/Linux)
source .venv/bin/activate

# Activate (Windows)
.venv\Scripts\activate

# Install amina-cli in the environment
pip install amina-cli
```

**Why this matters:**
- Many protein engineering workflows require Python scripts for data processing
- BioPython, NumPy, and other scientific packages may be needed
- Keeping dependencies isolated prevents conflicts

**If the user already has an environment:**
- Use their existing activation command
- Install additional packages as needed into that environment

## Quick Reference

```bash
# Authentication
amina auth set-key              # Set API key
amina auth status               # Check authentication

# Tool discovery
amina tools list                # List all available tools
amina tools <tool>              # Show tool description and outputs
amina run <tool> --help         # Show tool parameters and usage

# Running tools
amina run <tool> [params]       # Run a tool
amina run <tool> -o <dir>       # Specify output directory
amina run <tool> --background   # Run in background (async)
amina run <tool> --job-name <n> # Custom job name

# Job management
amina jobs list                 # List recent jobs
amina jobs status <job-id>      # Check job status
amina jobs wait <job-id>        # Wait for job completion
amina jobs download <job-id>    # Download job results
```

## Core Workflows

### 1. Before Running Any Tool

**Always discover parameters first:**
```bash
amina run <tool> --help
```

This returns exact parameter names, types, and descriptions. Never guess parameters.

### 2. Running Tools

**Synchronous (short-running):**
```bash
amina run pdb-cleaner --input protein.pdb -o ./results/
```

**Asynchronous (long-running):**
```bash
amina run esmfold --sequence "MKTV..." --background --job-name fold_run1
# Returns job ID immediately
amina jobs wait <job-id>
amina jobs download <job-id> -o ./results/
```

**When to use `--background`:**
- Structure prediction (ESMFold, Boltz)
- RFDiffusion (all modes)
- Molecular dynamics simulations
- Any job expected to run >30 seconds

### 3. Job Management

**Check running jobs:**
```bash
amina jobs list --status running
```

**Monitor a specific job:**
```bash
amina jobs status <job-id>
# Shows: status, progress, runtime, output files when complete
```

**Wait and download:**
```bash
amina jobs wait <job-id>           # Blocks until complete
amina jobs download <job-id> -o ./ # Download all outputs
```

## File Organization

Adapt organization to task complexity:

**Single tool run:**
```
./results/<tool>_<jobname>/
```

**Multi-step workflow:**
```
./project/
├── 01_structure_prediction/
├── 02_binder_design/
├── 03_sequence_design/
└── 04_validation/
```

**Iterative design:**
```
./project/
├── round_1/
├── round_2/
└── round_3/
```

Rename cryptic output files for clarity when needed.

## Error Handling

**Authentication errors:**
```bash
amina auth status               # Verify key is set
amina auth set-key              # Re-authenticate if needed
```

**Job failures:**
```bash
amina jobs status <job-id>      # Check error message
amina jobs logs <job-id>        # View detailed logs
```

**Common issues:**
- Invalid input file format → Check tool's expected format with `amina tools <tool>`
- Job timeout → Use smaller input or increase resources if available
- Parameter errors → Verify exact parameter names via `amina tools <tool>`
