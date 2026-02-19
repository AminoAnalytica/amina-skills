# Amina Skills

Claude Code skills for protein engineering using the [Amina CLI](https://aminoanalytica.com).

## Installation

### Option 1: Using Amina CLI (Recommended)

```bash
pip install amina-cli
cd /your/project
amina init claude-skills
```

### Option 2: Manual Installation

```bash
cd /your/project
mkdir -p .claude/skills
curl -L https://github.com/AminoAnalytica/amina-skills/archive/main.zip -o skills.zip
unzip skills.zip "amina-skills-main/claude-skills/*" -d .
mv amina-skills-main/claude-skills/* .claude/skills/
rm -rf amina-skills-main skills.zip
```

## Usage

After installing skills, run Claude Code from your project directory:

```bash
claude
```

Then use the `/amina-init` command to initialize Claude Code for protein engineering tasks.

## Skill Types

Skills are organized into categories:

| Type | Purpose | Examples |
|------|---------|----------|
| **Foundation** | Core CLI usage patterns | amina-init |
| **Tool Skills** | Best practices for specific tools | rfdiffusion, proteinmpnn, boltz, diffdock |
| **Workflow Skills** | End-to-end pipelines | binder-design, enzyme-engineering |

## What Skills Teach Claude

With these skills installed, Claude Code can:

- Run protein structure prediction, sequence design, and docking tools
- Manage background jobs and download results
- Organize project files for multi-step workflows
- Interpret results and suggest next steps
- Handle errors and troubleshoot common issues

## Links

- [Amina CLI Documentation](https://app.aminoanalytica.com/docs/cli)
- [Tool Catalog](https://app.aminoanalytica.com/docs/tools)
- [Get API Key](https://app.aminoanalytica.com/settings/api)

## License

Apache 2.0
