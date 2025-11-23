# GitHub Actions Workflows

This directory contains automated workflows that validate and showcase the project.

## Active Workflows

### 1. Documentation Quality (`documentation.yml`)

**Trigger**: Push to main/master, Pull requests
**Purpose**: Validates documentation quality and completeness

**Checks**:
- ✅ Markdown formatting (via markdownlint)
- ✅ Link validation (checks all URLs work)
- ✅ Required sections in README (Overview, Installation, Usage, etc.)
- ✅ CITATION.cff validation
- ✅ Presence of required files (LICENSE, CONTRIBUTING.md, etc.)

### 2. MATLAB Code Validation (`matlab-validation.yml`)

**Trigger**: Push to main/master, Pull requests
**Purpose**: Validates MATLAB project structure and code quality

**Checks**:
- ✅ Project structure (main functions, models, utilities)
- ✅ Basic syntax validation (balanced parentheses, function definitions)
- ✅ Code quality indicators (documentation, end statements)
- 📊 Generates project statistics (LOC, file counts)

**Note**: This workflow does not require MATLAB to be installed. It performs structural validation only.

### 3. Archive Status Check (`archive-status.yml`)

**Trigger**: Monthly (1st of each month), Manual dispatch
**Purpose**: Monitors the health of the research archive

**Checks**:
- ✅ Archive notice presence
- ✅ Documentation freshness
- ✅ Citation information accuracy
- 📊 Generates archive status report

### 4. Project Statistics (`project-stats.yml`)

**Trigger**: Push to main/master, Manual dispatch
**Purpose**: Generates comprehensive project statistics

**Reports**:
- 📊 Code metrics (files, LOC, comments, documentation ratio)
- 📊 Model implementations (lists all network models)
- 📊 Repository statistics (commits, dates)
- 📚 Research publication information

## Configuration Files

- `.markdownlint.json` - Markdown linting configuration
- `.markdown-link-check.json` - Link validation configuration

## Why These Workflows?

As an archived research project, these workflows serve multiple purposes:

1. **Quality Assurance**: Ensures documentation remains accurate and accessible
2. **Portfolio Showcase**: Demonstrates CI/CD and DevOps practices
3. **Research Integrity**: Validates citation and publication information
4. **Visitor Experience**: Ensures links and documentation work correctly

## Running Workflows Manually

Most workflows can be triggered manually via GitHub Actions tab:

1. Go to repository → Actions
2. Select workflow from left sidebar
3. Click "Run workflow"
4. Select branch and click "Run workflow"

## Workflow Status

Check workflow status badges in the main README.md or visit the [Actions tab](https://github.com/domenicodigangi/MaximumEntropyNetworksErgmMatlab/actions).

## For Maintainers

These workflows are designed to be lightweight and not require:
- MATLAB license or installation
- External paid services
- Complex infrastructure

All workflows use free GitHub Actions runners and open-source tools.
