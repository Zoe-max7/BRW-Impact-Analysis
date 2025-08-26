# 🧬 BRW-Impact-Analysis: Biological Random Walks Impact Assessment

[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Streamlit](https://img.shields.io/badge/Streamlit-Web%20App-red.svg)](https://streamlit.io/)
[![GitHub](https://img.shields.io/badge/GitHub-Repository-black.svg)](https://github.com/Zoe-max7/BRW-Impact-Analysis)

> **Comprehensive framework for evaluating the impact of biological information sources in Biological Random Walks algorithm**

## Table of Contents

- [Overview](#overview)
- [Research Objective](#research-objective)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage Methods](#usage-methods)
- [Data Requirements](#data-requirements)
- [Results Interpretation](#results-interpretation)
- [Project Structure](#project-structure)
- [Algorithm Details](#algorithm-details)
- [Contributing](#contributing)
- [License](#license)
- [Citation](#citation)
- [Contact](#contact)

## Overview

**BRW-Impact-Analysis** is a comprehensive research framework designed to systematically evaluate the contribution and impact of different biological information sources within the Biological Random Walks (BRW) algorithm for disease gene prioritization.

This project extends the original BRW implementation by providing:
- **Interactive GUI** for parameter exploration
- **Automated batch processing** for systematic analysis
- **Comprehensive ablation studies** to assess data source contributions
- **Advanced visualization** for result interpretation
- **Reproducible research** workflows

## Research Objective

The primary goal is to **quantify and compare the impact** of different biological data sources on gene prioritization accuracy:

### **Data Sources Under Investigation:**
1. **Protein-Protein Interaction (PPI) Networks** - Core network structure
2. **Co-expression Networks** - Gene expression correlations
3. **Differentially Expressed (DE) Genes** - Disease-specific expression changes
4. **Ontology Annotations** - Functional and pathway information

### **Research Questions:**
- Which data source contributes most to prediction accuracy?
- How do different combinations affect results?
- What are the optimal parameters for each data type?
- How does the impact vary across different cancer types?

## Features

### ** Core Algorithm**
- **Biological Random Walks** implementation with restart mechanism
- **Multi-source data integration** (PPI, co-expression, DE genes, ontologies)
- **Flexible parameter tuning** (alpha, beta, restart probability)
- **Cancer-specific analysis** (7 cancer types supported)

### ** Interactive GUI**
- **Streamlit-based web interface** - No command-line knowledge required
- **Real-time parameter adjustment** - Interactive exploration
- **Multiple alpha-beta configurations** - Systematic testing
- **Live result visualization** - Charts, heatmaps, statistics
- **Export functionality** - CSV downloads and publication-ready figures

### ** Batch Processing**
- **Automated execution** across all cancer types
- **Comprehensive parameter grids** - Systematic exploration
- **Progress tracking** - Real-time monitoring
- **Result aggregation** - Summary reports and statistics

### ** Analysis Tools**
- **Ablation studies** - Remove individual data sources
- **Parameter optimization** - Find optimal configurations
- **OncoKB validation** - Compare against known cancer genes
- **Statistical analysis** - Hit rates, averages, comparisons

## Installation

### **Prerequisites**
- Python 3.7 or higher
- pip package manager
- Git (for cloning)

### **Step 1: Clone Repository**
```bash
git clone https://github.com/Zoe-max7/BRW-Impact-Analysis.git
cd BRW-Impact-Analysis
```

### **Step 2: Install Dependencies**
```bash
pip install -r requirements.txt
```

### **Step 3: Prepare Data**
- Place your biological data files in the `data_set/` directory
- Ensure `Dataset OncoKB.xlsx` is in the root directory
- See [Data Requirements](#data-requirements) for detailed file specifications

## Quick Start

### **Option 1: Interactive GUI (Recommended)**
```bash
streamlit run app_gui.py
```
- Opens web interface at `http://localhost:8501`
- Select cancer type and parameters
- Run analysis with real-time results
- Export findings for further analysis

### **Option 2: Command Line**
```bash
python main.py -p <ppi_network> -s <seed_set> -a <ontology> -do <disease_ontology> -o <output>
```

### **Option 3: Batch Processing**
```bash
python run_brw_batch.py
```
- Automatically processes all cancer types
- Tests multiple parameter combinations
- Generates comprehensive summary reports

## Usage Methods

### **GUI Workflow**
1. **Select Cancer Type** - Choose from 7 available types
2. **Configure Parameters** - Set ablation options and weights
3. **Run Analysis** - Execute BRW algorithm
4. **View Results** - Interactive charts and tables
5. **Export Data** - Download results for publication

### **Ablation Study Options**
- **FULL**: All data sources included
- **noCOEXP**: Exclude co-expression networks
- **noDE**: Exclude differentially expressed genes
- **noONTO**: Exclude ontology data
- **PPI_only**: Only protein-protein interactions

### **Parameter Tuning**
- **Alpha (α)**: Weight for PPI network (0.0 - 1.0)
- **Beta (β)**: Weight for co-expression network (0.0 - 1.0)
- **Restart Probability**: Random walk restart parameter (0.1 - 0.99)

## Data Requirements

### **Required Files Structure**
```
data_set/
├── ppi_network/
│   └── HIPPIE.tsv                    # Protein-protein interactions
├── co-expression_networks/
│   ├── TCGA-BRCA__co_expression__t_70%.tsv
│   ├── TCGA-COAD__co_expression__t_70%.tsv
│   └── ...                           # Other cancer types
├── differentially_expressed_genes/
│   ├── TCGA-BRCA_de_genes.tsv
│   ├── TCGA-COAD_de_genes.tsv
│   └── ...                           # Other cancer types
├── ontology/
│   └── ontology_graph.txt            # Gene-ontology associations
├── disease_specific_ontologies/
│   ├── TCGA-BRCA_disease_ontologies.txt
│   ├── TCGA-COAD_disease_ontologies.txt
│   └── ...                           # Other cancer types
└── seed_set/
    ├── TCGA-BRCA_seed.txt
    ├── TCGA-COAD_seed.txt
    └── ...                           # Other cancer types
```

### **File Formats**
- **PPI Networks**: TSV with `ensembl_id_1 \t ensembl_id_2`
- **Co-expression**: TSV with `ensembl_id_1 \t ensembl_id_2 \t correlation_score`
- **DE Genes**: TSV with one `ensembl_id` per row
- **Ontology**: TSV with `ensembl_id \t annotation_id \t dataset_name`
- **Seed Sets**: TXT with one `ensembl_id` per row

### **OncoKB Database**
- **File**: `Dataset OncoKB.xlsx` (root directory)
- **Purpose**: Validate gene prioritization results
- **Format**: Excel file with cancer gene annotations

## Project Structure

```
BRW-Impact-Analysis/
├──  app_gui.py                      # Streamlit GUI application
├──  main.py                         # Core BRW execution script
├──  run_brw_batch.py                # Batch processing automation
├──  check_genes.py                  # OncoKB validation script
├──  biological_random_walks/        # Core algorithm implementation
│   ├── core/                          # PageRank and core algorithms
│   ├── loader/                        # Data loading utilities
│   ├── matrix_creation/               # Network combination methods
│   └── personalization_vector_creation/ # Biological teleporting
├──  data_preprocessing/             # Data preparation scripts
├──  data_set/                       # Input data directory
│   ├── ppi_network/                   # Protein interaction data
│   ├── co-expression_networks/        # Gene expression correlations
│   ├── differentially_expressed_genes/ # Disease-specific genes
│   ├── ontology/                      # Functional annotations
│   └── seed_set/                      # Known disease genes
├──  outputs/                        # Results and analysis
├──  requirements.txt                # Python dependencies
├──  README.md                       # Detailed documentation
├──  LICENSE                         # MIT License
└──  README_GITHUB.md               # This file
```

### **Development Setup**
```bash
git clone https://github.com/Zoe-max7/BRW-Impact-Analysis.git
cd BRW-Impact-Analysis
pip install -r requirements.txt
pip install -e .  # For development mode
```

## Getting Started Checklist

- [ ] **Clone repository** - `git clone https://github.com/Zoe-max7/BRW-Impact-Analysis.git`
- [ ] **Install dependencies** - `pip install -r requirements.txt`
- [ ] **Prepare data files** - Place in `data_set/` directory
- [ ] **Launch GUI** - `streamlit run app_gui.py`
- [ ] **Run first analysis** - Test with default parameters
- [ ] **Explore results** - Use visualization tools
- [ ] **Customize parameters** - Try different configurations
- [ ] **Export findings** - Download results for publication

---

---

*Last updated: January 2025*
*Version: 1.0.0*
