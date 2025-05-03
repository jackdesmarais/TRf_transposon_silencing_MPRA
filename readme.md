# TRf transposon silencing MPRA

This repository contains the code for analyzing the TRf transposon PBS MPRA data and generating the figures for the paper.

[![Documentation Status](https://github.com/jackdesmarais/TRf_transposon_silencing_MPRA/actions/workflows/docs.yml/badge.svg)](https://jackdesmarais.github.io/TRf_transposon_silencing_MPRA/)

## Installation

### Using pip
```bash
# Install base package
pip install -e .

# Install with dependencies for compiling the docs
pip install -e ".[docs]"

# Install with dependencies for running graphing notebooks
pip install -e ".[notebooks]"

# Install all dependencies
pip install -e ".[all]"
```

### Using conda
```bash
# Create and activate environment
git clone git@github.com:jackdesmarais/TRf_transposon_silencing_MPRA.git
cd TRf_transposon_silencing_MPRA
conda create -f ./graph_env.yml
conda activate TRf_transposon_silencing_MPRA
```

## Usage

The package provides tools for analyzing MPRA library data and generating figures. See the [full documentation](https://jackdesmarais.github.io/TRf_transposon_silencing_MPRA/) for detailed API reference and examples.

### Generate a library from reads
```python
from library_analyzer import Library

replicates = {'pool1_RNA1500':['pool1_DNA', 'pool1_RNA1500'],
              'pool2_RNA1500':['pool2_DNA','pool2_RNA1500'],
              'pool3_RNA1500':['pool3_DNA','pool3_RNA1500'],

              'pool1_RNA500':['pool1_DNA', 'pool1_RNA500'],
              'pool2_RNA500':['pool2_DNA','pool2_RNA500'],
              'pool3_RNA500':['pool3_DNA','pool3_RNA500'],

              'pool1_RNA250':['pool1_DNA', 'pool1_RNA250'],
              'pool2_RNA250':['pool2_DNA','pool2_RNA250'],
              'pool3_RNA250':['pool3_DNA','pool3_RNA250'],
            }

id_cols = ['category']
group_cols = ['sequence']
rate_method = 'l10fc'
alphabet = ['A','C','G','T']
WT_seq = 'TGGCGCCCGAACAGGGACCTGA'
process_call = 'category'
name_regex = '(?P<sequence>[ATCG]*)_(?:fillmd_)?'
data_file = '../data/original/allNOindel_1PBS_seq_MPRA.count.tsv'

lib = Library.build_from_reads(data_file,replicates, id_cols, group_cols, 
                               rate_method, alphabet, WT_seq, process_call, name_regex, 
                               sheet_name=None)

lib.save('../data/processed/allNOindel_1PBS_seq_MPRA.lib.pkl')
lib.data_df.to_csv('../data/processed/allNOindel_1PBS_seq_MPRA.processed.csv')
```

### Example Notebooks
The `graphing_notebooks` directory contains Jupyter notebooks demonstrating the analysis workflow:
- `process_data_to_library.ipynb`: Data processing pipeline
- `make_figure_panels.ipynb`: Main figure generation
- `make_supplemental_figures.ipynb`: Supplementary figure generation

## Documentation

Full documentation, including API reference and example notebooks, is available at:
https://jackdesmarais.github.io/TRf_transposon_silencing_MPRA/
