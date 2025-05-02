# Library Analyzer

A Python package for analyzing MPRA (Massively Parallel Reporter Assay) library data.

## Installation

```bash
pip install .
```

## Usage

```python
from library_analyzer import Library
# Create a library object
lib = Library.build_from_reads(
data_file="your_data.csv",
replicates={"rep1": ["t0", "t1"]},
id_cols=["variant_id"],
group_cols=["group"],
rate_method="l2fc",
alphabet=["A", "C", "G", "T"],
WT_seq="ACGT",
process_call="process",
name_regex="(.)"
)
# Analyze library skew
fig = lib.check_skew()
```
## Features

- Library data analysis
- Statistical calculations
- Visualization tools
- Replicate analysis