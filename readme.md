# TRf transposon silencing MPRA

This repository contains the code for analyzing the TRf transposon PBS MPRA data and generating the figures for the paper.

[![Documentation Status](https://github.com/jackdesmarais/TRf_transposon_silencing_MPRA/actions/workflows/docs.yml/badge.svg)](https://jackdesmarais.github.io/TRf_transposon_silencing_MPRA/)

## Installation

### Using pip
```bash
# Clone repo
git clone git@github.com:jackdesmarais/TRf_transposon_silencing_MPRA.git
cd TRf_transposon_silencing_MPRA

# Install base package
pip install -e .

# Install with dependencies for compiling the docs
pip install -e ".[docs]"

# Install with dependencies for running graphing notebooks
pip install -e ".[notebooks]"

# Install all dependencies
pip install -e ".[all]"
```

Time for pip installation (with notebooks dependencies): 43 seconds
This estimate may be a little high due to pip caching speeding up the installation.
If you already have cached versions of the dependencies, expect 15-30 second 
speed ups.

### Using conda
```bash
# Create and activate environment
git clone git@github.com:jackdesmarais/TRf_transposon_silencing_MPRA.git
cd TRf_transposon_silencing_MPRA
conda create -f ./graph_env.yml
conda activate TRf_transposon_silencing_MPRA
```
Time for CONDA-based environment creation: 36 seconds
This estimate may be a little low due to pip caching speeding up the installation.
If you do not already have cached versions of the dependencies, expect an extra
15-30 seconds.

### Instalation timing
```bash
# Instalation timings were produced with the following command.

./create_and_time_env.sh | tee install_timing_results.txt
```

### Example Notebooks
The `graphing_notebooks` directory contains Jupyter notebooks demonstrating the analysis workflow:
- `process_data_to_library.ipynb`: Data processing pipeline
  - Expected runtime: 0.09 minutes
- `make_figure_panels.ipynb`: Main figure generation
  - Expected runtime if MAVENN model is not retrained: 0.33 minutes
  - Expected runtime if MAVENN model is retrained: 0.37 minutes
- `make_supplemental_figures.ipynb`: Supplementary figure generation
  - Expected runtime: 0.31 minutes

These notebooks provide a demonstration on applying this code to our data as well as reproducing all of the 
figures that derive from this code with the expected output on our provided dataset.

These notebooks can also be found in the [documentation](https://jackdesmarais.github.io/TRf_transposon_silencing_MPRA/).

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


## Documentation

Full documentation, including API reference and example notebooks, is available at:
https://jackdesmarais.github.io/TRf_transposon_silencing_MPRA/


## Environment
This code was run on a 2020 13-inch M1 MacBook Pro running MacOS 14.0
This code was run in an environment with the following package versions:
# packages in environment at TRf_transposon_silencing_MPRA:
#
# Name                    Version                   Build  Channel
absl-py                   2.2.2                    pypi_0    pypi
appnope                   0.1.4                    pypi_0    pypi
asttokens                 3.0.0                    pypi_0    pypi
astunparse                1.6.3                    pypi_0    pypi
ca-certificates           2025.4.26            hbd8a1cb_0    conda-forge
cachetools                5.5.2                    pypi_0    pypi
certifi                   2025.4.26                pypi_0    pypi
charset-normalizer        3.4.2                    pypi_0    pypi
comm                      0.2.2                    pypi_0    pypi
contourpy                 1.3.0                    pypi_0    pypi
cycler                    0.12.1                   pypi_0    pypi
debugpy                   1.8.14                   pypi_0    pypi
decorator                 5.2.1                    pypi_0    pypi
docutils                  0.18.1                   pypi_0    pypi
exceptiongroup            1.2.2                    pypi_0    pypi
executing                 2.2.0                    pypi_0    pypi
flatbuffers               25.2.10                  pypi_0    pypi
fonttools                 4.57.0                   pypi_0    pypi
gast                      0.4.0                    pypi_0    pypi
google-auth               2.39.0                   pypi_0    pypi
google-auth-oauthlib      1.0.0                    pypi_0    pypi
google-pasta              0.2.0                    pypi_0    pypi
grpcio                    1.71.0                   pypi_0    pypi
h5py                      3.13.0                   pypi_0    pypi
idna                      3.10                     pypi_0    pypi
importlib-metadata        8.7.0                    pypi_0    pypi
ipykernel                 6.29.5                   pypi_0    pypi
ipython                   8.18.1                   pypi_0    pypi
jedi                      0.19.2                   pypi_0    pypi
joblib                    1.4.2                    pypi_0    pypi
jupyter-client            8.6.3                    pypi_0    pypi
jupyter-core              5.7.2                    pypi_0    pypi
keras                     2.15.0                   pypi_0    pypi
kiwisolver                1.4.7                    pypi_0    pypi
libclang                  18.1.1                   pypi_0    pypi
libcxx                    20.1.4               ha82da77_0    conda-forge
libffi                    3.4.6                h1da3d7d_1    conda-forge
library-analyzer          0.1.0                     dev_0    <develop>
libzlib                   1.2.13               hfb2fe0b_6    conda-forge
logomaker                 0.8.7                    pypi_0    pypi
markdown                  3.8                      pypi_0    pypi
markdown-it-py            2.2.0                    pypi_0    pypi
markupsafe                3.0.2                    pypi_0    pypi
matplotlib                3.6.0                    pypi_0    pypi
matplotlib-inline         0.1.7                    pypi_0    pypi
mavenn                    1.0.2                    pypi_0    pypi
mdit-py-plugins           0.3.5                    pypi_0    pypi
mdurl                     0.1.2                    pypi_0    pypi
ml-dtypes                 0.2.0                    pypi_0    pypi
myst-parser               1.0.0                    pypi_0    pypi
namex                     0.0.9                    pypi_0    pypi
ncurses                   6.5                  h5e97a16_3    conda-forge
nest-asyncio              1.6.0                    pypi_0    pypi
numpy                     1.23.5                   pypi_0    pypi
oauthlib                  3.2.2                    pypi_0    pypi
openssl                   3.5.0                h81ee809_0    conda-forge
opt-einsum                3.4.0                    pypi_0    pypi
optree                    0.15.0                   pypi_0    pypi
packaging                 25.0                     pypi_0    pypi
pandas                    1.5.3                    pypi_0    pypi
parso                     0.8.4                    pypi_0    pypi
patsy                     1.0.1                    pypi_0    pypi
pexpect                   4.9.0                    pypi_0    pypi
pillow                    11.2.1                   pypi_0    pypi
pip                       25.1               pyhc872135_2  
platformdirs              4.3.7                    pypi_0    pypi
prompt-toolkit            3.0.51                   pypi_0    pypi
protobuf                  4.25.7                   pypi_0    pypi
psutil                    7.0.0                    pypi_0    pypi
ptyprocess                0.7.0                    pypi_0    pypi
pure-eval                 0.2.3                    pypi_0    pypi
pyasn1                    0.6.1                    pypi_0    pypi
pyasn1-modules            0.4.2                    pypi_0    pypi
pygments                  2.19.1                   pypi_0    pypi
pyparsing                 3.2.3                    pypi_0    pypi
python                    3.9.21               hb885b13_1  
python-dateutil           2.9.0.post0              pypi_0    pypi
pytz                      2025.2                   pypi_0    pypi
pyyaml                    6.0.2                    pypi_0    pypi
pyzmq                     26.4.0                   pypi_0    pypi
readline                  8.2                  h1d1bf99_2    conda-forge
requests                  2.32.3                   pypi_0    pypi
requests-oauthlib         2.0.0                    pypi_0    pypi
rich                      14.0.0                   pypi_0    pypi
rsa                       4.9.1                    pypi_0    pypi
scikit-learn              1.6.1                    pypi_0    pypi
scipy                     1.12.0                   pypi_0    pypi
seaborn                   0.13.1                   pypi_0    pypi
setuptools                78.1.1           py39hca03da5_0  
six                       1.17.0                   pypi_0    pypi
sphinx                    5.0.2                    pypi_0    pypi
sphinx-rtd-theme          2.0.0                    pypi_0    pypi
sphinxcontrib-jquery      4.1                      pypi_0    pypi
sqlite                    3.45.3               h80987f9_0  
stack-data                0.6.3                    pypi_0    pypi
statsmodels               0.14.4                   pypi_0    pypi
tensorboard               2.15.2                   pypi_0    pypi
tensorboard-data-server   0.7.2                    pypi_0    pypi
tensorflow                2.15.0                   pypi_0    pypi
tensorflow-estimator      2.15.0                   pypi_0    pypi
tensorflow-io-gcs-filesystem 0.37.1                   pypi_0    pypi
tensorflow-macos          2.15.0                   pypi_0    pypi
termcolor                 3.1.0                    pypi_0    pypi
threadpoolctl             3.6.0                    pypi_0    pypi
tk                        8.6.14               h6ba3021_0  
tornado                   6.4.2                    pypi_0    pypi
traitlets                 5.14.3                   pypi_0    pypi
typing-extensions         4.5.0                    pypi_0    pypi
tzdata                    2025b                h78e105d_0    conda-forge
urllib3                   2.4.0                    pypi_0    pypi
wcwidth                   0.2.13                   pypi_0    pypi
werkzeug                  3.1.3                    pypi_0    pypi
wheel                     0.45.1           py39hca03da5_0  
wrapt                     1.14.1                   pypi_0    pypi
xz                        5.6.4                h80987f9_1  
zipp                      3.21.0                   pypi_0    pypi
zlib                      1.2.13 