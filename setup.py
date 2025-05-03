from setuptools import setup, find_packages

# Core requirements that are always installed
install_requires = [
    "pandas==1.5.3",
    "numpy",
    "matplotlib==3.6",
    "seaborn==0.13.1",  # Pinned to match environment version
    "scipy==1.12",
    "statsmodels",
]

# Optional requirements for documentation
docs_require = [
    "sphinx==5.0.2",
    "sphinx_rtd_theme",
    "myst-parser",
    "nbsphinx",  # For rendering Jupyter notebooks
    "ipython",   # Required by nbsphinx
]

# Optional requirements for development
notebook_require = [
    "ipykernel",
    "scikit-learn",
    "tensorflow==2.15.0",  # Last version with integrated keras
    "h5py",  # Required for model saving/loading
    "mavenn==1.0.2",  # From bioconda channel
]

setup(
    name="library_analyzer",
    version="0.1.0",
    packages=find_packages(),
    install_requires=install_requires,
    extras_require={
        'docs': docs_require,  # Install with: pip install -e ".[docs]"
        'notebooks': notebook_require,    # Install with: pip install -e ".[notebooks]"
        'all': docs_require + notebook_require,  # Install with: pip install -e ".[all]"
    },
    author="John J Desmarais",
    author_email="john.j.desmarais@gmail.com",
    description="A package for analyzing MPRA library data",
    long_description="A package for analyzing MPRA library data",
    long_description_content_type="text/markdown",
    url="https://github.com/jackdesmarais/TRf_transposon_silencing_MPRA",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    # Prevent setuptools from trying to fetch remote data
    zip_safe=False,
    # Don't try to include non-Python files
    include_package_data=False,
) 