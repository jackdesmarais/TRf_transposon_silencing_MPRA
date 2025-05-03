from setuptools import setup, find_packages

setup(
    name="library_analyzer",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas==1.5.3",
        "numpy==1.24",
        "matplotlib==3.6",
        "seaborn",
        "scipy==1.12",
        "statsmodels",
        # Documentation dependencies
        "sphinx==5.0.2",
        "sphinx_rtd_theme",
        "myst-parser",
    ],
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