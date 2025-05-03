from setuptools import setup, find_packages
import os

# Get the absolute path to the directory containing setup.py
here = os.path.abspath(os.path.dirname(__file__))
# Get the path to the readme file
readme_path = os.path.join(here, "readme.md")

setup(
    name="library_analyzer",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas>=1.5.3",
        "numpy>=1.24",
        "matplotlib>=3.6",
        "seaborn",
        "scipy>=1.12",
        "statsmodels",
    ],
    author="John J Desmarais",
    author_email="john.j.desmarais@gmail.com",
    description="A package for analyzing MPRA library data",
    long_description=open(readme_path).read(),
    long_description_content_type="text/markdown",
    url="https://github.com/jackdesmarais/TRf_transposon_silencing_MPRA",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
) 