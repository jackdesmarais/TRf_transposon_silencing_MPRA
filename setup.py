from setuptools import setup, find_packages

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
    author="Your Name",
    author_email="your.email@example.com",
    description="A package for analyzing MPRA library data",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/TRf_transposon_silencing_MPRA",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
) 