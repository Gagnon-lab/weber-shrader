from setuptools import setup, find_packages

setup(
    name='testis_model',
    version='0.1.0',
    packages=find_packages(include=['testis_model', 'testis_model.*']),
    python_requires=">=3.11",
    install_requires=[
        "joblib==1.4.2",
        "matplotlib==3.10.8",
        "networkx==3.4.2",
        "numba==0.61.2",
        "numpy==2.2.5",
        "pandas==3.0.0",
        "pyfastx==2.3.0",
        "python-Levenshtein==0.27.3",
        "scipy==1.17.0",
        "seaborn==0.13.2",
        "tqdm==4.67.1",
    ],
)
