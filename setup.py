from setuptools import setup, find_packages

README = """
# PetroSim
PetroSim is a collection of petrological simulation tools.

## Installation
PetroSim can be installed from the source or with `pip` as follows:

    pip install petrosim

## EC-AFC
The EC-AFC<sup>1</sup> model is a tool for simulating the evolution of the geochemical composition of a magma body as it intrudes into country rock, assimilating wallrock into the melt and fractionally crystallizing out cumulates.

### Usage
An EC-AFC simulation is run with the standalone driver script using the command below, taking the initial parameters via a YAML formatted `input.in` file and outputting the results in the `output.csv` CSV file:

    python3 -m petrosim.models.ecafc.driver input.in output.csv

Please visit the [source](https://github.com/rynecjohnston/petrosim) for further reading.

## References
1. Spera, F. J.; Bohrson, W. A. Energy-Constrained Open-System Magmatic Processes I. General Model and Energy-Constrained Assimilation and Fractional Crystallization (EC-AFC) Formulation *J. Petrol.*, **2001**, *42*, 999-1018, DOI: [10.1093/petrology/42.5.999](https://doi.org/10.1093/petrology/42.5.999)
"""

setup(
    name='petrosim',
    version='0.0.3',
    author='Ryne C. Johnston',
    author_email='johnston.ryne@gmail.com',
    description='A collection of petrological simulation tools.',
    long_description=README,
    long_description_content_type="text/markdown",
    license='MIT',
    project_urls={
        "Source": "https://github.com/rynecjohnston/petrosim/",
        "Tracker": "https://github.com/rynecjohnston/petrosim/issues"
    },
    packages=[
        'petrosim',
        'petrosim.models',
        'petrosim.models.ecafc',
        'petrosim.models.ecafc.utils'
    ],
    package_data={
        '': [
            'LICENSE',
            'README.md',
        ],
        'petrosim': [
            'models/ecafc/example.in',
            'models/ecafc/test/*'
        ]
    },
    include_package_data=True,
    install_requires = [
        "setuptools>=61.0",
        "ruamel.yaml>=0.15",
        "voluptuous>=0.13",
        "numpy>=1.19",
        "matplotlib>=3.6.2",
        "pytest>=6.2"
    ],
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8'
)