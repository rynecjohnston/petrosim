from setuptools import setup, find_packages

README = """
i# PetroSim
PetroSim is a collection of petrological simulation tools.

## Installation
PetroSim can be installed from the source or with `pip` as follows:

    pip install petrosim

## EC-AFC
The EC-AFC<sup>1</sup> model is a tool for simulating the evolution of the geochemical composition of a magma body as it intrudes into country rock, assimilating wallrock into the melt and fractionally crystallizing out cumulates.

## EC-RAFC
The EC-RAFC<sup>2</sup> model is a tool for simulating the evolution of the geochemical composition of a magma body as it and subsequent recharges intrude into country rock, assimilating wallrock into the melt and fractionally crystallizing out cumulates.

### Usage
An EC-AFC simulation is run with the standalone driver script using the command below, taking the initial parameters via a YAML formatted `input.in` file and outputting the results in the `output.csv` CSV file:

    python3 -m petrosim.models.ecafc.driver input.in output.csv


Likewise, an EC-RAFC simulation can be run by replacing `petrosim.models.ecafc.driver` with `petrosim.models.ecrafc.driver` in the command above.

Example input files are located at `petrosim/models/<model_type>/example.in`. This file may be generated at any location by running the utility script from the desired location in the following manner by specifyting the model type (either "ecafc" or "ecrafc"):

    python3 -m petrosim.utils.generate_example_input <model_type>

Running the simulation is by default silent and only outputs the CSV file, but this may be changed to also print to screen the first and last *N* lines of the results by adding the `-print <N>` flag in the command.

Basic visualization of the simulation results may be achieved with the included plotter, where one passes the resulting CSV file, the x-axis column name and the y-axis column name. An example of how to use the plotter on the resulting CSV file is below.

    python3 -m petrosim.results.plotter output.csv -x Tm -y "Sr Cm"

This will both show a figure on screen and write a PNG file named `Tm-Sr_Cm.png`. To turn off saving to file, pass the `-nosave` flag, and to turn off showing the figure on screen, pass the `-noshow` flag in the command.

## References
1. Spera, F. J.; Bohrson, W. A. Energy-Constrained Open-System Magmatic Processes I. General Model and Energy-Constrained Assimilation and Fractional Crystallization (EC-AFC) Formulation *J. Petrol.*, **2001**, *42*, 999-1018, DOI: [10.1093/petrology/42.5.999](https://doi.org/10.1093/petrology/42.5.999)
2. Spera, F. J.; Bohrson, W. A. Energy-Constrained Open-System Magmatic Processes 3. Energy-Constrained Recharge, Assimilation and Fractional Crystallization (EC-RAFC) *Geochem. Geophys. Geosyst.*, **2002**, *3*, 8001, DOI: [10.1029/2007GC001781](https://doi.org/10.1029/2007GC001781)
"""

setup(
    name='petrosim',
    version='0.0.4',
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
        'petrosim.models.ecrafc',
        'petrosim.results',
        'petrosim.utils',
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
