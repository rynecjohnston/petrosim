# PetroSim
PetroSim is a collection of petrological simulation tools.

## Installation
PetroSim can be installed from the source or with `pip` as follows:

    pip install petrosim

## EC-AFC
The EC-AFC<sup>1</sup> model is a tool for simulating the evolution of the geochemical composition of a a magma body as it intrudes into country rock, assimilating wallrock into the melt and fractionally crystallizing out cumulates.

### Usage
An EC-AFC simulation is run with the standalone driver script using the command below, taking the initial parameters via a YAML formatted `input.in` file and outputting the results in the `output.csv` CSV file:

    python3 petrosim/models/ecafc/driver.py input.in output.csv

An example input file is located at `petrosim/models/ecafc/example.in`, and this file may be generated at any location by running the `petrosim/models/ecafc/utils/generate_example_input.py` utility script from the desired location. Running the simulation by default prints to screen the first and last three lines of the results, but this may be changed to a different number of *N* lines or turned off entirely with *N*=0 with the `-print <N>` flag in the command.

## References
1. Spera, F. J.; Bohrson, W. A. Energy-Constrained Open-System Magmatic Processes I. General Model and Energy-Constrained Assimilation and Fractional Crystallization (EC-AFC) Formulation *J. Petrol.*, **2001**, *42*, 999-1018, DOI: [10.1093/petrology/42.5.999](https://doi.org/10.1093/petrology/42.5.999)