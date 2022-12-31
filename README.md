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

An example input file is located at `petrosim/models/ecafc/example.in`. This file may be generated at any location by running the utility script from the desired location in the following manner:

    python3 -m petrosim.models.ecafc.utils.generate_example_input

Running the simulation is by default silent and only outputs the CSV file, but this may be changed to also print to screen the first and last *N* lines of the results by adding the `-print <N>` flag in the command.

Basic visualization of the simulation results may be achieved with the included plotter, where one passes the resulting CSV file, the x-axis column name and the y-axis column name. An example of how to use the plotter on the resulting CSV file is below.

    python3 -m petrosim.models.ecafc.plotter output.csv -x Tm -y "Sr Cm"

This will both show a figure on screen and write a PNG file named `Tm-Sr_Cm.png`. To turn off saving to file, pass the `-nosave` flag, and to turn off showing the figure on screen, pass the `-noshow` flag in the command.

## References
1. Spera, F. J.; Bohrson, W. A. Energy-Constrained Open-System Magmatic Processes I. General Model and Energy-Constrained Assimilation and Fractional Crystallization (EC-AFC) Formulation *J. Petrol.*, **2001**, *42*, 999-1018, DOI: [10.1093/petrology/42.5.999](https://doi.org/10.1093/petrology/42.5.999)
