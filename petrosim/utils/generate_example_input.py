import argparse

EXAMPLE_IN = 'example.in'
CONTENTS_ECAFC = """# This is the input file for the EC-AFC simulation in YAML format.
# YAML allows for comments, ignoring text following the # symbol, like this and
# the parameters' descriptions. While the BULK parameters are required, the
# TRACE parameters optional.

BULK:
    MAGMA:
        M0:     1       # Normalized mass of melt in magma body
        T0:     1320    # Initial standing melt temperature (ºC)
        Tl:     1320    # Pristine standing melt liquidus temperature (ºC)
        cp:     1484    # Specific heat capacity of magma (J K-1 kg-1)
        dh:     3.96e5  # Enthalpy of crystallization (dh_cry) (J kg-1)
        loga:   450     # Logistic parameter a
        logb:   -10.5   # Logisitc parameter b
        D0:     1.5     # Bulk distribution coefficient between cumulate and melt
        dH:     0       # Enthalpy of cumulate/standing melt partition coefficient (J mol-1)
        d0:     1       # Oxygen isotopic composition of standing magma

    ASSIMILANT:
        M0:     0       # Normalized mass of anatectic melt
        T0:     600     # Initial wall rock temperature (ºC)
        Tl:     1100    # Wall rock liquidus temperature (ºC)
        cp:     1388    # Specific heat capacity of assimilant (J K-1 kg-1)
        dh:     3.54e5  # Enthalpy of fusion (dh_fus) (J kg-1)
        loga:   400     # Logistic parameter a
        logb:   -11     # Logisitc parameter b
        D0:     1.5     # Bulk distribution coefficient between wall rock and anatectic melt
        dH:     0       # Enthalpy of wall rock/anatectic melt partition coefficient (J mol-1)
        d0:     1       # Oxygen isotopic composition of assimilant

    OTHER:
        Ts:       950     # Country rock solidus (ºC)
        dT:       0.001   # Temperature step (K)
        Teq:              # Equilibration temperature (ºC)
        Teq_norm: 0.827   # Normalized equilibration temperature (dimensionless)
        koxy:     1       # Ratio of mass fraction of oxygen in bulk assimilant and pristine magma

TRACE:
    - elem:       Sr         # Element
      isoratio:   87Sr/86Sr  # Isotopic ratio
      MAGMA:
          C0:     700      # Concentration of trace element in standing melt (ppm)
          e0:     0.7035   # Isotopic ratio in standing melt
      ASSIMILANT:
          C0:     350      # Concentration of trace element in assimilant (ppm)
          e0:     0.722    # Isotopic ratio in assimilant

    - elem:       Nd          # Element
      isoratio:   143Nd/144Nd # Isotopic ratio
      MAGMA:
          C0:     35       # Concentration of trace element in standing melt (ppm)
          e0:     0.513    # Isotopic ratio in standing melt
      ASSIMILANT:
          C0:     26       # Concentration of trace element in assimilant (ppm)
          e0:     0.5118   # Isotopic ratio in assimilant
"""

CONTENTS_ECRAFC = """# This is the input file for the EC-RAFC simulation in YAML format.
# YAML allows for comments, ignoring text following the # symbol, like this and
# the parameters' descriptions. While the BULK parameters are required, the
# TRACE parameters optional.

BULK:
    MAGMA:
        M0:     1       # Normalized mass of melt in magma body
        T0:     1300    # Initial standing melt temperature (ºC)
        Tl:     1300    # Pristine standing melt liquidus temperature (ºC)
        cp:     1484    # Specific heat capacity of magma (J K-1 kg-1)
        dh:     3.96e5  # Enthalpy of crystallization (dh_cry) (J kg-1)
        loga:   450     # Logistic parameter a
        logb:   -11     # Logisitc parameter b
        D0:     1.5     # Bulk distribution coefficient between cumulate and melt
        dH:     0       # Enthalpy of cumulate/standing melt partition coefficient (J mol-1)
        d0:     6.2     # Oxygen isotopic composition of standing magma

    ASSIMILANT:
        M0:     0       # Normalized mass of anatectic melt
        T0:     600     # Initial wall rock temperature (ºC)
        Tl:     1150    # Wall rock liquidus temperature (ºC)
        cp:     1370    # Specific heat capacity of assimilant (J K-1 kg-1)
        dh:     2.70e5  # Enthalpy of fusion (dh_fus) (J kg-1)
        loga:   400     # Logistic parameter a
        logb:   -11     # Logisitc parameter b
        D0:     0.05    # Bulk distribution coefficient between wall rock and anatectic melt
        dH:     0       # Enthalpy of wall rock/anatectic melt partition coefficient (J mol-1)
        d0:     12      # Oxygen isotopic composition of assimilant

    RECHARGE:
        M0:     0.3     # Normalized mass of recharge magma
        T0:     1320    # Initial recharge magma temperature (ºC)
        Tl:     1320    # Recharge magma liquidus temperature (ºC)
        cp:     1484    # Specific heat capacity of recharge magma (J K-1 kg-1)
        dh:     3.96e5  # Enthalpy of crystallization(dh_cry) (J kg-1)
        loga:   450     # Logistic parameter a
        logb:   -11     # Logisitc parameter b
        D0:     1.5     # Bulk distribution coefficient between enclave and recharge melt
        dH:     0       # Enthalpy of enclave/recharge melt partition coefficient (J mol-1)
        d0:     7       # Oxygen isotopic composition of recharge magma

    OTHER:
        Ts:       900   # Country rock solidus (ºC)
        dT:       0.0005 # Temperature step (K)
        Teq:            # Equilibration temperature (ºC)
        Teq_norm: 0.904 # Normalized equilibration temperature (dimensionless)
        koxy:     1     # Ratio of mass fraction of oxygen in bulk assimilant and pristine magma

RECHARGE_PARAMS:
    Nrechg:     3       # Number of recharges
    linear:     False   # Whether recharge event is linear (True) or episodic (False)
    # The lists below are only for episodic recharge events and are ignored if linear recharge is set to True
    masses:             # masses of each recharge event
        - 0.1
        - 0.1
        - 0.1
    temps:              # masses of each recharge event
        - 1250
        - 1175
        - 1155
    widths:             # widths of the recharge function for each recharge event
        - 0.3
        - 0.3
        - 0.3
    slopes:             # slopes of the recharge function for each recharge event
        - 2
        - 2
        - 2

TRACE:
    - elem:       Sr         # Element
      isoratio:   87Sr/86Sr  # Isotopic ratio
      MAGMA:
          C0:     700      # Concentration of trace element in standing melt (ppm)
          e0:     0.705    # Isotopic ratio in standing melt
      ASSIMILANT:
          C0:     230      # Concentration of trace element in assimilant (ppm)
          e0:     0.720    # Isotopic ratio in assimilant
      RECHARGE:
          C0:     1000     # Concentration of trace element in recharge magma (ppm)
          e0:     0.703    # Isotopic ratio in recharge magma
"""


def parse_args():
    p = argparse.ArgumentParser(description="Generate example input files for EC models.")
    p.add_argument(
        'model_type',
        choices=['ecafc', 'ecrafc'],
        help='Which EC model input file is desired?'
    )
    return p.parse_args()


def write(model_type):
    contents = {'ecafc': CONTENTS_ECAFC, 'ecrafc': CONTENTS_ECRAFC}
    with open(EXAMPLE_IN, 'w') as file_handle:
        file_handle.write(contents[model_type])


if __name__ == '__main__':
    args = parse_args()
    write(args.model_type)