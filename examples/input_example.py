# Storage directory
directory = 'FLAMES'

# Cantera mechanism
mechanism = 'kin20285-chem.cti'

# Species to exclude (NOx chemistry for example)
species_exclude = ()

# Major species, not in the optimization loop
species_major = ('N2', 'AR', 'CH3OCH3', 'CO2', 'H2O',
                 'O2', 'H2', 'CO', 'OH', 'HO2')

# Fuel and N2-O2 ratio
fuel = 'CH3OCH3'
n2_o2_ratio = 3.76

# Threshold value for Zero
threshold = 1e-3

# Ai cases
# Pressure
P_ai = [1e5]
# Temperature
T_ai = [1200.0, 1400.0, 1600.0]
# Equivalence ratio
phi_ai = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
          1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
# Error tolerance
tolerance_ai = 0.20

# Flame cases
# Pressure
P_fl = [1e5]
# Temperature
T_fl = [300.0]
# Equivalence ratio
phi_fl = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
          1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
# Error tolerance
tolerance_fl = 0.05

# Initial beta vector
x0 = [9.99999944e-01, 9.99999878e-01, 9.99999904e-01, 9.99999960e-01
, 5.37093170e-01, 5.08976098e-01, 6.01165069e-08, 6.86705244e-08
, 6.47897171e-08, 6.17028363e-02, 9.99999914e-01, 3.20745868e-08
, 9.99999951e-01, 1.02802446e-03, 1.36235472e-07, 9.99999694e-01
, 8.76082885e-08, 2.75813114e-02, 5.54320308e-08, 1.07509611e-07
, 7.68159552e-08, 9.66355456e-08, 5.96294581e-08, 1.07365982e-07
, 8.39643376e-08, 9.60614566e-08, 1.03017202e-07, 9.21567832e-08
, 9.99999793e-01, 7.08789395e-08, 1.06592199e-07, 4.66776158e-08
, 9.47730344e-08, 1.47151161e-07, 1.10629874e-07, 1.07826961e-07
, 7.89858822e-08, 1.05792538e-07, 1.01698189e-07, 9.23493802e-08
, 9.84981048e-08, 1.01632182e-07, 8.47171651e-08, 1.06656556e-07
, 1.05745840e-07]
