verbosity = 'DEBUG'
type_calc = 'SA'
# Storage directory
directory = 'FLAMES'

# Cantera mechanism
mechanism = 'gri30.xml'

# Species to exclude (NOx chemistry for example)
species_exclude_init = ()
species_exclude_zero = ()
# Major species, not in the optimization loop
species_major = ()

# Fuel and N2-O2 ratio
fuel = 'CH4'
n2_o2_ratio = 3.76

# Threshold value for Zero
threshold = 1e-3

# Ai cases
# Pressure
P_ai = [1e5]
# Temperature
T_ai = [1200.0, 1400.0, 1600.0, 1800.0]
# Equivalence ratio
phi_ai = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
          1.3, 1.4]
# Error tolerance
tolerance_ai = 0.05

# Flame cases
# Pressure
P_fl = [1e5]
# Temperature
T_fl = [300.0]
# Equivalence ratio
phi_fl = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
          1.3, 1.4]
# Error tolerance
tolerance_fl = 0.05

# Initial beta vector
# x0 = [ 1.0 for i in range(43) ]
