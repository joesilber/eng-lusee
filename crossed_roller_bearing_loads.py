'''Calculation of central bearing loads for LuSEE-Night lunar instrument's
rotation mechanism. Joe Silber, Lawrence Berkeley National Lab, 2023.'''
print(__doc__, '\n')

import numpy as np

# peak loads from random vibration FEA
# these are all 5-sigma (99.99994% probability envelope)
# FEA date 2023-06-15
nsigma = '5\u03C3'  # load probabilty envelope
Fr = 7291  # N, 5-sigma radial load
Fa = 8235  # N, 5-sigma axial load
Fm = 588  # N*m, 5-sigma moment load

# # units and such 
N_per_lbf = 4.448
# Nm_per_ozin = 1 / 141.61193227806
# Nm_per_inlb = Nm_per_ozin * 16
tab = f'{" ":4}'
# earth_1G = 9.81  # m/s^2

# # component sizes
# z_cg = 0.0003  # m,  payload center of gravity
# payload_mass = 5.13  # kg
# single_bearing_diam = 0.083  # pitch circle diameter of single bearing
# duplex_bearing_sep = 0.040  # assumed axial separation distance between duplexed bearings (where applicable)

# bearing_supported_load = payload_mass * earth_1G
# bearing_supported_moment = payload_mass * earth_1G * z_cg
# bearing_loads_stationary = {'thrust': bearing_supported_load,  # nominally upright
#                             'radial': bearing_supported_load,  # when tipped 90 deg
#                             'moment': bearing_supported_moment,  # when tipped 90 deg
#                            }
# bearing_independent_loads = {'rms': {name: load * Grms_design for name, load in bearing_loads_stationary.items()},
#                              'peak': {name: load * Gpeak_design for name, load in bearing_loads_stationary.items()},
#                             }

# Fm_single = bearing_independent_loads['peak']['moment'] / (single_bearing_diam/2)
# Fm_duplex = bearing_independent_loads['peak']['moment'] / duplex_bearing_sep
typ_relative_axial_load = Fa / (24 * 5.6**2)
print(f'Note relative axial load = Fa / (Z * D^2) ~ {typ_relative_axial_load:.1f} N/mm^2 = {typ_relative_axial_load/N_per_lbf*25.4**2:.1f} lbf/in^2'
      f'\nHere using approximations appropriate to this level (i.e. at 750-1000 lbf/in^2):\n')
bearing_factors = {'X1': 1.0, 'Y1': 0.44, 'e': 1.5, 'X2': 0.67, 'Y2': 0.67}  # per IKO catalog
possible_loads = {}
for r in np.linspace(1e-9, 1-1e-9, 10):  # radial fraction of vector acceleration components
    a = (1 - r**2)**0.5  # axial fraction of vector acceleration components
    load_ratio = (a*Fa) / (r*Fr + r*Fm)
    x = bearing_factors['X1'] if load_ratio > bearing_factors['e'] else bearing_factors['X2']
    y = bearing_factors['Y1'] if load_ratio > bearing_factors['e'] else bearing_factors['Y2']
    radial_term = r * Fr * x
    axial_term = a * Fa * y
    moment_term = r * Fm
    # note by convention these load terms add as a straight sum, *not* in quadrature
    possible_loads[load_ratio] = radial_term + axial_term + moment_term
max_equiv_load = max(possible_loads.values())
ratio_for_max_equiv_load = list(possible_loads.keys())[list(possible_loads.values()).index(max_equiv_load)]
print(f'{tab} ... X1 = {bearing_factors["X1"]:.2f}, Y1 = {bearing_factors["Y1"]:.2f}, e = {bearing_factors["e"]:.2f}, X2 = {bearing_factors["X2"]:.2f}, Y2 = {bearing_factors["Y2"]:.2f}')
print(f'{tab} ... max equiv load, {max_equiv_load:5.0f} N = {max_equiv_load/N_per_lbf:4.0f} lbf, occurs when accel vector is weighted like axial : radial ~ {ratio_for_max_equiv_load:.3f} : 1\n')
print('')

loads_str = f'{nsigma} loads:'
            f'\n{tab}Fa = {Fa} N (axial)'
            f'\n{tab}Fr = {Fr} N (radial)'
            f'\n{tab}Fm = {Fm} N*m (moment)'
            f'\nyield the following bearing design loads.'
bearing_loads_str += f' These include all safety factors.'
bearing_loads_str += f'\nIn calculating equivalent loads, the moment component is estimated such that single bearings (1x) are'
bearing_loads_str += f'\nassumed to have a diameter ~ {single_bearing_diam*1000:.0f} mm, while duplex configurations (2x bearings) are assumed to have'
bearing_loads_str += f'\nan axial separation distance of ~ {duplex_bearing_sep*1000:.0f} mm.\n'

for accel_type, loads in bearing_independent_loads.items():
    bearing_loads_str += f'\n{tab}{accel_type} radial (independent) = {loads["radial"]:.0f} N = {loads["radial"]/N_per_lbf:.0f} lbf'
    bearing_loads_str += f'\n{tab}{accel_type} thrust (independent) = {loads["thrust"]:.0f} N = {loads["thrust"]/N_per_lbf:.0f} lbf'
    bearing_loads_str += f'\n{tab}{accel_type} moment (independent) = {loads["moment"]:.1f} N*m = {loads["moment"]/Nm_per_inlb:.1f} in*lb'
bearing_loads_str += '\n'
for bearing_type, load in bearing_equiv_loads.items():
    bearing_loads_str += f'\n{tab}{bearing_type:<18} ({n_bearings[bearing_type]}x) ... peak equivalent load = {load:.0f} N = {load/N_per_lbf:.0f} lbf'
print(bearing_loads_str)
