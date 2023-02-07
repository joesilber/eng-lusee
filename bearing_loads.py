'''Calculation of central bearing loads for LuSEE-Night lunar instrument's
rotation mechanism. Joe Silber, Lawrence Berkeley National Lab, 2023.'''
print(__doc__, '\n')

import numpy as np

# design factors
Grms_qual = 14.1  # per GEVs Table 2.4-3
proof_factor = 1.25  # per GEVs section 2.4.1.4.1
design_factor = 1.2  # covers uncertainties in matl props and test methods
peak_accel_factor = 5.0  # meant to represent ~ 5-sigma events
Grms_design = Grms_qual * proof_factor * design_factor
Gpeak_design = Grms_design * peak_accel_factor
total_safety_factor = Grms_design / Grms_qual
print(f'Design mean load level = {Grms_design:.1f} G(rms)')
print(f'Design peak load level = {Gpeak_design:.1f} G')
print(f'Total safety factor applied = {total_safety_factor:.2f}')

# units and such 
N_per_lbf = 4.448
Nm_per_ozin = 1 / 141.61193227806
Nm_per_inlb = Nm_per_ozin * 16
tab = f'{" ":4}'
earth_1G = 9.81  # m/s^2

# component sizes
z_cg = 0.011  # mm
payload_mass = 5.21  # kg
single_bearing_diam = 0.083  # pitch circle diameter of single bearing
duplex_bearing_sep = 0.040  # assumed axial separation distance between duplexed bearings (where applicable)

bearing_supported_load = payload_mass * earth_1G
bearing_supported_moment = bearing_supported_load * z_cg
bearing_loads_stationary = {'thrust': bearing_supported_load,  # nominally upright
                            'radial': bearing_supported_load,  # when tipped 90 deg
                            'moment': bearing_supported_moment,  # when tipped 90 deg
                           }
bearing_independent_loads = {'rms': {name: load * Grms_design for name, load in bearing_loads_stationary.items()},
                             'peak': {name: load * Gpeak_design for name, load in bearing_loads_stationary.items()},
                            }
Fr = bearing_independent_loads['peak']['radial']
Fa = bearing_independent_loads['peak']['thrust']
Fm_single = bearing_independent_loads['peak']['moment'] / (single_bearing_diam/2)
Fm_duplex = bearing_independent_loads['peak']['moment'] / duplex_bearing_sep
typ_relative_axial_load = Fa / (24 * 5.6**2)
print(f'Note relative axial load = Fa / (Z * D^2) ~ {typ_relative_axial_load:.1f} N/mm^2 = {typ_relative_axial_load/N_per_lbf*25.4**2:.1f} lbf/in^2'
      f'\nHere using approximations appropriate to this level (i.e. at 750-1000 lbf/in^2):\n')
bearing_factors = {'crossed roller': {'X1': 1.0, 'Y1': 0.44, 'e': 1.5, 'X2': 0.67, 'Y2': 0.67}, # per IKO catalog
                   'four-point contact': {'X1': 0.60, 'Y1': 1.07, 'e': 1.27, 'X2': 1.0, 'Y2': 0.0},  # X1, Y1 per SKF (contact angle 35°), https://evolution.skf.com/us/four-point-contact-ball-bearings-two-in-one
                   'radial deep groove': {'X1': 0.56, 'Y1': 1.0, 'e': 0.44, 'X2': 1.0, 'Y2': 0.0}, # per Minnebea catalog
                   'angular contact': {'X1': 0.44, 'Y1': 1.0, 'e': 0.56, 'X2': 1.0, 'Y2': 0.0},  # per Minnebea catalog (contact angle 15°)
                  }
n_bearings = {'crossed roller': 1,
              'four-point contact': 1,
              'radial deep groove': 2,
              'angular contact': 2,
             }
bearing_equiv_loads = {}
for key, factors in bearing_factors.items():
    if n_bearings[key] == 1:
        Fm = Fm_single
        external_preload = 0.
    elif n_bearings[key] == 2:
        Fm = Fm_duplex
        external_preload = Fa / design_factor
    else:
        assert False, f'error with bearing key {key}, single vs duplex unknown'
    possible_loads = {}
    for r in np.linspace(1e-9, 1-1e-9, 10):  # radial fraction of vector acceleration components
        a = (1 - r**2)**0.5  # axial fraction of vector acceleration components
        n = n_bearings[key]
        load_ratio = (a*Fa) / (r*Fr/n + r*Fm)
        x = factors['X1'] if load_ratio > factors['e'] else factors['X2']
        y = factors['Y1'] if load_ratio > factors['e'] else factors['Y2']
        radial_term = r * Fr * x / n
        axial_term = a * (Fa + external_preload) * y
        moment_term = r * Fm
        # note by convention these load terms add as a straight sum, *not* in quadrature
        possible_loads[load_ratio] = radial_term + axial_term + moment_term
    max_load = max(possible_loads.values())
    ratio_for_max_load = list(possible_loads.keys())[list(possible_loads.values()).index(max_load)]
    bearing_equiv_loads[key] = max_load
    print(f'{tab}{key:<18} ... X1 = {factors["X1"]:.2f}, Y1 = {factors["Y1"]:.2f}, e = {factors["e"]:.2f}, X2 = {factors["X2"]:.2f}, Y2 = {factors["Y2"]:.2f}')
    print(f'{tab}{tab}... max equiv load, {max_load:5.0f} N = {max_load/N_per_lbf:4.0f} lbf, occurs when accel vector is weighted like axial : radial ~ {ratio_for_max_load:.3f} : 1\n')
print('')

bearing_loads_str = f'Payload mass of {payload_mass:.2f} kg yields the following bearing design loads.'
bearing_loads_str += f' These include all safety factors.'
bearing_loads_str += f'\nIn calculating equivalent loads, the moment component is estimated such that single bearings (1x) are'
bearing_loads_str += f'\nassumed to have a diameter ~ {single_bearing_diam*1000:.0f} mm, while duplex configurations (2x bearings) are assumed to have'
bearing_loads_str += f'\nan axial separation distance of ~ {duplex_bearing_sep*1000:.0f} mm.\n'

for accel_type, loads in bearing_independent_loads.items():
    bearing_loads_str += f'\n{tab}{accel_type} radial (independent) = {loads["radial"]:.0f} N = {loads["radial"]/N_per_lbf:.0f} lbf'
    bearing_loads_str += f'\n{tab}{accel_type} thrust (independent) = {loads["thrust"]:.0f} N = {loads["thrust"]/N_per_lbf:.0f} lbf'
    bearing_loads_str += f'\n{tab}{accel_type} moment (independent) = {loads["moment"]:.0f} N*m = {loads["moment"]/Nm_per_inlb:.0f} in*lb'
bearing_loads_str += '\n'
for bearing_type, load in bearing_equiv_loads.items():
    bearing_loads_str += f'\n{tab}{bearing_type:<18} ({n_bearings[bearing_type]}x) ... peak equivalent load = {load:.0f} N = {load/N_per_lbf:.0f} lbf'
print(bearing_loads_str)
