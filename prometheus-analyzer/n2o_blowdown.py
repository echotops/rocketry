import numpy as np
from scipy.optimize import minimize_scalar
import CoolProp.CoolProp as CP

PI = np.pi
DISCHARGE_COEFFICIENT = 0.75

"""
Iteratively solve for NOX mass flow rate by matching injector pressure drop
"""
def get_n2o_mdot(N2O_ID, N2O_OD, n2o_pressure, chamber_pressure):
    area = PI * ((N2O_OD / 2)**2 - (N2O_ID / 2)**2)
    return mdot(area, n2o_pressure, chamber_pressure)

def get_n2o_injector_density(N2O_ID, N2O_OD, n2o_pressure, chamber_pressure):
    area = PI * ((N2O_OD / 2)**2 - (N2O_ID / 2)**2)
    upstream_entropy = CP.PropsSI('S', 'P', n2o_pressure, 'Q', 0, 'N2O')
    critical_chamber_pressure = mdot_hem_crit(area, n2o_pressure)
    effective_downstream_pressure = max(chamber_pressure, critical_chamber_pressure)  # Clamp to critical pressure
    return CP.PropsSI('D', 'P', effective_downstream_pressure, 'S', upstream_entropy, 'N2O')

# Single-phase incompressible mass flow rate
def mdot_spi(area, upstream_pressure, downstream_pressure):
    pressure_difference = upstream_pressure - downstream_pressure
    rho_liq = CP.PropsSI('D', 'P', upstream_pressure, 'Q', 0, 'N2O')
    return DISCHARGE_COEFFICIENT * area * np.sqrt(2 * rho_liq * pressure_difference)

# Homogenous equilibrium model mass flow rate
def mdot_hem(area, upstream_pressure, downstream_pressure):
    upstream_enthalpy = CP.PropsSI('H', 'P', upstream_pressure, 'Q', 0, 'N2O')
    upstream_entropy = CP.PropsSI('S', 'P', upstream_pressure, 'Q', 0, 'N2O')
    downstream_enthalpy = CP.PropsSI('H', 'P', downstream_pressure, 'S', upstream_entropy, 'N2O')
    rho_chamber = CP.PropsSI('D', 'P', downstream_pressure, 'S', upstream_entropy, 'N2O')
    return DISCHARGE_COEFFICIENT * area * rho_chamber * np.sqrt(2 * (upstream_enthalpy - downstream_enthalpy))

# Find maximum of HEM mass flow rate
def mdot_hem_crit(area, upstream_pressure, downstream_lower_bound=None, downstream_upper_bound=None):
    if downstream_upper_bound is None:
        downstream_upper_bound = 0.999 * upstream_pressure
    if downstream_lower_bound is None:
        downstream_lower_bound = 0.01 * upstream_pressure

    result = minimize_scalar(
        lambda downstream_pressure: -mdot_hem(area, upstream_pressure, downstream_pressure),
        bounds=(downstream_lower_bound, downstream_upper_bound),
        method='bounded'
    )

    return result.x

# Dyer non-homogenous non-equilibrium mass flow rate
def mdot(area, upstream_pressure, chamber_pressure):
    k_dyer = np.sqrt((upstream_pressure - chamber_pressure) / (upstream_pressure - chamber_pressure)) # Just 1 for self-pressurized N2O
    spi_weight = k_dyer / (1 + k_dyer)
    hem_weight = 1 / (1 + k_dyer)
    
    # Check choked condition
    critical_chamber_pressure = mdot_hem_crit(area, upstream_pressure)    
    if chamber_pressure < critical_chamber_pressure:
        chamber_pressure = critical_chamber_pressure

    mass_flow = spi_weight * mdot_spi(area, upstream_pressure, chamber_pressure) + hem_weight * mdot_hem(area, upstream_pressure, chamber_pressure)
    return mass_flow