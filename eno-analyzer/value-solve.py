"""
Solve for roughness parameters based on experimental data

This function is currently a work in progress

We are solving for 6 parameters that will reproduce the correct mass flow
rates and pressure differences across the flowmeter and injector. Thus we
have 4 metrics we can calculate from the 6 parameters. We also have 4
ideal metrics that we are iterating towards.

The solving algorithm will work as follows

There are three key functions
    1. Error
        Takes in the 6 parameters and returns the difference in our 4 metrics
        from the ideal
    2. Partial derivative
        Calculates the partial derivative of the 4 metrics with respect to 1
        of the 6 parameters
    3. Iteration
        Uses a local linear approximation to decrease the error

The algorithm steps are
    1. Make an initial guess of the 6 parameters
    2. Calculate the error in each metric
    3. Find partial derivatives
        a. We can organize these as one 4-coordinate vector for each parameter
        b. Let that matrix be A
    4. Find a vector x from a linear combiation of the partial derivatives
        a. The vector x should have the property that e - Ax is as close to 0 as
        possible, where e is the error vector
        b. Note that we have 6 parameters but only 4 metrics. Assuming each
        partial derivative is linearly independent, we should use the free variables
        to bound our variables to a reasonable range (between 0 and 1)
    5. Add this x vector to our 6 parameters
    6. Iterate starting from step 2 until e - x is reasonably small

* Authors: Ethan Rosenfeld, Michael Danley
* Date: 11/14/2025
"""

from conversions import *
from objects import *
import numpy as np

braided_roughness = 0.0000005
stainless_roughness = 0.000002

def get_kerosene_line(inj_cd, fm_cd, fmpt_cd, cv_cd):
    """FluidLine object from 2 initial parameters"""
    kerosene_line = [
        Pipe("Tank to Panel", braided_roughness, inches_to_meters(0.1875), inches_to_meters(16)),
        Pipe("Panel to FM", stainless_roughness, inches_to_meters(0.18), inches_to_meters(24)),
        Fitting("FM Upstream PT", fmpt_cd, inches_to_meters(0.18)),
        Flowmeter("Kerosene FM", inches_to_meters(0.18), inches_to_meters(0.12), fm_cd),
        Fitting("FM Downstream PT", fmpt_cd, inches_to_meters(0.18)),
        Pipe("FM to Panel", stainless_roughness, inches_to_meters(0.18), inches_to_meters(16)),
        Pipe("Panel to Stand", braided_roughness, inches_to_meters(0.1875), inches_to_meters(16)),
        Fitting("Check Valve", cv_cd, inches_to_meters(0.18)),
        Pipe("Servo to Stand", stainless_roughness, inches_to_meters(0.18), inches_to_meters(22)),
        Pipe("Stand to Regen", braided_roughness, inches_to_meters(0.1875), inches_to_meters(16)),
        Pipe("Regen Channels", stainless_roughness, inches_to_meters(0.3768), inches_to_meters(6))
    ]
    KeroseneLine = FluidLine(
        "Kerosene",
        Propellant("Water", 998, 0.5, 0.00089),
        kerosene_line,
        Injector("Kerosene Injector", 0.001, 12, inj_cd), 
        psi_to_pascal(700), 
        psi_to_pascal(14.7)
    )
    return KeroseneLine

def get_nitrous_line(inj_cd, fm_cd):
    """FluidLine object from 2 initial parameters"""
    nitrous_line = [
        Pipe("N2O P1", braided_roughness, inches_to_meters(0.3125), inches_to_meters(16)),
        Pipe("N2O P2", stainless_roughness, inches_to_meters(0.43), inches_to_meters(5)),
        Fitting("N2O F1", 0.7, inches_to_meters(0.43)),
        Pipe("N2O P3", stainless_roughness, inches_to_meters(0.43), inches_to_meters(11)),
        Flowmeter("N2O FM", inches_to_meters(0.43), inches_to_meters(0.24), fm_cd),
        Pipe("N2O P4", stainless_roughness, inches_to_meters(0.43), inches_to_meters(6)),
        Pipe("N2O P5", braided_roughness, inches_to_meters(0.43), inches_to_meters(14)),
        Fitting("N2O F2", 0.7, inches_to_meters(0.43)),
        Pipe("N2O P6", stainless_roughness, inches_to_meters(0.43), inches_to_meters(12)),
        Pipe("N2O P7", braided_roughness, inches_to_meters(0.3125), inches_to_meters(36))
    ]
    NitrousLine = FluidLine(
        "N2O",
        Propellant("Water", 998, 0.5, 0.00089), 
        nitrous_line, 
        Injector("N2O Injector", 0.0015, 12, inj_cd), 
        psi_to_pascal(523), 
        psi_to_pascal(14.7)
)
    return NitrousLine

def partial(get_line, index, params):
    """Partial derivative of the metrics with respect to a parameter given by index"""
    delta = 1e-6
    fm_mdot_initial, fm_delta_initial, inj_mdot_initial, inj_delta_initial, fm_upstream_initial = get_line(*params).solve(1000, False, 1e-10)
    params_final = [p for p in parameters]
    params_final[index] += delta
    fm_mdot_final, fm_delta_final, inj_mdot_final, inj_delta_final, fm_upstream_final = get_line(*params_final).solve(1000, False, 1e-10)
    fm_mdot_derivative = (fm_mdot_final - fm_mdot_initial) / delta
    fm_delta_derivative = (fm_delta_final - fm_delta_initial) / delta
    inj_mdot_derivative = (inj_mdot_final - inj_mdot_initial) / delta
    inj_delta_derivative = (inj_delta_final - inj_delta_initial) / delta
    fm_upstream_derivative = (fm_upstream_final - fm_upstream_initial) / delta
    return fm_mdot_derivative, fm_delta_derivative, inj_mdot_derivative, inj_delta_derivative, fm_upstream_derivative

def error(get_line, params, display=False):
    """Error generated by a set of parameters"""
    fm_mdot, fm_delta, inj_mdot, inj_delta, fm_upstream = get_line(*params).solve(1000, False)
    fm_mdot_error = flowmeter_mass_flow_rate - fm_mdot
    fm_delta_error = flowmeter_delta - fm_delta
    inj_mdot_error = injector_mass_flow_rate - inj_mdot
    inj_delta_error = injector_delta - inj_delta
    fm_upstream_error = flowmeter_upstream - fm_upstream
    if display:
        print(f"Flowmeter mdot (error): {fm_mdot:.3f} + ({fm_mdot_error:.3f}) kg/s")
        print(f"Flowmeter delta (error): {pascal_to_psi(fm_delta):.3f} + ({pascal_to_psi(fm_delta_error):.3f}) psi")
        print(f"Injector mdot (error): {inj_mdot:.3f} + ({inj_mdot_error:.3f}) kg/s")
        print(f"Injector delta (error): {pascal_to_psi(inj_delta):.3f} + ({pascal_to_psi(inj_delta_error):.3f}) psi")
        print(f"Flowmeter upstream (error): {pascal_to_psi(fm_upstream):.3f} + ({pascal_to_psi(fm_upstream_error):.3f}) psi")
    return fm_mdot_error, fm_delta_error, inj_mdot_error, inj_delta_error, fm_upstream_error

def iterate(get_line, params, epsilon=1e-1):
    """Local linear approximation to decrease error"""
    # Get partial derivative with respect to each parameter
    partials = []
    for i in range(0, len(params)):
        partials.append(np.array(partial(get_line, i, params)))
    # Jacobian Matrix (columns of partial derivatives)
    J = np.column_stack(partials)
    # Error vector
    e = np.array(error(get_line, params))
    # Scale values to the same range
    scales = np.maximum(np.abs(e), 1e-12)
    J_scaled = J / scales[:, None]
    e_scaled = e / scales
    # Solve the linear system Ax = e (least squares)
    x_hat, residuals, rank, s = np.linalg.lstsq(J_scaled, e_scaled, rcond=None)
    # Scale down (low confidence in Jacobian approximation)
    x_hat *= epsilon
    for i in range(0, len(params)):
        params[i] += x_hat[i]
        if params[i] < 0:
            params[i] = 0
    

def squared_error(error_terms):
    """Sum the squares of the error"""
    total_squared_error = 0
    for term in error_terms:
        total_squared_error += pow(term, 2)
    return total_squared_error

def grad_descent(get_line, parameters, n, display=False):
    # print("----- Running Gradient Descent -----")
    if display:
        print(f"Initial parameters: {parameters}")
        initial_error_terms = error(get_line, parameters, True)
        initial_error = squared_error(initial_error_terms)
        print(f"Initial error: {initial_error}")
    for i in range(0, n):
        if display and (100*i/n) % 5 == 0:
            print(f'\r Progress: {100*i/n}%', end="\r", flush=True)
        iterate(get_line, parameters)
    if display:
        print("----- Results -----")
        print(f"Final parameters: {parameters}")
        final_error_terms = error(get_line, parameters, True)
    final_error = squared_error(final_error_terms)
    if display:
        print(f"Final error: {final_error}")
        print(f"Error change: {100 * (final_error - initial_error) / initial_error:.10f}%")
    return final_error

# Experimental values
# Kerosene
flowmeter_upstream = psi_to_pascal(543)
flowmeter_delta = psi_to_pascal(113)
flowmeter_mass_flow_rate = 0.33
injector_delta = psi_to_pascal(195)
injector_mass_flow_rate = 0.33

# Nitrous
# flowmeter_mass_flow_rate = 1
# flowmeter_delta = psi_to_pascal(37)
# injector_mass_flow_rate = 1
# injector_delta = psi_to_pascal(246)
# flowmeter_upstream = psi_to_pascal(372)

# Initial guess
# Kerosene
parameters = [
    0.7247694833872729, # Injector Cd
    0.7428607385478915, # Flowmeter Cd
    0.7,
    0.7
]
# Nitrous
# parameters = [
#     0.7183760151020631, # Injector Cd
#     1.4428221117425961, # Flowmeter Cd
# ]

grad_descent(get_kerosene_line, parameters, 250, True)