"""
Objects used in simulating fluid flow

Contains a Propellant class to store propellant properties, Pipe, Fitting,
Flowmeter, and Injector classes to measure mass flow rate and pressure
differences, and a FluidLine class to contain elements and iteratively solve to
maintain continuity.

* Authors: Arthur Gwozdz, Ethan Rosenfeld
* Date: 11/11/2025
"""

from conversions import *
import math

class Propellant:
    """
    Propellant to flow through elements
    
    Functions:
        velocity (float -> float): Propellant flow velocity
        reynoldsNumber (float -> float): Flow Reynolds number
    """
    def __init__(self, name, density, mass_flow_rate, dynamic_viscosity):
        self.name = name
        self.density = density
        self.mass_flow_rate = mass_flow_rate
        self.dynamic_viscosity = dynamic_viscosity
        
    def velocity(self, diameter):
        """Calculate velocity from pipe diameter based on continuity"""
        area = math.pi * pow(diameter / 2, 2)
        return self.mass_flow_rate / (self.density * area)
    
    def reynoldsNumber(self, diameter):
        """Calculate Reynolds number from pipe diameter"""
        velocity = self.velocity(diameter)
        return self.density * velocity * diameter / self.dynamic_viscosity

class Pipe:
    """
    Pipe through which propellant can flow
    
    Functions:
        pressureLoss (Propellant -> float): Pressure lost by the propellant
        frictionFactor (Propellant, int -> float): The Dacy-Weisbach friction factor
    """
    def __init__(self, name, roughness, diameter, length):
        self.name = name
        self.roughness = roughness
        self.diameter = diameter
        self.length = length

    def pressureLoss(self, propellant):
        """Pressure lost in the pipe based on the propellant"""
        velocity = propellant.velocity(self.diameter)
        return self.frictionFactor(propellant, 100) * self.length * propellant.density * pow(velocity, 2) / (2 * self.diameter)

    def frictionFactor(self, propellant, n):
        """Iteratively solve for the Darcy-Weishbach friction factor"""
        reynolds_number = propellant.reynoldsNumber(self.diameter)
         # Check for laminar approimation
        if reynolds_number <= 2300:
            return 64 / reynolds_number 
        # Initial guess via Swamee-Jain
        friction_factor = 0.25 / pow(math.log10((self.roughness / (3.7 * self.diameter)) + (5.74 / pow(reynolds_number, 0.9))), 2)
        # Iteration using Colebrook-White
        for _ in range(0, n):
            friction_factor = pow(-2 * math.log10((self.roughness / (3.7 * self.diameter)) + (2.51 / (reynolds_number * math.sqrt(friction_factor)))), -2)
        return friction_factor

class Fitting:
    """
    Fitting through which propellant can flow
    
    Functions:
        pressureLoss (Propellant -> float): Pressure lost by the propellant
    """
    def __init__(self, name, k_value, diameter):
        self.name = name
        self.k_value = k_value
        self.diameter = diameter

    def pressureLoss(self, propellant):
        """Pressure lost in the fitting based on the propellant"""
        density = propellant.density
        velocity = propellant.velocity(self.diameter)
        return self.k_value * density * pow(velocity, 2) / 2

class Flowmeter:
    """
    Flowmeter through which propellant can flow
    
    Functions:
        pressureLoss (Propellant -> float): Pressure lost by the propellant
        pressureLossIndicated (Propellant -> float): Indicated pressure lost by the propellant
        massFlow (Propellant, float -> float): Propellant mass flow
        massFlowIndicated (Propellant, float -> float): Propellant mass flow from indicated pressure loss
    """
    def __init__(self, name, diameter, orifice_diameter, discharge_coef):
        self.name = name
        self.diameter = diameter
        self.orifice_diameter = orifice_diameter
        self.discharge_coef = discharge_coef

    def pressureLossIndicated(self, propellant):
        """Indicated pressure lost in the flowmeter based on the propellant"""
        mdot = propellant.mass_flow_rate
        rho = propellant.density
        beta = self.orifice_diameter / self.diameter
        expansibility = 1
        # Mass flow and pressure relation across an orifice
        pressure_delta = (8 * pow(mdot, 2) * (1 - pow(beta, 4))) / (rho * pow(self.discharge_coef * expansibility * math.pi, 2) * pow(self.orifice_diameter, 4))
        return pressure_delta

    def pressureLoss(self, propellant):
        """Pressure lost in the fitting based on the propellant"""
        beta = self.orifice_diameter / self.diameter
        # Correction for overall pressure difference
        return self.pressureLossIndicated(propellant) * (1 - pow(beta, 1.9))

    def massFlow(self, propellant, pressure_delta):
        """Propellant mass flow based on the propellant and pressure loss"""
        rho = propellant.density
        beta = self.orifice_diameter / self.diameter
        expansibility = 1
        # Mass flow and pressure relation across an orifice
        return (self.discharge_coef * expansibility * math.pi * pow(self.orifice_diameter, 2) * math.sqrt(2 * rho * pressure_delta)) / (4 * math.sqrt(1 - pow(beta, 4)))

    def massFlowIndicated(self, propellant, pressure_delta_indicated):
        """Propellant mass flow based on the propellant and indicated pressure loss"""
        beta = self.orifice_diameter / self.diameter
        # Correction for mesured pressure difference
        return self.massFlow(propellant, pressure_delta_indicated * (1 - pow(beta, 1.9)))

class Injector:
    """
    Injector through which propellant can flow
    
    Functions:
        massFlow (Propellant, float -> float): Propellant mass flow
    """
    def __init__(self, name, diameter, count, discharge_coef):
        self.name = name
        self.diameter = diameter
        self.count = count
        self.discharge_coef = discharge_coef

    def massFlow(self, propellant, pressure_delta):
        """Propellant mass flow based on the propellant and pressure loss"""
        area = math.pi * pow(self.diameter / 2, 2)
        # Mass flow and pressure relation through an injector
        return self.discharge_coef * area * self.count * math.sqrt(2 * pressure_delta * propellant.density)

class FluidLine:
    """
    FluidLine to contain elements and iteratively solve for continuity
    
    Functions:
        solve (int -> None): Iteratively solve to match mass
    """
    def __init__(self, name, propellant, line, injector, tank_pressure, chamber_pressure):
        self.name = name
        self.propellant = propellant
        self.line = line
        self.injector = injector
        self.tank_pressure = tank_pressure
        self.chamber_pressure = chamber_pressure

    def solve(self, iterations, display=True, threshold=1e-6):
        """Iteratively solve for continuity and print results"""
        if display:
            print("----- Fluid Line Solve (" + self.name + ") -----")
            print("Parameters")
            print(f"\t Tank pressure: {pascal_to_psi(self.tank_pressure)} psi")
            print(f"\t Chamber pressure: {pascal_to_psi(self.chamber_pressure)} psi")
            print(f"\t Iterations: {iterations}")
        error = float('inf')
        # Iterate until reaching threshold or hitting specified iteration count
        for i in range(0, iterations):
            converged = False # Solution convergence
            # Check if solution is converged
            if abs(error) <= threshold:
                converged = True
                if display:
                    print(f"Solution converged in {i} iterations")
                    print("Element Losses")
            flowmeter = None
            total_loss = 0
            flowmeter_upstream = self.tank_pressure
            # Loop through elements and sum pressure losses
            for element in self.line:
                # Keep track of flowmeters
                if isinstance(element, Flowmeter):
                    flowmeter = element
                loss = element.pressureLoss(self.propellant)
                total_loss += loss
                if flowmeter == None:
                    flowmeter_upstream -= loss
                # Print individual element losses on last run
                if i == iterations - 1 or converged:
                    if display:
                        print(f"\t {element.name}: {pascal_to_psi(loss):.3f} psi")
            # Injector pressure loss
            injector_delta = self.tank_pressure - total_loss - self.chamber_pressure
            if injector_delta < 0: injector_delta = 0 
            # Injector mass flow rate
            injector_mdot = self.injector.massFlow(self.propellant, injector_delta)
            # Flowmeter pressure loss
            flowmeter_delta = flowmeter.pressureLoss(self.propellant)
            if flowmeter_delta < 0: flowmeter_delta = 0
            # Flowmeter mass flow rate
            flowmeter_mdot = flowmeter.massFlow(self.propellant, flowmeter_delta)
            # Compare injector and flowmeter mass flow rates
            error = injector_mdot - flowmeter_mdot
            # Correct based on the difference between them
            self.propellant.mass_flow_rate += error/10
            # Stop looping if the solution has converged
            if converged:
                break
        if display:
            print("Results")
            print(f"\t Continuity residual: {error:.3e} kg/s")
            print(f"\t Head loss: {pascal_to_psi(total_loss):.3f} psi")
            print(f"\t Flowmeter upstream: {pascal_to_psi(flowmeter_upstream):.3f} psi")
            print(f"\t Flowmeter delta: {pascal_to_psi(flowmeter_delta):.3f} psi")
            print(f"\t Flowmeter mdot: {flowmeter_mdot:.3f} kg/s")
            print(f"\t Injector delta: {pascal_to_psi(injector_delta):.3f} psi")
            print(f"\t Injector mdot: {injector_mdot:.3f} kg/s")
        return flowmeter_mdot, flowmeter_delta, injector_mdot, injector_delta, flowmeter_upstream
