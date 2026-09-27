"""
Calculator to find fluid properties through piping

Iteratively solves for mass flow rate and pressure loss through various
elements of the liquid rocket test stand

* Authors: Arthur Gwozdz, Ethan Rosenfeld
* Date: 11/13/2025
"""

import numpy as np

from conversions import *
from objects import *

braided_roughness = 0.0000005
stainless_roughness = 0.000002

n2o_injector = Injector("N2O Injector", 0.0015, 12, 0.7)
n2o_line = [
    Pipe("N2O P1", 0.0000005, inches_to_meters(0.3125), inches_to_meters(16)),
    Pipe("N2O P2", 0.000002, inches_to_meters(0.43), inches_to_meters(5)),
    Fitting("N2O F1", 1.12, inches_to_meters(0.43)),
    Pipe("N2O P3", 0.000002, inches_to_meters(0.43), inches_to_meters(11)),
    Flowmeter("N2O FM", inches_to_meters(0.43), inches_to_meters(0.24), 0.7),
    Pipe("N2O P4", 0.000002, inches_to_meters(0.43), inches_to_meters(6)),
    Pipe("N2O P5", 0.0000005, inches_to_meters(0.43), inches_to_meters(14)),
    Fitting("N2O F2", 0.56, inches_to_meters(0.43)),
    Pipe("N2O P6", 0.000002, inches_to_meters(0.43), inches_to_meters(12)),
    Pipe("N2O P7", 0.0000005, inches_to_meters(0.3125), inches_to_meters(36))
]
N2OLine = FluidLine(
    "N2O",
    Propellant("Water", 998, 0.5, 0.00089), 
    # Propellant("Nitrous Oxide", 1220, 0.5, 0.003237),
    n2o_line, 
    n2o_injector, 
    psi_to_pascal(575), 
    psi_to_pascal(14.7)
)
# N2OLine.solve(1000)

kerosene_injector = Injector("Kerosene Injector", 0.001, 12, 0.46887)
kerosene_line = [
    Pipe("Kerosene P1", braided_roughness, inches_to_meters(0.1875), inches_to_meters(30)),
    Pipe("Kerosene P2", stainless_roughness, inches_to_meters(0.18), inches_to_meters(24)),
    Fitting("Kerosene F1", 0.7, inches_to_meters(0.18)),
    Flowmeter("Kerosene FM", inches_to_meters(0.18), inches_to_meters(0.12), 0.67396),
    Fitting("Kerosene F2", 0.7, inches_to_meters(0.18)),
    Pipe("Kerosene P3", stainless_roughness, inches_to_meters(0.18), inches_to_meters(6)),
    Pipe("Kerosene P4", braided_roughness, inches_to_meters(0.1875), inches_to_meters(14)),
    Fitting("Kerosene F3", 0.7, inches_to_meters(0.18)),
    Pipe("Kerosene P5", stainless_roughness, inches_to_meters(0.18), inches_to_meters(6)),
    Pipe("Kerosene P6", braided_roughness, inches_to_meters(0.1875), inches_to_meters(24)),
    Fitting("Kerosene F4", 0.7, inches_to_meters(0.18))
]
KeroseneLine = FluidLine(
    "Kerosene",
    Propellant("Water", 998, 0.26, 0.00089),
    # Propellant("Kerosene", 807, 0.26, 0.00164),
    kerosene_line, 
    kerosene_injector, 
    psi_to_pascal(700), 
    psi_to_pascal(0)
)
KeroseneLine.solve(1000)

# Kerosene 700 psi tank pressure 3/1/2026
# print("Flowmeter upstream: 543 psi")
# print("Flowmeter delta: 113 psi")
# print("Flowmeter mdot: 0.33 kg/s")
# print("Injector delta: 195 psi")
# print("Injector mdot: 0.33 kg/s")

# Nitrous 660 psi tank pressure 3/1/2026
# print("Flowmeter upstream: 455 psi")
# print("Flowmeter delta: 48 psi")
# print("Flowmeter mdot: 0.99 kg/s")
# print("Injector delta: 320 psi")
# print("Injector mdot: 0.99 kg/s")

# Nitrous tank
# Diameter: 5.747000 in
# Height: 32.484000 in

# Kerosene tank
# Diameter: 5.747000 in
# Height: 12.484000 in

def tank_mass(prop, diameter, height):
    volume = np.pi * np.pow(diameter/2, 2) * height
    mass = volume * prop.density
    return mass

kerosene_mass = tank_mass(KeroseneLine.propellant, inches_to_meters(5.747), inches_to_meters(12.484))
kerosene_time = 16.03


n2o_mass = tank_mass(N2OLine.propellant, inches_to_meters(5.747), inches_to_meters(32.484))
n2o_time = 13.95

# print(kerosene_mass / kerosene_time, n2o_mass / n2o_time)