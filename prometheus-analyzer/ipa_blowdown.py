import numpy as np
from scipy.optimize import newton

PI = np.pi
DENSITY = 786 # kg/m3
ROUGHNESS = 0.000045
DYNAMIC_VISCOSITY = 0.002381

"""
Iteratively solve for IPA mass flow rate by matching injector pressure drop
"""
def get_ipa_mdot(HEIGHT, WIDTH, LENGTH, N, ipa_pressure, chamber_pressure):
    area = (WIDTH * HEIGHT) + (PI * (WIDTH / 2)**2) # Rectangular section + 2 semicircles
    perimeter = (2 * HEIGHT) + (PI * WIDTH) # Rectangular section + 2 semicircles
    hydraulic_diameter = 4 * area / perimeter

    # Darcy-Weisbach frictional loss
    def dw_friction_loss(mdot):
        velocity = (mdot / N) / (DENSITY * area)
        reynolds = DENSITY * velocity * hydraulic_diameter / DYNAMIC_VISCOSITY
        friction_factor = 0.25 * np.log10((ROUGHNESS / (3.7 * hydraulic_diameter)) + (5.74 / reynolds**0.9))**-2
        loss_per_meter = friction_factor * DENSITY * velocity**2 / (2 * hydraulic_diameter)
        return loss_per_meter * LENGTH

    # Residual for Newton's Method
    def loss_residual(mdot):
        required_loss = ipa_pressure - chamber_pressure
        return dw_friction_loss(mdot) - required_loss

    mdot_solution = newton(loss_residual, x0=0.480) # Iteratively solve for mdot

    return mdot_solution
  