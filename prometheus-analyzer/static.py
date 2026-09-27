from ipa_blowdown import get_ipa_mdot
from n2o_blowdown import get_n2o_mdot
from scipy.optimize import newton, brentq
import numpy as np

IN = 0.0254
PSI = 6894.76
R_UNIV = 8314.46
PI = np.pi
g = 9.8

# ----- NOX Injector -----
NOX_ID = 0.82 * IN
NOX_OD = NOX_ID + 2 * (0.027 * IN)
NOX_LENGTH = 1 * IN

# ----- IPA Injector -----
IPA_HEIGHT = 0.3 * IN
IPA_WIDTH = 0.006 * IN
IPA_LENGTH = 0.2 * IN
IPA_N = 24

# ----- NOX Tank -----
nox_mass = 9.46
nox_volume = 14.40
nox_pressure = 900 * PSI

# ----- IPA Tank -----
ipa_mass = 2.85
ipa_pressure = 460 * PSI
nox_pressure = 900 * PSI

chamber_pressure = 400 * PSI

AMBIENT_PRESSURE = 101325
FLAME_TEMP = 3000
GAMMA = 1.3
MOLECULAR_MASS = 26
C_STAR = 1 * np.sqrt(GAMMA * (R_UNIV / MOLECULAR_MASS) * FLAME_TEMP) / (GAMMA * np.sqrt(np.pow(2 / (GAMMA + 1), (GAMMA + 1) / (GAMMA - 1))))
THROAT_AREA = PI * (1.298 * IN / 2)**2
EXPANSION_RATIO = 4.2

def total_mdot(ipa_pressure, nox_pressure, chamber_pressure):
    ipa_mdot = get_ipa_mdot(IPA_HEIGHT, IPA_WIDTH, IPA_LENGTH, IPA_N, ipa_pressure, chamber_pressure)
    nox_mdot = get_n2o_mdot(NOX_ID, NOX_OD, nox_pressure, chamber_pressure)
    return ipa_mdot + nox_mdot

def chamber_pressure_residual(ipa_pressure, nox_pressure, chamber_pressure):
    mdot = total_mdot(ipa_pressure, nox_pressure, chamber_pressure)
    predicted_chamber_pressure = C_STAR * mdot / THROAT_AREA
    return predicted_chamber_pressure - chamber_pressure

actual_chamber_pressure = newton(
    lambda chamber_pressure : chamber_pressure_residual(ipa_pressure, nox_pressure, chamber_pressure),
    400 * PSI
)

ipa_mdot = get_ipa_mdot(IPA_HEIGHT, IPA_WIDTH, IPA_LENGTH, IPA_N, ipa_pressure, actual_chamber_pressure)
nox_mdot = get_n2o_mdot(NOX_ID, NOX_OD, nox_pressure, actual_chamber_pressure)
total_mdot = ipa_mdot + nox_mdot

print(f"Equilibrium Chamber Pressure (PSI): {actual_chamber_pressure / PSI}")
print(f"IPA Mass Flow (kg/s), NOX Mass Flow (kg/s): {ipa_mdot}, {nox_mdot}")
print(f"Total Mass Flow (kg/s): {total_mdot}")

def area_mach_residual(mach, area_ratio):
    term = (2 / (GAMMA + 1)) * (1 + (GAMMA - 1) / 2 * mach**2)
    return (1 / mach) * term**((GAMMA + 1)/(2 * (GAMMA - 1))) - area_ratio

def get_exit_mach(area_ratio):
    return brentq(area_mach_residual, 1.0001, 20, args=(area_ratio))  # supersonic root

def get_exit_pressure(chamber_pressure, exit_mach):
    return chamber_pressure / (1 + (GAMMA - 1) / 2 * exit_mach**2)**(GAMMA / (GAMMA - 1))

def get_exit_temperature(exit_mach):
    return FLAME_TEMP / (1 + (GAMMA - 1) / 2 * exit_mach**2)

exit_mach = get_exit_mach(EXPANSION_RATIO)
exit_pressure = get_exit_pressure(actual_chamber_pressure, exit_mach)
exit_temperature = get_exit_temperature(exit_mach)

exit_velocity = exit_mach * np.sqrt(GAMMA * (R_UNIV / MOLECULAR_MASS) * exit_temperature)
print(f"Exit velocity (m/s): {exit_velocity}")
print(f"Exit pressure (PSI): {exit_pressure / PSI}")
thrust = total_mdot * exit_velocity + (exit_pressure - AMBIENT_PRESSURE) * (EXPANSION_RATIO * THROAT_AREA)
specific_impulse = thrust / (total_mdot * g)
print(f"Thrust (N), ISP (s): {thrust}, {specific_impulse}")