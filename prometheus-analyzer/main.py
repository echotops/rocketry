from ipa_blowdown import get_ipa_mdot
from n2o_blowdown import get_n2o_mdot, get_n2o_injector_density
from scipy.optimize import brentq
import CoolProp.CoolProp as CP
from CoolProp import AbstractState
import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj

# ----- Physics Constants -----
IN = 0.0254
PSI = 6894.76
R_UNIV = 8314.46
PI = np.pi
g = 9.8
CEA = CEA_Obj(oxName='N2O', fuelName='Isopropanol')

# ----- NOX Injector -----
N2O_ID = 0.82 * IN
N2O_OD = N2O_ID + 2 * (0.0227 * IN)
N2O_LENGTH = 1 * IN

# ----- IPA Injector -----
IPA_HEIGHT = 0.3 * IN
IPA_WIDTH = 0.006 * IN
IPA_LENGTH = 0.2 * IN
IPA_N = 24

# ----- NOX Tank -----
n2o_mass = 9.46
n2o_volume = 14.40
n2o_pressure = 900 * PSI

# ----- IPA Tank -----
ipa_mass = 2.85
ipa_pressure = 460 * PSI

# ----- Combustion Properties -----
AMBIENT_PRESSURE = 101325
THROAT_AREA = PI * (1.298 * IN / 2)**2
EXPANSION_RATIO = 4.32

# ----- Initial NOX Tank State -----
rho_liq_init = CP.PropsSI('D', 'P', n2o_pressure, 'Q', 0, 'N2O')
rho_vap_init = CP.PropsSI('D', 'P', n2o_pressure, 'Q', 1, 'N2O')
u_liq_init = CP.PropsSI('U', 'P', n2o_pressure, 'Q', 0, 'N2O')
u_vap_init = CP.PropsSI('U', 'P', n2o_pressure, 'Q', 1, 'N2O')

rho_init = 1e3 * n2o_mass / n2o_volume
quality_init = (rho_init - rho_vap_init) / (rho_liq_init - rho_vap_init)
launch_liq_n2o = quality_init * n2o_mass
launch_vap_n2o = (1 - quality_init) * n2o_mass

m_ipa = ipa_mass
m_n2o = launch_liq_n2o + launch_vap_n2o
U_n2o = launch_liq_n2o * u_liq_init + launch_vap_n2o * u_vap_init

def get_combustion_properties(chamber_pressure, mixture_ratio):
    molecular_weight, gamma = CEA.get_Chamber_MolWt_gamma(Pc=chamber_pressure/PSI, MR=mixture_ratio, eps=EXPANSION_RATIO)
    chamber_temp = CEA.get_Tcomb(Pc=chamber_pressure/PSI, MR=mixture_ratio) * 5/9  # CEA returns Rankine
    c_star = np.sqrt(gamma * (R_UNIV / molecular_weight) * chamber_temp) / (gamma * np.sqrt(np.pow(2 / (gamma + 1), (gamma + 1) / (gamma - 1))))
    return c_star, gamma, chamber_temp, molecular_weight

def get_total_mdot(ipa_pressure, n2o_pressure, chamber_pressure):
    ipa_mdot = get_ipa_mdot(IPA_HEIGHT, IPA_WIDTH, IPA_LENGTH, IPA_N, ipa_pressure, chamber_pressure)
    n2o_mdot = get_n2o_mdot(N2O_ID, N2O_OD, n2o_pressure, chamber_pressure)
    return ipa_mdot + n2o_mdot

def chamber_pressure_residual(ipa_pressure, n2o_pressure, chamber_pressure):
    ipa_mdot = get_ipa_mdot(IPA_HEIGHT, IPA_WIDTH, IPA_LENGTH, IPA_N, ipa_pressure, chamber_pressure)
    n2o_mdot = get_n2o_mdot(N2O_ID, N2O_OD, n2o_pressure, chamber_pressure)
    mixture_ratio = n2o_mdot / ipa_mdot
    c_star, _, _, _ = get_combustion_properties(chamber_pressure, mixture_ratio)
    predicted_chamber_pressure = c_star * (ipa_mdot + n2o_mdot) / THROAT_AREA
    return predicted_chamber_pressure - chamber_pressure

def area_mach_residual(mach, area_ratio, gamma):
    term = (2 / (gamma + 1)) * (1 + (gamma - 1) / 2 * mach**2)
    return (1 / mach) * term**((gamma + 1)/(2 * (gamma - 1))) - area_ratio

def get_exit_mach(area_ratio, gamma):
    return brentq(area_mach_residual, 1.0001, 20, args=(area_ratio, gamma))  # supersonic root

def get_exit_pressure(chamber_pressure, exit_mach, gamma):
    return chamber_pressure / (1 + (gamma - 1) / 2 * exit_mach**2)**(gamma / (gamma - 1))

def get_exit_temp(chamber_temp, exit_mach, gamma):
    return chamber_temp / (1 + (gamma - 1) / 2 * exit_mach**2)

t = 0
dt = 0.05
performance_data = []
n2o_data = []
ipa_data = []
AS = AbstractState("HEOS", "N2O")

while m_n2o > 0 and m_ipa > 0:
    rho_avg = 1e3 * m_n2o / n2o_volume
    U_avg = U_n2o / m_n2o
    AS.update(CP.DmassUmass_INPUTS, rho_avg, U_avg)
    T = AS.T()
    P = AS.p()
    Q = AS.Q()

    if Q >= 0.999 or Q < 0:
        break
    
    chamber_pressure = brentq(
        lambda chamber_pressure: chamber_pressure_residual(ipa_pressure, P, chamber_pressure),
        50 * PSI,
        0.99 * min(ipa_pressure, P),
    )

    ipa_mdot = get_ipa_mdot(IPA_HEIGHT, IPA_WIDTH, IPA_LENGTH, IPA_N, ipa_pressure, chamber_pressure)
    ipa_mdot = 0.480 # Just assume the pump's gonna be insanely accurate
    ipa_area = (IPA_WIDTH * IPA_HEIGHT) + (PI * (IPA_WIDTH / 2)**2)
    ipa_velocity = (ipa_mdot / IPA_N) / (786 * ipa_area)

    n2o_mdot = get_n2o_mdot(N2O_ID, N2O_OD, P, chamber_pressure)
    n2o_rho_injector = get_n2o_injector_density(N2O_ID, N2O_OD, n2o_pressure, chamber_pressure)
    n2o_area = PI * ((N2O_OD/2)**2 - (N2O_ID/2)**2)
    n2o_velocity = n2o_mdot / (n2o_rho_injector * n2o_area)

    total_mdot = ipa_mdot + n2o_mdot
    mixture_ratio = n2o_mdot / ipa_mdot
    c_star, gamma, chamber_temp, molecular_weight = get_combustion_properties(chamber_pressure, mixture_ratio)

    h_liq_out = CP.PropsSI('H', 'P', P, 'Q', 0, 'N2O') if Q < 1 else CP.PropsSI('H', 'T', T, 'Q', 1, 'N2O')
    
    exit_mach = get_exit_mach(EXPANSION_RATIO, gamma)
    exit_pressure = get_exit_pressure(chamber_pressure, exit_mach, gamma)
    exit_temp = get_exit_temp(chamber_temp, exit_mach, gamma)
    exit_velocity = exit_mach * np.sqrt(gamma * (R_UNIV / molecular_weight) * exit_temp)

    thrust = total_mdot * exit_velocity + (exit_pressure - AMBIENT_PRESSURE) * (EXPANSION_RATIO * THROAT_AREA)
    specific_impulse = thrust / (total_mdot * g)

    dm_ipa = ipa_mdot * dt
    dm_n2o = n2o_mdot * dt
    dU_n2o = n2o_mdot * h_liq_out * dt

    m_ipa -= dm_ipa
    m_n2o -= dm_n2o
    U_n2o -= dU_n2o
    t += dt

    performance_data.append((t, chamber_pressure/PSI, exit_pressure/PSI, exit_velocity, thrust, mixture_ratio))
    ipa_data.append((m_ipa, ipa_mdot, ipa_velocity))
    n2o_data.append((m_n2o, n2o_mdot, n2o_velocity, P/PSI))

# ----- Results -----
times, chamber_pressures, exit_pressures, exit_velocities, thrusts, mixture_ratios = np.array(performance_data).T
ipa_masses, ipa_mdots, ipa_velocities = np.array(ipa_data).T
n2o_masses, n2o_mdots, n2o_velocities, n2o_pressures = np.array(n2o_data).T

print(f"Burn ended at t = {times[-1]:.2f} s")
print(f"IPA Mass flow: {ipa_mdots[0]:.3f} -> {ipa_mdots[-1]:.3f} kg/s")
print(f"NOX Mass flow: {n2o_mdots[0]:.3f} -> {n2o_mdots[-1]:.3f} kg/s")
print(f"Thrust: {thrusts[0]:.1f} -> {thrusts[-1]:.1f} N")

# ----- Plotting -----
COLORS = [
    '#D55E00', # vermillion
    '#E69F00',  # orange
    '#009E73', # bluish green
    '#56B4E9', # sky blue
    '#0072B2', # blue
    '#CC79A7', # reddish purple
    '#000000', # black
]

# IPA metrics
fig1 = [
    [ipa_masses, (0, 3), "IPA Mass (kg)"],
    [ipa_mdots, (0, 1), "IPA Mass Flow (kg/s)"],
    [ipa_velocities, (0, 50), "IPA Velocities (m/s)"],

]

# NOX metrics
fig2 = [
    [n2o_masses, (0, 10), "NOX Mass (kg)"],
    [n2o_mdots, (0, 2), "NOX Mass Flow (kg/s)"],
    [n2o_velocities, (0, 100), "NOX Velocities (m/s)"],
    [n2o_pressures, (0, 1000), "NOX Pressure (PSI)"],
]

# Performance metrics
fig3 = [
    [chamber_pressures, (0, 600), "Chamber Pressure (PSI)"],
    [exit_velocities, (0, 3000), "Exit Velocity (m/s)"],
    [ipa_mdots + n2o_mdots, (0, 3), "Total Mass Flow (kg/s)"],
    [thrusts, (0, 5000), "Thrust (N)"],
    [mixture_ratios, (0, 4), "O/F Ratio"],
]

for title, plots in [("IPA Metrics", fig1), ("NOX Metrics", fig2), ("Performance Metrics", fig3)]:
    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax1.set_title(title)
    num_plots = len(plots)
    ax1.set_xlabel('Time (s)')
    for i in range(len(plots)):
        ax = ax1
        if i != 0:
            ax = ax1.twinx()
            ax.spines['right'].set_position(('outward', (i - 1) * 60))
        plot = plots[i]
        ax.set_ylabel(plot[2], color=COLORS[i])
        ax.plot(times, plot[0], color=COLORS[i])
        ax.tick_params(axis='y', labelcolor=COLORS[i])
        ax.set_ylim(*plot[1])
    fig.tight_layout()

plt.show()
