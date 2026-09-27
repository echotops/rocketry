"""
Conversions between ISO and Imperial units

Contains a functions to convert lengths between inches and meters and pressures
between Pascals and PSI.

* Authors: Arthur Gwozdz, Ethan Rosenfeld
* Date: 11/11/2025
"""

def inches_to_meters(length):
    """Convert lengths from inches to meters"""
    return length / 39.3701
    
def meters_to_inches(length):
    """Convert lengths from meters to inches"""
    return length * 39.3701
    
def psi_to_pascal(pressure):
    """Convert pressures from PSI to Pascals"""
    return pressure * 6894.76
    
def pascal_to_psi(pressure):
    """Convert pressures from Pascals to PSI"""
    return pressure / 6894.76