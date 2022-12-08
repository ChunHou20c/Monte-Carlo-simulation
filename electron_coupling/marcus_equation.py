"""This module define Marcus's equation of electron transfer"""
from math import pi, exp

def transfer_rate(electron_coupling:float, reorganization_energy:float, temperature:float, Eij:float = 0):
    """
    This function define the calculation of marcus transfer rate
    all the value are in ev
    temperature in Kelvin
    """

    reduced_plack_constant = 6.582119569e-16
    boltzmann_constant = 8.617333262e-5
    
    term1 = 2*pi/reduced_plack_constant

    term2 = electron_coupling**2

    term3 = (1/(4*pi*boltzmann_constant*temperature*reorganization_energy))**0.5

    inside_exponent = -((Eij + reorganization_energy)**2)/(4*boltzmann_constant*temperature*reorganization_energy)
    
    term4 = exp(inside_exponent)

    return term1*term2*term3*term4
