import numpy as np

def switching_function(distance, threshold, power_num=8, power_den=12):
    """
    Calculates the fractional probability of a connection based on the 
    switching functions used in HbondWire.f90.
    """
    ratio = distance / threshold
    numerator = 1.0 - (ratio ** power_num)
    denominator = 1.0 - (ratio ** power_den)
    
    # Handle edge case where distance exactly equals threshold
    if denominator == 0:
        return 0.0
        
    return numerator / denominator

def calculate_hbond_probability(mod_rOO, mod_rOiH, mod_rOjH):
    """
    Evaluates the continuous hydrogen bond weight between two water molecules.
    Based on Fortran logic: 
    Hb(i,j) = (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**8) / (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**12)
    * (1 - ((mod_rOO - 2.7)/0.5)**10) / (1 - ((mod_rOO - 2.7)/0.5)**20)
    """
    # Angle/Hydrogen placement component
    h_dist_factor = mod_rOiH + mod_rOjH - mod_rOO
    p_angle = switching_function(h_dist_factor, threshold=0.6, power_num=8, power_den=12)
    
    # Oxygen-Oxygen distance component
    oo_dist_factor = mod_rOO - 2.7
    p_dist = switching_function(oo_dist_factor, threshold=0.5, power_num=10, power_den=20)
    
    # Total fractional probability
    return p_angle * p_dist
