import numpy as np

def switching_function(distance, threshold, power_num=8, power_den=12):
    """
    Calculates the fractional probability of a connection based on the
    switching functions used in HbondWire.f90.
    """
    ratio = distance / threshold
    if ratio >= 1.0:
        return 0.0
    if ratio <= 0.0:
        return 1.0

    numerator = 1.0 - (ratio ** power_num)
    denominator = 1.0 - (ratio ** power_den)

    # Handle edge case where distance exactly equals threshold
    if denominator == 0:
        return float(power_num) / float(power_den)

    return numerator / denominator

def calculate_hbond_probability(mod_rOO, mod_rOiH, mod_rOjH, r0_oo=2.80, r0_threshold=0.45):
    """
    Evaluates the continuous hydrogen bond weight between two molecules.
    Based on Fortran logic:
    Hb(i,j) = (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**8) / (1 - ((mod_rOiH+mod_rOjH-mod_rOO)/0.6)**12)
    * (1 - ((mod_rOO - r0_oo)/r0_threshold)**10) / (1 - ((mod_rOO - r0_oo)/r0_threshold)**20)
    """
    # Angle/Hydrogen placement component
    h_dist_factor = mod_rOiH + mod_rOjH - mod_rOO
    p_angle = switching_function(h_dist_factor, threshold=0.6, power_num=8, power_den=12)

    # Heavy atom distance component
    oo_dist_factor = mod_rOO - r0_oo
    p_dist = switching_function(oo_dist_factor, threshold=r0_threshold, power_num=10, power_den=20)

    # Total fractional probability
    return p_angle * p_dist
