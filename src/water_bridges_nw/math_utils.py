import numpy as np

def switching_function(distance, threshold, power_num=6, power_den=12):
    ratio = distance / threshold
    if ratio >= 1.0:
        return 0.0
    if ratio <= 0.0:
        return 1.0

    return (1.0 - (ratio ** power_num)) / (1.0 - (ratio ** power_den))

def calculate_hbond_probability(mod_rOO, mod_rOiH, mod_rOjH, r0_oo=2.80, r0_threshold=0.45, ignore_angle=False):
    """
    Evaluates the continuous hydrogen bond weight between two molecules.
    """
    # Angle/Hydrogen placement component
    if ignore_angle:
        p_angle = 1.0
    else:
        h_dist_factor = mod_rOiH + mod_rOjH - mod_rOO
        p_angle = switching_function(h_dist_factor, threshold=0.6)

    # Heavy atom distance component
    oo_dist_factor = mod_rOO - r0_oo
    p_dist = switching_function(oo_dist_factor, threshold=r0_threshold)

    # Symmetric Gaussian well to penalize steric clashes
    if mod_rOO < 2.80:
        sigma = 0.15
        p_steric = np.exp(-((mod_rOO - 2.80)**2) / (2 * sigma**2))
    else:
        p_steric = 1.0

    # Total fractional probability
    return p_angle * p_dist * p_steric
