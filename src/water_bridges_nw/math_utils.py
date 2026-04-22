import numpy as np

def switching_function(distance, threshold, k=15):
    """
    Calculates the fractional probability of a connection using a
    Fermi-Dirac logistic function for continuous asymptotic decay.
    P(r) = 1 / (1 + exp(k * (r - r_c)))
    """
    # Protect against overflow in exp
    exponent = k * (distance - threshold)
    if exponent > 100:
        return 0.0
    elif exponent < -100:
        return 1.0

    return 1.0 / (1.0 + np.exp(exponent))

def calculate_hbond_probability(mod_rOO, mod_rOiH, mod_rOjH, r0_oo=2.80, r0_threshold=0.45):
    """
    Evaluates the continuous hydrogen bond weight between two molecules.
    """
    # Angle/Hydrogen placement component
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
