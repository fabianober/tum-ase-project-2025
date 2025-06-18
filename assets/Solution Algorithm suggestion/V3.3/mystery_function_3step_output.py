import numpy as np
from time import sleep

"""
! This function is just a placeholder for the moment
"""

def approximate_fem_behavior(t_skin, t_stringer, w_stringer):
    # Constants
    N_elements = 66
    N_panels = 10
    N_stiffeners = 9
    LC = 3  # Load cases

    # Physical constants
    skin_area = 10000  # mm^2
    stringer_length = 1000  # mm
    density = 2.7e-3  # g/mm^3

    # --- Weight model ---
    skin_weight = density * skin_area * t_skin
    stringer_weight = density * N_stiffeners * 3 * w_stringer * t_stringer * stringer_length
    total_weight = skin_weight + stringer_weight

    # Introduce some mild nonlinear coupling in performance
    thickness_ratio = t_skin / max(t_stringer, 1e-3)
    geometry_factor = np.sin(0.1 * w_stringer) + np.log1p(t_stringer)

    # --- Strength RFs (element-level) ---
    base_strength_rf = 0.8 + 4.0 * (t_skin + 0.3 * t_stringer) + 0.5 * np.cos(thickness_ratio)
    strength_noise = np.random.normal(0, 0.6, N_elements * LC)
    rf_strength = np.clip(base_strength_rf + strength_noise, 0.2, 25).tolist()

    # --- Stability RFs (panel-level) ---
    base_stability_rf = 0.5 + 2.2 * t_skin - 0.3 * np.sqrt(w_stringer) + 0.2 * np.sin(t_stringer)
    stability_noise = np.random.normal(0, 0.08, N_panels * LC)
    rf_stability = np.clip(base_stability_rf + stability_noise, 0.25, 2.0).tolist()

    # --- Buckling RFs (stiffener-level) ---
    interaction_term = (t_stringer * w_stringer) ** 0.8
    base_buckling_rf = 0.4 + 1.5 * interaction_term - 0.1 * thickness_ratio + 0.1 * geometry_factor
    buckling_noise = np.random.normal(0, 0.12, N_stiffeners * LC)
    rf_buckling = np.clip(base_buckling_rf + buckling_noise, 0.2, 1.8).tolist()

    # Simulate latency
    sleep(0.01)

    return float(total_weight), rf_strength, rf_stability, rf_buckling

def run_fem_simulation(x):
    return approximate_fem_behavior(x[0], x[1], x[2])