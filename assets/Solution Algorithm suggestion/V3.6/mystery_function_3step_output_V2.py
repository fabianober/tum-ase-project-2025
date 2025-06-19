import numpy as np

def expand_9_sym(vec5):
    """Expand [a,b,c,d,e] -> [a,b,c,d,e,d,c,b,a] for 9 stringers."""
    return np.array(list(vec5[:4]) + [vec5[4]] + list(reversed(vec5[:4])))

def expand_10_sym(vec5):
    """Expand [a,b,c,d,e] -> [a,b,c,d,e,d,c,b,a,a] for 10 panels (mirrored, then repeat outer)."""
    return np.array(list(vec5[:4]) + [vec5[4]] + list(reversed(vec5[:4])) + [vec5[0]])

def linspace_index_fractions(length, src_len=9):
    """Return index mapping to spread RFs evenly over expanded stringers/panels."""
    return np.linspace(0, src_len-1, length)

def approximate_fem_behavior(x: np.ndarray,
                             rf_lengths=(150, 100, 50),
                             rf_threshold=1.05,
                             seed_offset=0,
                             noise_amplitude=0.0001):
    """
    Mock FEM with explicit expansion to 9 stringers and 10 panels.
    No (or minimal) noise. Fully symmetric modeling. 
    """

    x = np.asarray(x)
    assert len(x) == 25, "Expected 25 input variables."

    # -- Unpack inputs --
    web_height        = expand_9_sym(x[0:5])
    flange_width      = expand_9_sym(x[5:10])
    lip_height        = expand_9_sym(x[10:15])
    stringer_thickness= expand_9_sym(x[15:20])
    skin_thickness    = expand_10_sym(x[20:25])   # 10 panels

    # -- Weight calculation (physical-ish) --
    panel_areas = np.array([1.15, 1.1, 1.0, 0.8, 0.7, 0.8, 1.0, 1.1, 1.15, 1.15])  # Outer panels slightly bigger
    stringer_areas = np.array([1.1, 1.0, 0.9, 0.85, 0.8, 0.85, 0.9, 1.0, 1.1])      # Middle stringer is smallest

    # Use mean thickness for simplicity, but apply area weighting
    density = 2.75
    weight = (np.sum(stringer_thickness * stringer_areas) * 0.32 +
              np.sum(web_height * stringer_areas) * 0.12 +
              np.sum(flange_width * stringer_areas) * 0.09 +
              np.sum(lip_height * stringer_areas) * 0.04 +
              np.sum(skin_thickness * panel_areas) * density * 1.0)

    # Penalty for lack of symmetry (stddev of each group)
    penalty = (
        np.std(stringer_thickness) +
        np.std(web_height) +
        np.std(flange_width) +
        np.std(lip_height) +
        np.std(skin_thickness)
    )
    weight += 0.7 * penalty

    # -- RF calculation helpers --
    def rf_array(source_vec, rf_len, trend_factor=0.12, min_rf=0.65, max_rf=1.22, center_penalty=0.10, buckling=False):
        # Map RFs to indices in the source_vec
        idx = linspace_index_fractions(rf_len, src_len=len(source_vec))
        # "Central" index is most critical
        mid = (len(source_vec)-1) / 2
        trend = 1 - trend_factor * np.abs(idx - mid) / mid   # RF trend: lowest at center, highest at outer
        # Base is proportional to source thickness (or area)
        base = np.interp(idx, np.arange(len(source_vec)), source_vec)
        # For buckling: outboard panels are penalized more for thinness
        if buckling:
            buck_penalty = 0.12 * (1.5 - base) * (np.abs(idx - mid) / mid)
            buck_penalty = np.clip(buck_penalty, 0, 0.17)
        else:
            buck_penalty = 0.0
        rf = base / 2.2 + 1.01
        rf = rf * trend - buck_penalty
        # Extra penalty for center (mimic load peaking in the center)
        rf -= center_penalty * np.exp(-0.5 * (idx - mid) ** 2 / (0.9 * mid) ** 2)
        # Add *minimal* noise
        rng = np.random.default_rng(int(np.sum(source_vec) * 1e5 + rf_len + seed_offset))
        rf += rng.normal(0, noise_amplitude, rf_len)
        # Clip to realistic range
        rf = np.clip(rf, min_rf, max_rf)
        return rf

    # -- RF_strength: dominated by stringer_thickness (and to lesser extent, web_height) --
    rf_strength = rf_array(stringer_thickness, rf_lengths[0], trend_factor=0.10, min_rf=0.80, max_rf=1.20, center_penalty=0.09)

    # -- RF_stability: influenced by web_height and flange_width, penalize high asymmetry --
    rf_stability = rf_array(web_height*0.7 + flange_width*0.3, rf_lengths[1], trend_factor=0.11, min_rf=0.77, max_rf=1.18, center_penalty=0.10)
    rf_stability -= 0.04 * np.std(web_height)  # Explicit symmetry penalty

    # -- RF_buckling: sensitive to skin_thickness, especially outboard panels --
    rf_buckling = rf_array(skin_thickness, rf_lengths[2], trend_factor=0.14, min_rf=0.60, max_rf=1.16, center_penalty=0.04, buckling=True)

    # -- "Trouble spot": minimum RFs in center, one randomly slightly lower --
    mid_idx_strength = rf_lengths[0] // 2
    mid_idx_stability = rf_lengths[1] // 2
    mid_idx_buckling = rf_lengths[2] // 2
    rf_strength[mid_idx_strength] *= 0.95
    rf_stability[mid_idx_stability] *= 0.94
    rf_buckling[mid_idx_buckling] *= 0.93

    # -- If any variable is very low, some RFs get much worse (emulate "failure") --
    if np.any(x < 0.32):
        n_low = int(0.08 * rf_lengths[0])
        rf_strength[:n_low] = np.minimum(rf_strength[:n_low], 0.82)
        n_low_b = int(0.08 * rf_lengths[2])
        rf_buckling[-n_low_b:] = np.minimum(rf_buckling[-n_low_b:], 0.70)

    # -- Encourage symmetry (add bonus to all if very symmetric) --
    for group in [stringer_thickness, web_height, flange_width, skin_thickness]:
        group_sym = np.std(group)
        bonus = 0.009 / (1.0 + 10.0 * group_sym)
        rf_strength += bonus
        rf_stability += bonus
        rf_buckling += bonus

    # -- Final clipping --
    rf_strength = np.clip(rf_strength, 0.60, 1.21)
    rf_stability = np.clip(rf_stability, 0.60, 1.20)
    rf_buckling = np.clip(rf_buckling, 0.57, 1.17)

    return float(weight), list(rf_strength), list(rf_stability), list(rf_buckling)

# ------------- Example Usage -------------
if __name__ == "__main__":
    # Reasonable starting vector (center values a bit higher, outboard thinner)
    x = np.array([
        2.0, 1.9, 1.8, 1.7, 2.2,        # web_height (center is thickest)
        1.2, 1.1, 1.1, 1.1, 1.2,        # flange_width
        0.7, 0.6, 0.6, 0.6, 0.7,        # lip_height
        0.8, 0.7, 0.7, 0.7, 0.9,        # stringer_thickness
        1.3, 1.1, 1.0, 0.9, 1.5         # skin_thickness (thicker center)
    ])
    w, rf_s, rf_stab, rf_buck = approximate_fem_behavior(x)
    print(f"Weight: {w:.2f}")
    print(f"RF_strength min: {min(rf_s):.4f}, RF_stability min: {min(rf_stab):.4f}, RF_buckling min: {min(rf_buck):.4f}")
