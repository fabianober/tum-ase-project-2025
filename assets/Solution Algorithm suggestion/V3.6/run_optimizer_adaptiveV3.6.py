import os
import sys
import time
import csv
import glob
import pickle
import numpy as np
import random
from datetime import datetime
from tqdm import tqdm, trange
from collections import defaultdict

"""
This algorithm gives a decent result on each run (each run: 5000 FEM calls or ~17h)
The fun part about this one is, that it gets better each time you run it. It's taking past values into consideration when analysing, where to focus on, meaning usualy you'll 
always get that last, extra % out, if you just run it again. So, even if it's off on its first runthrough it'll usually correct on its second runthrough.
Also, this refines all 25 variables. Could obviously be expanded to run pultiple times by default to iron out some of the inconsistency.
"""

# --- USER INPUTS SECTION (CONFIG) ---
USER_INPUTS = {
    # Name and order of variables (later, user can edit this directly)
    'variable_names': [
        # 20 stringer variables: (web_height_1, ..., web_height_5, flange_width_1, ..., thickness_5)
        'web_height_1', 'web_height_2', 'web_height_3', 'web_height_4', 'web_height_5',
        'flange_width_1', 'flange_width_2', 'flange_width_3', 'flange_width_4', 'flange_width_5',
        'lip_height_1', 'lip_height_2', 'lip_height_3', 'lip_height_4', 'lip_height_5',
        'thickness_1', 'thickness_2', 'thickness_3', 'thickness_4', 'thickness_5',
        # 5 skin panel thicknesses
        'skin_thickness_1', 'skin_thickness_2', 'skin_thickness_3', 'skin_thickness_4', 'skin_thickness_5'
    ],
    'bounds': [(0.3, 15.0)] * 25,
    'target_rf': 1.05,
    'fem_call_budget': 5000,  # Can be increased
    'log_folder': "optimizer_logs",
    'run_id_prefix': "wing_optim",
    'random_seed': None,
    'sampling_fraction': 0.15,  # Fraction of budget for initial random sampling
    'clustering_fraction': 0.25,  # Fraction of budget for region exploration
    'ga_fraction': 0.5,  # Fraction of budget for GAs in clustered regions
    'refine_fraction': 0.10,  # Remaining for final local refinement
    'cluster_n': 4,  # Number of clusters/regions (can tune)
    'report_every': 500,
    'resume': True,  # Attempt to read and use previous logs
    'max_ga_generations': 24,
    'ga_pop_size': 16,
    'display_top_n': 3,
}

# --- END USER INPUTS ---

# Set seed for reproducibility
random.seed(USER_INPUTS['random_seed'])
np.random.seed(USER_INPUTS['random_seed'])

# Import the FEM call
try:
    from mystery_function_3step_output_V2 import approximate_fem_behavior
except ImportError:
    print("ERROR: Could not import approximate_fem_behavior.")
    sys.exit(1)

# --- Utility Functions ---

def make_log_folder(log_folder):
    if not os.path.exists(log_folder):
        os.makedirs(log_folder)

def unique_run_name(prefix, log_folder):
    date = datetime.now().strftime("%Y%m%d_%H%M%S")
    base = f"{prefix}_{date}"
    # If already exists, increment number
    idx = 1
    while os.path.exists(os.path.join(log_folder, f"{base}_calls.csv")):
        base = f"{prefix}_{date}_{idx}"
        idx += 1
    return base

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')

def flatten_rf_list(rf):
    # Converts list of floats to string, e.g. "1.22;1.01;1.37"
    return ";".join([f"{x:.5f}" for x in rf])

def parse_rf_list(rf_str):
    # Back to list of floats
    return [float(x) for x in rf_str.split(';')] if rf_str else []

def log_headers(n_strength, n_stability, n_buckling):
    return [
        "timestamp", "stage", "region", "fem_call_idx",
        *USER_INPUTS['variable_names'],
        "weight",
        *(f"rf_strength_{i}" for i in range(n_strength)),
        *(f"rf_stability_{i}" for i in range(n_stability)),
        *(f"rf_buckling_{i}" for i in range(n_buckling)),
    ]

def read_previous_logs(log_folder, run_prefix):
    # Find all logs matching prefix
    files = glob.glob(os.path.join(log_folder, f"{run_prefix}_*_calls.csv"))
    if not files:
        return []
    all_rows = []
    for file in files:
        try:
            with open(file, 'r', newline='') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    all_rows.append(row)
        except Exception as e:
            print(f"Warning: Failed to read {file} ({e})")
    return all_rows

def save_csv_row(file_path, headers, row):
    file_exists = os.path.isfile(file_path)
    with open(file_path, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)

def save_summary(log_folder, run_base, summary, header, append=False):
    file_path = os.path.join(log_folder, f"{run_base}_report.csv")
    file_exists = os.path.isfile(file_path)
    with open(file_path, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header)
        if not file_exists:
            writer.writeheader()
        writer.writerow(summary)

def get_rf_lengths(y):
    # y = (weight, rf_strength, rf_stability, rf_buckling)
    return len(y[1]), len(y[2]), len(y[3])

def print_summary(best, near_misses, stage, region, calls_used, calls_total):
    clear_console()
    print(f"=== Optimization Progress (Stage: {stage} | Region: {region}) ===")
    print(f"Total FEM calls used: {calls_used} / {calls_total}\n")
    print(f"Current Best (All RFs ≥ {USER_INPUTS['target_rf']}):")
    if best:
        print(f"  Mass: {best['weight']:.3f} | Inputs: {np.round(best['x'],3)}")
    else:
        print("  No valid solution found yet.")
    print("\nNear Misses (Lowest weights with few violations):")
    for i, nm in enumerate(near_misses):
        print(f"  [{i+1}] Mass: {nm['weight']:.3f} | Violations: {nm['n_bad']} | Inputs: {np.round(nm['x'],3)}")
    print("-" * 65)

# --- Core Optimization Functions ---

def sample_random_points(n, bounds, reduce_to_5d=False):
    # If reduce_to_5d, only sample the 5 meta-variables and spread across 25 as per pattern
    points = []
    for _ in range(n):
        if reduce_to_5d:
            # 5 meta-vars: [web_height, flange_width, lip_height, thickness, skin_thickness]
            metas = [random.uniform(0.01, 30.0) for _ in range(5)]
            # Expand across the 25D vector: (assume equal or simple pattern)
            # web_height_1-5, flange_width_1-5, ...
            x = []
            for i in range(5):  # For each type (stringer/skintype)
                x.append(metas[0])  # web_height
            for i in range(5):
                x.append(metas[1])  # flange_width
            for i in range(5):
                x.append(metas[2])  # lip_height
            for i in range(5):
                x.append(metas[3])  # thickness
            for i in range(5):
                x.append(metas[4] * (1 - i * 0.12))  # skin_thickness, linearly decrease
            # Clip all to bounds
            x = [max(bounds[i][0], min(bounds[i][1], v)) for i, v in enumerate(x)]
        else:
            x = [random.uniform(low, high) for (low, high) in bounds]
        points.append(x)
    return points

def fem_eval(x):
    try:
        y = approximate_fem_behavior(np.array(x))
        # y: (weight, rf_strength, rf_stability, rf_buckling)
        return y
    except Exception as e:
        print(f"FEM call failed for x={x}: {e}")
        return None

def evaluate_points(xs, stage, region, start_fem_idx, log_path, headers, n_strength, n_stability, n_buckling, progress_bar=None):
    results = []
    for i, x in enumerate(xs):
        t0 = time.time()
        y = fem_eval(x)
        t1 = time.time()
        if y is None:
            continue
        weight, rf_strength, rf_stability, rf_buckling = y
        row = {
            "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "stage": stage,
            "region": region,
            "fem_call_idx": start_fem_idx + i,
        }
        for idx, name in enumerate(USER_INPUTS['variable_names']):
            row[name] = x[idx]
        row["weight"] = weight
        # Flatten RFs for logging
        for idx in range(n_strength):
            row[f"rf_strength_{idx}"] = rf_strength[idx] if idx < len(rf_strength) else ""
        for idx in range(n_stability):
            row[f"rf_stability_{idx}"] = rf_stability[idx] if idx < len(rf_stability) else ""
        for idx in range(n_buckling):
            row[f"rf_buckling_{idx}"] = rf_buckling[idx] if idx < len(rf_buckling) else ""
        save_csv_row(log_path, headers, row)
        results.append({
            "x": np.array(x),
            "weight": weight,
            "rf_strength": np.array(rf_strength),
            "rf_stability": np.array(rf_stability),
            "rf_buckling": np.array(rf_buckling),
            "call_idx": start_fem_idx + i,
            "runtime": t1-t0,
        })
        if progress_bar:
            progress_bar.update(1)
    return results

def cluster_points(points, n_clusters):
    # points: list of dicts with x
    from sklearn.cluster import KMeans
    X = np.array([p['x'] for p in points])
    # Use only the 5 "meta"-features: mean web_height, mean flange_width, mean lip_height, mean thickness, mean skin_thickness
    # Assuming order as per config
    X5 = np.stack([
        np.mean(X[:, 0:5], axis=1),  # web_height
        np.mean(X[:, 5:10], axis=1), # flange_width
        np.mean(X[:, 10:15], axis=1), # lip_height
        np.mean(X[:, 15:20], axis=1), # thickness
        np.mean(X[:, 20:25], axis=1), # skin_thickness
    ], axis=1)
    # Remove points with NaN or Inf
    mask = np.all(np.isfinite(X5), axis=1)
    X5 = X5[mask]
    clean_points = [p for i,p in enumerate(points) if mask[i]]
    if len(X5) < n_clusters:
        n_clusters = max(2, len(X5))
    kmeans = KMeans(n_clusters=n_clusters, random_state=USER_INPUTS['random_seed'])
    labels = kmeans.fit_predict(X5)
    regions = [[] for _ in range(n_clusters)]
    for i, lbl in enumerate(labels):
        regions[lbl].append(clean_points[i])
    centers = kmeans.cluster_centers_
    return regions, centers

def genetic_algorithm(region_points, n_strength, n_stability, n_buckling, log_path, headers, fem_call_start, max_calls, stage, region_idx):
    # Simple (μ+λ) GA, population = region_points + randoms
    pop_size = USER_INPUTS['ga_pop_size']
    generations = USER_INPUTS['max_ga_generations']
    bounds = USER_INPUTS['bounds']
    # Select top N solutions so far as seeds
    sorted_pts = sorted(region_points, key=lambda p: (np.min(np.concatenate([p['rf_strength'], p['rf_stability'], p['rf_buckling']])) < USER_INPUTS['target_rf'], p['weight']))
    seeds = [p['x'] for p in sorted_pts[:max(2, pop_size//2)]]
    # Fill with randoms
    while len(seeds) < pop_size:
        seeds.append(np.array([random.uniform(l, h) for l, h in bounds]))
    population = [np.copy(x) for x in seeds]
    results = []
    fem_call_idx = fem_call_start
    for gen in range(generations):
        new_pop = []
        # Selection: Tournament on best (by valid weight, then weight)
        pop_metrics = []
        for x in population:
            y = fem_eval(x)
            if y is None:
                continue
            weight, rf_strength, rf_stability, rf_buckling = y
            rf_min = min(np.min(rf_strength), np.min(rf_stability), np.min(rf_buckling))
            valid = rf_min >= USER_INPUTS['target_rf']
            pop_metrics.append((valid, weight, x, y))
            # Log
            row = {
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "stage": stage,
                "region": region_idx,
                "fem_call_idx": fem_call_idx,
            }
            for idx, name in enumerate(USER_INPUTS['variable_names']):
                row[name] = x[idx]
            row["weight"] = weight
            for idx in range(n_strength):
                row[f"rf_strength_{idx}"] = rf_strength[idx] if idx < len(rf_strength) else ""
            for idx in range(n_stability):
                row[f"rf_stability_{idx}"] = rf_stability[idx] if idx < len(rf_stability) else ""
            for idx in range(n_buckling):
                row[f"rf_buckling_{idx}"] = rf_buckling[idx] if idx < len(rf_buckling) else ""
            save_csv_row(log_path, headers, row)
            fem_call_idx += 1
            results.append({
                "x": np.copy(x),
                "weight": weight,
                "rf_strength": np.array(rf_strength),
                "rf_stability": np.array(rf_stability),
                "rf_buckling": np.array(rf_buckling),
                "call_idx": fem_call_idx,
                "runtime": None,
            })
            if fem_call_idx >= fem_call_start + max_calls:
                return results
        # Elitism: top 2 retained
        pop_metrics = sorted(pop_metrics, key=lambda t: (not t[0], t[1]))
        elites = [np.copy(t[2]) for t in pop_metrics[:2]]
        # Crossover and Mutation to fill population
        while len(new_pop) < pop_size:
            p1 = random.choice(population)
            p2 = random.choice(population)
            cross = np.where(np.random.rand(len(p1)) < 0.5, p1, p2)
            # Small mutation: add Gaussian noise to each variable
            mutation = np.random.normal(0, 0.04 * (np.array([h-l for l,h in bounds])), len(p1))
            child = cross + mutation
            # Clip to bounds
            child = np.clip(child, [l for l,_ in bounds], [h for _,h in bounds])
            new_pop.append(child)
        population = elites + new_pop[:pop_size-2]
        if fem_call_idx >= fem_call_start + max_calls:
            break
    return results

def find_best_and_near_misses(results, target_rf, top_n=3):
    # Only count as valid if *all* RFs ≥ target_rf
    valids = []
    near_misses = []
    for r in results:
        rf_all = np.concatenate([r['rf_strength'], r['rf_stability'], r['rf_buckling']])
        rf_min = np.min(rf_all)
        n_bad = np.sum(rf_all < target_rf)
        if n_bad == 0:
            valids.append(r)
        elif n_bad < 5:
            near_misses.append({**r, 'n_bad': n_bad})
    if valids:
        best = min(valids, key=lambda r: r['weight'])
    else:
        best = None
    near_misses = sorted(near_misses, key=lambda r: (r['n_bad'], r['weight']))[:top_n]
    return best, near_misses

# --- Main Program ---

def main():
    make_log_folder(USER_INPUTS['log_folder'])
    # Unique run name
    run_base = unique_run_name(USER_INPUTS['run_id_prefix'], USER_INPUTS['log_folder'])
    log_path = os.path.join(USER_INPUTS['log_folder'], f"{run_base}_calls.csv")

    # Load any previous logs
    all_results = []
    n_strength = n_stability = n_buckling = None
    if USER_INPUTS['resume']:
        prev_rows = read_previous_logs(USER_INPUTS['log_folder'], USER_INPUTS['run_id_prefix'])
        if prev_rows:
            print(f"Found {len(prev_rows)} rows from previous runs. Incorporating for initial analysis.")
            for row in prev_rows:
                try:
                    x = np.array([float(row[v]) for v in USER_INPUTS['variable_names']])
                    weight = float(row['weight'])
                    rf_strength = [float(row[f"rf_strength_{i}"]) for i in range(1000) if f"rf_strength_{i}" in row and row[f"rf_strength_{i}"] != ""]
                    rf_stability = [float(row[f"rf_stability_{i}"]) for i in range(1000) if f"rf_stability_{i}" in row and row[f"rf_stability_{i}"] != ""]
                    rf_buckling = [float(row[f"rf_buckling_{i}"]) for i in range(1000) if f"rf_buckling_{i}" in row and row[f"rf_buckling_{i}"] != ""]
                    all_results.append({
                        "x": x, "weight": weight,
                        "rf_strength": np.array(rf_strength),
                        "rf_stability": np.array(rf_stability),
                        "rf_buckling": np.array(rf_buckling),
                        "call_idx": int(row.get("fem_call_idx", -1))
                    })
                    if n_strength is None: n_strength = len(rf_strength)
                    if n_stability is None: n_stability = len(rf_stability)
                    if n_buckling is None: n_buckling = len(rf_buckling)
                except Exception as e:
                    print(f"Error parsing previous row: {e}")
    fem_calls_used = len(all_results)
    # If no prior data, need to run 1 to get RF lengths
    if n_strength is None or n_stability is None or n_buckling is None:
        print("Performing an initial FEM call to determine RF output lengths...")
        test_x = [random.uniform(l, h) for l, h in USER_INPUTS['bounds']]
        y = fem_eval(test_x)
        if y is None:
            print("ERROR: Initial FEM call failed, cannot proceed.")
            return
        n_strength, n_stability, n_buckling = get_rf_lengths(y)

    headers = log_headers(n_strength, n_stability, n_buckling)

    # --- STAGE 1: Initial Random Sampling (using reduced 5D) ---
    print("Stage 1: Initial broad random sampling (reduced 5D -> 25D)...")
    n_sample = int(USER_INPUTS['fem_call_budget'] * USER_INPUTS['sampling_fraction'])
    random_points = sample_random_points(n_sample, USER_INPUTS['bounds'], reduce_to_5d=True)
    pbar1 = tqdm(total=n_sample, desc="Random Sampling", ncols=70)
    results1 = evaluate_points(random_points, "sampling", -1, fem_calls_used, log_path, headers, n_strength, n_stability, n_buckling, pbar1)
    pbar1.close()
    all_results.extend(results1)
    fem_calls_used += len(results1)
    best, near_misses = find_best_and_near_misses(all_results, USER_INPUTS['target_rf'], USER_INPUTS['display_top_n'])
    print_summary(best, near_misses, "sampling", "-", fem_calls_used, USER_INPUTS['fem_call_budget'])
    save_summary(USER_INPUTS['log_folder'], run_base, {
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "stage": "sampling",
        "calls_used": fem_calls_used,
        "best_weight": best['weight'] if best else "",
        "n_valid": len([r for r in all_results if np.min(np.concatenate([r['rf_strength'], r['rf_stability'], r['rf_buckling']])) >= USER_INPUTS['target_rf']])
    }, ["timestamp","stage","calls_used","best_weight","n_valid"], append=True)

    # --- STAGE 2: Clustering/Region Detection ---
    print("Stage 2: Clustering and region selection...")
    n_clusters = USER_INPUTS['cluster_n']
    regions, centers = cluster_points(all_results, n_clusters)
    print(f"Identified {n_clusters} regions for further search.")
    region_budgets = [int(USER_INPUTS['fem_call_budget'] * USER_INPUTS['clustering_fraction'] / n_clusters)] * n_clusters
    for region_idx, reg_points in enumerate(regions):
        print(f"\nExploring Region {region_idx+1}/{n_clusters} ({len(reg_points)} seeds)...")
        # Sample additional points near center, add small noise to seeds
        center = centers[region_idx]
        # Reconstruct 25D vectors from center
        x_seed = []
        for i in range(5): x_seed.append(center[0])
        for i in range(5): x_seed.append(center[1])
        for i in range(5): x_seed.append(center[2])
        for i in range(5): x_seed.append(center[3])
        for i in range(5): x_seed.append(center[4] * (1 - i*0.12))
        n_extra = max(0, region_budgets[region_idx] - len(reg_points))
        nearby_points = [np.clip(x_seed + np.random.normal(0, 0.08, 25), [l for l,_ in USER_INPUTS['bounds']], [h for _,h in USER_INPUTS['bounds']]).tolist() for _ in range(n_extra)]
        pbar2 = tqdm(total=n_extra, desc=f"Region {region_idx+1} Sampling", ncols=70)
        results2 = evaluate_points(nearby_points, "region_sampling", region_idx, fem_calls_used, log_path, headers, n_strength, n_stability, n_buckling, pbar2)
        pbar2.close()
        all_results.extend(results2)
        fem_calls_used += len(results2)
        if fem_calls_used % USER_INPUTS['report_every'] < len(results2):
            best, near_misses = find_best_and_near_misses(all_results, USER_INPUTS['target_rf'], USER_INPUTS['display_top_n'])
            print_summary(best, near_misses, "region_sampling", region_idx, fem_calls_used, USER_INPUTS['fem_call_budget'])
            save_summary(USER_INPUTS['log_folder'], run_base, {
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "stage": f"region_{region_idx}",
                "calls_used": fem_calls_used,
                "best_weight": best['weight'] if best else "",
                "n_valid": len([r for r in all_results if np.min(np.concatenate([r['rf_strength'], r['rf_stability'], r['rf_buckling']])) >= USER_INPUTS['target_rf']])
            }, ["timestamp","stage","calls_used","best_weight","n_valid"], append=True)
        if fem_calls_used >= USER_INPUTS['fem_call_budget']:
            break

    # --- STAGE 3: Local Genetic Algorithms in Top Regions ---
    print("\nStage 3: Local Genetic Algorithm optimization in each region...")
    # Choose best regions (by valid count or best near-threshold)
    region_best_weights = [min((r['weight'] for r in reg if np.min(np.concatenate([r['rf_strength'], r['rf_stability'], r['rf_buckling']])) >= USER_INPUTS['target_rf']), default=np.inf) for reg in regions]
    region_order = np.argsort(region_best_weights)
    ga_budget_per_region = int(USER_INPUTS['fem_call_budget'] * USER_INPUTS['ga_fraction'] / n_clusters)
    for idx in region_order:
        reg = regions[idx]
        print(f"\nRunning GA in Region {idx+1}/{n_clusters} ...")
        results_ga = genetic_algorithm(reg, n_strength, n_stability, n_buckling, log_path, headers, fem_calls_used, ga_budget_per_region, "ga", idx)
        all_results.extend(results_ga)
        fem_calls_used += len(results_ga)
        if fem_calls_used % USER_INPUTS['report_every'] < len(results_ga):
            best, near_misses = find_best_and_near_misses(all_results, USER_INPUTS['target_rf'], USER_INPUTS['display_top_n'])
            print_summary(best, near_misses, "GA", idx, fem_calls_used, USER_INPUTS['fem_call_budget'])
            save_summary(USER_INPUTS['log_folder'], run_base, {
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "stage": f"GA_{idx}",
                "calls_used": fem_calls_used,
                "best_weight": best['weight'] if best else "",
                "n_valid": len([r for r in all_results if np.min(np.concatenate([r['rf_strength'], r['rf_stability'], r['rf_buckling']])) >= USER_INPUTS['target_rf']])
            }, ["timestamp","stage","calls_used","best_weight","n_valid"], append=True)
        if fem_calls_used >= USER_INPUTS['fem_call_budget']:
            break

    # --- STAGE 4: Final Refinement (Pattern search/Local) ---
    print("\nStage 4: Final refinement in best region...")
    best, _ = find_best_and_near_misses(all_results, USER_INPUTS['target_rf'], USER_INPUTS['display_top_n'])
    if best is not None:
        from scipy.optimize import minimize
        def penalty_obj(x):
            y = fem_eval(x)
            if y is None:
                return 1e8
            weight, rf_strength, rf_stability, rf_buckling = y
            rfs = np.concatenate([rf_strength, rf_stability, rf_buckling])
            rf_min = np.min(rfs)
            penalty = 0
            if rf_min < USER_INPUTS['target_rf']:
                penalty = 1e5 * (USER_INPUTS['target_rf'] - rf_min)
            return weight + penalty
        bounds = USER_INPUTS['bounds']
        refine_calls = int(USER_INPUTS['fem_call_budget'] * USER_INPUTS['refine_fraction'])
        x0 = best['x']
        options = {'maxiter': refine_calls, 'disp': False}
        res = minimize(penalty_obj, x0, method='Nelder-Mead', bounds=bounds, options=options)
        # Final eval of result
        y = fem_eval(res.x)
        if y is not None:
            weight, rf_strength, rf_stability, rf_buckling = y
            row = {
                "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "stage": "refine",
                "region": "best",
                "fem_call_idx": fem_calls_used,
            }
            for idx, name in enumerate(USER_INPUTS['variable_names']):
                row[name] = res.x[idx]
            row["weight"] = weight
            for idx in range(n_strength):
                row[f"rf_strength_{idx}"] = rf_strength[idx] if idx < len(rf_strength) else ""
            for idx in range(n_stability):
                row[f"rf_stability_{idx}"] = rf_stability[idx] if idx < len(rf_stability) else ""
            for idx in range(n_buckling):
                row[f"rf_buckling_{idx}"] = rf_buckling[idx] if idx < len(rf_buckling) else ""
            save_csv_row(log_path, headers, row)
            all_results.append({
                "x": np.copy(res.x),
                "weight": weight,
                "rf_strength": np.array(rf_strength),
                "rf_stability": np.array(rf_stability),
                "rf_buckling": np.array(rf_buckling),
                "call_idx": fem_calls_used,
                "runtime": None,
            })
            fem_calls_used += 1

    # --- Final Summary ---
    best, near_misses = find_best_and_near_misses(all_results, USER_INPUTS['target_rf'], USER_INPUTS['display_top_n'])
    print_summary(best, near_misses, "FINISHED", "-", fem_calls_used, USER_INPUTS['fem_call_budget'])
    print("\n=== MIDPOINTS OF FOUND REGIONS ===")
    for i, ctr in enumerate(centers):
        print(f"Region {i+1}: {np.round(ctr,3)}")
    print("\nBest found solution (input vector, mass, min RF):")
    if best:
        print(f"Input vector: {np.round(best['x'],3)}\nWeight: {best['weight']:.3f}\nMin RF: {min(np.min(best['rf_strength']), np.min(best['rf_stability']), np.min(best['rf_buckling'])):.4f}")
    else:
        print("No valid solution found.")
    print("\nAll logs are in:", USER_INPUTS['log_folder'])
    print("Run complete.")

if __name__ == '__main__':
    main()
