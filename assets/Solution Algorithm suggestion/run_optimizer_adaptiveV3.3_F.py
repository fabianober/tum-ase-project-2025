import numpy as np
from mystery_function_3step_output import run_fem_simulation
import csv
import os
import pickle
from sklearn.cluster import KMeans
import logging
import traceback
from typing import List, Tuple, Any
import sys
import time

"""
Algorithm takes around 5-7 hours to finish when approximating a single FEM call time with a time of 5 seconds per call.
Algorithm is optimized to reduce the number of calles to the FEM function (max. 4900 per run)
Initially this was designed for 10 input variables, so it might be a tad over-engineered. 
All expensive calls to the FEM function are being logged in the CSV file, which will be created as the ptogream runs.
"""

# === USER SETTINGS: SET THESE! ===
NUM_VARS = 3  # <<<--- Set this to match your mystery function!
EXPECTED_RESERVE_FACTOR_COUNT = 255  # <<<--- Set to correct output size (number of RFs returned)

BOUNDS = [
    (0.1, 4.5),    # skin thickness
    (0.1, 5.0),    # stringer thickness
    (0.1, 20.0)   # stringer width
]  # <<<--- Must have NUM_VARS tuples!

"""
The clearer the bounds are set the better the result. Large bounds will lead to slightly longer computation
If your final result variables are very close to the set bounds, rethink your bounds. You might be limiting a better solution
"""


# === DO NOT EDIT BELOW (unless you know what you're doing) ===

SAMPLES = 1200
NUM_CLUSTERS = 6
GENS_PER_CLUSTER = 12
POP_CLUSTER = 20
REFINE_GENS = 25
REFINE_POP = 30
RF_THRESHOLD = 1.01
MUTATION_RATE = 0.25
ELITE_PER_CLUSTER = 2
MAX_MEMORY_SIZE = 20
DIVERSITY_THRESHOLD = 0.05

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MEMORY_FILE = os.path.join(BASE_DIR, "elite_memory.pkl")
LOG_FILE = os.path.join(BASE_DIR, "smart_ga_log.csv")
ERROR_LOG_FILE = os.path.join(BASE_DIR, "smart_ga_errors.log")

GLOBAL_RANDOM_SEED = 42
#np.random.seed(GLOBAL_RANDOM_SEED)  # <--- This will give repeatable results if necessary.
np.random.seed(None) # <--- This gives random seeds and a fresh run each time.

def cursor_up(n=1): sys.stdout.write(f"\033[{n}A")
def clear_line(): sys.stdout.write('\033[2K')

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')

def normalize_vector(v):
    lows = np.array([b[0] for b in BOUNDS])
    highs = np.array([b[1] for b in BOUNDS])
    return (v - lows) / (highs - lows + 1e-12)

def clip_to_bounds(ind):
    for i in range(NUM_VARS):
        ind[i] = np.clip(ind[i], BOUNDS[i][0], BOUNDS[i][1])
    return ind

def safe_file_remove(filename):
    try:
        if os.path.exists(filename): os.remove(filename)
    except Exception as e: logging.error(f"Failed to remove file {filename}: {e}")

def write_log_header():
    param_headers = [f"x{i}" for i in range(NUM_VARS)]
    rf_headers = [f"RF_{i}" for i in range(EXPECTED_RESERVE_FACTOR_COUNT)]
    with open(LOG_FILE, mode="w", newline="") as file:
        csv.writer(file).writerow(
            ["phase", "step", "call_id", "mass", "feasible", "violating_RF_count"] +
            param_headers + rf_headers
        )

def safe_csv_write(phase, step, call_id, mass, feasible, violating_RF_count, params, RFs):
    param_list = list(params) if isinstance(params, (np.ndarray, list)) else [params]
    rf_list = list(RFs) if isinstance(RFs, (np.ndarray, list)) else [RFs]
    row = [phase, step, call_id, mass, feasible, violating_RF_count] + \
        [float(x) for x in param_list] + [float(x) for x in rf_list]
    with open(LOG_FILE, mode="a", newline="") as file:
        csv.writer(file).writerow(row)

def fitness(ind, penalty_scale=1000, log_label=None, log_step=None, call_counter=[0], progress=None):
    """Compute penalized mass, log all evaluations, and update progress."""
    try:
        mass, RFs = safe_run_fem(ind)
        total_violation = sum(max(0, RF_THRESHOLD - rf) for rf in RFs)
        penalty = penalty_scale * total_violation ** 2
        feasible = all(rf >= RF_THRESHOLD for rf in RFs)
        violating_count = sum(rf < RF_THRESHOLD for rf in RFs)
        call_counter[0] += 1
        if progress is not None:
            progress['calls_done'] += 1
            progress['last_update_time'] = time.time()
        safe_csv_write(log_label or "", log_step or "", call_counter[0], mass, feasible, violating_count, ind, RFs)
        return (mass + penalty, feasible, violating_count, mass, RFs)
    except Exception as e:
        logging.error(f"Fitness evaluation failed: {e}")
        logging.error(traceback.format_exc())
        return float('inf'), False, EXPECTED_RESERVE_FACTOR_COUNT, float('inf'), [0.0]*EXPECTED_RESERVE_FACTOR_COUNT

def safe_run_fem(ind, fail_limit=3):
    for attempt in range(fail_limit):
        try:
            result = run_fem_simulation(ind)
            if not isinstance(result, tuple) or len(result) != 4: raise ValueError("Blackbox output must be a tuple (mass, RFs)")
            mass, rf_strength, rf_stability, rf_buckling = result
            RFs = list(rf_strength) + list(rf_stability) + list(rf_buckling)
            if not isinstance(mass, float) or not all(isinstance(r, (list, np.ndarray)) for r in [rf_strength, rf_stability, rf_buckling]): raise ValueError("Unexpected types from blackbox function.")
            if len(RFs) != EXPECTED_RESERVE_FACTOR_COUNT: raise ValueError(f"Expected {EXPECTED_RESERVE_FACTOR_COUNT} reserve factors, got {len(RFs)}.")
            return mass, list(RFs)
        except Exception as e:
            logging.error(f"safe_run_fem failed for input {ind}: {e}")
            logging.error(traceback.format_exc())
    return float('inf'), [0.0] * EXPECTED_RESERVE_FACTOR_COUNT

def generate_population(n): return [np.array([np.random.uniform(low, high) for (low, high) in BOUNDS]) for _ in range(n)]
def crossover(p1, p2): alpha = np.random.rand(); child = alpha * p1 + (1 - alpha) * p2; return clip_to_bounds(child)
def mutate(ind, rate): 
    for i in range(NUM_VARS):
        if np.random.rand() < rate:
            rng = BOUNDS[i][1] - BOUNDS[i][0]
            ind[i] += np.random.normal(0, 0.1 * rng)
    return clip_to_bounds(ind)
def is_unique(new_ind, others): new_norm = normalize_vector(new_ind); return all(np.linalg.norm(new_norm - normalize_vector(o)) > DIVERSITY_THRESHOLD for o in others)

def render_bar(name, value, total, width=40, show_perc=True, show_counter=True, extra=""):
    perc = value / total if total > 0 else 1.0
    filled = int(width*perc)
    bar = "#"*filled + "-"*(width-filled)
    pc = f"{perc*100:5.1f}%" if show_perc else ""
    cnt = f"{value}/{total}" if show_counter else ""
    s = f"{name:10s} [{bar}] {pc} {cnt}"
    if extra: s += "  " + extra
    return s

def eta_str(start_time, done, total):
    now = time.time()
    elapsed = now - start_time
    if done==0: return "--:--:--"
    remain = (elapsed/done)*(total-done)
    return time.strftime('%H:%M:%S', time.gmtime(remain))

def refresh_display(lines):
    cursor_up(len(lines))
    for l in lines: clear_line(); print(l)
    sys.stdout.flush()

def main():
    safe_file_remove(LOG_FILE)
    safe_file_remove(ERROR_LOG_FILE)
    write_log_header()

    main_total = SAMPLES + (NUM_CLUSTERS * GENS_PER_CLUSTER * POP_CLUSTER) + (3 * REFINE_GENS * REFINE_POP)
    progress = dict(calls_done=0, total_calls=main_total, start_time=time.time(), last_update_time=time.time())

    def update_all_bars(main_txt, sub_txt="", subsub_txt="", subsubsub_txt="", subsubsubsub_txt=""):
        eta = eta_str(progress['start_time'], progress['calls_done'], progress['total_calls'])
        lines = [
            render_bar("TOTAL", progress['calls_done'], progress['total_calls'], show_perc=True, show_counter=False, extra=f"ETA: {eta}"),
            main_txt,
            sub_txt,
            subsub_txt,
            subsubsub_txt,
            ""
        ]
        refresh_display(lines)

    # ==== PHASE 1 ====
    clear_console()
    print("="*65)
    print("PHASE 1: GLOBAL SAMPLING")
    print("="*65)
    print("Sampling random designs to explore the search space...\n")
    print('\n' * 5)  # Reserve lines for bars

    sample_data = []
    feasible_results = []
    all_vectors = generate_population(SAMPLES)

    for i, vec in enumerate(all_vectors):
        main_txt = render_bar("Sampling", i+1, SAMPLES, extra="")
        update_all_bars(main_txt)
        fit, feasible, violating, mass, RFs = fitness(vec, log_label="Sampling", log_step=i, progress=progress)
        if feasible:
            sample_data.append(vec)
        feasible_results.append((fit, vec))
    update_all_bars(render_bar("Sampling", SAMPLES, SAMPLES))
    print(f"[INFO] Sampling complete. Feasible designs found: {len(sample_data)} out of {SAMPLES}\n")
    time.sleep(5.2)  # brief pause so you can see it finished

    # ==== PHASE 2 ====
    clear_console()
    print("="*65)
    print("PHASE 2: CLUSTERING AND REGIONAL GENETIC ALGORITHMS")
    print("="*65)
    print("Grouping the best designs into clusters, then running GAs in each region.\n")
    print('\n'*5)

    # --- Filtered sample data for clustering ---
    filtered_sample_data = [
        vec for (fit, vec) in sorted(feasible_results, key=lambda x: x[0])[:max(100, NUM_CLUSTERS * 10)]
    ]
    kmeans = KMeans(n_clusters=NUM_CLUSTERS, n_init=10, random_state=GLOBAL_RANDOM_SEED).fit(
        [np.array(x) for x in filtered_sample_data]
    )
    cluster_centers = kmeans.cluster_centers_
    best_cluster_results = []
    call_counter = [progress['calls_done']]

    # ==== Robust elite memory loading (auto-delete if wrong shape) ====
    past_elites = []
    elite_memory_corrupt = False
    if os.path.exists(MEMORY_FILE):
        try:
            with open(MEMORY_FILE, "rb") as f:
                loaded_elites = pickle.load(f)
            for ind in loaded_elites:
                if len(ind) != NUM_VARS:
                    elite_memory_corrupt = True
                    break
            if not elite_memory_corrupt:
                past_elites = loaded_elites
            else:
                print("[WARNING] Detected old elite_memory.pkl with mismatched dimensionality. Deleting it for safety.")
                os.remove(MEMORY_FILE)
        except Exception as e:
            print(f"[WARNING] Could not load elite memory: {e}")
            past_elites = []

    for c_id, center in enumerate(cluster_centers):
        cluster_txt = render_bar("Cluster", c_id+1, NUM_CLUSTERS, extra=f"(C {c_id+1}/{NUM_CLUSTERS})")
        for gen in range(GENS_PER_CLUSTER):
            gen_txt = render_bar("Gen", gen+1, GENS_PER_CLUSTER, extra=f"Cluster {c_id+1}")
            main_txt = render_bar("Clustering", c_id*GENS_PER_CLUSTER*POP_CLUSTER+gen*POP_CLUSTER, NUM_CLUSTERS*GENS_PER_CLUSTER*POP_CLUSTER)
            subsub_txt = render_bar("Indiv", 0, POP_CLUSTER, show_perc=False, extra=f"Gen {gen+1}")
            subsubsub_txt = render_bar("Eval", 0, POP_CLUSTER, show_perc=True, show_counter=True, extra=f"Indivs in Gen")
            update_all_bars(main_txt, cluster_txt, gen_txt, subsub_txt, subsubsub_txt)
            population = [mutate(center.copy(), 1.0) for _ in range(POP_CLUSTER)] if gen==0 else population
            best_local = (float("inf"), None)
            stagnant = 0
            last_best = float("inf")
            results = []
            for i, ind in enumerate(population):
                subsub_txt = render_bar("Indiv", i+1, POP_CLUSTER, show_perc=False, extra=f"Gen {gen+1}")
                subsubsub_txt = render_bar("Eval", i+1, POP_CLUSTER, show_perc=True, show_counter=True, extra=f"Indivs in Gen")
                update_all_bars(main_txt, cluster_txt, gen_txt, subsub_txt, subsubsub_txt)
                fit, feasible, violating, mass, RFs = fitness(
                    ind, log_label=f"Cluster{c_id+1}_Gen{gen+1}", log_step=i, call_counter=call_counter, progress=progress
                )
                results.append((fit, ind, mass, feasible, violating))
                if feasible and fit < best_local[0]: best_local = (fit, ind)
            results.sort(key=lambda x: x[0])
            if abs(results[0][0] - last_best) < 0.001: stagnant += 1
            else: stagnant = 0
            last_best = results[0][0]
            if stagnant >= 3: break
            new_population = [results[0][1]]
            parents = [ind for (_, ind, *_ ) in results[:POP_CLUSTER//2]]
            while len(new_population) < POP_CLUSTER:
                idxs = np.random.choice(len(parents), 2, replace=False)
                p1, p2 = parents[idxs[0]], parents[idxs[1]]
                child = crossover(p1, p2)
                child = mutate(child, MUTATION_RATE)
                new_population.append(child)
            population = new_population
        best_cluster_results.extend(sorted(results, key=lambda x: x[0])[:ELITE_PER_CLUSTER])
    update_all_bars(render_bar("Clustering", NUM_CLUSTERS*GENS_PER_CLUSTER*POP_CLUSTER, NUM_CLUSTERS*GENS_PER_CLUSTER*POP_CLUSTER))
    time.sleep(1.2)

    new_elites = [ind for (_, ind, feasible, *_ ) in best_cluster_results if feasible]
    unique_new_elites = []
    for ind in new_elites:
        if is_unique(ind, past_elites + unique_new_elites): unique_new_elites.append(ind)
    past_elites.extend(unique_new_elites)
    past_elites = sorted(past_elites, key=lambda ind: fitness(ind, log_label="MemorySort", progress=progress)[0])[:MAX_MEMORY_SIZE]
    try:
        with open(MEMORY_FILE, "wb") as f: pickle.dump(past_elites, f)
    except Exception as e: logging.error(f"Failed to write elite memory: {e}")

    # ==== PHASE 3 ====
    clear_console()
    print("="*65)
    print("PHASE 3: MICRO-ENSEMBLE REFINEMENT")
    print("="*65)
    print("Refining top designs with micro-GAs for best and most robust solution.\n")
    print('\n'*5)

    def refine_phase(starting_population, tag, ensemble_num):
        population = starting_population.copy()
        best_local = (float("inf"), None)
        stagnant = 0
        last_best = float("inf")
        for gen in range(REFINE_GENS):
            gen_txt = render_bar("Gen", gen+1, REFINE_GENS, extra=f"Ensemble {ensemble_num}")
            main_txt = render_bar("Refining", (ensemble_num-1)*REFINE_GENS*REFINE_POP+gen*REFINE_POP, 3*REFINE_GENS*REFINE_POP)
            subsub_txt = render_bar("Indiv", 0, REFINE_POP, show_perc=False, extra=f"Gen {gen+1}")
            subsubsub_txt = render_bar("Eval", 0, REFINE_POP, show_perc=True, show_counter=True, extra=f"Indivs in Gen")
            update_all_bars(main_txt, render_bar("Ensemble", ensemble_num, 3), gen_txt, subsub_txt, subsubsub_txt)
            results = []
            for i, ind in enumerate(population):
                subsub_txt = render_bar("Indiv", i+1, REFINE_POP, show_perc=False, extra=f"Gen {gen+1}")
                subsubsub_txt = render_bar("Eval", i+1, REFINE_POP, show_perc=True, show_counter=True, extra=f"Indivs in Gen")
                update_all_bars(main_txt, render_bar("Ensemble", ensemble_num, 3), gen_txt, subsub_txt, subsubsub_txt)
                fit, feasible, violating, mass, RFs = fitness(
                    ind, penalty_scale=1000 + 50 * gen, log_label=tag, log_step=gen, call_counter=call_counter, progress=progress
                )
                results.append((fit, ind, mass, feasible, violating))
                if feasible and fit < best_local[0]: best_local = (fit, ind)
            results.sort(key=lambda x: x[0])
            if abs(results[0][0] - last_best) < 1e-4: stagnant += 1
            else: stagnant = 0
            last_best = results[0][0]
            if stagnant >= 5: break
            new_population = [results[0][1]]
            parents = [ind for (_, ind, *_ ) in results[:len(results)//2]]
            for i in range(1, REFINE_POP):
                idxs = np.random.choice(len(parents), 2, replace=False)
                p1, p2 = parents[idxs[0]], parents[idxs[1]]
                child = crossover(p1, p2)
                adaptive_rate = MUTATION_RATE * (i / REFINE_POP)
                child = mutate(child, adaptive_rate)
                new_population.append(child)
            population = new_population
        return best_local

    refine_pool = [ind for (_, ind, feasible, *_ ) in sorted(best_cluster_results, key=lambda x: x[0]) if feasible][:REFINE_POP]
    memory_additions = 0
    for ind in past_elites:
        if is_unique(ind, refine_pool) and memory_additions < 5:
            refine_pool.append(mutate(ind.copy(), 0.1)); memory_additions += 1

    ensemble_results = []
    for i in range(3):
        ensemble_results.append(refine_phase(refine_pool[:REFINE_POP], f"Ensemble{i+1}", i+1))
    update_all_bars(render_bar("TOTAL", progress['total_calls'], progress['total_calls'], extra="ETA: 00:00:00"))
    time.sleep(1.2)

    best_global = sorted(ensemble_results, key=lambda x: x[0])[0]

    # Final output section: clear and print only final result
    clear_console()
    print("\n" + "="*65)
    print("[RESULT] FINAL BEST DESIGN - Standard variation of ~0.5 expected")
    print("="*65)
    print(f"Best found mass: {best_global[0]:.6f}")
    print(f"Design parameters (x): {np.round(best_global[1], 6)}\n")

if __name__ == "__main__":
    main()
