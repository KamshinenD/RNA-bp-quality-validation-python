#!/usr/bin/env python3
"""
Generate backbone torsion cluster parameters by edge type (Leontis-Westhof).

Fits GMM clusters for (alpha, beta, gamma) and (delta, epsilon, zeta) angle groups,
and computes chi angle ranges. Results are saved to data/backbone_clusters.json
for use by the scorer to penalize unusual backbone conformations.

Uses data from:
  - data/uniqueRNAS.csv (unique PDB IDs to process)
  - data/basepairs/ (base pair geometry + bp_type, lw)
  - data/torsions/ (per-residue torsion angles)

Output: data/backbone_clusters.json
"""

import json
import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
# scipy.spatial.distance.mahalanobis replaced by vectorized _vectorized_mahalanobis()
from g_quads import g_quads

# Paths
DATA_DIR = Path(__file__).resolve().parent / "data"
UNIQUE_RNAS_CSV = DATA_DIR / "uniqueRNAS.csv"
BASEPAIR_DIR = DATA_DIR / "basepairs"
TORSION_DIR = DATA_DIR / "torsions"
OUTPUT_JSON = DATA_DIR / "backbone_clusters.json"

# Thresholds
MIN_EDGE_COUNT = 1000
MIN_PUCKER_GROUP = 100
MIN_CHI_COUNT = 50
MAX_GMM_K = 12
MAHALANOBIS_PERCENTILE = 97.5
PERCENTILE_LO = 2.5
PERCENTILE_HI = 97.5
MIN_SAMPLE_FOR_BIMODAL = 50
# Merge GMM components whose means are this close in 6D sin/cos space (Euclidean).
# Ensures "same place" clusters are identified as one unique cluster.
MERGE_DISTANCE_THRESHOLD = 0.35


# ---------------------------------------------------------------------------
# Reused utility functions (copied from existing scripts per project pattern)
# ---------------------------------------------------------------------------

def normalize_bp_type(bp_type: str) -> str:
    """Normalize bp_type so A-U == U-A, G-C == C-G, etc."""
    if not bp_type or "-" not in bp_type:
        return bp_type or "unknown"
    parts = bp_type.split("-")
    if len(parts) != 2:
        return bp_type
    return "-".join(sorted(parts))


def load_unique_pdb_ids() -> set:
    """Load unique PDB IDs from data/uniqueRNAS.csv."""
    if not UNIQUE_RNAS_CSV.exists():
        raise FileNotFoundError(f"Unique RNAs CSV not found: {UNIQUE_RNAS_CSV}")
    df = pd.read_csv(UNIQUE_RNAS_CSV)
    col = "unique_pdb_id"
    if col not in df.columns:
        raise ValueError(f"Column '{col}' not found in {UNIQUE_RNAS_CSV}")
    pdb_ids = set(df[col].astype(str).str.strip().str.upper().dropna().unique())
    return pdb_ids


def get_excluded_pdb_ids() -> set:
    """PDB IDs to exclude (e.g., G-quadruplex structures)."""
    return {p.upper() for p in g_quads}


def load_basepairs(pdb_id: str) -> List[Dict]:
    """Load base pair data for a PDB."""
    fp = BASEPAIR_DIR / f"{pdb_id}.json"
    if not fp.exists():
        return []
    with open(fp, "r") as f:
        data = json.load(f)
    if isinstance(data, list):
        return data
    return data.get("base_pairs", [])


def load_torsions(pdb_id: str) -> Dict:
    """Load torsion data for a PDB. Keys are residue IDs (e.g. 'AN1-A-8-')."""
    fp = TORSION_DIR / f"{pdb_id}.json"
    if not fp.exists():
        return {}
    with open(fp, "r") as f:
        return json.load(f)


def is_adjacent_pair(res_1: str, res_2: str) -> bool:
    """
    Return True if the two residues are adjacent in sequence (same chain,
    residue numbers differ by 1).
    """
    if not res_1 or not res_2:
        return False
    parts1 = res_1.split('-')
    parts2 = res_2.split('-')
    if len(parts1) < 3 or len(parts2) < 3:
        return False
    try:
        chain1, num1 = parts1[0], int(parts1[2])
        chain2, num2 = parts2[0], int(parts2[2])
    except (ValueError, IndexError):
        return False
    if chain1 != chain2:
        return False
    return abs(num1 - num2) == 1


def get_sugar_pucker(delta: float) -> str:
    """Classify sugar pucker based on delta torsion angle."""
    if delta is None:
        return "unknown"
    if 60 <= delta <= 110:
        return "C3'-endo"
    elif 125 <= delta <= 165:
        return "C2'-endo"
    else:
        return "intermediate"


def _is_bimodal(arr: np.ndarray) -> bool:
    """Use GMM to test if distribution is bimodal. Returns True if bimodal."""
    if len(arr) < MIN_SAMPLE_FOR_BIMODAL:
        return False
    try:
        X = arr.reshape(-1, 1)
        gmm1 = GaussianMixture(n_components=1, random_state=42, max_iter=100)
        gmm2 = GaussianMixture(n_components=2, random_state=42, max_iter=100)
        gmm1.fit(X)
        gmm2.fit(X)
        bic1 = gmm1.bic(X)
        bic2 = gmm2.bic(X)
        if bic2 >= bic1:
            return False
        means = np.sort(gmm2.means_.ravel())
        sep = abs(means[1] - means[0])
        if sep > 180:
            sep = 360 - sep
        if sep < 60:
            return False
        weights = gmm2.weights_
        if min(weights) < 0.05:
            return False
        return True
    except Exception:
        return False


def compute_chi_stats(values: List[float]) -> Optional[Dict]:
    """
    Compute allowed ranges for chi: unimodal -> 2.5th-97.5th percentile;
    bimodal -> GMM clustering, 2.5th-97.5th within each cluster.
    """
    arr = np.array([float(x) for x in values if x is not None], dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) < 2:
        return None
    count = int(len(arr))
    bimodal = _is_bimodal(arr)
    ranges = []
    if bimodal:
        try:
            X = arr.reshape(-1, 1)
            gmm = GaussianMixture(n_components=2, random_state=42, max_iter=100)
            labels = gmm.fit_predict(X)
            for c in range(2):
                sub = arr[labels == c]
                if len(sub) >= 5:
                    lo = float(np.percentile(sub, PERCENTILE_LO))
                    hi = float(np.percentile(sub, PERCENTILE_HI))
                    ranges.append([round(lo, 2), round(hi, 2)])
            if len(ranges) < 2:
                bimodal = False
        except Exception:
            bimodal = False
    if not bimodal or not ranges:
        lo = float(np.percentile(arr, PERCENTILE_LO))
        hi = float(np.percentile(arr, PERCENTILE_HI))
        ranges = [[round(lo, 2), round(hi, 2)]]
        bimodal = False
    return {
        "ranges": ranges,
        "count": count,
        "modality": "bimodal" if bimodal else "unimodal",
    }


# ---------------------------------------------------------------------------
# New functions for 6D sin/cos GMM clustering
# ---------------------------------------------------------------------------

def angles_to_sincos(angle_tuples: np.ndarray) -> np.ndarray:
    """
    Convert Nx3 array of angle tuples (degrees) to Nx6 sin/cos space.
    Input columns: (a1, a2, a3)
    Output columns: (sin_a1, cos_a1, sin_a2, cos_a2, sin_a3, cos_a3)
    """
    rads = np.radians(angle_tuples)
    n = rads.shape[0]
    result = np.empty((n, 6), dtype=float)
    for i in range(3):
        result[:, 2 * i] = np.sin(rads[:, i])
        result[:, 2 * i + 1] = np.cos(rads[:, i])
    return result


BIC_SUBSAMPLE = 50_000  # subsample for BIC model selection (picking k)


def fit_gmm_with_bic(X: np.ndarray, max_k: int = MAX_GMM_K) -> GaussianMixture:
    """
    Fit GMM with k=1..max_k, select best k by BIC (lower is better).
    BIC selection uses a subsample with n_init=1 (fast).
    Final refit uses ALL data with n_init=1 (stable since k is already chosen).
    """
    rng = np.random.RandomState(42)

    # Subsample for BIC selection
    if len(X) > BIC_SUBSAMPLE:
        X_bic = X[rng.choice(len(X), BIC_SUBSAMPLE, replace=False)]
    else:
        X_bic = X

    best_k = 1
    best_bic = np.inf
    for k in range(1, max_k + 1):
        try:
            gmm = GaussianMixture(
                n_components=k,
                covariance_type="full",
                random_state=42,
                max_iter=200,
                n_init=1,  # single init for BIC search (fast)
            )
            gmm.fit(X_bic)
            bic = gmm.bic(X_bic)
            if bic < best_bic:
                best_bic = bic
                best_k = k
        except Exception:
            continue

    # Refit winning k on ALL data
    try:
        gmm = GaussianMixture(
            n_components=best_k,
            covariance_type="full",
            random_state=42,
            max_iter=200,
            n_init=1,  # single init on full data (k already chosen)
        )
        gmm.fit(X)
        return gmm
    except Exception:
        return None


def _vectorized_mahalanobis(points: np.ndarray, mean: np.ndarray, cov_inv: np.ndarray) -> np.ndarray:
    """Compute Mahalanobis distances for all points at once (vectorized)."""
    diff = points - mean  # (N, D)
    # Mahalanobis: sqrt(diff @ cov_inv @ diff.T) per row
    left = diff @ cov_inv  # (N, D)
    dists_sq = np.sum(left * diff, axis=1)  # (N,)
    return np.sqrt(np.maximum(dists_sq, 0.0))


def compute_mahalanobis_threshold(
    gmm: GaussianMixture, X: np.ndarray
) -> List[float]:
    """
    For each GMM component, compute Mahalanobis distances of assigned points
    and return the 97.5th percentile as the threshold.
    """
    labels = gmm.predict(X)
    thresholds = []
    for k in range(gmm.n_components):
        mask = labels == k
        points = X[mask]
        mean_k = gmm.means_[k]
        cov_k = gmm.covariances_[k]
        try:
            cov_inv = np.linalg.inv(cov_k)
        except np.linalg.LinAlgError:
            cov_inv = np.linalg.pinv(cov_k)
        if len(points) > 0:
            dists = _vectorized_mahalanobis(points, mean_k, cov_inv)
            thr = float(np.percentile(dists, MAHALANOBIS_PERCENTILE))
        else:
            thr = 0.0
        thresholds.append(thr)
    return thresholds


def _merge_two_components(
    mean_i: np.ndarray,
    cov_i: np.ndarray,
    w_i: float,
    thr_i: float,
    mean_j: np.ndarray,
    cov_j: np.ndarray,
    w_j: float,
    thr_j: float,
) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """
    Merge two GMM components into one (weighted mean, merged covariance, max threshold).
    """
    w_tot = w_i + w_j
    mean_new = (w_i * mean_i + w_j * mean_j) / w_tot
    diff_i = mean_i - mean_new
    diff_j = mean_j - mean_new
    cov_new = (
        w_i * (cov_i + np.outer(diff_i, diff_i))
        + w_j * (cov_j + np.outer(diff_j, diff_j))
    ) / w_tot
    thr_new = max(thr_i, thr_j)
    return mean_new, cov_new, w_tot, thr_new


def _merge_duplicate_clusters(
    gmm: GaussianMixture,
    thresholds: List[float],
    distance_threshold: float = MERGE_DISTANCE_THRESHOLD,
) -> Tuple[GaussianMixture, List[float]]:
    """
    Merge GMM components whose means are within distance_threshold in sin/cos space,
    so that clusters in the "same place" are identified as one unique cluster.
    Returns a new GMM and updated thresholds (same length as n_components).
    """
    n = gmm.n_components
    means = [gmm.means_[k].copy() for k in range(n)]
    covs = [gmm.covariances_[k].copy() for k in range(n)]
    weights = list(gmm.weights_.copy())
    thrs = list(thresholds)

    while True:
        best_i, best_j, best_d = None, None, float("inf")
        for i in range(len(means)):
            for j in range(i + 1, len(means)):
                d = float(np.linalg.norm(means[i] - means[j]))
                if d < distance_threshold and d < best_d:
                    best_d, best_i, best_j = d, i, j
        if best_i is None:
            break

        mean_new, cov_new, w_new, thr_new = _merge_two_components(
            means[best_i], covs[best_i], weights[best_i], thrs[best_i],
            means[best_j], covs[best_j], weights[best_j], thrs[best_j],
        )
        means[best_i] = mean_new
        covs[best_i] = cov_new
        weights[best_i] = w_new
        thrs[best_i] = thr_new
        # Remove component best_j (higher index first so indices stay valid)
        del means[best_j], covs[best_j], weights[best_j], thrs[best_j]

    if len(means) == 0:
        return gmm, thresholds

    K = len(means)
    weights_arr = np.array(weights, dtype=float)
    weights_arr /= weights_arr.sum()

    new_gmm = GaussianMixture(
        n_components=K,
        covariance_type="full",
        random_state=42,
    )
    new_gmm.means_ = np.array(means)
    new_gmm.covariances_ = np.array(covs)
    new_gmm.weights_ = weights_arr
    # precisions_cholesky_ required for predict(); Cholesky of precision = inv(cov).
    # Robustify with jitter in case merged covariances are near-singular.
    precisions_chol = []
    n_features = new_gmm.means_.shape[1]
    eye = np.eye(n_features)
    for c in covs:
        c = np.asarray(c, dtype=float)
        # try increasing jitter until Cholesky succeeds
        jitter = 1e-8
        L = None
        for _ in range(8):
            try:
                cov_j = c + jitter * eye
                try:
                    p = np.linalg.inv(cov_j)
                except np.linalg.LinAlgError:
                    p = np.linalg.pinv(cov_j)
                L = np.linalg.cholesky(p)
                break
            except np.linalg.LinAlgError:
                jitter *= 10.0
        if L is None:
            # Last-resort fallback: identity precision (very permissive / avoids crash)
            L = np.linalg.cholesky(eye)
        precisions_chol.append(L)
    new_gmm.precisions_cholesky_ = np.array(precisions_chol)
    # Mark as "fitted enough" for sklearn checks
    new_gmm.converged_ = True
    new_gmm.n_iter_ = 0
    new_gmm.lower_bound_ = -np.inf
    new_gmm.n_features_in_ = n_features
    return new_gmm, thrs


def _serialize_global_clusters(
    gmm: GaussianMixture, thresholds: List[float], n_points: int
) -> Dict:
    """Convert GMM parameters to JSON-serializable dict (no weights — those are edge-specific)."""
    clusters = []
    for k in range(gmm.n_components):
        clusters.append({
            "mean": [round(float(v), 6) for v in gmm.means_[k]],
            "covariance": [
                [round(float(v), 8) for v in row]
                for row in gmm.covariances_[k]
            ],
            "mahalanobis_threshold": round(thresholds[k], 4),
        })
    return {
        "n_components": gmm.n_components,
        "n_points": n_points,
        "clusters": clusters,
    }


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def collect_backbone_data(
    pdb_ids: set, exclude: set
) -> Tuple[Dict[str, List], Dict[str, int]]:
    """
    Main data collection loop. For each non-adjacent base-pair, collect
    torsion angles grouped by edge type (lw).

    Returns:
        data: dict mapping edge_type -> list of per-residue torsion dicts
              Each dict has keys: alpha, beta, gamma, delta, epsilon, zeta, chi
        counts: dict mapping edge_type -> total observations (for _OTHER logic)
    """
    data = defaultdict(list)
    counts = defaultdict(int)

    working_pdb_ids = sorted(pdb_ids - exclude)
    processed = 0
    skipped_no_bp = 0
    skipped_no_torsion = 0

    for pdb_id in working_pdb_ids:
        basepairs = load_basepairs(pdb_id)
        torsions = load_torsions(pdb_id)

        if not basepairs:
            skipped_no_bp += 1
            continue
        if not torsions:
            skipped_no_torsion += 1
            continue

        for bp in basepairs:
            res_1 = bp.get("res_1")
            res_2 = bp.get("res_2")
            lw = bp.get("lw", "unknown")

            if not res_1 or not res_2 or not lw or lw == "unknown":
                continue
            if is_adjacent_pair(res_1, res_2):
                continue

            counts[lw] += 1

            for res_id in (res_1, res_2):
                t = torsions.get(res_id, {})
                if not t:
                    continue

                entry = {}
                has_any = False
                for angle in ("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"):
                    val = t.get(angle)
                    if val is not None:
                        try:
                            v = float(val)
                            if np.isfinite(v):
                                entry[angle] = v
                                has_any = True
                        except (ValueError, TypeError):
                            pass

                if has_any:
                    data[lw].append(entry)

        processed += 1
        if processed % 500 == 0:
            print(f"  Processed {processed}/{len(working_pdb_ids)} PDBs...")

    print(f"Processed {processed} PDBs with both basepairs and torsions")
    print(f"Skipped {skipped_no_bp} (no basepairs), {skipped_no_torsion} (no torsions)")
    return dict(data), dict(counts)


# ---------------------------------------------------------------------------
# Cluster fitting orchestration
# ---------------------------------------------------------------------------

def _pool_angle_data(
    all_data: Dict[str, List], angle_group: str
) -> Dict[str, np.ndarray]:
    """
    Pool angle triples across ALL edge types for a given angle group.

    Args:
        all_data: dict mapping edge_type -> list of torsion dicts
        angle_group: "abg" or "dez"

    Returns:
        For "abg": {"c3_endo": Nx3 array, "c2_endo": Nx3 array}
        For "dez": {"combined": Nx3 array}
    """
    if angle_group == "abg":
        required = ("alpha", "beta", "gamma")
        c3_rows, c2_rows = [], []
        for entries in all_data.values():
            for e in entries:
                if all(k in e for k in required):
                    pucker = get_sugar_pucker(e.get("delta"))
                    triple = [e["alpha"], e["beta"], e["gamma"]]
                    if pucker == "C3'-endo":
                        c3_rows.append(triple)
                    elif pucker == "C2'-endo":
                        c2_rows.append(triple)
        result = {}
        if c3_rows:
            result["c3_endo"] = np.array(c3_rows)
        if c2_rows:
            result["c2_endo"] = np.array(c2_rows)
        return result
    else:  # dez
        required = ("delta", "epsilon", "zeta")
        rows = []
        for entries in all_data.values():
            for e in entries:
                if all(k in e for k in required):
                    rows.append([e["delta"], e["epsilon"], e["zeta"]])
        if rows:
            return {"combined": np.array(rows)}
        return {}


def _fit_pooled_clusters(
    pooled_angles: np.ndarray,
) -> Tuple[Optional[GaussianMixture], Optional[np.ndarray], Optional[List[float]]]:
    """
    Fit a single GMM on pooled angle data (across all edge types).
    Merges components whose means fall in the same place (see _merge_duplicate_clusters)
    so that unique conformational clusters are identified.

    Args:
        pooled_angles: Nx3 angle array

    Returns:
        (fitted GMM, sin/cos data X, global Mahalanobis thresholds)
        or (None, None, None) if fitting fails
    """
    X = angles_to_sincos(pooled_angles)
    gmm = fit_gmm_with_bic(X)
    if gmm is None:
        return None, None, None
    thresholds = compute_mahalanobis_threshold(gmm, X)
    k_before = gmm.n_components
    gmm, thresholds = _merge_duplicate_clusters(gmm, thresholds, MERGE_DISTANCE_THRESHOLD)
    k_after = gmm.n_components
    if k_after < k_before:
        print(f"    Merged duplicate clusters: k {k_before} -> {k_after} (unique)")
    return gmm, X, thresholds


def _compute_edge_projection(
    gmm: GaussianMixture,
    edge_entries: List[Dict],
    angle_group: str,
) -> Optional[Dict]:
    """
    Project an edge type's data onto the global GMM clusters.

    Computes weight vector (fraction of this edge's data per cluster)
    and edge-specific Mahalanobis thresholds.

    Args:
        gmm: fitted global GMM
        edge_entries: list of torsion dicts for this edge type
        angle_group: one of "abg_c3_endo", "abg_c2_endo", "dez"

    Returns:
        {"n_points": N, "weights": [...], "mahalanobis_thresholds": [...]}
        or None if no valid data
    """
    # Extract the right angle triple based on angle_group
    if angle_group.startswith("abg"):
        required = ("alpha", "beta", "gamma")
        target_pucker = "C3'-endo" if "c3" in angle_group else "C2'-endo"
        valid = [
            e for e in edge_entries
            if all(k in e for k in required)
            and get_sugar_pucker(e.get("delta")) == target_pucker
        ]
        if not valid:
            return None
        angles = np.array([[e["alpha"], e["beta"], e["gamma"]] for e in valid])
    else:  # dez
        required = ("delta", "epsilon", "zeta")
        valid = [e for e in edge_entries if all(k in e for k in required)]
        if not valid:
            return None
        angles = np.array([[e["delta"], e["epsilon"], e["zeta"]] for e in valid])

    X = angles_to_sincos(angles)
    labels = gmm.predict(X)
    n_components = gmm.n_components
    n_points = len(valid)

    # Weight vector: fraction of points in each cluster
    weights = []
    for k in range(n_components):
        weights.append(round(float(np.sum(labels == k)) / n_points, 6))

    # Edge-specific Mahalanobis thresholds (vectorized)
    edge_thresholds = []
    for k in range(n_components):
        mask = labels == k
        points = X[mask]
        if len(points) == 0:
            edge_thresholds.append(0.0)
            continue
        mean_k = gmm.means_[k]
        cov_k = gmm.covariances_[k]
        try:
            cov_inv = np.linalg.inv(cov_k)
        except np.linalg.LinAlgError:
            cov_inv = np.linalg.pinv(cov_k)
        dists = _vectorized_mahalanobis(points, mean_k, cov_inv)
        thr = float(np.percentile(dists, MAHALANOBIS_PERCENTILE))
        edge_thresholds.append(round(thr, 4))

    return {
        "n_points": n_points,
        "weights": weights,
        "mahalanobis_thresholds": edge_thresholds,
    }


def _fit_chi(entries: List[Dict]) -> Optional[Dict]:
    """Compute chi angle ranges using bimodal GMM approach."""
    chi_vals = [e["chi"] for e in entries if "chi" in e]
    if len(chi_vals) < MIN_CHI_COUNT:
        return None
    return compute_chi_stats(chi_vals)


def build_cluster_results(
    data: Dict[str, List], counts: Dict[str, int]
) -> Dict:
    """
    Orchestrator: pool all data across edge types, fit global GMMs,
    then compute per-edge-type weight profiles.
    """
    # --- Step 1: Pool angle data across ALL edge types ---
    print("  Pooling angle data across all edge types...")
    abg_pooled = _pool_angle_data(data, "abg")
    dez_pooled = _pool_angle_data(data, "dez")

    # --- Step 2: Fit global GMMs ---
    global_clusters = {}
    global_gmms = {}  # keep fitted GMMs for projection step

    for label, angles in [
        ("abg_c3_endo", abg_pooled.get("c3_endo")),
        ("abg_c2_endo", abg_pooled.get("c2_endo")),
        ("dez", dez_pooled.get("combined")),
    ]:
        if angles is None or len(angles) < MIN_PUCKER_GROUP:
            print(f"  Skipping {label}: insufficient data")
            continue
        print(f"  Fitting pooled GMM for {label} (n={len(angles):,})...")
        gmm, X, thresholds = _fit_pooled_clusters(angles)
        if gmm is not None:
            global_clusters[label] = _serialize_global_clusters(gmm, thresholds, len(angles))
            global_gmms[label] = gmm
            print(f"    -> k={gmm.n_components}")

    # --- Step 3: Determine qualifying edge types ---
    qualifying = []
    other_entries = []

    for lw in sorted(data.keys()):
        edge_count = counts.get(lw, 0)
        if edge_count < MIN_EDGE_COUNT or "." in lw or lw == "--":
            other_entries.extend(data[lw])
        else:
            qualifying.append(lw)

    # --- Step 4: Compute per-edge-type projections ---
    edge_results = {}

    for lw in qualifying:
        print(f"  Computing projections for {lw}...")
        entries = data[lw]
        entry_result = {}

        for group_label in ("abg_c3_endo", "abg_c2_endo", "dez"):
            if group_label not in global_gmms:
                continue
            proj = _compute_edge_projection(
                global_gmms[group_label], entries, group_label
            )
            if proj is not None:
                entry_result[group_label] = proj

        chi = _fit_chi(entries)
        if chi is not None:
            entry_result["chi"] = chi

        if entry_result:
            edge_results[lw] = entry_result

    # --- Step 5: Handle _OTHER bucket ---
    if other_entries:
        print(f"  Computing projections for _OTHER (n={len(other_entries):,})...")
        other_result = {}

        for group_label in ("abg_c3_endo", "abg_c2_endo", "dez"):
            if group_label not in global_gmms:
                continue
            proj = _compute_edge_projection(
                global_gmms[group_label], other_entries, group_label
            )
            if proj is not None:
                other_result[group_label] = proj

        chi = _fit_chi(other_entries)
        if chi is not None:
            other_result["chi"] = chi

        if other_result:
            edge_results["_OTHER"] = other_result

    return {"global_clusters": global_clusters, "edge_types": edge_results}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("BACKBONE TORSION CLUSTER GENERATION (POOLED)")
    print("=" * 70)

    print("\nLoading unique PDB IDs from data/uniqueRNAS.csv...")
    unique_pdb_ids = load_unique_pdb_ids()
    exclude = get_excluded_pdb_ids()
    print(f"  {len(unique_pdb_ids)} unique PDB IDs")
    print(f"  Excluding {len(exclude)} G-quad PDBs from g_quads.py")
    print(f"  Processing {len(unique_pdb_ids - exclude)} PDBs")

    print("\nCollecting backbone torsion data by edge type...")
    data, counts = collect_backbone_data(unique_pdb_ids, exclude)

    print(f"\nEdge type counts:")
    for lw in sorted(counts.keys(), key=lambda x: counts[x], reverse=True):
        marker = "" if counts[lw] >= MIN_EDGE_COUNT else " -> _OTHER"
        print(f"  {lw:6s}: {counts[lw]:>8,} base pairs{marker}")

    print("\nFitting pooled GMM clusters...")
    results = build_cluster_results(data, counts)

    # Build output
    output = {
        "metadata": {
            "description": "Pooled backbone torsion clusters with per-edge-type weight profiles",
            "sin_cos_encoding": "Angles transformed to [sin,cos,...] before clustering",
            "mahalanobis_percentile": MAHALANOBIS_PERCENTILE,
            "min_edge_count": MIN_EDGE_COUNT,
            "min_cluster_group": MIN_PUCKER_GROUP,
            "max_gmm_k": MAX_GMM_K,
            "merge_distance_threshold": MERGE_DISTANCE_THRESHOLD,
        },
        "global_clusters": results["global_clusters"],
        "edge_types": results["edge_types"],
    }

    print(f"\nWriting results to {OUTPUT_JSON}...")
    with open(OUTPUT_JSON, "w") as f:
        json.dump(output, f, indent=2)

    # Summary
    print(f"\nGlobal clusters:")
    for label, gc in results["global_clusters"].items():
        print(f"  {label}: k={gc['n_components']}, n={gc['n_points']:,}")

    n_edges = len(results["edge_types"])
    print(f"\n{n_edges} edge type entries (including _OTHER):")
    for lw in sorted(results["edge_types"].keys()):
        et = results["edge_types"][lw]
        parts = []
        for key in ("abg_c3_endo", "abg_c2_endo", "dez"):
            if key in et:
                n = et[key]["n_points"]
                top_w = max(et[key]["weights"])
                parts.append(f"{key}(n={n:,}, top_w={top_w:.2f})")
        if "chi" in et:
            parts.append(f"chi({et['chi']['modality']})")
        print(f"  {lw}: {', '.join(parts)}")
    print(f"\nOutput: {OUTPUT_JSON}")


if __name__ == "__main__":
    main()


# python3 generate_backbone_clusters.py
