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
from scipy.spatial.distance import mahalanobis
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
MAX_GMM_K = 6
MAHALANOBIS_PERCENTILE = 97.5
PERCENTILE_LO = 2.5
PERCENTILE_HI = 97.5
MIN_SAMPLE_FOR_BIMODAL = 50


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
    if 55 <= delta <= 110:
        return "C3'-endo"
    elif 120 <= delta <= 165:
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


def fit_gmm_with_bic(X: np.ndarray, max_k: int = MAX_GMM_K) -> GaussianMixture:
    """
    Fit GMM with k=1..max_k, select best k by BIC (lower is better).
    Returns fitted GaussianMixture model.
    """
    best_gmm = None
    best_bic = np.inf
    for k in range(1, max_k + 1):
        try:
            gmm = GaussianMixture(
                n_components=k,
                covariance_type="full",
                random_state=42,
                max_iter=200,
                n_init=3,
            )
            gmm.fit(X)
            bic = gmm.bic(X)
            if bic < best_bic:
                best_bic = bic
                best_gmm = gmm
        except Exception:
            continue
    return best_gmm


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
        dists = np.array([
            mahalanobis(p, mean_k, cov_inv) for p in points
        ])
        if len(dists) > 0:
            thr = float(np.percentile(dists, MAHALANOBIS_PERCENTILE))
        else:
            thr = 0.0
        thresholds.append(thr)
    return thresholds


def serialize_gmm_clusters(
    gmm: GaussianMixture, thresholds: List[float], n_points: int
) -> Dict:
    """Convert GMM parameters to JSON-serializable dict."""
    clusters = []
    for k in range(gmm.n_components):
        clusters.append({
            "weight": round(float(gmm.weights_[k]), 6),
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

def _fit_abg_clusters(entries: List[Dict], pucker_split: bool = True) -> Dict:
    """
    Fit GMM clusters for (alpha, beta, gamma) angles.
    If pucker_split is True, split by sugar pucker first.
    """
    result = {}

    # Collect entries that have all three angles
    valid = [
        e for e in entries
        if all(k in e for k in ("alpha", "beta", "gamma"))
    ]

    if not valid:
        return result

    if pucker_split:
        c3_endo = [e for e in valid if get_sugar_pucker(e.get("delta")) == "C3'-endo"]
        c2_endo = [e for e in valid if get_sugar_pucker(e.get("delta")) == "C2'-endo"]

        for label, subset in [("abg_c3_endo", c3_endo), ("abg_c2_endo", c2_endo)]:
            if len(subset) >= MIN_PUCKER_GROUP:
                angles = np.array([[e["alpha"], e["beta"], e["gamma"]] for e in subset])
                X = angles_to_sincos(angles)
                gmm = fit_gmm_with_bic(X)
                if gmm is not None:
                    thresholds = compute_mahalanobis_threshold(gmm, X)
                    result[label] = serialize_gmm_clusters(gmm, thresholds, len(subset))

    # If either pucker group was too small or pucker_split is off, do combined
    if "abg_c3_endo" not in result or "abg_c2_endo" not in result:
        angles = np.array([[e["alpha"], e["beta"], e["gamma"]] for e in valid])
        X = angles_to_sincos(angles)
        gmm = fit_gmm_with_bic(X)
        if gmm is not None:
            thresholds = compute_mahalanobis_threshold(gmm, X)
            result["abg_combined"] = serialize_gmm_clusters(gmm, thresholds, len(valid))

    return result


def _fit_dez_clusters(entries: List[Dict]) -> Optional[Dict]:
    """
    Fit GMM clusters for (delta, epsilon, zeta) angles.
    No pucker split — delta IS the pucker, clusters emerge naturally.
    """
    valid = [
        e for e in entries
        if all(k in e for k in ("delta", "epsilon", "zeta"))
    ]
    if not valid:
        return None

    angles = np.array([[e["delta"], e["epsilon"], e["zeta"]] for e in valid])
    X = angles_to_sincos(angles)
    gmm = fit_gmm_with_bic(X)
    if gmm is None:
        return None
    thresholds = compute_mahalanobis_threshold(gmm, X)
    return serialize_gmm_clusters(gmm, thresholds, len(valid))


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
    Orchestrator: fit clusters per edge type, handle _OTHER merging for
    edge types with fewer than MIN_EDGE_COUNT observations.
    """
    edge_results = {}
    other_entries = []

    # Sort edge types for deterministic output
    for lw in sorted(data.keys()):
        edge_count = counts.get(lw, 0)

        if edge_count < MIN_EDGE_COUNT or "." in lw or lw == "--":
            # Merge into _OTHER: low-count edges, partial annotations (.), and unclassified (--)
            other_entries.extend(data[lw])
            continue

        print(f"  Fitting clusters for {lw} (n={edge_count:,})...")
        entries = data[lw]
        entry_result = {}

        abg = _fit_abg_clusters(entries, pucker_split=True)
        entry_result.update(abg)

        dez = _fit_dez_clusters(entries)
        if dez is not None:
            entry_result["dez"] = dez

        chi = _fit_chi(entries)
        if chi is not None:
            entry_result["chi"] = chi

        if entry_result:
            edge_results[lw] = entry_result

    # Handle _OTHER bucket
    if other_entries:
        print(f"  Fitting clusters for _OTHER (n={len(other_entries):,} residue observations)...")
        other_result = {}

        abg = _fit_abg_clusters(other_entries, pucker_split=True)
        other_result.update(abg)

        dez = _fit_dez_clusters(other_entries)
        if dez is not None:
            other_result["dez"] = dez

        chi = _fit_chi(other_entries)
        if chi is not None:
            other_result["chi"] = chi

        if other_result:
            edge_results["_OTHER"] = other_result

    return edge_results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("BACKBONE TORSION CLUSTER GENERATION")
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

    print("\nFitting GMM clusters...")
    edge_results = build_cluster_results(data, counts)

    # Build output
    output = {
        "metadata": {
            "description": "Backbone torsion cluster parameters by edge type",
            "sin_cos_encoding": "Angles transformed to [sin,cos,...] before clustering",
            "mahalanobis_percentile": MAHALANOBIS_PERCENTILE,
            "min_edge_count": MIN_EDGE_COUNT,
            "min_cluster_group": MIN_PUCKER_GROUP,
        },
        "edge_types": edge_results,
    }

    print(f"\nWriting results to {OUTPUT_JSON}...")
    with open(OUTPUT_JSON, "w") as f:
        json.dump(output, f, indent=2)

    # Summary
    n_edges = len(edge_results)
    print(f"\nDone. {n_edges} edge type entries (including _OTHER).")
    for lw, result in sorted(edge_results.items()):
        parts = []
        for key in ("abg_c3_endo", "abg_c2_endo", "abg_combined", "dez", "chi"):
            if key in result:
                if key == "chi":
                    parts.append(f"chi({result[key]['modality']})")
                elif key.startswith("abg"):
                    parts.append(f"{key}(k={result[key]['n_components']})")
                else:
                    parts.append(f"{key}(k={result[key]['n_components']})")
        print(f"  {lw}: {', '.join(parts)}")
    print(f"\nOutput: {OUTPUT_JSON}")


if __name__ == "__main__":
    main()


# python3 generate_backbone_clusters.py
