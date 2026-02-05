#!/usr/bin/env python3
"""
Analyze backbone torsion angles by base-pair type and edge type (Leontis-Westhof).

Uses data from:
  - data/uniqueRNAS.csv (unique PDB IDs to process)
  - data/basepairs/ (base pair geometry + bp_type, lw)
  - data/torsions/ (per-residue torsion angles: alpha, beta, gamma, delta, epsilon, zeta, eta, chi, v1-v4)

Outputs: torsion_thresholds_analysis.json with allowed ranges per (bp_type, edge) per angle.
Unimodal: 2.5th-97.5th percentile. Bimodal: GMM clustering, 2.5th-97.5th within each cluster.
"""

import json
import numpy as np
import pandas as pd

try:
    from sklearn.mixture import GaussianMixture
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False
from pathlib import Path
from typing import Dict, List, Optional
from collections import defaultdict
from g_quads import g_quads

# Paths
DATA_DIR = Path(__file__).resolve().parent / "data"
UNIQUE_RNAS_CSV = DATA_DIR / "uniqueRNAS.csv"
BASEPAIR_DIR = DATA_DIR / "basepairs"
TORSION_DIR = DATA_DIR / "torsions"
OUTPUT_JSON = Path(__file__).resolve().parent / "torsion_thresholds_analysis.json"

# Torsion angles to analyze
TORSION_ANGLES = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "chi", "v1", "v2", "v3", "v4"]
MIN_EDGE_COUNT = 1000


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


def collect_torsions_by_bp_and_edge(
    unique_pdb_ids: set,
    exclude_pdb_ids: set,
) -> Dict[str, Dict[str, Dict[str, List[float]]]]:
    """
    For each (bp_type_norm, edge), collect lists of torsion values from base-pair residues.
    grouped[bp_type][edge][angle] = [values...]
    """
    grouped = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    counts = defaultdict(lambda: defaultdict(int))
    totals = defaultdict(int)

    # Exclude G-quadruplex structures (from g_quads.py)
    pdb_ids = unique_pdb_ids - exclude_pdb_ids
    pdb_ids = sorted(pdb_ids)

    processed = 0
    skipped_no_bp = 0
    skipped_no_torsion = 0
    total_bps = 0

    for pdb_id in pdb_ids:
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
            bp_type = bp.get("bp_type", "unknown")
            lw = bp.get("lw", "unknown")

            if not res_1 or not res_2 or not bp_type or not lw or lw == "unknown":
                continue

            bp_type_norm = normalize_bp_type(bp_type)
            totals[bp_type_norm] += 1
            counts[bp_type_norm][lw] += 1

            t1 = torsions.get(res_1, {}) if isinstance(torsions, dict) else {}
            t2 = torsions.get(res_2, {}) if isinstance(torsions, dict) else {}

            for angle in TORSION_ANGLES:
                for t in (t1, t2):
                    val = t.get(angle)
                    if val is not None:
                        try:
                            v = float(val)
                            if np.isfinite(v):
                                grouped[bp_type_norm][lw][angle].append(v)
                        except (ValueError, TypeError):
                            pass

        processed += 1
        total_bps += len(basepairs)
        if processed % 500 == 0:
            print(f"  Processed {processed} PDBs, {total_bps:,} base pairs...")

    print(f"Processed {processed} PDBs with both basepairs and torsions")
    print(f"Skipped {skipped_no_bp} PDBs (no basepairs), {skipped_no_torsion} PDBs (no torsions)")
    return dict(grouped), dict(counts), dict(totals)


MIN_SAMPLE_FOR_BIMODAL = 50
PERCENTILE_LO = 2.5
PERCENTILE_HI = 97.5


def _is_bimodal(arr: np.ndarray) -> bool:
    """Use GMM to test if distribution is bimodal. Returns True if bimodal."""
    if not HAS_SKLEARN or len(arr) < MIN_SAMPLE_FOR_BIMODAL:
        return False
    try:
        X = arr.reshape(-1, 1)
        gmm1 = GaussianMixture(n_components=1, random_state=42, max_iter=100)
        gmm2 = GaussianMixture(n_components=2, random_state=42, max_iter=100)
        gmm1.fit(X)
        gmm2.fit(X)
        # BIC: lower is better. If 2-component improves fit enough, treat as bimodal
        bic1 = gmm1.bic(X)
        bic2 = gmm2.bic(X)
        if bic2 >= bic1:
            return False
        # Check separation: means should be > 60° apart for angular data
        means = np.sort(gmm2.means_.ravel())
        sep = abs(means[1] - means[0])
        if sep > 180:
            sep = 360 - sep
        if sep < 60:
            return False
        # Both components should have non-negligible weight
        weights = gmm2.weights_
        if min(weights) < 0.05:
            return False
        return True
    except Exception:
        return False


def compute_stats(values: List[float]) -> Optional[Dict]:
    """
    Compute allowed ranges: unimodal -> 2.5th-97.5th percentile;
    bimodal -> GMM clustering, 2.5th-97.5th within each cluster.
    Returns dict with "ranges" (list of (min, max) tuples), "count", "modality".
    """
    arr = np.array([float(x) for x in values if x is not None], dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) < 2:
        return None
    count = int(len(arr))
    bimodal = _is_bimodal(arr)
    ranges = []
    if bimodal and HAS_SKLEARN:
        try:
            X = arr.reshape(-1, 1)
            gmm = GaussianMixture(n_components=2, random_state=42, max_iter=100)
            labels = gmm.fit_predict(X)
            for c in range(2):
                sub = arr[labels == c]
                if len(sub) >= 5:
                    lo = float(np.percentile(sub, PERCENTILE_LO))
                    hi = float(np.percentile(sub, PERCENTILE_HI))
                    ranges.append((round(lo, 2), round(hi, 2)))
            if len(ranges) < 2:
                bimodal = False
        except Exception:
            bimodal = False
    if not bimodal or not ranges:
        lo = float(np.percentile(arr, PERCENTILE_LO))
        hi = float(np.percentile(arr, PERCENTILE_HI))
        ranges = [(round(lo, 2), round(hi, 2))]
        bimodal = False
    return {
        "ranges": ranges,
        "count": count,
        "modality": "bimodal" if bimodal else "unimodal",
    }


def build_results(
    grouped: Dict,
    counts: Dict,
    totals: Dict,
) -> Dict:
    """Build output structure: edge -> bp_type -> angle -> stats (edge first, then base-pair)."""
    # First build bp_type -> edge -> stats (same logic as before)
    by_bp_then_edge = {}

    for bp_type, edges in grouped.items():
        by_bp_then_edge[bp_type] = {}
        total_bp = totals.get(bp_type, 0)
        other_bucket = defaultdict(list)

        for lw, angle_lists in edges.items():
            edge_count = counts.get(bp_type, {}).get(lw, 0)
            use_other = edge_count < MIN_EDGE_COUNT and total_bp >= MIN_EDGE_COUNT

            if use_other:
                for angle, vals in angle_lists.items():
                    other_bucket[angle].extend(vals)
            else:
                edge_stats = {}
                for angle, vals in angle_lists.items():
                    st = compute_stats(vals)
                    if st:
                        edge_stats[angle] = st
                if edge_stats:
                    by_bp_then_edge[bp_type][lw] = edge_stats

        if other_bucket:
            other_stats = {}
            for angle, vals in other_bucket.items():
                st = compute_stats(vals)
                if st:
                    other_stats[angle] = st
            if other_stats:
                by_bp_then_edge[bp_type]["_OTHER"] = other_stats

    # Restructure: edge -> bp_type -> angle -> stats (edge first, then base-pair)
    results = {}
    for bp_type, edges in by_bp_then_edge.items():
        for edge, angle_stats in edges.items():
            if edge not in results:
                results[edge] = {}
            results[edge][bp_type] = angle_stats

    # Sort for consistent output: edges alphabetically, _OTHER last; bp_types alphabetically within each edge
    def edge_sort_key(item):
        edge = item[0]
        return (edge == "_OTHER", edge)  # _OTHER sorts last
    sorted_results = {
        edge: dict(sorted(bps.items(), key=lambda x: x[0]))
        for edge, bps in sorted(results.items(), key=edge_sort_key)
    }
    return sorted_results


def main():
    print("Loading unique PDB IDs from data/uniqueRNAS.csv...")
    unique_pdb_ids = load_unique_pdb_ids()
    exclude = get_excluded_pdb_ids()
    print(f"  {len(unique_pdb_ids)} unique PDB IDs")
    print(f"  Excluding {len(exclude)} G-quad PDBs from g_quads.py")
    print(f"  Processing {len(unique_pdb_ids - exclude)} PDBs")

    print("\nCollecting torsion angles by base-pair type and edge...")
    grouped, counts, totals = collect_torsions_by_bp_and_edge(unique_pdb_ids, exclude)

    print("\nComputing allowed ranges (unimodal: 2.5-97.5%%; bimodal: GMM clusters)...")
    results = build_results(grouped, counts, totals)

    total_base_pairs = sum(totals.values())
    output = {"total_base_pairs": total_base_pairs, **results}

    print(f"\nWriting results to {OUTPUT_JSON}...")
    with open(OUTPUT_JSON, "w") as f:
        json.dump(output, f, indent=2)

    n_edges = len(results)
    n_bp_types = sum(len(bps) for bps in results.values())
    print(f"Done. Total base pairs: {total_base_pairs:,}. {n_edges} edge types, {n_bp_types} base-pair categories.")
    print(f"Output: {OUTPUT_JSON}")


if __name__ == "__main__":
    main()
    
    
#python3 analyze_torsions_by_edge_type.py