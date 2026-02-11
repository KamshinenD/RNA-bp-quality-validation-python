#!/usr/bin/env python3
"""
Export per-base-pair data for all motifs into a single CSV.

Relies on cached metadata (metadata_cache.json) and existing basepair/H-bond files.
"""

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd

from config import Config
from utils.data_loader import DataLoader
from scorer2 import Scorer
from app import filter_motif_data  # reuse existing motif filtering logic
from analyze_by_edge_type import is_adjacent_pair  # exclude adjacent base pairs (res diff=1)


def load_metadata_cache(cache_file: str = "metadata_cache.json") -> dict:
    if Path(cache_file).exists():
        try:
            with open(cache_file, "r") as f:
                return json.load(f)
        except Exception:
            return {}
    return {}


def get_pdb_metadata(pdb_id: str, cache: dict) -> dict:
    metrics = cache.get("validation_metrics", {}).get(pdb_id, {}) or {}
    deposition_date = metrics.get("deposition_date", "") or metrics.get("Deposition_Date", "")
    deposition_year = ""
    if deposition_date:
        m = re.match(r"(\\d{4})", deposition_date)
        if m:
            deposition_year = m.group(1)

    resolution = (
        metrics.get("refinement_resolution")
        or metrics.get("em_resolution")
        or metrics.get("em_diffraction_resolution")
        or ""
    )
    method = (
        metrics.get("structure_determination_method")
        or metrics.get("experimental_method")
        or metrics.get("Experimental_method")
        or ""
    )

    return {
        "resolution": resolution,
        "method": method,
        "deposition_year": deposition_year,
    }


def parse_motif_cif(motif_path: Path):
    """
    Returns (pdb_id, chain, residues_set, motif_range_start, motif_range_end).
    """
    pdb_match = re.search(r"([0-9][A-Z0-9]{3})", motif_path.stem)
    pdb_id = pdb_match.group(1).upper() if pdb_match else None
    chain = None
    residues = []
    with open(motif_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                parts = line.split()
                if len(parts) >= 6:
                    if chain is None:
                        chain = parts[5]
                    res_num = parts[4]
                    base = parts[3]
                    if res_num.strip() == "":
                        continue
                    residues.append(f"{chain}-{base}-{res_num}-")
    if not residues:
        return pdb_id, chain, set(), None, None
    res_nums = []
    for r in residues:
        try:
            res_nums.append(int(r.split("-")[2]))
        except (IndexError, ValueError):
            continue
    if not res_nums:
        return pdb_id, chain, set(residues), None, None
    res_nums_sorted = sorted(set(res_nums))
    return pdb_id, chain, set(residues), res_nums_sorted[0], res_nums_sorted[-1]


def has_protein_binding(all_hbond_df: pd.DataFrame, res1: str, res2: str) -> bool:
    if all_hbond_df is None or all_hbond_df.empty:
        return False
    mask = (
        ((all_hbond_df["res_1"] == res1) & (all_hbond_df["res_type_2"] != "RNA"))
        | ((all_hbond_df["res_2"] == res1) & (all_hbond_df["res_type_1"] != "RNA"))
        | ((all_hbond_df["res_1"] == res2) & (all_hbond_df["res_type_2"] != "RNA"))
        | ((all_hbond_df["res_2"] == res2) & (all_hbond_df["res_type_1"] != "RNA"))
    )
    return bool(mask.any())


def summarize_hbonds(hbond_rows: pd.DataFrame) -> dict:
    if hbond_rows is None or hbond_rows.empty:
        return {
            "distance": "",
            "angle1": "",
            "angle2": "",
            "dihedral": "",
            "hbond_quality": "",
            "num_hbonds": 0,
        }
    distance = hbond_rows["distance"].mean() if "distance" in hbond_rows else ""
    angle1 = hbond_rows["angle_1"].mean() if "angle_1" in hbond_rows else ""
    angle2 = hbond_rows["angle_2"].mean() if "angle_2" in hbond_rows else ""
    dihedral = hbond_rows["dihedral_angle"].mean() if "dihedral_angle" in hbond_rows else ""
    quality = ""
    if "quality" in hbond_rows:
        quality_counts = Counter(hbond_rows["quality"].dropna().tolist())
        if quality_counts:
            quality = quality_counts.most_common(1)[0][0]
    return {
        "distance": round(distance, 2) if distance != "" else "",
        "angle1": round(angle1, 2) if angle1 != "" else "",
        "angle2": round(angle2, 2) if angle2 != "" else "",
        "dihedral": round(dihedral, 2) if dihedral != "" else "",
        "hbond_quality": quality,
        "num_hbonds": len(hbond_rows),
    }


def export_motif_basepairs(
    motifs_dir: Path,
    output_csv: Path,
    cache: dict,
    shard_by_pdb: bool = False,
    shard_dir: Path = None,
    pdb_filters=None,
):
    config = Config()
    data_loader = DataLoader(config)
    scorer = Scorer(config)
    pdb_filters = {p.upper() for p in pdb_filters} if pdb_filters else None

    motif_files = sorted(motifs_dir.glob("*.cif"))
    if not motif_files:
        print(f"No motif CIF files found in {motifs_dir}")
        return

    # Keyed storage to allow overwrite: (pdb_id, motif_name, res1, res2)
    rows_by_key = {}
    if not shard_by_pdb and output_csv.exists():
        try:
            with open(output_csv, "r", newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    res1_existing = row.get("res1") or (row.get("base_pair", "").split("-")[0] if row.get("base_pair") else "")
                    res2_existing = row.get("res2") or (row.get("base_pair", "").split("-")[-1] if row.get("base_pair") else "")
                    key = (row.get("pdb_id", ""), row.get("motif", ""), res1_existing, res2_existing)
                    rows_by_key[key] = row
        except Exception:
            rows_by_key = {}

    fieldnames = [
        "pdb_id",
        "resolution",
        "method",
        "deposition_year",
        "motif",
        "motif_type",
        "motif_chain",
        "motif_range",
        "res1",
        "res2",
        "base_pair",
        "lw_notation",
        "bp_type",
        "basepair_score",
        "isPoor",
        "shear",
        "stretch",
        "stagger",
        "buckle",
        "propeller",
        "opening",
        "dihedral_angle",
        "distance",
        "angle1",
        "angle2",
        "hbond_quality",
        "hbond_score",
        "number_of_hbonds",
        "HasProtein_binding",
        "geometry_penalty",
        "hbond_penalty",
        "avg_suiteness",
        "res1_conformer",
        "res1_suiteness",
        "res2_conformer",
        "res2_suiteness",
        "backbone_outlier",
        "chi_outlier",
        "res1_chi_conf",
        "res2_chi_conf",
        "issues",
    ]
    total_basepairs_written = 0

    for idx, motif_path in enumerate(motif_files, 1):
        try:
            motif_name = motif_path.stem
            motif_type = motif_name.split("-")[0] if "-" in motif_name else motif_name
            pdb_id, chain, motif_residues, start_res, end_res = parse_motif_cif(motif_path)

            if pdb_filters and pdb_id and pdb_id.upper() not in pdb_filters:
                continue

            if not pdb_id or not motif_residues:
                print(f"[{idx}/{len(motif_files)}] Skipping {motif_name}: could not parse residues/PDB ID")
                continue

            # No per-motif print to reduce noise; progress will be per 100 base pairs

            # Load data
            basepairs = data_loader.load_basepairs(pdb_id, quiet=True)
            hbonds = data_loader.load_hbonds(pdb_id, quiet=True)
            all_hbonds = data_loader.load_all_hbonds(pdb_id, quiet=True)
            torsion_data = data_loader.load_torsions(pdb_id, quiet=True)

            if basepairs is None or hbonds is None:
                print(f"  ✗ Missing data for {pdb_id}, skipping")
                continue

            # Filter to motif residues
            motif_bps, motif_hbonds = filter_motif_data(
                basepairs,
                hbonds,
                motif_residues=motif_residues,
                start_res=start_res,
                end_res=end_res,
                chain=chain,
            )

            # Exclude adjacent base pairs (same chain, residue numbers differ by 1)
            motif_bps = [bp for bp in motif_bps if not is_adjacent_pair(bp.get('res_1', ''), bp.get('res_2', ''))]

            pdb_meta = get_pdb_metadata(pdb_id, cache)

            # Load existing rows for this motif (for dedup/overwrite)
            existing_rows = {}
            if shard_by_pdb:
                motif_csv = (shard_dir / pdb_id / f"{motif_name}.csv")
                if motif_csv.exists():
                    try:
                        with open(motif_csv, "r", newline="") as f:
                            reader = csv.DictReader(f)
                            for row in reader:
                                key = (
                                    row.get("res1", ""),
                                    row.get("res2", ""),
                                    row.get("bp_type", ""),
                                    row.get("lw_notation", ""),
                                )
                                existing_rows[key] = row
                    except Exception:
                        existing_rows = {}

            for bp in motif_bps:
                try:
                    bp_score = scorer._score_base_pair(bp, motif_hbonds, torsion_data)  # reuse scoring logic
                except Exception:
                    continue  # skip problematic base pair
                bp_info = bp_score.get("bp_info", {})
                res1 = bp_info.get("res_1", "")
                res2 = bp_info.get("res_2", "")
                lw = bp_info.get("lw") or bp.get("lw", "")
                bp_type = bp_info.get("bp_type") or bp.get("bp_type", "")

                # Get H-bonds associated with this base pair
                hbond_rows = scorer._get_basepair_hbonds(res1, res2, motif_hbonds)
                hbond_summary = summarize_hbonds(hbond_rows)

                issues = []
                for issue, present in bp_score.get("geometry_issues", {}).items():
                    if present:
                        issues.append(f"geom_{issue}")
                for issue, present in bp_score.get("hbond_issues", {}).items():
                    if present:
                        issues.append(f"hbond_{issue}")
                # Append backbone suiteness details
                backbone = bp_score.get("backbone", [])
                for bd in backbone:
                    if bd.get('is_outlier', False):
                        issues.append(f"backbone_outlier({bd.get('residue','?')})")
                    elif bd.get('suiteness', 1.0) < 0.5:
                        issues.append(f"low_suiteness({bd.get('residue','?')},s={bd.get('suiteness',0):.2f})")

                has_binding = has_protein_binding(all_hbonds, res1, res2)

                row = {
                    "pdb_id": pdb_id,
                    "resolution": pdb_meta["resolution"],
                    "method": pdb_meta["method"],
                    "deposition_year": pdb_meta["deposition_year"],
                    "motif": motif_name,
                    "motif_type": motif_type,
                    "motif_chain": chain or "",
                    "motif_range": f"{start_res}-{end_res}" if start_res and end_res else "",
                    "res1": res1,
                    "res2": res2,
                    "base_pair": f"{res1}-{res2}",
                    "lw_notation": lw,
                    "bp_type": bp_type,
                    "basepair_score": bp_score.get("score", ""),
                    "isPoor": bp_score.get("score", 100) < config.BASELINE,
                    "shear": bp.get("shear", ""),
                    "stretch": bp.get("stretch", ""),
                    "stagger": bp.get("stagger", ""),
                    "buckle": bp.get("buckle", ""),
                    "propeller": bp.get("propeller", ""),
                    "opening": bp.get("opening", ""),
                    "dihedral_angle": hbond_summary["dihedral"],
                    "distance": hbond_summary["distance"],
                    "angle1": hbond_summary["angle1"],
                    "angle2": hbond_summary["angle2"],
                    "hbond_quality": hbond_summary["hbond_quality"],
                    "hbond_score": bp.get("hbond_score", ""),
                    "number_of_hbonds": hbond_summary["num_hbonds"],
                    "HasProtein_binding": has_binding,
                    "geometry_penalty": bp_score.get("geometry_penalty", ""),
                    "hbond_penalty": bp_score.get("hbond_penalty", ""),
                    "avg_suiteness": "",
                    "res1_conformer": "",
                    "res1_suiteness": "",
                    "res2_conformer": "",
                    "res2_suiteness": "",
                    "backbone_outlier": False,
                    "chi_outlier": False,
                    "res1_chi_conf": "",
                    "res2_chi_conf": "",
                    "issues": ",".join(issues) if issues else "",
                }

                # Populate backbone suiteness columns
                if backbone:
                    suiteness_vals = [bd.get('suiteness', 0.0) for bd in backbone]
                    row["avg_suiteness"] = round(sum(suiteness_vals) / len(suiteness_vals), 3) if suiteness_vals else ""
                    row["backbone_outlier"] = any(bd.get('is_outlier', False) for bd in backbone)
                    if len(backbone) >= 1:
                        row["res1_conformer"] = backbone[0].get('conformer', '')
                        row["res1_suiteness"] = backbone[0].get('suiteness', '')
                    if len(backbone) >= 2:
                        row["res2_conformer"] = backbone[1].get('conformer', '')
                        row["res2_suiteness"] = backbone[1].get('suiteness', '')

                # Populate chi conformation columns
                chi_details = bp_score.get('chi_details', {})
                if chi_details:
                    row["chi_outlier"] = chi_details.get('chi_outlier', False)
                    row["res1_chi_conf"] = chi_details.get('res1_chi_conf', '')
                    row["res2_chi_conf"] = chi_details.get('res2_chi_conf', '')

                if shard_by_pdb:
                    key = (row["res1"], row["res2"], row["bp_type"], row["lw_notation"])
                    existing_rows[key] = row  # overwrite if already present
                else:
                    key = (pdb_id, motif_name, res1, res2)
                    rows_by_key[key] = row  # overwrite if key already present

                total_basepairs_written += 1
                if total_basepairs_written % 100 == 0:
                    print(f"Progress: {total_basepairs_written} base pairs written...")

            if shard_by_pdb:
                shard_dir.mkdir(parents=True, exist_ok=True)
                pdb_dir = shard_dir / pdb_id
                pdb_dir.mkdir(parents=True, exist_ok=True)
                motif_csv = pdb_dir / f"{motif_name}.csv"
                with open(motif_csv, "w", newline="") as f:
                    writer = csv.DictWriter(f, fieldnames=fieldnames)
                    writer.writeheader()
                    writer.writerows(existing_rows.values())
        except Exception:
            # Skip any motif that raises unexpected errors
            continue

    if shard_by_pdb:
        print(f"\nExport complete. Base pairs written: {total_basepairs_written}. Sharded CSVs saved under: {shard_dir}")
        return

    # Write all rows back (existing + new), ensuring overwrite behavior
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_by_key.values())

    print(f"\nExport complete. Base pairs written: {len(rows_by_key)}. CSV saved to: {output_csv}")


def merge_sharded_csvs(shard_dir: Path, output_csv: Path):
    fieldnames = [
        "pdb_id",
        "resolution",
        "method",
        "deposition_year",
        "motif",
        "motif_type",
        "motif_chain",
        "motif_range",
        "motif_length",
        "res1",
        "res1_chain",
        "res2",
        "res2_chain",
        "base_pair",
        "lw_notation",
        "bp_type",
        "bp_stability_class",
        "is_interchain",
        "basepair_score",
        "isPoor",
        "shear",
        "stretch",
        "stagger",
        "buckle",
        "propeller",
        "opening",
        "dihedral_angle",
        "distance",
        "angle1",
        "angle2",
        "hbond_quality",
        "hbond_score",
        "number_of_hbonds",
        "HasProtein_binding",
        "Protein_Binding_Details",
        "Protein_AA_Types",
        "Protein_HBond_Count",
        "HasLigand_binding",
        "Ligand_Binding_Details",
        "Ligand_Types",
        "Ligand_HBond_Count",
        "geometry_penalty",
        "hbond_penalty",
        "avg_suiteness",
        "res1_conformer",
        "res1_suiteness",
        "res2_conformer",
        "res2_suiteness",
        "backbone_outlier",
        "chi_outlier",
        "res1_chi_conf",
        "res2_chi_conf",
        "issues",
    ]

    rows_by_key = {}
    for motif_csv in shard_dir.glob("**/*.csv"):
        try:
            with open(motif_csv, "r", newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    res1_existing = row.get("res1") or (row.get("base_pair", "").split("-")[0] if row.get("base_pair") else "")
                    res2_existing = row.get("res2") or (row.get("base_pair", "").split("-")[-1] if row.get("base_pair") else "")
                    key = (row.get("pdb_id", ""), row.get("motif", ""), res1_existing, res2_existing)
                    rows_by_key[key] = row
        except Exception:
            continue

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_by_key.values())

    print(f"\nMerge complete. CSV saved to: {output_csv} (rows: {len(rows_by_key)})")


def main():
    parser = argparse.ArgumentParser(description="Export per-base-pair motif data to CSV")
    parser.add_argument("--motifs-dir", default="motifs", help="Directory containing motif CIF files")
    parser.add_argument("--output", default="motif_basepairs.csv", help="Output CSV path (used when not sharding)")
    parser.add_argument("--cache", default="metadata_cache.json", help="Metadata cache JSON file")
    parser.add_argument(
        "--shard-by-pdb",
        action="store_true",
        help="If set, writes one CSV per motif under shard-dir/pdb_id/ to avoid race conditions",
    )
    parser.add_argument(
        "--shard-dir",
        default="motif_base_pair",
        help="Directory to write sharded motif CSVs when --shard-by-pdb is set (default: motif_base_pair)",
    )
    parser.add_argument(
        "--merge-shards",
        action="store_true",
        help="Merge all sharded CSVs under shard-dir into the --output file and exit",
    )
    parser.add_argument(
        "--pdb",
        action="append",
        dest="pdb_filters",
        help="Limit export to one or more PDB IDs (repeatable); default is all",
    )
    args = parser.parse_args()

    motifs_dir = Path(args.motifs_dir)
    output_csv = Path(args.output)
    shard_dir = Path(args.shard_dir)
    cache = load_metadata_cache(args.cache)

    if args.merge_shards:
        merge_sharded_csvs(shard_dir, output_csv)
        return

    export_motif_basepairs(
        motifs_dir,
        output_csv,
        cache,
        shard_by_pdb=args.shard_by_pdb,
        shard_dir=shard_dir,
        pdb_filters=args.pdb_filters,
    )


if __name__ == "__main__":
    main()

