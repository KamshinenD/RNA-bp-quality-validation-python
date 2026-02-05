#!/usr/bin/env python3
"""
Generate torsion thresholds for config from torsion_thresholds_analysis.json.

Rules:
- Threshold = mean ± 1SD (outside range = outlier)
- Only edge types with >1k total data
- Only base-pair types with >1k data (within each edge)
- Exclude edge types containing a dot (.)
- Edge type "--" goes to global _OTHER
- _OTHER bin within each qualifying edge for bp_types with <1k
- Global _OTHER for edges that don't qualify (., --, or <1k total)

Output: TORSION_THRESHOLDS_BY_EDGE_BY_BASE_PAIR written to config.
"""

import json
from pathlib import Path
from collections import defaultdict

SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_JSON = SCRIPT_DIR / "torsion_thresholds_analysis.json"
CONFIG_PY = SCRIPT_DIR / "config.py"
MIN_COUNT = 1000
TORSION_ANGLES = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "chi", "v1", "v2", "v3", "v4"]


def get_bp_count(bp_data: dict) -> int:
    """Get representative count for (edge, bp_type) - use min across angles."""
    counts = [a.get("count", 0) for a in bp_data.values() if isinstance(a, dict) and "count" in a]
    return min(counts) if counts else 0


def get_edge_total_count(edge_data: dict) -> int:
    """Total count for an edge (sum of max per bp_type)."""
    total = 0
    for bp_type, bp_data in edge_data.items():
        if bp_type == "total_base_pairs":
            continue
        counts = [a.get("count", 0) for a in bp_data.values() if isinstance(a, dict) and "count" in a]
        if counts:
            total += max(counts)
    return total


def edge_qualifies(edge_name: str, edge_data: dict) -> bool:
    """Edge qualifies if: no dot, not '--', total count >= MIN_COUNT."""
    if "." in edge_name or edge_name == "--":
        return False
    total = get_edge_total_count(edge_data)
    return total >= MIN_COUNT


def build_thresholds(data: dict) -> dict:
    """
    Build TORSION_THRESHOLDS_BY_EDGE_BY_BASE_PAIR.
    Structure: edge -> bp_type -> {angle: (min, max)} where min/max = mean ± 1SD
    """
    result = {}
    global_other_bucket = defaultdict(lambda: defaultdict(dict))  # angle -> agg values for stats

    for edge_name, edge_data in data.items():
        if edge_name == "total_base_pairs":
            continue

        if not edge_qualifies(edge_name, edge_data):
            # Add all bp_types from this edge to global _OTHER
            for bp_type, bp_data in edge_data.items():
                if not isinstance(bp_data, dict):
                    continue
                for angle, stats in bp_data.items():
                    if angle not in TORSION_ANGLES or not isinstance(stats, dict):
                        continue
                    mn = stats.get("mean_minus_1sd")
                    mx = stats.get("mean_plus_1sd")
                    if mn is not None and mx is not None:
                        if angle not in global_other_bucket[bp_type]:
                            global_other_bucket[bp_type][angle] = {"mins": [], "maxs": []}
                        global_other_bucket[bp_type][angle]["mins"].append(mn)
                        global_other_bucket[bp_type][angle]["maxs"].append(mx)
            continue

        # Edge qualifies - build edge entry with _OTHER for rare bp_types
        edge_result = {}
        other_bucket = defaultdict(lambda: {"mins": [], "maxs": [], "counts": []})

        for bp_type, bp_data in edge_data.items():
            if not isinstance(bp_data, dict):
                continue
            count = get_bp_count(bp_data)
            if count < MIN_COUNT:
                for angle, stats in bp_data.items():
                    if angle not in TORSION_ANGLES or not isinstance(stats, dict):
                        continue
                    mn = stats.get("mean_minus_1sd")
                    mx = stats.get("mean_plus_1sd")
                    if mn is not None and mx is not None:
                        other_bucket[angle]["mins"].append(mn)
                        other_bucket[angle]["maxs"].append(mx)
                        other_bucket[angle]["counts"].append(stats.get("count", 0))
            else:
                bp_thresholds = {}
                for angle, stats in bp_data.items():
                    if angle not in TORSION_ANGLES or not isinstance(stats, dict):
                        continue
                    mn = stats.get("mean_minus_1sd")
                    mx = stats.get("mean_plus_1sd")
                    if mn is not None and mx is not None:
                        bp_thresholds[angle] = (round(mn, 2), round(mx, 2))
                if bp_thresholds:
                    edge_result[bp_type] = bp_thresholds

        if other_bucket:
            other_thresholds = {}
            for angle, agg in other_bucket.items():
                if not agg["mins"]:
                    continue
                # Use weighted average of ranges, or min of mins and max of maxs
                min_val = min(agg["mins"])
                max_val = max(agg["maxs"])
                other_thresholds[angle] = (round(min_val, 2), round(max_val, 2))
            if other_thresholds:
                edge_result["_OTHER"] = other_thresholds

        if edge_result:
            result[edge_name] = edge_result

    # Build global _OTHER: just angles, no bp_type breakdown (aggregate across all "other" bp_types)
    if global_other_bucket:
        angle_aggs = defaultdict(lambda: {"mins": [], "maxs": []})
        for bp_type, angles in global_other_bucket.items():
            for angle, agg in angles.items():
                if not agg.get("mins"):
                    continue
                angle_aggs[angle]["mins"].extend(agg["mins"])
                angle_aggs[angle]["maxs"].extend(agg["maxs"])
        other_thresholds = {}
        for angle, agg in angle_aggs.items():
            if agg["mins"]:
                min_val = min(agg["mins"])
                max_val = max(agg["maxs"])
                other_thresholds[angle] = (round(min_val, 2), round(max_val, 2))
        if other_thresholds:
            result["_OTHER"] = other_thresholds

    return result


def main():
    if not INPUT_JSON.exists():
        print(f"Error: {INPUT_JSON} not found. Run analyze_torsions_by_edge_type.py first.")
        return 1

    with open(INPUT_JSON) as f:
        data = json.load(f)

    thresholds = build_thresholds(data)
    n_edges = len([k for k in thresholds if k != "_OTHER"])
    n_bp_total = sum(len(bps) for bps in thresholds.values())

    # Generate Python config block
    block_lines = [
        "",
        "    # ===== BACKBONE TORSION THRESHOLDS (mean ± 1SD, outside = outlier) =====",
        "    # Structure: edge -> bp_type -> angle -> (min, max)",
        "    # Only edge/bp_type with >=1000 data; edges with '.', '--', or <1k go to _OTHER",
        f"    # Generated from torsion_thresholds_analysis.json",
        "    TORSION_THRESHOLDS_BY_EDGE_BY_BASE_PAIR = {",
    ]

    for edge in sorted(thresholds.keys(), key=lambda x: (x == "_OTHER", x)):
        bps = thresholds[edge]
        block_lines.append(f'        "{edge}": {{')
        if edge == "_OTHER":
            # _OTHER: just angles, no bp_type breakdown
            ang_str = ", ".join(f'"{a}": ({lo}, {hi})' for a, (lo, hi) in sorted(bps.items()))
            block_lines.append(f'            {ang_str},')
        else:
            for bp_type in sorted(bps.keys(), key=lambda x: (x == "_OTHER", x)):
                angles = bps[bp_type]
                ang_str = ", ".join(f'"{a}": ({lo}, {hi})' for a, (lo, hi) in sorted(angles.items()))
                block_lines.append(f'            "{bp_type}": {{{ang_str}}},')
        block_lines.append("        },")
    block_lines.append("    }")

    config_block = "\n".join(block_lines)

    config_content = CONFIG_PY.read_text()
    if "TORSION_THRESHOLDS_BY_EDGE_BY_BASE_PAIR" in config_content:
        # Replace existing block - find start and end by bracket matching
        start = config_content.find("    # ===== BACKBONE TORSION THRESHOLDS")
        if start >= 0:
            # Find the dict start
            dict_start = config_content.find("TORSION_THRESHOLDS_BY_EDGE_BY_BASE_PAIR = {", start)
            if dict_start >= 0:
                depth = 0
                i = config_content.find("{", dict_start)
                for j, c in enumerate(config_content[i:], i):
                    if c == "{":
                        depth += 1
                    elif c == "}":
                        depth -= 1
                    if depth == 0:
                        end = j + 1
                        break
                else:
                    end = len(config_content)
                config_content = config_content[:start] + config_block + "\n" + config_content[end:]
                CONFIG_PY.write_text(config_content)
                print(f"Replaced TORSION_THRESHOLDS_BY_EDGE_BY_BASE_PAIR in {CONFIG_PY}")
            else:
                print(config_block)
        else:
            print(config_block)
    else:
        # Insert before MAX_HOTSPOT_SCORE line (at end of Config class)
        insert_marker = "    MAX_HOTSPOT_SCORE = 70.0"
        if insert_marker in config_content:
            idx = config_content.find(insert_marker) + len(insert_marker)
            config_content = config_content[:idx] + "\n\n" + config_block + config_content[idx:]
            CONFIG_PY.write_text(config_content)
            print(f"Added TORSION_THRESHOLDS_BY_EDGE_BY_BASE_PAIR to {CONFIG_PY}")
        else:
            print("Config block to add manually:\n")
            print(config_block)

    n_bp = sum(len(bps) for k, bps in thresholds.items() if k != "_OTHER")
    print(f"\nSummary: {n_edges} edge types, {n_bp} base-pair categories; _OTHER = angles only (no bp breakdown)")
    return 0


if __name__ == "__main__":
    exit(main())
    
    
    
#python3 generate_torsion_config.py
