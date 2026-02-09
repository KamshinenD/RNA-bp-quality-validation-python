#!/usr/bin/env python3
"""
Generate torsion thresholds for config from torsion_thresholds_analysis.json.

Rules:
- Unimodal: 2.5th-97.5th percentile. Bimodal: GMM clusters, 2.5th-97.5th per cluster.
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
TORSION_ANGLES = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "chi"]


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


def _ranges_to_config_val(ranges: list) -> tuple:
    """Convert ranges list to config value: (min,max) if single, [(m1,M1),(m2,M2)] if multi."""
    if not ranges:
        return None
    if len(ranges) == 1:
        return tuple(ranges[0])
    return [tuple(r) for r in ranges]


def _extract_ranges(stats: dict, angle: str) -> list:
    """Get ranges from stats. Handles old format (mean_minus_2sd) and new (ranges)."""
    if "ranges" in stats and stats["ranges"]:
        return stats["ranges"]
    mn = stats.get("mean_minus_2sd")
    mx = stats.get("mean_plus_2sd")
    if mn is not None and mx is not None:
        return [(mn, mx)]
    return []


def build_thresholds(data: dict) -> dict:
    """
    Build TORSION_THRESHOLDS_BY_EDGE_BY_BASE_PAIR.
    Structure: edge -> bp_type -> {angle: (min, max) or [(m1,M1),(m2,M2)]}
    """
    result = {}
    global_other_bucket = defaultdict(lambda: defaultdict(lambda: {"mins": [], "maxs": []}))

    for edge_name, edge_data in data.items():
        if edge_name == "total_base_pairs":
            continue

        if not edge_qualifies(edge_name, edge_data):
            for bp_type, bp_data in edge_data.items():
                if not isinstance(bp_data, dict):
                    continue
                for angle, stats in bp_data.items():
                    if angle not in TORSION_ANGLES or not isinstance(stats, dict):
                        continue
                    for lo, hi in _extract_ranges(stats, angle):
                        global_other_bucket[bp_type][angle]["mins"].append(lo)
                        global_other_bucket[bp_type][angle]["maxs"].append(hi)
            continue

        edge_result = {}
        other_bucket = defaultdict(lambda: {"mins": [], "maxs": []})

        for bp_type, bp_data in edge_data.items():
            if not isinstance(bp_data, dict):
                continue
            count = get_bp_count(bp_data)
            if count < MIN_COUNT:
                for angle, stats in bp_data.items():
                    if angle not in TORSION_ANGLES or not isinstance(stats, dict):
                        continue
                    for lo, hi in _extract_ranges(stats, angle):
                        other_bucket[angle]["mins"].append(lo)
                        other_bucket[angle]["maxs"].append(hi)
            else:
                bp_thresholds = {}
                for angle, stats in bp_data.items():
                    if angle not in TORSION_ANGLES or not isinstance(stats, dict):
                        continue
                    ranges = _extract_ranges(stats, angle)
                    if ranges:
                        val = _ranges_to_config_val(ranges)
                        if val is not None:
                            bp_thresholds[angle] = val
                if bp_thresholds:
                    edge_result[bp_type] = bp_thresholds

        if other_bucket:
            other_thresholds = {}
            for angle, agg in other_bucket.items():
                if agg["mins"]:
                    min_val = min(agg["mins"])
                    max_val = max(agg["maxs"])
                    other_thresholds[angle] = (round(min_val, 2), round(max_val, 2))
            if other_thresholds:
                edge_result["_OTHER"] = other_thresholds

        if edge_result:
            result[edge_name] = edge_result

    if global_other_bucket:
        angle_aggs = defaultdict(lambda: {"mins": [], "maxs": []})
        for bp_type, angles in global_other_bucket.items():
            for angle, agg in angles.items():
                if agg.get("mins"):
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
        "    # ===== BACKBONE TORSION THRESHOLDS (2.5-97.5 pct; unimodal or bimodal GMM) =====",
        "    # angle -> (min, max) or [(m1,M1),(m2,M2)] for bimodal",
        "    # Only edge/bp_type with >=1000 data; edges with '.', '--', or <1k go to _OTHER",
        f"    # Generated from torsion_thresholds_analysis.json",
        "    TORSION_THRESHOLDS_BY_EDGE_BY_BASE_PAIR = {",
    ]

    def _fmt_val(v):
        if isinstance(v, tuple):
            return f"({v[0]}, {v[1]})"
        return "[" + ", ".join(f"({a}, {b})" for a, b in v) + "]"

    for edge in sorted(thresholds.keys(), key=lambda x: (x == "_OTHER", x)):
        bps = thresholds[edge]
        block_lines.append(f'        "{edge}": {{')
        if edge == "_OTHER":
            ang_str = ", ".join(f'"{a}": {_fmt_val(val)}' for a, val in sorted(bps.items()))
            block_lines.append(f'            {ang_str},')
        else:
            for bp_type in sorted(bps.keys(), key=lambda x: (x == "_OTHER", x)):
                angles = bps[bp_type]
                ang_str = ", ".join(f'"{a}": {_fmt_val(val)}' for a, val in sorted(angles.items()))
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
