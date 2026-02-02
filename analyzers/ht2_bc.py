"""Optimized hotspot analyzer with pre-computed lookups for speed."""

import pandas as pd
from typing import Dict, List, Set, Tuple
from dataclasses import dataclass
from collections import defaultdict
from .analyzers_utils import BasePairScoring, HBondScoring, ScoringUtils

@dataclass
class Hotspot:
    """Represents a poorly refined region."""
    chain: str
    start_res: int
    end_res: int
    residue_count: int
    score: float
    severity: str
    issue_density: float
    dominant_issues: List[str]
    base_pairs_affected: int
    hbonds_affected: int
    details: Dict

class HotspotAnalyzer:
    """Finds hotspots using connectivity-based approach with hierarchical classification."""
    
    def __init__(self, config, bp_analyzer, hb_analyzer):
        self.config = config
        self.bp_analyzer = bp_analyzer
        self.hb_analyzer = hb_analyzer
    
    def find_hotspots(self, basepair_data: list, hbond_data: pd.DataFrame) -> List[Hotspot]:
        """Find hotspots using hierarchical geometry + H-bond analysis."""
        if not basepair_data:
            return []
        
        self._basepair_cache = basepair_data
        self._hbond_cache = hbond_data
        
        # ===== CRITICAL OPTIMIZATION: Pre-compute all H-bond lookups ONCE =====
        print("  Pre-computing H-bond lookups...")
        self._precompute_hbond_data()
        print(f"  Found {len(self._pair_to_hbonds)} base pairs with H-bonds")
        
        # Step 1: Score all residues
        print("  Scoring residues...")
        residue_scores = self._score_all_residues(basepair_data, hbond_data)
        
        if not residue_scores:
            print("  ERROR: No residue scores calculated!")
            return []
        
        # Calculate distribution statistics
        scores = list(residue_scores.values())
        scores_sorted = sorted(scores)
        mean_score = sum(scores) / len(scores)
        std_dev = (sum((s - mean_score)**2 for s in scores) / len(scores))**0.5
        
        percentile_5 = scores_sorted[int(len(scores) * 0.05)]
        statistical = mean_score - 1.0 * std_dev
        threshold = max(percentile_5, statistical, self.config.DAMAGE_THRESHOLD_FOR_HOTSPOTS)
        
        print(f"  Damage threshold: {threshold:.1f}")
        
        # Find damaged residues
        damaged_residues = self._find_damaged_residues(residue_scores, threshold)
        total_damaged = sum(len(res_set) for res_set in damaged_residues.values())

        if total_damaged == 0:
            print("  No damaged residues")
            return []
        
        print(f"  Found {total_damaged} damaged residues")
        
        # Build connectivity and find components
        print("  Building connectivity graph...")
        connectivity = self._build_connectivity_graph(basepair_data)
        
        print("  Finding connected components...")
        components = self._find_connected_components(damaged_residues, connectivity)
        
        min_size = getattr(self.config, 'MIN_HOTSPOT_RESIDUES', 3)
        
        print("  Creating hotspots...")
        hotspots = []
        for chain, residue_sets in components.items():
            for residues in residue_sets:
                if len(residues) >= min_size:
                    hotspot = self._create_hotspot_from_residues(chain, residues)
                    if hotspot:
                        hotspots.append(hotspot)
        
        print(f"  Created {len(hotspots)} initial hotspots")
        
        # Merge and filter
        print("  Merging overlapping hotspots...")
        hotspots = self._merge_overlapping_hotspots(hotspots)
        
        print("  Stitching hotspot chains...")
        hotspots = self._stitch_hotspot_chains(hotspots)
        
        print("  Filtering by context...")
        hotspots = self._filter_by_context(hotspots)
        
        hotspots.sort(key=lambda h: h.score)
        
        if not hotspots:
            print("\n  === NO HOTSPOTS - Overall structure quality is good ===")
        else:
            print(f"\n  === FOUND {len(hotspots)} HOTSPOTS ===")
        
        return hotspots
    
    def _precompute_hbond_data(self):
        """
        Pre-compute H-bond lookups ONCE to avoid repeated scans.
        This is the KEY optimization that speeds up everything.
        """
        # Initialize lookup dictionaries
        self._pair_to_hbonds = defaultdict(list)  # pair_key -> list of H-bonds
        self._pair_hbond_count = defaultdict(int)  # pair_key -> count
        self._pair_has_issues = set()  # set of pair_keys with H-bond issues
        
        if self._hbond_cache is None or self._hbond_cache.empty:
            return
        
        # Build base pair set for filtering
        base_pair_set = ScoringUtils.build_base_pair_set(self._basepair_cache)
        
        # Single pass through H-bond data
        for _, hb in self._hbond_cache.iterrows():
            # Filter 1: Only base-base H-bonds
            if not ScoringUtils.is_base_base_hbond(hb['atom_1'], hb['atom_2']):
                continue
            
            # Filter 2: Skip adjacent pairs
            if BasePairScoring.check_adjacent_pairing(hb['res_1'], hb['res_2']):
                continue
            
            # Filter 3: Only H-bonds in base pairs
            pair_key = tuple(sorted([hb['res_1'], hb['res_2']]))
            if pair_key not in base_pair_set:
                continue
            
            # Store H-bond for this pair
            self._pair_to_hbonds[pair_key].append(hb)
            self._pair_hbond_count[pair_key] += 1
            
            # Check if this H-bond has issues
            hb_issues = HBondScoring.score_hbond(hb, self.config)
            if any(hb_issues.values()):
                self._pair_has_issues.add(pair_key)
    def _classify_geometry_severity(self, bp: dict) -> str:
        """
        Classify geometry severity - MAXIMUM SENSITIVITY.
        
        Returns: 'SEVERE', 'MODERATE', 'MINOR', or 'NONE'
        """
        shear = abs(bp.get('shear', 0))
        buckle = abs(bp.get('buckle', 0))
        opening = abs(bp.get('opening', 0))
        propeller = bp.get('propeller', 0)
        stretch = bp.get('stretch', 0)
        stagger = abs(bp.get('stagger', 0))
        hbond_score = bp.get('hbond_score', 0)
        
        # SEVERE: Major structural distortion
        if (shear > 1.0 or              # More sensitive
            buckle > 22 or              # More sensitive
            opening > 16 or             # More sensitive
            hbond_score == 0.0 or
            stagger > 0.7):             # More sensitive
            return 'SEVERE'
        
        # MODERATE: Significant deviation
        if (shear > 0.5 or              # More sensitive
            buckle > 15 or              # More sensitive
            not (-30 <= propeller <= 10) or  # Wider range
            not (-0.7 <= stretch <= 0.3) or  # Wider range
            stagger > 0.4 or            # More sensitive
            (0 < hbond_score < 2.2)):   # More sensitive
            return 'MODERATE'
        
        # MINOR: Slight deviation
        if (shear > 0.3 or              # Very sensitive!
            buckle > 10 or              # Very sensitive!
            opening > 6 or              # Very sensitive!
            stagger > 0.25 or           # Very sensitive!
            (hbond_score > 0 and hbond_score < 2.8)):  # Very sensitive!
            return 'MINOR'
        
        return 'NONE'
    
    def _is_problematic_base_pair(self, bp: dict) -> Tuple[bool, str]:
        """
        Hierarchical classification - MORE SENSITIVE - OPTIMIZED.
        
        Returns: (is_problematic, reason)
        """
        geometry_severity = self._classify_geometry_severity(bp)
        has_hbond_issues = self._bp_has_hbond_issues_fast(bp)  # ← FAST VERSION
        
        # TIER 1: Severe geometry = automatic problem
        if geometry_severity == 'SEVERE':
            return True, 'severe_geometry'
        
        # TIER 2: Moderate geometry + H-bond issues = compounding
        if geometry_severity == 'MODERATE' and has_hbond_issues:
            return True, 'geometry_plus_hbond'
        
        # TIER 3: Moderate geometry alone
        if geometry_severity == 'MODERATE':
            return True, 'moderate_geometry'
        
        # TIER 4: Minor geometry + H-bond issues
        if geometry_severity == 'MINOR' and has_hbond_issues:
            return True, 'minor_geometry_with_hbond'
        
        # TIER 5: Minor geometry alone = NOW PROBLEMATIC
        if geometry_severity == 'MINOR':
            return True, 'minor_geometry'
        
        # TIER 6: Good geometry but H-bond issues = artifact
        if geometry_severity == 'NONE' and has_hbond_issues:
            return False, 'hbond_artifact'
        
        return False, 'acceptable'
    
    def _bp_has_hbond_issues_fast(self, bp: dict) -> bool:
        """
        OPTIMIZED: Check if base pair has H-bond issues using pre-computed data.
        This is O(1) instead of O(n) where n = number of H-bonds.
        """
        pair_key = tuple(sorted([bp['res_1'], bp['res_2']]))
        
        # Check geometry issues (pre-computed)
        if pair_key in self._pair_has_issues:
            return True
        
        # Check count issues (pre-computed)
        actual_count = self._pair_hbond_count.get(pair_key, 0)
        hbond_score = bp.get('hbond_score', 0)
        
        if hbond_score == 0.0:
            return False
        
        is_incorrect, _ = HBondScoring.check_hbond_count(
            bp.get('bp_type', 'unknown'), actual_count, bp.get('lw', ''), self.config
        )
        
        return is_incorrect

    def _merge_overlapping_hotspots(self, hotspots: List[Hotspot]) -> List[Hotspot]:
        """Merge hotspots that overlap significantly on the same chain."""
        if not hotspots:
            return []
        
        valid_hotspots = [h for h in hotspots if h is not None]
        if not valid_hotspots:
            return []
        
        by_chain = defaultdict(list)
        for h in valid_hotspots:
            by_chain[h.chain].append(h)
        
        merged = []
        
        for chain, chain_hotspots in by_chain.items():
            chain_hotspots.sort(key=lambda h: h.start_res)
            current = chain_hotspots[0]
            
            for next_hotspot in chain_hotspots[1:]:
                overlap_threshold = 0.3
                
                overlap_start = max(current.start_res, next_hotspot.start_res)
                overlap_end = min(current.end_res, next_hotspot.end_res)
                overlap_size = max(0, overlap_end - overlap_start + 1)
                
                min_size = min(
                    current.end_res - current.start_res + 1,
                    next_hotspot.end_res - next_hotspot.start_res + 1
                )
                
                gap = next_hotspot.start_res - current.end_res
                max_gap = 2 if current.severity in ["CRITICAL", "SEVERE"] else 1
                
                if (overlap_size / min_size > overlap_threshold) or (0 < gap <= max_gap):
                    merged_result = self._merge_two_hotspots(current, next_hotspot)
                    if merged_result is not None:
                        current = merged_result
                    else:
                        merged.append(current)
                        current = next_hotspot
                else:
                    merged.append(current)
                    current = next_hotspot
            
            if current is not None:
                merged.append(current)
        
        return merged
    
    def _merge_two_hotspots(self, h1: Hotspot, h2: Hotspot) -> Hotspot:
        """Merge two overlapping hotspots into one."""
        all_residues = set(range(h1.start_res, h1.end_res + 1))
        all_residues.update(range(h2.start_res, h2.end_res + 1))
        return self._create_hotspot_from_residues(h1.chain, all_residues)
    
    def _filter_by_context(self, hotspots: List[Hotspot]) -> List[Hotspot]:
        """Filter hotspots based on structure quality context."""
        if not hotspots:
            return []
        
        filtered = [h for h in hotspots if h is not None and h.severity != 'MINOR']
        return filtered
    
    def _score_all_residues(self, basepair_data: list, hbond_data: pd.DataFrame) -> Dict[Tuple[str, int], float]:
        """Score each residue using weighted penalties."""
        residue_penalties = defaultdict(lambda: {'penalty_sum': 0.0, 'bp_count': 0})

        # BASE-PAIR GEOMETRY PENALTIES
        for bp in basepair_data:
            bp_issues = self._count_bp_issues(bp)
            bp_penalty = sum(self.config.PENALTY_WEIGHTS.get(issue, 0) for issue in bp_issues)

            for res_id in [bp['res_1'], bp['res_2']]:
                chain, res_num = self._parse_residue(res_id)
                key = (chain, res_num)
                residue_penalties[key]['penalty_sum'] += bp_penalty
                residue_penalties[key]['bp_count'] += 1

        base_pair_set = ScoringUtils.build_base_pair_set(basepair_data)

        # HYDROGEN-BOND GEOMETRY PENALTIES
        if hbond_data is not None and not hbond_data.empty:
            for _, row in hbond_data.iterrows():
                if not ScoringUtils.is_base_base_hbond(row['atom_1'], row['atom_2']):
                    continue
                if BasePairScoring.check_adjacent_pairing(row['res_1'], row['res_2']):
                    continue
                
                pair_key = tuple(sorted([row['res_1'], row['res_2']]))
                if pair_key not in base_pair_set:
                    continue

                hb_issues = self._count_hb_issues(row)
                hb_penalty = sum(self.config.PENALTY_WEIGHTS.get(issue, 0) for issue in hb_issues)

                for res_id in [row['res_1'], row['res_2']]:
                    chain, res_num = self._parse_residue(res_id)
                    key = (chain, res_num)
                    residue_penalties[key]['penalty_sum'] += hb_penalty
                    residue_penalties[key]['bp_count'] += 1

        # H-BOND COUNT MISMATCH PENALTIES - OPTIMIZED
        if hbond_data is not None and not hbond_data.empty:
            for bp in basepair_data:
                res_1, res_2 = bp['res_1'], bp['res_2']
                pair_key = tuple(sorted([res_1, res_2]))
                actual_count = self._pair_hbond_count.get(pair_key, 0)  # ← FAST LOOKUP
                
                hbond_score = bp.get('hbond_score', 0)
                if hbond_score == 0.0:
                    continue

                count_issues = self._detect_hbond_count_issue(bp, actual_count)
                if count_issues:
                    count_penalty = sum(self.config.PENALTY_WEIGHTS.get(issue, 0) for issue in count_issues)
                    for res_id in [res_1, res_2]:
                        chain, res_num = self._parse_residue(res_id)
                        key = (chain, res_num)
                        residue_penalties[key]['penalty_sum'] += count_penalty
                        residue_penalties[key]['bp_count'] += 1

        # COMPUTE FINAL WEIGHTED SCORE PER RESIDUE
        residue_scores = {}
        for key, data in residue_penalties.items():
            avg_penalty = data['penalty_sum'] / max(data['bp_count'], 1)
            score = max(0, self.config.BASE_SCORE - avg_penalty)
            residue_scores[key] = score

        return residue_scores

    def _count_bp_issues(self, bp: dict) -> List[str]:
        """Return list of specific issue types found in a base pair."""
        issues = []
        
        bp_issues_dict = BasePairScoring.score_base_pair(bp, self.config)
        
        if bp_issues_dict is None:
            return []
        
        if bp_issues_dict['misaligned']:
            issues.append('misaligned_pairs')
        if bp_issues_dict['twisted']:
            issues.append('twisted_pairs')
        if bp_issues_dict['non_coplanar']:
            issues.append('non_coplanar_pairs')
        if bp_issues_dict['poor_hbond']:
            issues.append('poor_hbond_pairs')
        if bp_issues_dict['zero_hbond']:
            issues.append('zero_hbond_pairs')
        if bp_issues_dict['self_pairing']:
            issues.append('self_pairing')
        
        return issues
    
    def _count_hb_issues(self, hb: dict) -> list:
        """Identify specific hydrogen-bond issues."""
        issues = []
        hb_series = pd.Series(hb)
        hb_issues_dict = HBondScoring.score_hbond(hb_series, self.config)
        
        if hb_issues_dict['bad_distance']:
            issues.append('bad_hbond_distance')
        if hb_issues_dict['bad_angles']:
            issues.append('bad_hbond_angles')
        if hb_issues_dict['bad_dihedral']:
            issues.append('bad_hbond_dihedrals')
        if hb_issues_dict['weak_quality']:
            issues.append('weak_hbond_quality')
        
        return issues

    def _detect_hbond_count_issue(self, bp: dict, actual_count: int) -> list:
        """Check whether a base pair has an incorrect number of H-bonds."""
        issues = []
        bp_type = bp.get('bp_type', 'unknown')
        lw = bp.get('lw', '')
        
        is_incorrect, _ = HBondScoring.check_hbond_count(bp_type, actual_count, lw, self.config)
        
        if is_incorrect:
            issues.append('incorrect_hbond_count')
        
        return issues

    def _find_damaged_residues(self, residue_scores: Dict[Tuple[str, int], float], 
                              threshold: float) -> Dict[str, Set[int]]:
        """Find all residues with score below threshold."""
        damaged = defaultdict(set)
        
        for (chain, res_num), score in residue_scores.items():
            if score < threshold:
                damaged[chain].add(res_num)
        
        return damaged
    
    def _build_connectivity_graph(self, basepair_data: list) -> Dict[str, Dict[int, Set[int]]]:
        """Build connectivity graph for residues."""
        graph = defaultdict(lambda: defaultdict(set))
        
        for bp in basepair_data:
            chain1, res1 = self._parse_residue(bp['res_1'])
            chain2, res2 = self._parse_residue(bp['res_2'])
            
            graph[chain1][res1].add(res2)
            graph[chain2][res2].add(res1)
            graph[chain1][res1]
            graph[chain2][res2]
        
        neighbor_range = self.config.SEQUENTIAL_NEIGHBOR_RANGE
        
        for chain in list(graph.keys()):
            all_residues = sorted(graph[chain].keys())
            
            for i, res in enumerate(all_residues):
                for j in range(max(0, i - neighbor_range), min(len(all_residues), i + neighbor_range + 1)):
                    if i != j:
                        neighbor = all_residues[j]
                        graph[chain][res].add(neighbor)
                        graph[chain][neighbor].add(res)
        
        return graph
    
    def _find_connected_components(self, damaged_residues: Dict[str, Set[int]], 
                                  connectivity: Dict[str, Dict[int, Set[int]]]) -> Dict[str, List[Set[int]]]:
        """Find connected components of damaged residues."""
        components = defaultdict(list)
        
        for chain, damaged in damaged_residues.items():
            if chain not in connectivity:
                continue
            
            visited = set()
            for start_res in damaged:
                if start_res in visited:
                    continue
                
                component = self._flood_fill(chain, start_res, damaged, connectivity, visited)
                if component:
                    components[chain].append(component)
        
        return components
    
    def _flood_fill(self, chain: str, start_res: int, 
                   damaged_set: Set[int], 
                   connectivity: Dict[str, Dict[int, Set[int]]], 
                   visited: Set[int]) -> Set[int]:
        """Flood fill to find all damaged residues connected to start_res."""
        if start_res in visited or start_res not in damaged_set:
            return set()
        
        component = set()
        queue = [start_res]
        
        while queue:
            current = queue.pop(0)
            
            if current in visited:
                continue
            
            visited.add(current)
            component.add(current)
            
            if current in connectivity[chain]:
                for neighbor in connectivity[chain][current]:
                    if neighbor in damaged_set and neighbor not in visited:
                        queue.append(neighbor)
        
        return component
    
    def _create_hotspot_from_residues(self, chain: str, residues: Set[int]) -> Hotspot:
        """Create a hotspot from a set of connected damaged residues."""
        residue_list = sorted(residues)
        start_res = min(residue_list)
        end_res = max(residue_list)
        
        # Check for huge gaps
        gaps = [residue_list[i+1] - residue_list[i] for i in range(len(residue_list)-1)]
        if gaps and max(gaps) > 50:
            return None
        
        span = end_res - start_res + 1
        density = len(residue_list) / span
        
        if density < 0.02 and span > 100:
            return None

        # GET ALL BASE PAIRS IN REGION
        all_region_bps = []
        for bp in self._basepair_cache:
            chain1, res1 = self._parse_residue(bp['res_1'])
            chain2, res2 = self._parse_residue(bp['res_2'])
            
            res1_in_range = (chain == chain1 and start_res <= res1 <= end_res)
            res2_in_range = (chain == chain2 and start_res <= res2 <= end_res)
            
            if res1_in_range or res2_in_range:
                all_region_bps.append(bp)
        
        if not all_region_bps:
            return None
        
        # HIERARCHICAL CLASSIFICATION
        problematic_bps = []
        problem_types = defaultdict(int)
        
        for bp in all_region_bps:
            is_bad, reason = self._is_problematic_base_pair(bp)
            
            if is_bad:
                problematic_bps.append(bp)
                problem_types[reason] += 1
        
        # CALCULATE METRICS
        total_bps = len(all_region_bps)
        num_problematic = len(problematic_bps)
        fraction_bad = num_problematic / total_bps if total_bps > 0 else 0
        
        # VALIDATION FILTERS
        if num_problematic < 1:
            return None
        
        if fraction_bad < 0.15:
            return None
        
        # ANALYZE REGION - USE all_region_bps (not problematic_bps)
        region_hbs = self._filter_hbs_for_region(all_region_bps)
        bp_results = self.bp_analyzer.analyze(all_region_bps)
        hb_results = self.hb_analyzer.analyze(region_hbs, all_region_bps)
        
        region_score = (
            self.config.BASE_SCORE - 
            bp_results['penalty'] - 
            hb_results['penalty']
        )
        
        issue_density = self.calculate_issue_density(residues, all_region_bps, region_hbs)
        severity = self._classify_severity(region_score)
        dominant = self._identify_dominant_issues(bp_results, hb_results)
        detailed_issues = self._collect_all_issues_for_hotspot(all_region_bps, region_hbs)
        
        # CRITICAL FIX: Store ALL base pairs (not just {res_1, res_2})
        all_base_pairs_data = [
            {'res_1': bp['res_1'], 'res_2': bp['res_2']}
            for bp in all_region_bps  # ← ALL base pairs in region
        ]
        
        print(f"DEBUG: Creating hotspot with {len(all_base_pairs_data)} all_base_pairs and {num_problematic} problematic")

        
        return Hotspot(
            chain=chain,
            start_res=start_res,
            end_res=end_res,
            residue_count=len(residue_list),
            score=round(region_score, 1),
            severity=severity,
            issue_density=round(issue_density, 2),
            dominant_issues=dominant,
            base_pairs_affected=total_bps,
            hbonds_affected=len(region_hbs) if region_hbs is not None else 0,
            details={
                'base_pairs': bp_results,
                'hbonds': hb_results,
                'residues': residue_list,
                'detailed_issues': detailed_issues,
                'all_base_pairs': all_base_pairs_data,  # ← This is the key!
                'problematic_base_pairs_count': num_problematic,  # ← And this!
                'fraction_problematic': round(fraction_bad, 2),
                'problem_type_breakdown': dict(problem_types),
                'issue_statistics': {
                    'total_bases': len(residues),
                    'problematic_bp_count': num_problematic,
                    'problematic_bp_fraction': round(fraction_bad, 2),
                }
            }
        )
    def _filter_hbs_for_region(self, region_bps: list) -> pd.DataFrame:
        """Filter H-bonds for base pairs in region - OPTIMIZED."""
        if not self._pair_to_hbonds:
            return pd.DataFrame()
        
        # Collect all H-bonds for base pairs in this region
        region_hbonds = []
        for bp in region_bps:
            pair_key = tuple(sorted([bp['res_1'], bp['res_2']]))
            if pair_key in self._pair_to_hbonds:
                region_hbonds.extend(self._pair_to_hbonds[pair_key])
        
        if region_hbonds:
            return pd.DataFrame(region_hbonds)
        return pd.DataFrame()
    
    def calculate_issue_density(self, hotspot_residues: Set[int], region_bps: list, region_hbs: pd.DataFrame) -> float:
        """Calculate fraction of bases with problems in the hotspot."""
        problematic_residues = set()
        
        problematic_residues.update(self._find_residues_with_bp_geometry_issues(region_bps))
        problematic_residues.update(self._find_residues_with_hbond_geometry_issues(region_hbs, region_bps))
        problematic_residues.update(self._find_residues_with_hbond_count_issues(region_bps, region_hbs))
        problematic_residues.update(self._find_residues_with_zero_hbond_issues(region_bps))
        
        problematic_in_hotspot = problematic_residues.intersection(hotspot_residues)
        
        total_bases = len(hotspot_residues)
        if total_bases == 0:
            return 0.0
        
        bases_with_issues = len(problematic_in_hotspot)
        issue_density = bases_with_issues / total_bases
        
        return round(issue_density, 3)
    
    def _find_residues_with_bp_geometry_issues(self, region_bps: list) -> Set[int]:
        """Find residues involved in base pairs with geometry issues."""
        problematic_residues = set()
        
        for bp in region_bps:
            bp_issues = BasePairScoring.score_base_pair(bp, self.config)
            
            if bp_issues and any(bp_issues.values()):
                chain1, res1 = self._parse_residue(bp['res_1'])
                chain2, res2 = self._parse_residue(bp['res_2'])
                
                problematic_residues.add(res1)
                problematic_residues.add(res2)
        
        return problematic_residues

    def _find_residues_with_hbond_geometry_issues(self, region_hbs: pd.DataFrame, region_bps: list) -> Set[int]:
        """Find residues involved in H-bonds with geometry issues."""
        problematic_residues = set()
        
        if region_hbs is None or region_hbs.empty:
            return problematic_residues
        
        base_pair_set = ScoringUtils.build_base_pair_set(region_bps)
        
        for _, hb in region_hbs.iterrows():
            hb_issues = HBondScoring.score_hbond(hb, self.config)
            
            if any(hb_issues.values()):
                pair_key = tuple(sorted([hb['res_1'], hb['res_2']]))
                if pair_key in base_pair_set:
                    chain1, res1 = self._parse_residue(hb['res_1'])
                    chain2, res2 = self._parse_residue(hb['res_2'])
                    
                    problematic_residues.add(res1)
                    problematic_residues.add(res2)
        
        return problematic_residues

    def _find_residues_with_hbond_count_issues(self, region_bps: list, region_hbs: pd.DataFrame) -> Set[int]:
        """Find residues involved in base pairs with incorrect H-bond counts - OPTIMIZED."""
        problematic_residues = set()
        
        if region_hbs is None or region_hbs.empty:
            return problematic_residues
        
        for bp in region_bps:
            res_1, res_2 = bp['res_1'], bp['res_2']
            pair_key = tuple(sorted([res_1, res_2]))
            actual_count = self._pair_hbond_count.get(pair_key, 0)  # ← FAST LOOKUP
            
            hbond_score = bp.get('hbond_score', 0)
            if hbond_score == 0.0:
                continue
            
            is_incorrect, _ = HBondScoring.check_hbond_count(
                bp.get('bp_type', 'unknown'), actual_count, bp.get('lw', ''), self.config
            )
            
            if is_incorrect:
                chain1, res1 = self._parse_residue(res_1)
                chain2, res2 = self._parse_residue(res_2)
                
                problematic_residues.add(res1)
                problematic_residues.add(res2)
        
        return problematic_residues

    def _find_residues_with_zero_hbond_issues(self, region_bps: list) -> Set[int]:
        """Find residues involved in base pairs with zero H-bonds."""
        problematic_residues = set()
        
        for bp in region_bps:
            hbond_score = bp.get('hbond_score', 0)
            
            if hbond_score == 0.0:
                chain1, res1 = self._parse_residue(bp['res_1'])
                chain2, res2 = self._parse_residue(bp['res_2'])
                
                problematic_residues.add(res1)
                problematic_residues.add(res2)
        
        return problematic_residues

    def _count_hbonds_per_pair_in_region(self, region_bps: list) -> Dict[Tuple[str, str], int]:
        """Count H-bonds for each base pair in the region - OPTIMIZED."""
        pair_hbond_counts = {}
        
        for bp in region_bps:
            pair_key = tuple(sorted([bp['res_1'], bp['res_2']]))
            pair_hbond_counts[pair_key] = self._pair_hbond_count.get(pair_key, 0)
        
        return pair_hbond_counts
    
    def _parse_residue(self, res_id: str) -> Tuple[str, int]:
        """Parse residue ID into (chain, res_num)."""
        parts = res_id.split('-')
        return parts[0], int(parts[2].rstrip('-'))
    
    def _classify_severity(self, score: float) -> str:
        """Classify severity based on score."""
        if score < self.config.SEVERITY_CRITICAL:
            return "CRITICAL"
        elif score < self.config.SEVERITY_SEVERE:
            return "SEVERE"
        elif score < self.config.SEVERITY_MODERATE:
            return "MODERATE"
        else:
            return "MINOR"
        
        
    def _identify_dominant_issues(self, bp_results: Dict, hb_results: Dict) -> List[str]:
        """Identify the most prevalent issue types - INCLUDING H-BOND ISSUES."""
        issues = []
        threshold = self.config.ISSUE_DENSITY_MEDIUM
        
        # BASE PAIR GEOMETRY ISSUES
        if bp_results['misaligned_frac'] > threshold:
            issues.append(f"Misaligned pairs ({bp_results['misaligned_frac']:.0%})")
        if bp_results['twisted_frac'] > threshold:
            issues.append(f"Twisted pairs ({bp_results['twisted_frac']:.0%})")
        if bp_results['non_coplanar_frac'] > threshold:
            issues.append(f"Non-coplanar pairs ({bp_results['non_coplanar_frac']:.0%})")
        if bp_results['poor_hbond_frac'] > threshold:
            issues.append(f"Low DSSR quality ({bp_results['poor_hbond_frac']:.0%})")
        if bp_results.get('zero_hbond_frac', 0) > threshold:
            issues.append(f"No H-bonds detected ({bp_results['zero_hbond_frac']:.0%})")
        if bp_results.get('self_pairing_frac', 0) > 0.0:
            issues.append(f"Self-pairing ({bp_results['self_pairing_frac']:.0%})")
        
        # H-BOND GEOMETRY ISSUES - NOW INCLUDED!
        if hb_results['bad_distance_frac'] > threshold:
            issues.append(f"Bad H-bond distances ({hb_results['bad_distance_frac']:.0%})")
        if hb_results['bad_angles_frac'] > threshold:
            issues.append(f"Bad H-bond angles ({hb_results['bad_angles_frac']:.0%})")
        if hb_results['bad_dihedral_frac'] > threshold:
            issues.append(f"Poor H-bond dihedrals ({hb_results['bad_dihedral_frac']:.0%})")
        if hb_results['weak_quality_frac'] > threshold:
            issues.append(f"Weak H-bond quality ({hb_results['weak_quality_frac']:.0%})")
        if hb_results['incorrect_hbond_count_frac'] > threshold:
            issues.append(f"Wrong H-bond count ({hb_results['incorrect_hbond_count_frac']:.0%})")
        
        return issues if issues else ["Multiple minor issues"]
    
    
    
    def _stitch_hotspot_chains(self, hotspots: List[Hotspot]) -> List[Hotspot]:
        """Merge hotspots that are close together."""
        if not hotspots:
            return []

        stitched = []
        by_chain = defaultdict(list)
        
        for h in hotspots:
            if h is not None:
                by_chain[h.chain].append(h)
        
        for chain, chain_hotspots in by_chain.items():
            if not chain_hotspots:
                continue
                
            chain_hotspots.sort(key=lambda h: h.start_res)
            current = chain_hotspots[0]
            
            for nxt in chain_hotspots[1:]:
                gap = nxt.start_res - current.end_res
                allowed_gap = 1 if current.severity in ["CRITICAL", "SEVERE"] else 1
                    
                if gap <= allowed_gap:
                    merged_hotspot = self._merge_two_hotspots(current, nxt)
                    if merged_hotspot is not None:
                        current = merged_hotspot
                    else:
                        stitched.append(current)
                        current = nxt
                else:
                    stitched.append(current)
                    current = nxt
            
            if current is not None:
                stitched.append(current)
        
        return stitched
        
    def _collect_all_issues_for_hotspot(self, base_pairs: list, hbonds: pd.DataFrame) -> List[Dict]:
        """Collect ALL types of issues for a hotspot region."""
        combined_issues = {}
        
        self._add_basepair_geometry_to_combined(base_pairs, combined_issues)
        self._add_hbond_geometry_to_combined(hbonds, combined_issues)
        self._add_hbond_count_to_combined(base_pairs, hbonds, combined_issues)
        self._finalize_combined_issues(combined_issues)
        
        return list(combined_issues.values())
    
    def _add_basepair_geometry_to_combined(self, base_pairs: list, combined_issues: dict):
        """Add base pair geometry issues to combined_issues dictionary."""
        for bp in base_pairs:
            bp_issues = BasePairScoring.score_base_pair(bp, self.config)
            
            if bp_issues is None:
                continue

            if not any(bp_issues.values()):
                continue
            
            res_1 = bp['res_1']
            res_2 = bp['res_2']
            pair_key = tuple(sorted([res_1, res_2]))
            
            if pair_key not in combined_issues:
                combined_issues[pair_key] = self._create_empty_issue_entry(res_1, res_2, bp)
            
            entry = combined_issues[pair_key]
            
            if 'base_pair_geometry' not in entry['issue_types']:
                entry['issue_types'].append('base_pair_geometry')
            
            entry['geometry_parameters'] = {
                'shear': round(bp.get('shear', 0), 2),
                'stretch': round(bp.get('stretch', 0), 2),
                'stagger': round(bp.get('stagger', 0), 2),
                'buckle': round(bp.get('buckle', 0), 2),
                'propeller': round(bp.get('propeller', 0), 2),
                'opening': round(bp.get('opening', 0), 2)
            }
            
            hbond_score = bp.get('hbond_score', 0)
            entry['dssr_hbond_score'] = round(hbond_score, 2)
            
            entry['geometry_severity'] = self._classify_geometry_severity(bp)
            
            self._add_bp_specific_issues(bp_issues, entry)
                
    def _add_hbond_geometry_to_combined(self, hbonds: pd.DataFrame, combined_issues: dict):
        """Add H-bond geometry issues to the unified issue list."""
        if hbonds is None or hbonds.empty:
            return
        
        for _, hb in hbonds.iterrows():
            hb_issues = HBondScoring.score_hbond(hb, self.config)
            
            if not any(hb_issues.values()):
                continue
            
            res_1 = hb['res_1']
            res_2 = hb['res_2']
            pair_key = tuple(sorted([res_1, res_2]))
            
            if pair_key not in combined_issues:
                combined_issues[pair_key] = self._create_empty_issue_entry(res_1, res_2)
            
            entry = combined_issues[pair_key]
            
            if 'hbond_geometry' not in entry['issue_types']:
                entry['issue_types'].append('hbond_geometry')
            
            hb_issue_detail = self._create_hbond_issue_detail(hb, hb_issues)
            
            for issue_name in hb_issue_detail['issues']:
                if issue_name not in entry['specific_issues']:
                    entry['specific_issues'].append(issue_name)
            
            entry['issue_details'].append({
                'issue_category': 'hbond_geometry',
                'atoms': hb_issue_detail['atoms'],
                'issues': hb_issue_detail['issues'],
                'parameters': hb_issue_detail['parameters']
            })

    def _add_hbond_count_to_combined(self, base_pairs: list, hbonds: pd.DataFrame, combined_issues: dict):
        """Add H-bond count issues to the unified issue list - OPTIMIZED."""
        if hbonds is None or hbonds.empty:
            return
        
        for bp in base_pairs:
            res_1, res_2 = bp['res_1'], bp['res_2']
            pair_key = tuple(sorted([res_1, res_2]))
            actual_count = self._pair_hbond_count.get(pair_key, 0)  # ← FAST LOOKUP
            
            hbond_score = bp.get('hbond_score', 0)
            
            if hbond_score == 0.0 or actual_count == 0:
                continue
            
            is_incorrect, issue_desc = HBondScoring.check_hbond_count(
                bp.get('bp_type', 'unknown'), actual_count, bp.get('lw', ''), self.config
            )
            
            if is_incorrect and pair_key in combined_issues:
                entry = combined_issues[pair_key]
                
                if 'hbond_count' not in entry['issue_types']:
                    entry['issue_types'].append('hbond_count')
                
                if 'incorrect_hbond_count' not in entry['specific_issues']:
                    entry['specific_issues'].append('incorrect_hbond_count')
                
                entry['issue_details'].append({
                    'issue_category': 'hbond_count',
                    'actual_count': actual_count,
                    'expected_range': self.config.EXPECTED_HBOND_COUNTS.get(
                        bp.get('bp_type', 'unknown'), (0, 0, 0)
                    )[:2],
                    'description': issue_desc
                })
    
    def _create_empty_issue_entry(self, res_1: str, res_2: str, bp: dict = None) -> dict:
        """Create an empty issue entry for a residue pair."""
        res1_parts = res_1.split('-')
        res2_parts = res_2.split('-')
        
        return {
            'residues': f"{res_1} - {res_2}",
            'chain_1': res1_parts[0],
            'chain_2': res2_parts[0],
            'residue_1': int(res1_parts[2]),
            'residue_2': int(res2_parts[2]),
            'bp_type': bp.get('bp_type', 'unknown') if bp else 'unknown',
            'lw_notation': bp.get('lw', 'unknown') if bp else 'unknown',
            'issue_types': [],
            'specific_issues': [],
            'geometry_parameters': {},
            'geometry_severity': 'NONE',
            'issue_details': []
        }

    def _add_bp_specific_issues(self, bp_issues: dict, entry: dict):
        """Add specific base pair issues to entry."""
        if bp_issues['misaligned']:
            entry['specific_issues'].append('misaligned')
        if bp_issues['twisted']:
            entry['specific_issues'].append('twisted')
        if bp_issues['non_coplanar']:
            entry['specific_issues'].append('non_coplanar')
        if bp_issues['poor_hbond'] and not bp_issues['zero_hbond']:
            entry['specific_issues'].append('low_dssr_quality')
        if bp_issues['zero_hbond']:
            entry['specific_issues'].append('no_hbonds_detected')
        if bp_issues['self_pairing']:
            entry['specific_issues'].append('self_pairing')

    def _create_hbond_issue_detail(self, hb: pd.Series, hb_issues: dict) -> dict:
        """Create detailed H-bond issue entry."""
        issues = []
        
        if hb_issues['bad_distance']:
            issues.append('bad_hbond_distance')
        if hb_issues['bad_angles']:
            issues.append('bad_hbond_angles')
        if hb_issues['bad_dihedral']:
            issues.append('bad_hbond_dihedral')
        if hb_issues['weak_quality']:
            issues.append('weak_hbond_quality')
        
        return {
            'atoms': f"{hb['atom_1']} - {hb['atom_2']}",
            'issues': issues,
            'parameters': {
                'distance': round(hb['distance'], 2),
                'angle_1': round(hb['angle_1'], 1),
                'angle_2': round(hb['angle_2'], 1),
                'dihedral': round(hb['dihedral_angle'], 1),
                'quality_score': round(hb['score'], 2)
            }
        }
    
    def _count_hbonds_per_pair(self, hbonds: pd.DataFrame) -> dict:
        """Count H-bonds between each pair of residues."""
        pair_hbond_counts = defaultdict(int)
        for _, hb in hbonds.iterrows():
            pair_key = tuple(sorted([hb['res_1'], hb['res_2']]))
            pair_hbond_counts[pair_key] += 1
        return pair_hbond_counts

    def _finalize_combined_issues(self, combined_issues: dict):
        """Create combined issue_type string for each entry."""
        for entry in combined_issues.values():
            if len(entry['issue_types']) > 1:
                entry['issue_type'] = ' + '.join(entry['issue_types'])
            else:
                entry['issue_type'] = entry['issue_types'][0] if entry['issue_types'] else 'unknown'