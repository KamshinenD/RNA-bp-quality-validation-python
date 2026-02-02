"""Improved hotspot analyzer with better filtering and merging."""

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
    """Finds hotspots using connectivity-based approach with improved filtering."""
    
    def __init__(self, config, bp_analyzer, hb_analyzer):
        self.config = config
        self.bp_analyzer = bp_analyzer
        self.hb_analyzer = hb_analyzer
    
    def find_hotspots(self, basepair_data: list, hbond_data: pd.DataFrame) -> List[Hotspot]:
        """Find hotspots"""
        if not basepair_data:
            return []
        
        self._basepair_cache = basepair_data
        self._hbond_cache = hbond_data
        
        # Step 1: Score all residues
        residue_scores = self._score_all_residues(basepair_data, hbond_data)
        
        if not residue_scores:
            print("  ERROR: No residue scores calculated!")
            return []
        
        
        # I am using a percentile-based threshold. I am just checking and finding the worst 15% of residues (percentile-based)
        # Calculate distribution statistics
        scores = list(residue_scores.values())
        scores_sorted = sorted(scores)
        mean_score = sum(scores) / len(scores)
        std_dev = (sum((s - mean_score)**2 for s in scores) / len(scores))**0.5
        
        # Find the score below which worst 15% of residues fall(percentile-based)
        percentile_15 = scores_sorted[int(len(scores) * 0.05)]
        
        #Find the score below which residues are 1 standard deviation below the mean score
        statistical = mean_score - 1.0 * std_dev
        
        #Smart Threshold Selection - 3 way validation
        #out of the 3, Use the most conservative(max is the most strict) of score as threshold
        threshold = max(percentile_15, statistical, self.config.DAMAGE_THRESHOLD_FOR_HOTSPOTS)
        
        # Find damaged residues
        damaged_residues = self._find_damaged_residues(residue_scores, threshold)
        
        total_damaged = sum(len(res_set) for res_set in damaged_residues.values())

        if total_damaged == 0:
            print("  No damaged residues")
            return []
        
        # Continue with rest of pipeline
        connectivity = self._build_connectivity_graph(basepair_data)
        components = self._find_connected_components(damaged_residues, connectivity)
        
        min_size = getattr(self.config, 'MIN_HOTSPOT_RESIDUES', 5)
        
        
        hotspots = []
        for chain, residue_sets in components.items():
            for residues in residue_sets:
                print(f"    Component size: {len(residues)}")
                if len(residues) >= min_size:
                    hotspot = self._create_hotspot_from_residues(chain, residues)
                    if hotspot:
                        hotspots.append(hotspot)
                        print(f"      Created hotspot: {hotspot.chain}:{hotspot.start_res}-{hotspot.end_res}")
                    else:
                        print(f"      Failed to create hotspot from {len(residues)} residues")
        
        #print(f"\n  Hotspots before merging: {len(hotspots)}")
        
        # Rest of pipeline
        hotspots = self._merge_overlapping_hotspots(hotspots)
        #print(f"  After merging: {len(hotspots)}")
        
        hotspots = self._stitch_hotspot_chains(hotspots)
        #print(f"  After stitching: {len(hotspots)}")
        
        hotspots = self._filter_by_context(hotspots)
        
        hotspots.sort(key=lambda h: h.score)
        
        if not hotspots:
            print("\n  === NO HOTSPOTS - LIKELY CAUSES ===")
            print("    - Overall structure quality is good")
        return hotspots 

    # No overlap, save current and move to next
    #This allows long damaged regions with short clean patches (like 10–12 nucleotides) to remain one continuous hotspot, not artificially split.
    def _merge_overlapping_hotspots(self, hotspots: List[Hotspot]) -> List[Hotspot]:
        """Merge hotspots that overlap significantly on the same chain."""
        if not hotspots:
            return []
        
        # Group by chain
        by_chain = defaultdict(list)
        for h in hotspots:
            by_chain[h.chain].append(h)
        
        merged = []
        for chain, chain_hotspots in by_chain.items():
            # Sort by start position
            chain_hotspots.sort(key=lambda h: h.start_res)
            
            current = chain_hotspots[0]
            
            for next_hotspot in chain_hotspots[1:]:
                # Check if they overlap (allowing small gaps)
                overlap_threshold = 0.3  # 30% overlap
                
                overlap_start = max(current.start_res, next_hotspot.start_res)
                overlap_end = min(current.end_res, next_hotspot.end_res)
                overlap_size = max(0, overlap_end - overlap_start + 1)
                
                min_size = min(
                    current.end_res - current.start_res + 1,
                    next_hotspot.end_res - next_hotspot.start_res + 1
                )
                
                # Merge if significant overlap OR small gap between them
                # Adaptive gap tolerance based on severity and size
                gap = next_hotspot.start_res - current.end_res
                max_gap = 5

                if current.severity == "CRITICAL" or next_hotspot.severity == "CRITICAL":
                    max_gap = 2 #25
                elif current.severity == "SEVERE" or next_hotspot.severity == "SEVERE":
                    max_gap = 2 #20
                elif current.severity == "MODERATE" or next_hotspot.severity == "MODERATE":
                    max_gap = 1 #15
                else:
                    max_gap = 1 #10 
                
                if (overlap_size / min_size > overlap_threshold) or (0 < gap <= max_gap):
                    # Merge: expand current to include next_hotspot
                    current = self._merge_two_hotspots(current, next_hotspot)
                else:
                    # No overlap, save current and move to next
                    merged.append(current)
                    current = next_hotspot
            # Add the last one
            merged.append(current)
        
        return merged
    
    def _merge_two_hotspots(self, h1: Hotspot, h2: Hotspot) -> Hotspot:
        """Merge two overlapping hotspots into one."""
        # Combine residue ranges
        all_residues = set(range(h1.start_res, h1.end_res + 1))
        all_residues.update(range(h2.start_res, h2.end_res + 1))
        
        # Recreate hotspot from merged residue set
        return self._create_hotspot_from_residues(h1.chain, all_residues)
    
    def _filter_by_context(self, hotspots: List[Hotspot]) -> List[Hotspot]:
        """Filter hotspots based on structure quality context."""
        if not hotspots:
            return []
        
        # Only keep SEVERE and CRITICAL hotspots
        # (MODERATE and MINOR are not true "hotspots" in a good structure)
        #filtered = [h for h in hotspots if h.severity in ['SEVERE', 'CRITICAL']]
        filtered = [h for h in hotspots if h.severity != 'MINOR']
        
        return filtered
    
    def _score_all_residues(self, basepair_data: list, hbond_data: pd.DataFrame) -> Dict[Tuple[str, int], float]:
        """
        Score each residue using weighted penalties from:
        1. Base-pair geometry issues
        2. Hydrogen-bond geometry issues
        3. Incorrect H-bond counts per base pair
        """
        residue_penalties = defaultdict(lambda: {'penalty_sum': 0.0, 'bp_count': 0})

        # === BASE-PAIR GEOMETRY PENALTIES ===
        for bp in basepair_data:
            bp_issues = self._count_bp_issues(bp)
            bp_penalty = sum(self.config.PENALTY_WEIGHTS.get(issue, 0) for issue in bp_issues)

            for res_id in [bp['res_1'], bp['res_2']]:
                chain, res_num = self._parse_residue(res_id)
                key = (chain, res_num)
                residue_penalties[key]['penalty_sum'] += bp_penalty
                residue_penalties[key]['bp_count'] += 1

        # **BUILD BASE PAIR SET for H-bond filtering**
        base_pair_set = ScoringUtils.build_base_pair_set(basepair_data)

        # === HYDROGEN-BOND GEOMETRY PENALTIES ===
        if hbond_data is not None and not hbond_data.empty:
            for _, row in hbond_data.iterrows():
                
                # **FILTER 1: Skip non-base H-bonds**
                if not ScoringUtils.is_base_base_hbond(row['atom_1'], row['atom_2']):
                    continue
                
                # **FILTER 2: Skip adjacent pairs**
                if BasePairScoring.check_adjacent_pairing(row['res_1'], row['res_2']):
                    continue
                
                # **FILTER 3: Only H-bonds belonging to base pairs**
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

        # === H-BOND COUNT MISMATCH PENALTIES ===
        if hbond_data is not None and not hbond_data.empty:
            # Count all H-bonds between residue pairs
            pair_hbond_counts = defaultdict(int)
            for _, hb in hbond_data.iterrows():
                
                # **FILTER 1: Only count base-base H-bonds in base pairs**
                if not ScoringUtils.is_base_base_hbond(hb['atom_1'], hb['atom_2']):
                    continue
                
                # **FILTER 2: Skip adjacent pairs**
                if BasePairScoring.check_adjacent_pairing(hb['res_1'], hb['res_2']):
                    continue
                
                # **FILTER 3: Only in base pairs**
                pair_key = tuple(sorted([hb['res_1'], hb['res_2']]))
                if pair_key not in base_pair_set:
                    continue
                
                pair_hbond_counts[pair_key] += 1

            for bp in basepair_data:
                res_1, res_2 = bp['res_1'], bp['res_2']
                pair_key = tuple(sorted([res_1, res_2]))
                actual_count = pair_hbond_counts.get(pair_key, 0)
                
                # **SKIP if DSSR score is 0 (avoid redundancy)**
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

        # === COMPUTE FINAL WEIGHTED SCORE PER RESIDUE ===
        residue_scores = {}
        for key, data in residue_penalties.items():
            avg_penalty = data['penalty_sum'] / max(data['bp_count'], 1)
            score = max(0, self.config.BASE_SCORE - avg_penalty)
            residue_scores[key] = score

        return residue_scores


    
    def _count_bp_issues(self, bp: dict) -> List[str]:
        """Return list of specific issue types found in a base pair."""
        issues = []
        
        # Use shared scoring logic
        bp_issues_dict = BasePairScoring.score_base_pair(bp, self.config)
        
        if bp_issues_dict is None:
            return []
        
        # Convert boolean flags to issue keys
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
        # if bp_issues_dict['adjacent_pairing']:
        #     issues.append('adjacent_pairing')
        if bp_issues_dict['self_pairing']:
            issues.append('self_pairing')
        
        return issues
    
    
    def _count_hb_issues(self, hb: dict) -> list:
        """
        Identify specific hydrogen-bond issues in a single H-bond record.

        Returns:
            List[str]: list of issue keys corresponding to Config.PENALTY_WEIGHTS.
        """
        issues = []
        
        # Convert dict to pandas Series for compatibility with shared scoring
        hb_series = pd.Series(hb)
        
        # Use shared scoring logic
        hb_issues_dict = HBondScoring.score_hbond(hb_series, self.config)
        
        # Convert boolean flags to issue keys
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
        """
        Check whether a base pair has an incorrect number of H-bonds.
        Returns a list of issue keys (usually ['incorrect_hbond_count']) if mismatch detected.
        """
        issues = []
        bp_type = bp.get('bp_type', 'unknown')
        lw = bp.get('lw', '')
        
        # Use shared scoring logic
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
        """
        Build connectivity graph for residues.
        Connect residues that are:
        1. Base-paired together
        2. Sequential neighbors within SEQUENTIAL_NEIGHBOR_RANGE
        """
        graph = defaultdict(lambda: defaultdict(set))
        
        # Step 1: Add base pair connections
        for bp in basepair_data:
            chain1, res1 = self._parse_residue(bp['res_1'])
            chain2, res2 = self._parse_residue(bp['res_2'])
            
            # Connect base-paired residues (even across chains for now)
            graph[chain1][res1].add(res2)
            graph[chain2][res2].add(res1)
            
            # Also ensure the residues exist in their own chain's graph
            graph[chain1][res1]  # Initialize if not exists
            graph[chain2][res2]  # Initialize if not exists
        
        # Step 2: Add sequential connectivity MORE AGGRESSIVELY
        # This is critical to prevent fragmentation
        neighbor_range = self.config.SEQUENTIAL_NEIGHBOR_RANGE
        
        for chain in list(graph.keys()):  # Use list() to avoid dict modification during iteration
            all_residues = sorted(graph[chain].keys())
            
            # Connect each residue to its neighbors within ±neighbor_range
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
    
        # Identify contiguous segments within the hotspot
        segments = []
        current_segment = [residue_list[0]]
    
        for i in range(1, len(residue_list)):
            if residue_list[i] - residue_list[i-1] <= 3:  # Contiguous (allowing small gaps)
                current_segment.append(residue_list[i])
            else:
                # Gap found, save current segment and start new one
                segments.append((current_segment[0], current_segment[-1]))
                current_segment = [residue_list[i]]
    
        # Add the last segment
        segments.append((current_segment[0], current_segment[-1]))
    
        # CRITICAL FILTER: Check for huge gaps between residues
        # If any gap is > 30 positions, this is likely an artifact
        gaps = [residue_list[i+1] - residue_list[i] for i in range(len(residue_list)-1)]
        if gaps and max(gaps) > 50:
            # Residues are too sparse - connectivity artifact
            return None
    
        # Filter out very sparse hotspots
        span = end_res - start_res + 1
        density = len(residue_list) / span
        
        # Continuity boost: extend end points if nearby residues are also damaged
        if density < 0.25 and span > 40:
            # Expand range slightly if nearby residues are within 5 nt and moderately damaged
            expansion = 3
            start_res = max(1, start_res - expansion)
            end_res = end_res + expansion

    
        # If density < 5% AND span > 80, it's suspicious
        if density < 0.02 and span > 100:
            return None
    
    
        # === FIX: EXPAND TO INCLUDE BASE PAIR PARTNERS ===
        expanded_residues = set()
        
        for res_num in residues:
            expanded_residues.add((chain, res_num))


        
        # Find all base pairs involving these residues
        for bp in self._basepair_cache:
            chain1, res1 = self._parse_residue(bp['res_1'])
            chain2, res2 = self._parse_residue(bp['res_2'])
            
            # If either side is in hotspot, include BOTH sides
            if (chain1 == chain and res1 in residues) or (chain2 == chain and res2 in residues):
                expanded_residues.add((chain1, res1))
                expanded_residues.add((chain2, res2))
    
        #Use expanded residues for analysis
        region_bps = self._filter_bps_by_expanded_residues(expanded_residues)
        region_hbs = self._filter_hbs_by_expanded_residues(expanded_residues)

    
    
        # Get all base pairs involving these residues
        # region_bps = self._filter_bps_by_residues(chain, residues)
        # region_hbs = self._filter_hbs_by_residues(chain, residues)
    
        if not region_bps:
            return None
        # Analyze the region
        bp_results = self.bp_analyzer.analyze(region_bps)
        hb_results = self.hb_analyzer.analyze(region_hbs, region_bps)
        
        detailed_issues = self._collect_all_issues_for_hotspot(region_bps, region_hbs)

    
        region_score = (
            self.config.BASE_SCORE - 
            bp_results['penalty'] - 
            hb_results['penalty']
        )
    
        issue_density = self.calculate_issue_density(bp_results, hb_results, region_bps)
    
        severity = self._classify_severity(region_score)
        dominant = self._identify_dominant_issues(bp_results, hb_results)
        all_base_pairs_data = [
            {'res_1': bp['res_1'], 'res_2': bp['res_2']}
            for bp in region_bps
        ]
        
        return Hotspot(
            chain=chain,
            start_res=start_res,
            end_res=end_res,
            residue_count=len(residue_list),
            score=round(region_score, 1),
            severity=severity,
            issue_density=round(issue_density, 2),
            dominant_issues=dominant,
            base_pairs_affected=len(region_bps),
            hbonds_affected=len(region_hbs) if region_hbs is not None else 0,
            details={
                'base_pairs': bp_results,
                'hbonds': hb_results,
                'segments': segments,
                'residues': residue_list,
                'detailed_issues': detailed_issues,
                'all_base_pairs': all_base_pairs_data
            }
        )
    
    
    def calculate_weighted_issues(self, bp_results: Dict, hb_results: Dict) -> float:
        """Calculate weighted issue count with safe dictionary access."""
        
        weighted_total = 0.0
        
        # Safe access to base pair stats
        bp_stats = bp_results.get('stats', {})
        weighted_total += bp_stats.get('misaligned', 0) * self.config.PENALTY_WEIGHTS['misaligned_pairs']
        weighted_total += bp_stats.get('twisted', 0) * self.config.PENALTY_WEIGHTS['twisted_pairs']
        weighted_total += bp_stats.get('non_coplanar', 0) * self.config.PENALTY_WEIGHTS['non_coplanar_pairs']
        weighted_total += bp_stats.get('poor_hbond', 0) * self.config.PENALTY_WEIGHTS['poor_hbond_pairs']
        weighted_total += bp_stats.get('zero_hbond', 0) * self.config.PENALTY_WEIGHTS['zero_hbond_pairs']  # Safe access
        weighted_total += bp_stats.get('self_pairing', 0) * self.config.PENALTY_WEIGHTS['self_pairing']    # Safe access
        # weighted_total += bp_stats.get('adjacent_pairing', 0) * self.config.PENALTY_WEIGHTS['adjacent_pairing']  # Safe access
        
        # Safe access to H-bond stats
        hb_stats = hb_results.get('stats', {})
        weighted_total += hb_stats.get('incorrect_hbond_count', 0) * self.config.PENALTY_WEIGHTS['incorrect_hbond_count']
        weighted_total += hb_stats.get('bad_distance', 0) * self.config.PENALTY_WEIGHTS['bad_hbond_distance']
        weighted_total += hb_stats.get('bad_angles', 0) * self.config.PENALTY_WEIGHTS['bad_hbond_angles']
        weighted_total += hb_stats.get('bad_dihedral', 0) * self.config.PENALTY_WEIGHTS['bad_hbond_dihedrals']
        weighted_total += hb_stats.get('weak_quality', 0) * self.config.PENALTY_WEIGHTS['weak_hbond_quality']
        
        return weighted_total
    
    
    def calculate_issue_density(self, bp_results: Dict, hb_results: Dict, region_bps: list) -> float:
        """Calculate simple issue density (issues per base pair)."""
        
        total_issues = 0
        
        # Count base pair issues
        bp_stats = bp_results.get('stats', {})
        total_issues += bp_stats.get('misaligned', 0)
        total_issues += bp_stats.get('twisted', 0)
        total_issues += bp_stats.get('non_coplanar', 0)
        total_issues += bp_stats.get('poor_hbond', 0)
        total_issues += bp_stats.get('zero_hbond', 0)
        total_issues += bp_stats.get('self_pairing', 0)
        # total_issues += bp_stats.get('adjacent_pairing', 0)
        
        # Count H-bond issues (normalize since they're per H-bond, not per base pair)
        hb_stats = hb_results.get('stats', {})
        total_hb_issues = (
            hb_stats.get('incorrect_hbond_count', 0) +
            hb_stats.get('bad_distance', 0) +
            hb_stats.get('bad_angles', 0) +
            hb_stats.get('bad_dihedral', 0) +
            hb_stats.get('weak_quality', 0)
        )
        
        # Normalize H-bond issues to base pair scale
        # Each base pair typically has 2-3 H-bonds
        normalized_hb_issues = total_hb_issues / 2.5  # Adjust this factor as needed
        
        total_issues += normalized_hb_issues
        
        # Calculate density (issues per base pair)
        if region_bps:
            issue_density = total_issues / len(region_bps)
        else:
            issue_density = 0.0
        
        return round(issue_density, 2)

    def _filter_bps_by_residues(self, chain: str, residues: Set[int]) -> list:
        """Filter base pairs involving specific residues."""
        def in_region(res_id):
            c, r = self._parse_residue(res_id)
            return c == chain and r in residues
        
        return [bp for bp in self._basepair_cache 
                if in_region(bp['res_1']) or in_region(bp['res_2'])]
    
    def _filter_hbs_by_residues(self, chain: str, residues: Set[int]) -> pd.DataFrame:
        """Filter H-bonds involving specific residues."""
        if self._hbond_cache is None or len(self._hbond_cache) == 0:
            return pd.DataFrame()
        
        def in_region(res_id):
            c, r = self._parse_residue(res_id)
            return c == chain and r in residues
        
        mask = (self._hbond_cache['res_1'].apply(in_region) | 
                self._hbond_cache['res_2'].apply(in_region))
        return self._hbond_cache[mask]
    
    def _filter_bps_by_expanded_residues(self, residues: Set[Tuple[str, int]]) -> list:
        """Filter base pairs where EITHER residue is in the expanded set."""
        def in_region(res_id):
            chain, res_num = self._parse_residue(res_id)
            return (chain, res_num) in residues  # Check BOTH chain and number!
        
        return [bp for bp in self._basepair_cache 
                if in_region(bp['res_1']) or in_region(bp['res_2'])]



    def _filter_hbs_by_expanded_residues(self, residues: Set[Tuple[str, int]]) -> pd.DataFrame:
        """Filter H-bonds where EITHER residue is in the expanded set AND they form a base pair."""
        if self._hbond_cache is None or len(self._hbond_cache) == 0:
            return pd.DataFrame()
        
        # **Build base pair set**
        base_pair_set = ScoringUtils.build_base_pair_set(self._basepair_cache)

        def in_region(res_id):
            chain, res_num = self._parse_residue(res_id)
            return (chain, res_num) in residues  # Check BOTH chain and number!
        
        def is_base_pair_hbond(row):
            pair_key = tuple(sorted([row['res_1'], row['res_2']]))
            return pair_key in base_pair_set
        
        def is_not_adjacent(row):
            """Check that residues are NOT adjacent."""
            return not BasePairScoring.check_adjacent_pairing(row['res_1'], row['res_2'])


        
        # **Filter 1: Residues in region**
        mask = (self._hbond_cache['res_1'].apply(in_region) | 
                self._hbond_cache['res_2'].apply(in_region))
        filtered = self._hbond_cache[mask]

        # **FILTER 2: Only base-base H-bonds**
        if not filtered.empty:
            base_mask = filtered.apply(
                lambda row: ScoringUtils.is_base_base_hbond(row['atom_1'], row['atom_2']),
                axis=1
            )
            filtered = filtered[base_mask]
            
        # **Filter 3: NOT adjacent pairs**
        if not filtered.empty:
            adjacent_mask = filtered.apply(is_not_adjacent, axis=1)
            filtered = filtered[adjacent_mask]
        
        # **Filter 4: Only H-bonds that belong to base pairs**
        if not filtered.empty:
            bp_mask = filtered.apply(is_base_pair_hbond, axis=1)
            filtered = filtered[bp_mask]
            
    
        return filtered
    
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
        """Identify the most prevalent issue types."""
        issues = []
        threshold = self.config.ISSUE_DENSITY_MEDIUM
        
        if bp_results['misaligned_frac'] > threshold:
            issues.append(f"Misaligned pairs ({bp_results['misaligned_frac']:.0%})")
        if bp_results['twisted_frac'] > threshold:
            issues.append(f"Twisted pairs ({bp_results['twisted_frac']:.0%})")
        if bp_results['non_coplanar_frac'] > threshold:
            issues.append(f"Non-coplanar pairs ({bp_results['non_coplanar_frac']:.0%})")
        if bp_results['poor_hbond_frac'] > threshold:
            issues.append(f"Poor H-bonding ({bp_results['poor_hbond_frac']:.0%})")
        if bp_results.get('zero_hbond_frac', 0) > threshold:
            issues.append(f"Zero H-bonding ({bp_results['zero_hbond_frac']:.0%})")
        # if bp_results.get('adjacent_pairing_frac', 0) > threshold:
        #     issues.append(f"Adjacent pairing ({bp_results['adjacent_pairing_frac']:.0%})")
        if bp_results.get('self_pairing_frac', 0) > 0.0:
            issues.append(f"Self-pairing ({bp_results['self_pairing_frac']:.0%})")



        
        if hb_results['incorrect_hbond_count_frac'] > threshold:
            issues.append(f"Wrong H-bond count ({hb_results['incorrect_hbond_count_frac']:.0%})")
        if hb_results['bad_distance_frac'] > threshold:
            issues.append(f"Bad H-bond distances ({hb_results['bad_distance_frac']:.0%})")
        if hb_results['bad_angles_frac'] > threshold:
            issues.append(f"Bad H-bond angles ({hb_results['bad_angles_frac']:.0%})")
        
        return issues if issues else ["Multiple minor issues"]
    
    
    def _stitch_hotspot_chains(self, hotspots: List[Hotspot]) -> List[Hotspot]:
        """Merge hotspots that are close together (<20 nt apart) even if not overlapping."""
        if not hotspots:
            return []
    
        stitched = []
        by_chain = defaultdict(list)
        for h in hotspots:
            by_chain[h.chain].append(h)
    
        for chain, chain_hotspots in by_chain.items():
            chain_hotspots.sort(key=lambda h: h.start_res)
            current = chain_hotspots[0]
            for nxt in chain_hotspots[1:]:
                gap = nxt.start_res - current.end_res
                allowed_gap = 20  # Default gap tolerance
                
                if "CRITICAL" in (current.severity, nxt.severity):
                    allowed_gap = 1 #15
                elif "SEVERE" in (current.severity, nxt.severity):
                    allowed_gap = 1 #10
                else:
                    allowed_gap = 1 #6
                    
                if gap <= allowed_gap:  # allow small clean patches
                    current = self._merge_two_hotspots(current, nxt)
                else:
                    stitched.append(current)
                    current = nxt
            stitched.append(current)
        return stitched
    
    
    
    def _collect_all_issues_for_hotspot(self, base_pairs: list, hbonds: pd.DataFrame) -> List[Dict]:
        """
        Collect ALL types of issues for a hotspot region.
        COMBINES issues for the same base pair into single entries.
        
        Returns:
            List of dictionaries with combined issues per base pair.
        """
        combined_issues = {}
        
        # Step 1: Add base pair geometry issues
        self._add_basepair_geometry_to_combined(base_pairs, combined_issues)
        
        # Step 2: Add H-bond geometry issues
        self._add_hbond_geometry_to_combined(hbonds, combined_issues)
        
        # Step 3: Add H-bond count issues
        self._add_hbond_count_to_combined(base_pairs, hbonds, combined_issues)
        
        # Step 4: Finalize entries (create issue_type string)
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
            
            # Create entry if needed
            if pair_key not in combined_issues:
                combined_issues[pair_key] = self._create_empty_issue_entry(res_1, res_2, bp)
            
            entry = combined_issues[pair_key]
            
            # Add base pair geometry type
            entry['issue_types'].append('base_pair_geometry')
            
            # Add geometry parameters
            entry['geometry_parameters'] = {
                'shear': round(bp.get('shear', 0), 2),
                'stretch': round(bp.get('stretch', 0), 2),
                'stagger': round(bp.get('stagger', 0), 2),
                'buckle': round(bp.get('buckle', 0), 2),
                'propeller': round(bp.get('propeller', 0), 2),
                'opening': round(bp.get('opening', 0), 2)
            }
            
            # Add specific issues
            self._add_bp_specific_issues(bp_issues, entry)
            
            # Add DSSR H-bond score
            hbond_score = bp.get('hbond_score', 0)
            if hbond_score > 0:
                entry['dssr_hbond_score'] = round(hbond_score, 2)
                
    
    def _add_hbond_geometry_to_combined(self, hbonds: pd.DataFrame, combined_issues: dict):
        """Add H-bond geometry issues to combined_issues dictionary."""
        if hbonds is None or hbonds.empty:
            return
        
        for _, hb in hbonds.iterrows():
            
            #FILTER: Skip non-base H-bonds** THIS IS ALREADY HANDLED IN _fiter_hbs_by_expanded_residues
            # if not ScoringUtils.is_base_base_hbond(hb['atom_1'], hb['atom_2']):
            #     continue
            hb_issues = HBondScoring.score_hbond(hb, self.config)
            
            if not any(hb_issues.values()):
                continue
            
            res_1 = hb['res_1']
            res_2 = hb['res_2']
            pair_key = tuple(sorted([res_1, res_2]))
            
            # Create entry if needed (for H-bonds between non-paired residues)
            if pair_key not in combined_issues:
                combined_issues[pair_key] = self._create_empty_issue_entry(res_1, res_2)
            
            entry = combined_issues[pair_key]
            
            # Add H-bond geometry type
            if 'hbond_geometry' not in entry['issue_types']:
                entry['issue_types'].append('hbond_geometry')
            
            # Create H-bond issue detail
            hb_detail = self._create_hbond_issue_detail(hb, hb_issues)
            entry['hbond_issues'].append(hb_detail)


    def _add_hbond_count_to_combined(self, base_pairs: list, hbonds: pd.DataFrame, combined_issues: dict):
        """Add H-bond count issues to combined_issues dictionary."""
        if hbonds is None or hbonds.empty:
            return
        
        # Count H-bonds per base pair
        pair_hbond_counts = self._count_hbonds_per_pair(hbonds)
        
        # Check each base pair for incorrect count
        for bp in base_pairs:
            res_1, res_2 = bp['res_1'], bp['res_2']
            pair_key = tuple(sorted([res_1, res_2]))
            actual_count = pair_hbond_counts.get(pair_key, 0)
            
            # **SKIP if DSSR already flagged zero H-bonds to avoid redundancy**
            hbond_score = bp.get('hbond_score', 0)
            if hbond_score == 0.0:
                continue
            
            is_incorrect, issue_desc = HBondScoring.check_hbond_count(
                bp.get('bp_type', 'unknown'), actual_count, bp.get('lw', ''), self.config
            )
            
            if is_incorrect and pair_key in combined_issues:
                entry = combined_issues[pair_key]
                
                # Add H-bond count type
                if 'hbond_count' not in entry['issue_types']:
                    entry['issue_types'].append('hbond_count')
                
                entry['specific_issues'].append('incorrect_hbond_count')
                entry['hbond_count_details'] = {
                    'actual_count': actual_count,
                    'expected_range': self.config.EXPECTED_HBOND_COUNTS.get(
                        bp.get('bp_type', 'unknown'), (0, 0, 0)
                    )[:2],
                    'issue_description': issue_desc
                }
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
            'hbond_issues': [],
        }

    
    def _add_bp_specific_issues(self, bp_issues: dict, entry: dict):
        """Add specific base pair issues to entry."""
        if bp_issues['misaligned']:
            entry['specific_issues'].append('misaligned')
        if bp_issues['twisted']:
            entry['specific_issues'].append('twisted')
        if bp_issues['non_coplanar']:
            entry['specific_issues'].append('non_coplanar')
        if bp_issues['poor_hbond']:
            entry['specific_issues'].append('poor_hbond_geometry')
        if bp_issues['zero_hbond']:
            entry['specific_issues'].append('zero_hbond_geometry')
        # if bp_issues['adjacent_pairing']:
        #     entry['specific_issues'].append('adjacent_pairing')
        if bp_issues['self_pairing']:
            entry['specific_issues'].append('self_pairing')


    def _create_hbond_issue_detail(self, hb: pd.Series, hb_issues: dict) -> dict:
        """Create detailed H-bond issue entry."""
        hb_detail = {
            'atoms': f"{hb['atom_1']} - {hb['atom_2']}",
            'issues': [],
            'parameters': {
                'distance': round(hb['distance'], 2),
                'angle_1': round(hb['angle_1'], 1),
                'angle_2': round(hb['angle_2'], 1),
                'dihedral': round(hb['dihedral_angle'], 1),
                'quality_score': round(hb['score'], 2)
            }
        }
        
        # Add specific H-bond issues
        if hb_issues['bad_distance']:
            hb_detail['issues'].append('bad_hbond_distance')
        if hb_issues['bad_angles']:
            hb_detail['issues'].append('bad_hbond_angles')
        if hb_issues['bad_dihedral']:
            hb_detail['issues'].append('bad_hbond_dihedral')
        if hb_issues['weak_quality']:
            hb_detail['issues'].append('weak_hbond_quality')
        
        return hb_detail
    
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


    
    
    
    