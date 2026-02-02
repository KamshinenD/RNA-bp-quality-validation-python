#!/usr/bin/env python3
"""
G-Quadruplex Detection Script

Detects G-quadruplexes (G-quads) in RNA structures using multi-tier validation:
- Tier 1: Structural connectivity (4-cycles of G residues)
- Tier 2: Hoogsteen edge type validation
- Tier 3: Geometric planarity validation
- Tier 4: H-bond pattern validation (optional)
- Tier 5: Stacking detection (for multi-tetrad G-quads)
"""

import json
import csv
from pathlib import Path
from collections import defaultdict, Counter
from typing import List, Dict, Set, Tuple, Optional
import itertools


class GQuadDetector:
    """Detects G-quadruplexes in RNA structures."""
    
    # Hoogsteen edge types (Leontis-Westhof notation)
    # Includes: tHS, cHS (trans/cis Hoogsteen-Sugar), tSH, cSH (trans/cis Sugar-Hoogsteen)
    #           tHH, cHH (trans/cis Hoogsteen-Hoogsteen)
    #           tWH, cWH (trans/cis Watson-Hoogsteen), tHW, cHW (trans/cis Hoogsteen-Watson)
    HOOGSTEEN_EDGES = {'tHS', 'cHS', 'tSH', 'cSH', 'tHH', 'cHH', 'tWH', 'cWH', 'tHW', 'cHW'}
    
    # G-quad specific H-bond atom pairs
    GQUAD_HBOND_PAIRS = {
        ('O6', 'N1'), ('N1', 'O6'),
        ('O6', 'N2'), ('N2', 'O6'),
        ('N2', 'N7'), ('N7', 'N2')
    }
    
    # Geometric thresholds
    MAX_BUCKLE = 25.0  # degrees (planarity threshold)
    MAX_PROPELLER = 20.0  # degrees (planarity threshold)
    
    def __init__(self):
        self.basepairs_dir = Path('data/basepairs')
        self.hbonds_dir = Path('data/hbonds')
    
    def find_gg_pairs(self, basepairs: List[Dict]) -> List[Dict]:
        """Find all G-G base pairs in the structure."""
        gg_pairs = []
        for bp in basepairs:
            if bp.get('bp_type') == 'G-G':
                gg_pairs.append({
                    'res_1': bp.get('res_1', ''),
                    'res_2': bp.get('res_2', ''),
                    'lw': bp.get('lw', ''),
                    'hbond_score': bp.get('hbond_score', 0),
                    'buckle': bp.get('buckle', 0),
                    'propeller': bp.get('propeller', 0),
                    'shear': bp.get('shear', 0),
                    'stretch': bp.get('stretch', 0),
                    'stagger': bp.get('stagger', 0)
                })
        return gg_pairs
    
    def build_connectivity_graph(self, gg_pairs: List[Dict]) -> Dict[str, Set[str]]:
        """Build connectivity graph where nodes are G residues and edges are G-G pairs."""
        graph = defaultdict(set)
        for pair in gg_pairs:
            res_1 = pair['res_1']
            res_2 = pair['res_2']
            graph[res_1].add(res_2)
            graph[res_2].add(res_1)
        return dict(graph)
    
    def find_4_cycles(self, graph: Dict[str, Set[str]]) -> List[Tuple[str, str, str, str]]:
        """
        Find all 4-cycles in the graph (potential G-tetrads).
        Uses a more efficient approach: find nodes with 2+ connections, then search for cycles.
        """
        cycles = []
        
        # Filter to nodes with at least 2 connections (required for cycles)
        candidate_nodes = {node: conns for node, conns in graph.items() if len(conns) >= 2}
        
        if len(candidate_nodes) < 4:
            return []
        
        nodes = list(candidate_nodes.keys())
        
        # Try all combinations of 4 nodes
        for combo in itertools.combinations(nodes, 4):
            # Check if this forms a 4-cycle
            # Each node should connect to exactly 2 others in the cycle
            connections = {node: len(graph[node] & set(combo)) for node in combo}
            
            # Valid 4-cycle: each node connects to exactly 2 others
            if all(count == 2 for count in connections.values()):
                # Verify it's actually a cycle (not just 4 nodes with 2 connections each)
                # Build subgraph and check if it forms a cycle
                subgraph = {node: graph[node] & set(combo) for node in combo}
                if self._is_valid_4_cycle(subgraph, combo):
                    cycles.append(tuple(sorted(combo)))
        
        # Remove duplicates (same cycle, different order)
        return list(set(cycles))
    
    def find_g_clusters(self, graph: Dict[str, Set[str]], gg_pairs: List[Dict]) -> List[List[str]]:
        """
        Alternative approach: Find clusters of G residues that might form G-quads.
        Looks for groups of 4+ G residues with high connectivity.
        """
        clusters = []
        visited = set()
        
        # Find all nodes with 2+ connections
        candidate_nodes = {node for node, conns in graph.items() if len(conns) >= 2}
        
        for node in candidate_nodes:
            if node in visited:
                continue
            
            # BFS to find connected cluster
            cluster = []
            queue = [node]
            visited.add(node)
            
            while queue:
                current = queue.pop(0)
                cluster.append(current)
                
                for neighbor in graph.get(current, set()):
                    if neighbor not in visited and neighbor in candidate_nodes:
                        visited.add(neighbor)
                        queue.append(neighbor)
            
            # Only consider clusters with 4+ nodes
            if len(cluster) >= 4:
                clusters.append(cluster)
        
        return clusters
    
    def _is_valid_4_cycle(self, subgraph: Dict[str, Set[str]], nodes: Tuple[str, ...]) -> bool:
        """Check if the subgraph forms a valid 4-cycle."""
        # Try to traverse the cycle
        if len(nodes) != 4:
            return False
        
        # Start from first node, try to visit all 4 nodes and return to start
        start = nodes[0]
        visited = {start}
        current = start
        path = [start]
        
        for _ in range(3):  # Need to visit 3 more nodes
            neighbors = subgraph[current] - visited
            if not neighbors:
                return False
            current = next(iter(neighbors))
            visited.add(current)
            path.append(current)
        
        # Check if last node connects back to start
        if start in subgraph[current]:
            return True
        
        return False
    
    def validate_hoogsteen(self, cycle: Tuple[str, ...], gg_pairs: List[Dict]) -> Tuple[int, int]:
        """
        Validate Hoogsteen edge types in the cycle.
        Returns: (hoogsteen_count, total_edges)
        """
        cycle_set = set(cycle)
        cycle_edges = []
        
        for pair in gg_pairs:
            if pair['res_1'] in cycle_set and pair['res_2'] in cycle_set:
                cycle_edges.append(pair)
        
        hoogsteen_count = sum(1 for edge in cycle_edges if edge['lw'] in self.HOOGSTEEN_EDGES)
        return hoogsteen_count, len(cycle_edges)
    
    def validate_geometry(self, cycle: Tuple[str, ...], gg_pairs: List[Dict]) -> Tuple[bool, Dict]:
        """
        Validate geometric constraints (planarity).
        Returns: (is_planar, geometry_stats)
        """
        cycle_set = set(cycle)
        cycle_pairs = [p for p in gg_pairs 
                      if p['res_1'] in cycle_set and p['res_2'] in cycle_set]
        
        if not cycle_pairs:
            return False, {}
        
        avg_buckle = sum(abs(p['buckle']) for p in cycle_pairs) / len(cycle_pairs)
        avg_propeller = sum(abs(p['propeller']) for p in cycle_pairs) / len(cycle_pairs)
        
        is_planar = (avg_buckle < self.MAX_BUCKLE and 
                    avg_propeller < self.MAX_PROPELLER)
        
        return is_planar, {
            'avg_buckle': round(avg_buckle, 2),
            'avg_propeller': round(avg_propeller, 2),
            'max_buckle': round(max(abs(p['buckle']) for p in cycle_pairs), 2),
            'max_propeller': round(max(abs(p['propeller']) for p in cycle_pairs), 2)
        }
    
    def validate_hbonds(self, cycle: Tuple[str, ...], pdb_id: str) -> Tuple[int, int]:
        """
        Validate G-quad specific H-bond patterns.
        Returns: (gquad_hbonds, total_hbonds)
        """
        hbond_file = self.hbonds_dir / f"{pdb_id}.csv"
        if not hbond_file.exists():
            return 0, 0
        
        cycle_set = set(cycle)
        gquad_hbonds = 0
        total_hbonds = 0
        
        try:
            with open(hbond_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    res_1 = row.get('res_1', '')
                    res_2 = row.get('res_2', '')
                    atom_1 = row.get('atom_1', '')
                    atom_2 = row.get('atom_2', '')
                    
                    # Check if both residues are in the cycle
                    if res_1 in cycle_set and res_2 in cycle_set:
                        total_hbonds += 1
                        if (atom_1, atom_2) in self.GQUAD_HBOND_PAIRS:
                            gquad_hbonds += 1
        except Exception:
            pass
        
        return gquad_hbonds, total_hbonds
    
    def detect_stacking(self, cycles: List[Tuple[str, ...]]) -> Dict[Tuple[str, ...], int]:
        """
        Detect stacked G-tetrads.
        Two tetrads are stacked if they have the same residue numbers but different chains.
        Returns: dict mapping cycle to number of stacked tetrads
        """
        stacking = {}
        
        # Extract residue numbers from each cycle
        def get_residue_numbers(cycle):
            """Extract residue numbers from cycle (e.g., ['A-G-4-', 'B-G-4-'] -> {'4'} or {'4', '5'} if mixed)"""
            res_nums = set()
            for res in cycle:
                parts = res.split('-')
                if len(parts) >= 3:
                    res_nums.add(parts[2])
            return res_nums
        
        # Group cycles by residue numbers
        cycles_by_resnum = defaultdict(list)
        for cycle in cycles:
            res_nums = get_residue_numbers(cycle)
            # Use sorted tuple as key for grouping
            key = tuple(sorted(res_nums))
            cycles_by_resnum[key].append(cycle)
        
        # Count stacking: cycles with same residue numbers are stacked
        for cycle in cycles:
            res_nums = get_residue_numbers(cycle)
            key = tuple(sorted(res_nums))
            # Count how many cycles share the same residue numbers (stacked tetrads)
            stacking[cycle] = len(cycles_by_resnum[key])
        
        return stacking
    
    def calculate_confidence_score(self, hoogsteen_ratio: float, is_planar: bool,
                                   has_hbonds: bool, stacking_count: int) -> float:
        """Calculate confidence score based on validation tiers."""
        score = 1.0  # Base score for 4-cycle detection (Tier 1)
        
        # Tier 2: Hoogsteen validation
        if hoogsteen_ratio >= 0.67:  # At least 4/6 edges
            score += 0.3
        elif hoogsteen_ratio >= 0.5:  # At least 3/6 edges
            score += 0.15
        
        # Tier 3: Geometric validation
        if is_planar:
            score += 0.2
        
        # Tier 4: H-bond validation
        if has_hbonds:
            score += 0.3
        
        # Tier 5: Stacking detection
        if stacking_count > 1:
            score += 0.2
        
        return round(score, 2)
    
    def find_g_residues_from_hbonds(self, pdb_id: str) -> Set[str]:
        """Find all G residues from H-bond data."""
        hbond_file = self.hbonds_dir / f"{pdb_id}.csv"
        g_residues = set()
        
        if hbond_file.exists():
            try:
                with open(hbond_file, 'r') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        res_1 = row.get('res_1', '')
                        res_2 = row.get('res_2', '')
                        res_type_1 = row.get('res_type_1', '')
                        res_type_2 = row.get('res_type_2', '')
                        
                        if res_type_1 == 'RNA' and 'G-' in res_1:
                            g_residues.add(res_1)
                        if res_type_2 == 'RNA' and 'G-' in res_2:
                            g_residues.add(res_2)
            except Exception:
                pass
        
        return g_residues
    
    def detect_parallel_gquad(self, gg_pairs: List[Dict], pdb_id: str) -> List[Dict]:
        """
        Detect parallel G-quads where G-G pairs are loop connections, not tetrad connections.
        Looks for patterns where G residues at specific positions form tetrads.
        """
        # Extract all G residues and their positions
        g_residues = set()
        g_positions = {}  # residue -> position number
        
        for pair in gg_pairs:
            res_1 = pair['res_1']
            res_2 = pair['res_2']
            g_residues.add(res_1)
            g_residues.add(res_2)
            
            # Extract position numbers
            parts_1 = res_1.split('-')
            parts_2 = res_2.split('-')
            if len(parts_1) >= 3:
                try:
                    g_positions[res_1] = int(parts_1[2])
                except:
                    pass
            if len(parts_2) >= 3:
                try:
                    g_positions[res_2] = int(parts_2[2])
                except:
                    pass
        
        if len(g_residues) < 8:  # Need at least 8 G residues for 2 tetrads
            return []
        
        # Group G residues by chain
        chains = defaultdict(list)
        for g in g_residues:
            parts = g.split('-')
            if len(parts) >= 3:
                chain = parts[0]
                pos = g_positions.get(g)
                if pos is not None:
                    chains[chain].append((g, pos))
        
        detected_quads = []
        
        # For each chain with 8+ G residues, try to find tetrad pattern
        for chain, gs in chains.items():
            if len(gs) < 8:
                continue
            
            # Sort by position
            sorted_gs = sorted(gs, key=lambda x: x[1])
            positions = [g[1] for g in sorted_gs]
            
            # Look for pattern: parallel G-quad with two stacked tetrads
            # Pattern: positions like [1,2,9,10] and [6,7,12,13]
            # The G-G pairs connect positions that are NOT in the same tetrad
            # (e.g., 1-6, 2-7, 9-12, 10-13 are loop connections)
            
            # Strategy: Find two groups where G-G pairs connect between groups, not within groups
            # Try all ways to split into two groups of 4
            for combo1 in itertools.combinations(positions, 4):
                remaining = [p for p in positions if p not in combo1]
                if len(remaining) == 4:  # Exactly 4 remaining
                    combo2 = tuple(remaining)
                    tetrad1_positions = sorted(combo1)
                    tetrad2_positions = sorted(combo2)
                    
                    # Check if G-G pairs connect between groups (loop pattern)
                    # In parallel G-quad, G-G pairs should connect tetrad1[i] to tetrad2[i]
                    between_group_pairs = 0
                    within_group1_pairs = 0
                    within_group2_pairs = 0
                    
                    for pair in gg_pairs:
                        res_1 = pair['res_1']
                        res_2 = pair['res_2']
                        pos_1 = g_positions.get(res_1)
                        pos_2 = g_positions.get(res_2)
                        
                        if pos_1 is None or pos_2 is None:
                            continue
                        
                        in_group1_1 = pos_1 in tetrad1_positions
                        in_group1_2 = pos_2 in tetrad1_positions
                        in_group2_1 = pos_1 in tetrad2_positions
                        in_group2_2 = pos_2 in tetrad2_positions
                        
                        if (in_group1_1 and in_group2_2) or (in_group1_2 and in_group2_1):
                            between_group_pairs += 1
                        elif in_group1_1 and in_group1_2:
                            within_group1_pairs += 1
                        elif in_group2_1 and in_group2_2:
                            within_group2_pairs += 1
                    
                    # Parallel G-quad pattern: most pairs should be between groups (loop connections)
                    # and few/no pairs within groups
                    if between_group_pairs >= 3 and within_group1_pairs == 0 and within_group2_pairs == 0:
                        # Form tetrads
                        tetrad1 = sorted([g[0] for g in sorted_gs if g[1] in tetrad1_positions])
                        tetrad2 = sorted([g[0] for g in sorted_gs if g[1] in tetrad2_positions])
                        
                        if len(tetrad1) == 4 and len(tetrad2) == 4:
                            detected_quads.append({
                                'tetrad': tetrad1,
                                'confidence': 1.4,
                                'hoogsteen_edges': 'N/A (parallel G-quad)',
                                'hoogsteen_ratio': 0.0,
                                'is_planar': None,
                                'geometry': {},
                                'gquad_hbonds': 0,
                                'total_hbonds': 0,
                                'stacking_count': 2
                            })
                            detected_quads.append({
                                'tetrad': tetrad2,
                                'confidence': 1.4,
                                'hoogsteen_edges': 'N/A (parallel G-quad)',
                                'hoogsteen_ratio': 0.0,
                                'is_planar': None,
                                'geometry': {},
                                'gquad_hbonds': 0,
                                'total_hbonds': 0,
                                'stacking_count': 2
                            })
                            # Return first valid combination
                            return detected_quads
        
        return detected_quads
    
    def detect_g_quads_by_pattern(self, pdb_id: str, g_residues: Set[str]) -> List[Dict]:
        """
        Alternative detection: Look for G-quad pattern (4 chains, 4 G residues each, or 2 chains, 4 G residues each).
        This handles cases where G-G base pairs aren't in base pair data.
        """
        if len(g_residues) < 8:  # Need at least 8 G residues (2 chains x 4 G residues)
            return []
        
        # Group G residues by chain
        chains = defaultdict(list)
        for g in g_residues:
            parts = g.split('-')
            if len(parts) >= 3:
                chain = parts[0]
                res_num = parts[2]
                chains[chain].append((g, res_num))
        
        # Look for chains with 4 G residues (potential G-quad chains)
        candidate_chains = {ch: gs for ch, gs in chains.items() if len(gs) >= 4}
        
        # Handle both 4-chain and 2-chain G-quad patterns
        if len(candidate_chains) < 2:
            return []
        
        # If we have 2 chains with 4 G residues each, that's also a valid pattern
        min_chains_required = 2 if len(candidate_chains) == 2 else 4
        if len(candidate_chains) < min_chains_required:
            return []
        
        # Try to find 4 chains with G residues at similar positions
        # Group by relative position within chain
        detected_quads = []
        
        # For each set of chains (2 or 4), check if they have G residues that could form tetrads
        num_chains = len(candidate_chains)
        if num_chains >= 4:
            chain_combos = list(itertools.combinations(candidate_chains.keys(), 4))
        elif num_chains == 2:
            # 2-chain G-quad: each chain contributes 2 G residues to each tetrad
            chain_combos = [tuple(candidate_chains.keys())]
        else:
            return []
        
        for combo in chain_combos[:10]:  # Limit to avoid too many combinations
            chains_g_residues = {ch: sorted(candidate_chains[ch], key=lambda x: int(x[1]) if x[1].isdigit() else 999) 
                               for ch in combo}
            
            # Check if each chain has at least required G residues
            required_count = 4 if len(combo) == 4 else 4  # 4 G residues per chain
            all_have_required = all(len(gs) >= required_count for gs in chains_g_residues.values())
            
            if all_have_required:
                # Try to form tetrads
                tetrads = []
                min_g_count = min(len(gs) for gs in chains_g_residues.values())
                
                if len(combo) == 4:
                    # 4-chain G-quad: take G residues at same position from each chain
                    for i in range(min(min_g_count, 4)):  # Up to 4 tetrads
                        tetrad = tuple(sorted([gs[i][0] for gs in chains_g_residues.values()]))
                        if len(tetrad) == 4:
                            tetrads.append(tetrad)
                elif len(combo) == 2:
                    # 2-chain G-quad: each chain contributes 2 G residues to each tetrad
                    # Form tetrads by pairing G residues from both chains
                    chain_list = list(chains_g_residues.keys())
                    chain1_gs = chains_g_residues[chain_list[0]]
                    chain2_gs = chains_g_residues[chain_list[1]]
                    
                    # Sequential tetrads: chain1[i] + chain1[i+1] + chain2[i] + chain2[i+1]
                    # This forms stacked tetrads (most common pattern)
                    for i in range(min(len(chain1_gs), len(chain2_gs)) - 1):
                        if i + 1 < min(len(chain1_gs), len(chain2_gs)):
                            tetrad = tuple(sorted([
                                chain1_gs[i][0], chain1_gs[i+1][0],
                                chain2_gs[i][0], chain2_gs[i+1][0]
                            ]))
                            if len(tetrad) == 4:
                                tetrads.append(tetrad)
                
                if len(tetrads) >= 1:
                    # This looks like a G-quad pattern
                    for tetrad in tetrads:
                        detected_quads.append({
                            'tetrad': list(tetrad),
                            'confidence': 1.3,  # Lower confidence since we can't validate geometry
                            'hoogsteen_edges': 'N/A (pattern-based)',
                            'hoogsteen_ratio': 0.0,
                            'is_planar': None,
                            'geometry': {},
                            'gquad_hbonds': 0,
                            'total_hbonds': 0,
                            'stacking_count': len(tetrads)
                        })
        
        return detected_quads
    
    def detect_g_quads(self, pdb_id: str) -> Dict:
        """
        Main detection function for a single structure.
        Returns: dict with detection results
        """
        basepair_file = self.basepairs_dir / f"{pdb_id}.json"
        if not basepair_file.exists():
            return {
                'pdb_id': pdb_id,
                'g_quad_detected': False,
                'confidence_score': 0.0,
                'num_g_quads': 0,
                'g_quad_details': []
            }
        
        try:
            with open(basepair_file, 'r') as f:
                basepairs = json.load(f)
        except Exception:
            return {
                'pdb_id': pdb_id,
                'g_quad_detected': False,
                'confidence_score': 0.0,
                'num_g_quads': 0,
                'g_quad_details': []
            }
        
        # Tier 1: Find G-G pairs and build connectivity graph
        gg_pairs = self.find_gg_pairs(basepairs)
        
        # Special case: If we have G-G pairs but they don't form 4-cycles,
        # try to detect parallel G-quad pattern (single strand with loop connections)
        if len(gg_pairs) >= 4:
            graph = self.build_connectivity_graph(gg_pairs)
            cycles = self.find_4_cycles(graph)
            
            # If no 4-cycles but we have G-G pairs, try parallel G-quad detection
            if not cycles:
                parallel_quads = self.detect_parallel_gquad(gg_pairs, pdb_id)
                if parallel_quads:
                    return {
                        'pdb_id': pdb_id,
                        'g_quad_detected': True,
                        'confidence_score': max(q['confidence'] for q in parallel_quads),
                        'num_g_quads': len(parallel_quads),
                        'g_quad_details': parallel_quads
                    }
        
        # If no G-G pairs found, try pattern-based detection
        if len(gg_pairs) < 4:
            g_residues = self.find_g_residues_from_hbonds(pdb_id)
            # Check for G-quad patterns: 16+ G residues (4 chains) or 8+ G residues (2 chains)
            if len(g_residues) >= 8:  # Potential G-quad pattern (2 or 4 chains)
                pattern_quads = self.detect_g_quads_by_pattern(pdb_id, g_residues)
                if pattern_quads:
                    return {
                        'pdb_id': pdb_id,
                        'g_quad_detected': True,
                        'confidence_score': max(q['confidence'] for q in pattern_quads),
                        'num_g_quads': len(pattern_quads),
                        'g_quad_details': pattern_quads
                    }
            
            return {
                'pdb_id': pdb_id,
                'g_quad_detected': False,
                'confidence_score': 0.0,
                'num_g_quads': 0,
                'g_quad_details': []
            }
        
        graph = self.build_connectivity_graph(gg_pairs)
        cycles = self.find_4_cycles(graph)
        
        # If no perfect 4-cycles found, try cluster-based approach
        if not cycles:
            clusters = self.find_g_clusters(graph, gg_pairs)
            # Try to find 4-cycles within clusters
            for cluster in clusters:
                if len(cluster) >= 4:
                    # Try all 4-combinations within cluster
                    for combo in itertools.combinations(cluster, 4):
                        connections = {node: len(graph[node] & set(combo)) for node in combo}
                        if all(count >= 2 for count in connections.values()):
                            # Check if it forms a valid cycle
                            subgraph = {node: graph[node] & set(combo) for node in combo}
                            if self._is_valid_4_cycle(subgraph, combo):
                                cycles.append(tuple(sorted(combo)))
            cycles = list(set(cycles))
        
        if not cycles:
            return {
                'pdb_id': pdb_id,
                'g_quad_detected': False,
                'confidence_score': 0.0,
                'num_g_quads': 0,
                'g_quad_details': []
            }
        
        # Detect stacking
        stacking = self.detect_stacking(cycles)
        
        # Validate each cycle
        validated_quads = []
        for cycle in cycles:
            # Tier 2: Hoogsteen validation
            hoogsteen_count, total_edges = self.validate_hoogsteen(cycle, gg_pairs)
            hoogsteen_ratio = hoogsteen_count / total_edges if total_edges > 0 else 0
            
            # Tier 3: Geometric validation
            is_planar, geom_stats = self.validate_geometry(cycle, gg_pairs)
            
            # Tier 4: H-bond validation
            gquad_hbonds, total_hbonds = self.validate_hbonds(cycle, pdb_id)
            has_hbonds = gquad_hbonds >= 4
            
            # Calculate confidence
            stacking_count = stacking.get(cycle, 1)
            confidence = self.calculate_confidence_score(
                hoogsteen_ratio, is_planar, has_hbonds, stacking_count
            )
            
            # Only include if confidence >= 1.5 (high confidence threshold)
            if confidence >= 1.5:
                validated_quads.append({
                    'tetrad': list(cycle),
                    'confidence': confidence,
                    'hoogsteen_edges': f"{hoogsteen_count}/{total_edges}",
                    'hoogsteen_ratio': round(hoogsteen_ratio, 2),
                    'is_planar': is_planar,
                    'geometry': geom_stats,
                    'gquad_hbonds': gquad_hbonds,
                    'total_hbonds': total_hbonds,
                    'stacking_count': stacking_count
                })
        
        # Sort by confidence (highest first)
        validated_quads.sort(key=lambda x: x['confidence'], reverse=True)
        
        return {
            'pdb_id': pdb_id,
            'g_quad_detected': len(validated_quads) > 0,
            'confidence_score': validated_quads[0]['confidence'] if validated_quads else 0.0,
            'num_g_quads': len(validated_quads),
            'g_quad_details': validated_quads
        }
    
    def process_all_structures(self, output_csv: str = 'g_quads_detected.csv'):
        """Process all RNA structures and save results to CSV."""
        print("="*80)
        print("G-QUADRUPLEX DETECTION")
        print("="*80)
        
        # Get all PDB IDs
        basepair_files = list(self.basepairs_dir.glob('*.json'))
        pdb_ids = sorted([f.stem for f in basepair_files])
        
        print(f"\nFound {len(pdb_ids)} RNA structures to analyze")
        print("Starting detection...\n")
        
        results = []
        detected_count = 0
        
        for i, pdb_id in enumerate(pdb_ids, 1):
            if i % 100 == 0:
                print(f"Progress: {i}/{len(pdb_ids)} | Detected: {detected_count}")
            
            result = self.detect_g_quads(pdb_id)
            results.append(result)
            
            if result['g_quad_detected']:
                detected_count += 1
        
        # Save to CSV
        print(f"\nSaving results to {output_csv}...")
        with open(output_csv, 'w', newline='') as f:
            fieldnames = ['PDB_ID', 'G_Quad_Detected', 'Confidence_Score', 
                         'Num_G_Quads', 'G_Quad_Details']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            
            for result in results:
                writer.writerow({
                    'PDB_ID': result['pdb_id'],
                    'G_Quad_Detected': result['g_quad_detected'],
                    'Confidence_Score': result['confidence_score'],
                    'Num_G_Quads': result['num_g_quads'],
                    'G_Quad_Details': json.dumps(result['g_quad_details']) if result['g_quad_details'] else ''
                })
        
        # Summary
        print(f"\n{'='*80}")
        print("DETECTION COMPLETE")
        print(f"{'='*80}")
        print(f"Total structures analyzed: {len(pdb_ids)}")
        print(f"Structures with G-quads detected: {detected_count}")
        print(f"Detection rate: {100*detected_count/len(pdb_ids):.2f}%")
        print(f"\nResults saved to: {output_csv}")
        print(f"{'='*80}\n")


def main():
    detector = GQuadDetector()
    detector.process_all_structures()


if __name__ == '__main__':
    main()

