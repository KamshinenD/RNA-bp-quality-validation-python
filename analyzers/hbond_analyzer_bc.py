"""Analyzer for hydrogen bond quality."""

import pandas as pd
from typing import Dict
from .analyzers_utils import HBondScoring, ScoringUtils, BasePairScoring


class HBondAnalyzer:
    """Analyzes hydrogen bond quality."""
    
    name = "hbonds"
    
    def __init__(self, config):
        self.config = config
        
    def analyze(self, hbond_df: pd.DataFrame, basepair_data: list = None) -> Dict:
        if hbond_df is None or len(hbond_df) == 0:
            return self._empty_result()
        
         # **FILTER 1: Build base pair set for validation**
        base_pair_set = set()
        if basepair_data:
            base_pair_set = ScoringUtils.build_base_pair_set(basepair_data)
            
        # **FILTER 2: Only keep base-base H-bonds that belong to base pairs**
        valid_hbonds = []
        for _, row in hbond_df.iterrows():
            # Check if atoms are base atoms (not backbone)
            if not ScoringUtils.is_base_base_hbond(row['atom_1'], row['atom_2']):
                continue
            
            # **Skip adjacent pairs**
            if BasePairScoring.check_adjacent_pairing(row['res_1'], row['res_2']):
                continue
            
            # Check if this H-bond belongs to a recognized base pair
            pair_key = tuple(sorted([row['res_1'], row['res_2']]))
            if pair_key not in base_pair_set:
                continue
            
            valid_hbonds.append(row)
        
        if not valid_hbonds:
            return self._empty_result()
        
        # Convert back to DataFrame
        hbond_df = pd.DataFrame(valid_hbonds)
        
        # Continue with rest of analysis...
        issues = []
        detailed_issues = []
        stats = {
            'bad_distance': 0,
            'bad_angles': 0,
            'bad_dihedral': 0,
            'weak_quality': 0,
            'incorrect_hbond_count': 0
        }
        
         # Analyze individual H-bonds
        for index, row in hbond_df.iterrows():
            hbond_id = f"{row['res_1']}-{row['res_2']}"
            atoms = f"{row['atom_1']}-{row['atom_2']}"
            
            # Use shared scoring logic to evaluate the H-bond
            hb_issues = HBondScoring.score_hbond(row, self.config)
            
            # Create detailed issue information
            if any(hb_issues.values()):
                specific_issues = []
                for issue_type, is_present in hb_issues.items():
                    if is_present:
                        specific_issues.append(issue_type)
                
                # Extract residue numbers safely
                try:
                    chain_1 = row['res_1'].split('-')[0]
                    residue_1 = int(row['res_1'].split('-')[2])
                    chain_2 = row['res_2'].split('-')[0]
                    residue_2 = int(row['res_2'].split('-')[2])
                except (IndexError, ValueError):
                    # Fallback if residue format is unexpected
                    chain_1 = row['res_1']
                    residue_1 = 0
                    chain_2 = row['res_2']
                    residue_2 = 0
                
                detailed_issue = {
                    "type": "hydrogen_bond",
                    "residues": f"{row['res_1']} - {row['res_2']}",
                    "chain_1": chain_1,
                    "chain_2": chain_2,
                    "residue_1": residue_1,
                    "residue_2": residue_2,
                    "atoms": f"{row['atom_1']} - {row['atom_2']}",
                    "specific_issues": specific_issues,
                    "hbond_parameters": {
                        "distance": row['distance'],
                        "angle_1": row['angle_1'],
                        "angle_2": row['angle_2'],
                        "dihedral": row['dihedral_angle'],
                        "quality_score": row['score']
                    }
                }
                detailed_issues.append(detailed_issue)
    
           
            
            # Update stats and issues based on scoring results
            if hb_issues['bad_distance']:
                stats['bad_distance'] += 1
                issues.append(f"Bad distance: {hbond_id} ({atoms}): {row['distance']:.2f}Å")
            
            if hb_issues['bad_angles']:
                stats['bad_angles'] += 1
                issues.append(f"Bad angles: {hbond_id} ({atoms})")
            
            if hb_issues['bad_dihedral']:
                stats['bad_dihedral'] += 1
                # Get deviation for detailed reporting
                _, deviation = HBondScoring.check_dihedral(row['dihedral_angle'], self.config)
                issues.append(f"Bad dihedral: {hbond_id} ({atoms}): {row['dihedral_angle']:.1f}° (dev: {deviation:.1f}°)")
            
            if hb_issues['weak_quality']:
                stats['weak_quality'] += 1
                issues.append(f"Weak quality: {hbond_id} ({atoms}): {row['score']:.2f}")
        
        # Validate H-bond counts per base pair using shared logic
        total_pairs_checked = 0
            
        if basepair_data is not None:
            # Use shared utility to count H-bonds per pair
            pair_hbond_counts = ScoringUtils.count_hbonds_per_pair(hbond_df)
            
            for bp in basepair_data:
                res_1 = bp.get('res_1', '?')
                res_2 = bp.get('res_2', '?')
                bp_type = bp.get('bp_type', 'unknown')
                lw = bp.get('lw', '')
                
                # Create canonical key
                pair_key = tuple(sorted([res_1, res_2]))
                actual_hbonds = pair_hbond_counts.get(pair_key, 0)
                
                
                # Use shared scoring to check H-bond count
                is_incorrect, issue_desc = HBondScoring.check_hbond_count(
                    bp_type, actual_hbonds, lw, self.config
                )
                
                if bp_type in self.config.EXPECTED_HBOND_COUNTS:
                    total_pairs_checked += 1
                    
                    if is_incorrect:
                        stats['incorrect_hbond_count'] += 1
                        min_exp, max_exp, ideal = self.config.EXPECTED_HBOND_COUNTS[bp_type]
                        issues.append(f"{issue_desc} for {res_1}-{res_2} ({bp_type})")
        
        total_hbonds = len(hbond_df)
        
        # Calculate fractions for reporting
        bad_distance_frac = stats['bad_distance'] / total_hbonds
        bad_angles_frac = stats['bad_angles'] / total_hbonds
        bad_dihedral_frac = stats['bad_dihedral'] / total_hbonds
        weak_quality_frac = stats['weak_quality'] / total_hbonds
        incorrect_hbond_count_frac = stats['incorrect_hbond_count'] / total_pairs_checked if total_pairs_checked > 0 else 0
        
        # Calculate penalty using weights from config
        # penalty = (
        #     self.config.PENALTY_WEIGHTS['bad_hbond_distance'] * bad_distance_frac +
        #     self.config.PENALTY_WEIGHTS['bad_hbond_angles'] * bad_angles_frac +
        #     self.config.PENALTY_WEIGHTS['bad_hbond_dihedrals'] * bad_dihedral_frac +
        #     self.config.PENALTY_WEIGHTS['weak_hbond_quality'] * weak_quality_frac +
        #     self.config.PENALTY_WEIGHTS['incorrect_hbond_count'] * incorrect_hbond_count_frac
        # )
        
        # Calculate fractions for reporting
        bad_distance_frac = stats['bad_distance'] / total_hbonds
        bad_angles_frac = stats['bad_angles'] / total_hbonds
        bad_dihedral_frac = stats['bad_dihedral'] / total_hbonds
        weak_quality_frac = stats['weak_quality'] / total_hbonds
        incorrect_hbond_count_frac = (stats['incorrect_hbond_count'] / total_pairs_checked 
                                     if total_pairs_checked > 0 else 0)
        
        # Use shared scoring to calculate penalty
        penalty = HBondScoring.calculate_hbond_penalty(
            stats, total_hbonds, total_pairs_checked, self.config
        )


        return {
            'total_hbonds': total_hbonds,
            'bad_distance_frac': round(bad_distance_frac, 3),
            'bad_angles_frac': round(bad_angles_frac, 3),
            'bad_dihedral_frac': round(bad_dihedral_frac, 3),
            'weak_quality_frac': round(weak_quality_frac, 3),
            'incorrect_hbond_count_frac': round(incorrect_hbond_count_frac, 3),
            'penalty': round(penalty, 1),
            'stats': stats,
            'issues': issues[:self.config.MAX_ISSUES_DISPLAYED],
            'detailed_issues': detailed_issues
        }
    
    
    def _empty_result(self) -> Dict:
        """Return empty result structure when no valid H-bonds found."""
        return {
            'total_hbonds': 0,
            'bad_distance_frac': 0.0,
            'bad_angles_frac': 0.0,
            'bad_dihedral_frac': 0.0,
            'weak_quality_frac': 0.0,
            'incorrect_hbond_count_frac': 0.0,
            'penalty': 0.0,
            'stats': {},
            'issues': [],
            'detailed_issues': []
        }