"""Analyzer for base pair geometry using precomputed data."""

from typing import Dict, List
from .analyzers_utils import BasePairScoring

class BasePairAnalyzer:
    """Analyzes base pair geometry quality."""
    
    name = "base_pairs"
    
    def __init__(self, config):
        self.config = config
    
    def analyze(self, basepair_data: list) -> Dict:
        if not basepair_data:
            return {
                'total_pairs': 0,
                'misaligned_frac': 0.0,
                'twisted_frac': 0.0,
                'non_coplanar_frac': 0.0,
                'low_hbond_score_frac': 0.0,
                'zero_hbond_frac': 0.0,
                'penalty': 0.0,
                'stats': {},
                'issues': [],
                'detailed_issues': []
            }
        
        # **STEP 1: Filter out adjacent and self pairs**
        valid_pairs = []
        for bp in basepair_data:
            res_1 = bp.get('res_1', '')
            res_2 = bp.get('res_2', '')
            
            # Skip adjacent pairs (stacking interactions - not issues)
            if BasePairScoring.check_adjacent_pairing(res_1, res_2):
                continue
            
            # Skip self pairs (structural errors in data)
            if BasePairScoring.check_self_pairing(res_1, res_2):
                continue
            
            valid_pairs.append(bp)
        
        # If no valid pairs after filtering
        if not valid_pairs:
            return {
                'total_pairs': 0,
                'misaligned_frac': 0.0,
                'twisted_frac': 0.0,
                'non_coplanar_frac': 0.0,
                'low_hbond_score_frac': 0.0,
                'zero_hbond_frac': 0.0,
                'penalty': 0.0,
                'stats': {},
                'issues': [],
                'detailed_issues': []
            }
        
        detailed_issues = []
        issues = []
        stats = {
            'canonical_wc': 0,
            'non_canonical': 0,
            'misaligned': 0,
            'twisted': 0,
            'non_coplanar': 0,
            'low_hbond_score': 0,
            'zero_hbond': 0,
        }
        
        # **STEP 2: Analyze only valid (non-adjacent) pairs**
        for bp in valid_pairs:
            lw = bp.get('lw', '')
            res_1 = bp.get('res_1', '?')
            res_2 = bp.get('res_2', '?')
            bp_id = f"{res_1}-{res_2}"
            bp_type = bp.get('bp_type', 'unknown')
            
            if lw == 'cWW':
                stats['canonical_wc'] += 1
            else:
                stats['non_canonical'] += 1
                
            # Use shared scoring logic to evaluate the base pair
            bp_issues = BasePairScoring.score_base_pair(bp, self.config)
            
            # Skip if None (shouldn't happen since we filtered, but safe check)
            if bp_issues is None:
                continue
            
            # Collect detailed issue information
            if any(bp_issues.values()):
                specific_issues = []
                for issue_type, is_present in bp_issues.items():
                    if is_present:
                        specific_issues.append(issue_type)
                
                # Extract residue numbers safely
                try:
                    chain_1 = res_1.split('-')[0]
                    residue_1 = int(res_1.split('-')[2])
                    chain_2 = res_2.split('-')[0]
                    residue_2 = int(res_2.split('-')[2])
                except (IndexError, ValueError):
                    # Fallback if residue format is unexpected
                    chain_1 = res_1
                    residue_1 = 0
                    chain_2 = res_2
                    residue_2 = 0
                
                detailed_issue = {
                    "type": "base_pair",
                    "residues": f"{res_1} - {res_2}",
                    "chain_1": chain_1,
                    "chain_2": chain_2,
                    "residue_1": residue_1,
                    "residue_2": residue_2,
                    "bp_type": bp_type,
                    "lw_notation": lw,
                    "specific_issues": specific_issues,
                    "geometry_parameters": {
                        "shear": bp.get('shear', 0),
                        "stretch": bp.get('stretch', 0), 
                        "stagger": bp.get('stagger', 0),
                        "buckle": bp.get('buckle', 0),
                        "propeller": bp.get('propeller', 0),
                        "opening": bp.get('opening', 0),
                        "hbond_score": bp.get('hbond_score', 0)
                    }
                }
                detailed_issues.append(detailed_issue)
            
            # Update stats and issues based on scoring results
            if bp_issues['misaligned']:
                stats['misaligned'] += 1
                issues.append(f"Misaligned: {bp_id} ({bp_type})")
            
            if bp_issues['twisted']:
                stats['twisted'] += 1
                issues.append(f"Twisted: {bp_id} ({bp_type})")
            
            if bp_issues['non_coplanar']:
                stats['non_coplanar'] += 1
                issues.append(f"Non-coplanar: {bp_id} ({bp_type})")
            
            if bp_issues['low_hbond_score']:
                stats['low_hbond_score'] += 1
                issues.append(f"Poor H-bonding: {bp_id} ({bp_type})")
            
            if bp_issues['zero_hbond']:
                stats['zero_hbond'] += 1
                issues.append(f"Zero H-bonding in pair: {bp_id} ({bp_type})")
        
        # **STEP 3: Calculate fractions and penalties using only valid pairs**
        total_pairs = len(valid_pairs)
        
        # Calculate fractions
        misaligned_frac = stats['misaligned'] / total_pairs
        twisted_frac = stats['twisted'] / total_pairs
        non_coplanar_frac = stats['non_coplanar'] / total_pairs
        low_hbond_score_frac = stats['low_hbond_score'] / total_pairs
        zero_hbond_frac = stats['zero_hbond'] / total_pairs
        
        # Use shared scoring to calculate penalty
        penalty = BasePairScoring.calculate_basepair_penalty(stats, total_pairs, self.config)
        
        return {
            'total_pairs': total_pairs,
            'misaligned_frac': round(misaligned_frac, 3),
            'twisted_frac': round(twisted_frac, 3),
            'non_coplanar_frac': round(non_coplanar_frac, 3),
            'low_hbond_score_frac': round(low_hbond_score_frac, 3),
            'zero_hbond_frac': round(zero_hbond_frac, 3),
            'penalty': round(penalty, 2),
            'stats': stats,
            'issues': issues[:self.config.MAX_ISSUES_DISPLAYED],
            'detailed_issues': detailed_issues
        }