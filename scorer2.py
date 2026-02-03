"""Baseline RNA structure quality scorer - Base pairs only, no hotspots."""

import pandas as pd
from typing import Dict, List, Tuple
from dataclasses import dataclass


@dataclass
class BaselineResult:
    """Results from base pair scoring."""
    overall_score: float  # 0-100
    total_base_pairs: int
    
    # Per base pair metrics
    avg_basepair_score: float
    basepair_scores: List[Dict]  # Individual BP scores with details
    
    # Aggregate statistics
    geometry_issues: Dict[str, int]
    geometry_fractions: Dict[str, float]
    hbond_issues: Dict[str, int]
    hbond_fractions: Dict[str, float]
    
    # Detailed breakdown
    grade: str  # EXCELLENT, GOOD, FAIR, POOR
    summary: str
    detailed_issues: List[Dict]


class Scorer:
    """
    Score RNA structures based on base pair quality only.
    
    Scoring approach:
    1. Each base pair gets scored based on:
       - Geometry (alignment, twist, coplanarity)
       - Its hydrogen bonds (distance, angles, quality)
    2. Overall score = average of all base pair scores
    3. No hotspot detection - pure base pair quality assessment
    """
    
    def __init__(self, config, bp_analyzer=None, hb_analyzer=None):
        """
        Args:
            config: Configuration with thresholds and weights
            bp_analyzer: Not used in this version (kept for compatibility)
            hb_analyzer: Not used in this version (kept for compatibility)
        """
        self.config = config
    
    def score_structure(self, basepair_data: list, hbond_data: pd.DataFrame) -> BaselineResult:
        """
        Score the entire structure based on base pair quality.
        
        Args:
            basepair_data: List of base pair dictionaries
            hbond_data: DataFrame of hydrogen bonds
            
        Returns:
            Result with overall quality score and breakdown
        """
        print("\n" + "="*60)
        print("RNA STRUCTURE QUALITY ASSESSMENT")
        print("="*60)
        
        if not basepair_data:
            return self._empty_result()
        
        print(f"\nAnalyzing {len(basepair_data)} base pairs...")
        
        # Score each base pair individually
        basepair_scores = []
        geometry_stats = {
            'misaligned': 0, 'rotational_distortion': 0, 'non_coplanar': 0,
            'improper_opening': 0, 'zero_hbond': 0
        }
        hbond_stats = {
            'low_dssr_score': 0, 'bad_distance': 0, 'bad_angles': 0, 'bad_dihedral': 0,
            'weak_quality': 0, 'incorrect_count': 0
        }
        
        for bp in basepair_data:
            bp_score = self._score_base_pair(bp, hbond_data)
            basepair_scores.append(bp_score)
            
            # Accumulate statistics
            for issue in geometry_stats:
                if bp_score['geometry_issues'].get(issue, False):
                    geometry_stats[issue] += 1
            
            for issue in hbond_stats:
                if bp_score['hbond_issues'].get(issue, False):
                    hbond_stats[issue] += 1
        
        # Calculate overall score
        total_pairs = len(basepair_scores)
        avg_score = sum(bp['score'] for bp in basepair_scores) / total_pairs if total_pairs > 0 else 0
        
        print(f"\n{'='*60}")
        print(f"Average base pair score: {avg_score:.1f}/100")
        
        # Calculate fractions
        geometry_fractions = {k: v/total_pairs for k, v in geometry_stats.items()}
        hbond_fractions = {k: v/total_pairs for k, v in hbond_stats.items()}
        
        # Determine grade
        grade = self._assign_grade(avg_score)
        
        # Create summary
        summary = self._create_summary(
            avg_score, grade, total_pairs,
            geometry_stats, geometry_fractions,
            hbond_stats, hbond_fractions
        )
        
        # Collect detailed issues (only problematic base pairs)
        # Use BASELINE threshold (75) to match Detailed_Issues column in CSV
        detailed_issues = [
            bp for bp in basepair_scores 
            if bp['score'] < self.config.BASELINE
        ]
        
        result = BaselineResult(
            overall_score=round(avg_score, 1),
            total_base_pairs=total_pairs,
            avg_basepair_score=round(avg_score, 1),
            basepair_scores=basepair_scores,
            geometry_issues=geometry_stats,
            geometry_fractions=geometry_fractions,
            hbond_issues=hbond_stats,
            hbond_fractions=hbond_fractions,
            grade=grade,
            summary=summary,
            detailed_issues=detailed_issues
        )
        
        self._print_results(result)
        
        return result
    
    def _score_base_pair(self, bp: Dict, hbond_data: pd.DataFrame) -> Dict:
        """
        Score a single base pair based on geometry and its H-bonds.
        Uses edge-specific thresholds for geometry and H-bonds, global for dihedrals.
        
        Args:
            bp: Base pair dictionary with geometry info
            hbond_data: DataFrame of all hydrogen bonds
            
        Returns:
            Dictionary with:
                - score: 0-100 quality score
                - geometry_penalty: penalty from geometry issues
                - hbond_penalty: penalty from H-bond issues
                - geometry_issues: dict of boolean flags
                - hbond_issues: dict of boolean flags
                - bp_info: identifier info
        """
        nt1_id = bp.get('res_1', '')
        nt2_id = bp.get('res_2', '')
        
        # Resolve edge and base-pair type
        edge_type = bp.get('lw', '') or 'uncommon_edges'
        bp_type = bp.get('bp_type', '')

        # Get base-pair + edge thresholds for geometry (fallback to empty -> globals below)
        geo_by_bp = getattr(self.config, 'GEOMETRY_THRESHOLDS_BY_BASE_PAIR_BY_EDGE', {})
        geo_thresh = geo_by_bp.get(bp_type, {}).get(edge_type, {})

        # Get base-pair + edge thresholds for H-bonds (fallback to _OTHER within bp_type, else empty)
        hb_by_bp = getattr(self.config, 'HBOND_THRESHOLDS_BY_BASE_PAIR_BY_EDGE', {})
        hbond_thresh = hb_by_bp.get(bp_type, {}).get(edge_type) or hb_by_bp.get(bp_type, {}).get('_OTHER', {})
        
        # Get H-bonds for this base pair
        bp_hbonds = self._get_basepair_hbonds(nt1_id, nt2_id, hbond_data)
        
        # Analyze geometry using edge-specific thresholds
        geometry_issues = {}
        geometry_penalty = 0.0
        
        # Check alignment (shear and stretch)
        shear = abs(bp.get('shear', 0))
        stretch = bp.get('stretch', 0)
        shear_max = geo_thresh.get('SHEAR_MAX', 1.5)  # Fallback to old global if missing
        stretch_min = geo_thresh.get('STRETCH_MIN', -1.3)
        stretch_max = geo_thresh.get('STRETCH_MAX', 1.3)
        
        if shear > shear_max or stretch < stretch_min or stretch > stretch_max:
            geometry_issues['misaligned'] = True
            geometry_penalty += self.config.PENALTY_WEIGHTS['misaligned_pairs']
        
        # Check rotational distortion (buckle and propeller)
        buckle = abs(bp.get('buckle', 0))
        propeller = bp.get('propeller', 0)
        buckle_max = geo_thresh.get('BUCKLE_MAX', 16.1)
        propeller_min = geo_thresh.get('PROPELLER_MIN', -14.7)
        propeller_max = geo_thresh.get('PROPELLER_MAX', 14.7)

        if buckle > buckle_max or propeller < propeller_min or propeller > propeller_max:
            geometry_issues['rotational_distortion'] = True
            geometry_penalty += self.config.PENALTY_WEIGHTS['rotational_distortion_pairs']

        # Check coplanarity (stagger only)
        stagger = abs(bp.get('stagger', 0))
        stagger_max = geo_thresh.get('STAGGER_MAX', 0.4)

        if stagger > stagger_max:
            geometry_issues['non_coplanar'] = True
            geometry_penalty += self.config.PENALTY_WEIGHTS['non_coplanar_pairs']

        # Check opening angle separately
        opening = bp.get('opening', 0)
        opening_min = geo_thresh.get('OPENING_MIN', -30.8)
        opening_max = geo_thresh.get('OPENING_MAX', 30.8)

        if opening < opening_min or opening > opening_max:
            geometry_issues['improper_opening'] = True
            geometry_penalty += self.config.PENALTY_WEIGHTS['improper_pair_opening']
        
        # Check DSSR quality score (hbond_score in the JSON)
        dssr_quality = bp.get('hbond_score', 0.0)
        has_dssr_hbonds = dssr_quality > 0.0  # DSSR detected H-bonds if score > 0.0

        # Analyze H-bonds for this base pair
        hbond_issues = {}
        hbond_penalty = 0.0

        # Flag low DSSR hbond_score for reporting (not penalised)
        # If score == 0.0, zero_hbond penalty will be applied instead
        if dssr_quality > 0.0 and dssr_quality < self.config.HBOND_SCORE_MIN:
            hbond_issues['low_dssr_score'] = True  # (not penalised)
        
        # Reconcile DSSR hbond_score with CSV H-bond data
        # Only flag "zero_hbond" if BOTH sources agree there are no H-bonds
        if len(bp_hbonds) == 0:
            # If DSSR says no H-bonds (score == 0.0), then flag zero_hbond
            if not has_dssr_hbonds:
                geometry_issues['zero_hbond'] = True
                geometry_penalty += self.config.PENALTY_WEIGHTS['zero_hbond_pairs']
            # If DSSR says H-bonds exist but CSV doesn't, it's a data mismatch
            # Don't penalize - trust DSSR's detection (CSV might be incomplete or filtered)
        else:
            # Check expected H-bond count using range-based checking
            bp_type = bp.get('bp_type', '')
            actual_count = len(bp_hbonds)
            
            if bp_type in self.config.EXPECTED_HBOND_COUNTS:
                # Interpret as (min, max, ideal) to match analyzers_utils.py
                min_expected, max_expected, ideal = self.config.EXPECTED_HBOND_COUNTS[bp_type]
                lw = bp.get('lw', '')
                
                # Check if count is outside expected range
                is_outside_range = actual_count < min_expected or actual_count > max_expected
                
                # For canonical Watson-Crick pairs, also flag if not ideal (stricter check)
                is_suboptimal_cww = (lw == 'cWW' and actual_count != ideal)
                
                # Penalize if outside range OR suboptimal cWW
                if is_outside_range or is_suboptimal_cww:
                    hbond_issues['incorrect_count'] = True
                    hbond_penalty += self.config.PENALTY_WEIGHTS['incorrect_hbond_count']
            
            # Check each H-bond quality using edge-specific thresholds
            # Collect flags first
            has_bad_distance = False
            has_bad_angles = False
            has_bad_dihedral = False
            has_weak_quality = False
            
            # Get edge-specific H-bond thresholds
            distance_min = hbond_thresh.get('DIST_MIN', self.config.HBOND_DISTANCE_MIN)
            distance_max = hbond_thresh.get('DIST_MAX', self.config.HBOND_DISTANCE_MAX)
            angle_min = hbond_thresh.get('ANGLE_MIN', 80.0)
            quality_min = hbond_thresh.get('QUALITY_MIN', 0.70)
            
            # Dihedral uses GLOBAL thresholds (only flag forbidden zone)
            dihedral_cis_min = self.config.HBOND_DIHEDRAL_CIS_MIN
            dihedral_cis_max = self.config.HBOND_DIHEDRAL_CIS_MAX
            dihedral_trans_min = self.config.HBOND_DIHEDRAL_TRANS_MIN
            
            for _, hb in bp_hbonds.iterrows():
                # Distance check using edge-specific threshold
                distance = hb.get('distance', 0)
                if distance < distance_min or distance > distance_max:
                    has_bad_distance = True
                
                # Angle checks using edge-specific threshold
                angle_1 = hb.get('angle_1', 0)
                angle_2 = hb.get('angle_2', 0)
                if angle_1 < angle_min or angle_2 < angle_min:
                    has_bad_angles = True
                
                # Dihedral check - use GLOBAL thresholds (only flag forbidden zone)
                dihedral = hb.get('dihedral_angle', 0)
                is_cis = dihedral_cis_min <= dihedral <= dihedral_cis_max
                is_trans = dihedral >= dihedral_trans_min or dihedral <= -dihedral_trans_min
                if not (is_cis or is_trans):
                    has_bad_dihedral = True
                
                # Quality score check using edge-specific threshold
                quality = hb.get('score', 1.0)
                if quality < quality_min:
                    has_weak_quality = True
            
            # Apply penalties ONCE per issue type
            if has_bad_distance:
                hbond_issues['bad_distance'] = True
                hbond_penalty += self.config.PENALTY_WEIGHTS['bad_hbond_distance']
            
            if has_bad_angles:
                hbond_issues['bad_angles'] = True
                hbond_penalty += self.config.PENALTY_WEIGHTS['bad_hbond_angles']
            
            if has_bad_dihedral:
                hbond_issues['bad_dihedral'] = True
                hbond_penalty += self.config.PENALTY_WEIGHTS['bad_hbond_dihedrals']
            
            if has_weak_quality:
                hbond_issues['weak_quality'] = True
                hbond_penalty += self.config.PENALTY_WEIGHTS['weak_hbond_quality']
        
        # Calculate base pair score
        total_penalty = geometry_penalty + hbond_penalty
        bp_score = max(0, min(100, self.config.BASE_SCORE - total_penalty))
        
        return {
            'score': bp_score,
            'geometry_penalty': geometry_penalty,
            'hbond_penalty': hbond_penalty,
            'geometry_issues': geometry_issues,
            'hbond_issues': hbond_issues,
            'bp_info': {
                'res_1': nt1_id,
                'res_2': nt2_id,
                'bp_type': bp.get('bp_type', ''),
                'edge_type': edge_type,
                'num_hbonds': len(bp_hbonds),
                'dssr_score': bp.get('hbond_score', 0)
            }
        }
    
    def _get_basepair_hbonds(self, nt1_id: str, nt2_id: str, 
                            hbond_data: pd.DataFrame) -> pd.DataFrame:
        """
        Get hydrogen bonds between two nucleotides in a base pair.
        Only includes H-bonds between base atoms (not backbone/sugar).
        
        Args:
            nt1_id: First nucleotide identifier (e.g., "Q-C-0-")
            nt2_id: Second nucleotide identifier (e.g., "Q-G-22-")
            hbond_data: DataFrame of all hydrogen bonds
            
        Returns:
            DataFrame of base-base H-bonds for this base pair
        """
        if hbond_data.empty:
            return pd.DataFrame()
        
        # Find H-bonds where res_1/res_2 match this base pair (bidirectional)
        mask = (
            ((hbond_data['res_1'] == nt1_id) & (hbond_data['res_2'] == nt2_id)) |
            ((hbond_data['res_1'] == nt2_id) & (hbond_data['res_2'] == nt1_id))
        )
        
        bp_hbonds = hbond_data[mask]
        
        # Filter to only base-base H-bonds (exclude backbone/sugar)
        if not bp_hbonds.empty:
            base_base_mask = bp_hbonds.apply(
                lambda row: self._is_base_base_hbond(row['atom_1'], row['atom_2']),
                axis=1
            )
            bp_hbonds = bp_hbonds[base_base_mask]
        
        return bp_hbonds
    
    def _is_base_base_hbond(self, atom_1: str, atom_2: str) -> bool:
        """
        Check if H-bond is between two base atoms (not backbone/sugar).
        
        Args:
            atom_1: First atom name
            atom_2: Second atom name
            
        Returns:
            True if both atoms are base atoms
        """
        return (self._is_base_atom(atom_1) and self._is_base_atom(atom_2))
    
    def _is_base_atom(self, atom_name: str) -> bool:
        """
        Check if an atom is part of the base (not backbone/sugar).
        
        Args:
            atom_name: Atom name (e.g., 'N1', 'O2', 'C5')
            
        Returns:
            True if base atom, False if backbone/sugar atom
        """
        # Remove whitespace
        atom_name = atom_name.strip()
        
        # Pattern-based exclusion for robustness
        # Phosphate patterns
        if atom_name.startswith('P') and len(atom_name) <= 3:  # P, PA, PB, PG
            return False
        if atom_name.startswith('OP'):  # OP1, OP2, OP3
            return False
        if 'P' in atom_name and 'O' in atom_name:  # O1P, O2P, O3P
            return False
        
        # Sugar patterns (contains prime ' or asterisk *)
        if "'" in atom_name or '*' in atom_name:
            return False
        
        # Explicit exclusion list (backup)
        backbone_sugar_atoms = {
            # Phosphate
            'P', 'OP1', 'OP2', 'OP3', 'O1P', 'O2P', 'O3P',
            'PA', 'PB', 'PG',
            
            # Ribose (prime notation)
            "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
            
            # Ribose (asterisk notation)
            'O5*', 'C5*', 'C4*', 'O4*', 'C3*', 'O3*', 'C2*', 'O2*', 'C1*',
            
            # Hydrogens
            "HO5'", "HO2'", "HO3'", 'HO5*', 'HO2*', 'HO3*',
        }
        
        return atom_name not in backbone_sugar_atoms
    
    
    def _assign_grade(self, score: float) -> str:
        """Assign letter grade based on score."""
        if score >= self.config.GRADE_EXCELLENT:
            return "EXCELLENT"
        elif score >= self.config.GRADE_GOOD:
            return "GOOD"
        elif score >= self.config.GRADE_FAIR:
            return "FAIR"
        else:
            return "POOR"
    
    def _create_summary(self, score: float, grade: str, total_pairs: int,
                       geometry_stats: Dict, geometry_fractions: Dict,
                       hbond_stats: Dict, hbond_fractions: Dict) -> str:
        """Create human-readable summary of results."""
        lines = []
        lines.append(f"Overall Quality: {score:.1f}/100 ({grade})")
        lines.append(f"Base Pairs Analyzed: {total_pairs}")
        lines.append("")
        
        # Geometry issues
        lines.append("Geometry Issues:")
        has_geom_issues = False
        for issue, frac in geometry_fractions.items():
            if frac > 0.01:
                has_geom_issues = True
                lines.append(f"  - {frac*100:.1f}% {issue.replace('_', ' ')}")
        if not has_geom_issues:
            lines.append("  - No significant geometry issues")
        
        lines.append("")
        
        # H-bond issues
        lines.append("H-bond Issues:")
        has_hbond_issues = False
        for issue, frac in hbond_fractions.items():
            if frac > 0.01:
                has_hbond_issues = True
                lines.append(f"  - {frac*100:.1f}% {issue.replace('_', ' ')}")
        if not has_hbond_issues:
            lines.append("  - No significant H-bond issues")
        
        return "\n".join(lines)
    
    def _print_results(self, result: BaselineResult):
        """Print formatted results to console."""
        print("\n" + "="*60)
        print("RESULTS")
        print("="*60)
        print(f"\nOverall Score: {result.overall_score}/100 ({result.grade})")
        print(f"Base Pairs: {result.total_base_pairs}")
        
        print("\n" + "-"*60)
        print("GEOMETRY BREAKDOWN")
        print("-"*60)
        for issue, count in result.geometry_issues.items():
            if count > 0:
                frac = result.geometry_fractions[issue]
                print(f"{issue:20s}: {count:4d} ({frac*100:5.1f}%)")
        
        print("\n" + "-"*60)
        print("H-BOND BREAKDOWN")
        print("-"*60)
        for issue, count in result.hbond_issues.items():
            if count > 0:
                frac = result.hbond_fractions[issue]
                print(f"{issue:20s}: {count:4d} ({frac*100:5.1f}%)")
        
        print("\n" + "="*60)
        print("SUMMARY")
        print("="*60)
        print(result.summary)
        print("="*60 + "\n")
    
    def _empty_result(self) -> BaselineResult:
        """Return empty result when no base pairs found."""
        return BaselineResult(
            overall_score="N/A",
            total_base_pairs=0,
            avg_basepair_score=0.0,
            basepair_scores=[],
            geometry_issues={},
            geometry_fractions={},
            hbond_issues={},
            hbond_fractions={},
            grade="N/A",
            summary="No base pairs found for analysis",
            detailed_issues=[]
        )
    
    def export_to_dict(self, result: BaselineResult) -> Dict:
        """Export result to dictionary for JSON/CSV export."""
        return {
            'overall_score': result.overall_score,
            'grade': result.grade,
            'total_base_pairs': result.total_base_pairs,
            'avg_basepair_score': result.avg_basepair_score,
            'geometry_issues': result.geometry_issues,
            'geometry_fractions': result.geometry_fractions,
            'hbond_issues': result.hbond_issues,
            'hbond_fractions': result.hbond_fractions,
            'summary': result.summary,
            'num_problematic_bps': len(result.detailed_issues),
            'basepair_scores': result.basepair_scores
        }