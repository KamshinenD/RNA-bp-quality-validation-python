"""Baseline RNA structure quality scorer - Base pairs only, no hotspots."""

import json
import math
import pandas as pd
import numpy as np
from pathlib import Path
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
    avg_suiteness: float = None  # Average suiteness across all scored residues


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
    
    def __init__(self, config):
        """
        Args:
            config: Configuration with thresholds and weights
        """
        self.config = config
        self._load_richardson_suites()
        self._build_predecessor_map_cache = {}

    def _load_richardson_suites(self):
        """Load Richardson suite conformer definitions for backbone scoring."""
        suite_path = Path(__file__).parent / 'data' / 'richardson_suites.json'
        if not suite_path.exists():
            self.richardson_suites = None
            return

        with open(suite_path, 'r') as f:
            self.richardson_suites = json.load(f)

        # Pre-index conformers by bin for fast lookup
        self._conformers_by_bin = {}
        for conf in self.richardson_suites['conformers']:
            bin_key = conf['bin']
            if bin_key not in self._conformers_by_bin:
                self._conformers_by_bin[bin_key] = []
            self._conformers_by_bin[bin_key].append(conf)

        self._normal_widths = np.array(self.richardson_suites['widths']['normal'], dtype=float)

    def _find_predecessor(self, res_id: str, torsion_data: dict) -> str:
        """Find the preceding residue in the chain (same chain, resnum - 1).

        Returns the predecessor residue ID or None if not found.
        """
        parts = res_id.split('-')
        if len(parts) < 3:
            return None
        try:
            chain = parts[0]
            resnum = int(parts[2])
        except (ValueError, IndexError):
            return None

        # Search for chain-X-(resnum-1)- in torsion_data
        prev_num = resnum - 1
        for key in torsion_data:
            kparts = key.split('-')
            if len(kparts) >= 3 and kparts[0] == chain:
                try:
                    if int(kparts[2]) == prev_num:
                        return key
                except ValueError:
                    continue
        return None

    @staticmethod
    def _angle_diff(a: float, b: float) -> float:
        """Compute angular difference handling wraparound (0-360 space)."""
        d = abs(a - b)
        return min(d, 360.0 - d)

    def _l3_distance(self, suite_angles: np.ndarray, cluster_angles: np.ndarray,
                     widths: np.ndarray, indices: range) -> float:
        """Compute L3 (Minkowski p=3) hyperellipsoid distance.

        Args:
            suite_angles: 7-element array of suite angles [dM, eps, zeta, alpha, beta, gamma, delta]
            cluster_angles: 7-element cluster center
            widths: 7-element width array
            indices: range of indices to use (e.g., range(7) for 7D, range(1,5) for 4D eps-beta)
        """
        total = 0.0
        for k in indices:
            d = self._angle_diff(suite_angles[k], cluster_angles[k])
            total += (d / widths[k]) ** 3
        return total ** (1.0 / 3.0)

    def _classify_pucker(self, delta: float) -> int:
        """Classify sugar pucker: 3 = C3'-endo, 2 = C2'-endo, 0 = unclassified."""
        sieve = self.richardson_suites['sieve']
        c3 = sieve['delta']['C3_endo']
        c2 = sieve['delta']['C2_endo']
        if c3[0] <= delta <= c3[1]:
            return 3
        elif c2[0] <= delta <= c2[1]:
            return 2
        return 0

    def _classify_gamma(self, gamma: float) -> str:
        """Classify gamma angle: 'p', 't', 'm', or '' for unclassified."""
        sieve = self.richardson_suites['sieve']
        for name, rng in sieve['gamma'].items():
            if rng[0] <= gamma <= rng[1]:
                return name
        return ''

    def _score_single_suite(self, suite_angles: np.ndarray) -> dict:
        """Score a single suite (7 angles) against Richardson conformers.

        Args:
            suite_angles: [deltaMinus, epsilon, zeta, alpha, beta, gamma, delta] in 0-360 range

        Returns:
            dict with conformer, suiteness, bin, is_outlier, distance, issue
        """
        if self.richardson_suites is None:
            return None

        result = {
            'conformer': '!!', 'suiteness': 0.0, 'bin': None,
            'is_outlier': True, 'distance': None, 'issue': None
        }

        triage = self.richardson_suites['triage']

        # Triage checks
        eps = suite_angles[1]
        if not (triage['epsilon'][0] <= eps <= triage['epsilon'][1]):
            result['issue'] = f'epsilon={eps:.1f} outside [{triage["epsilon"][0]}, {triage["epsilon"][1]}]'
            return result

        alpha = suite_angles[3]
        if not (triage['alpha'][0] <= alpha <= triage['alpha'][1]):
            result['issue'] = f'alpha={alpha:.1f} outside [{triage["alpha"][0]}, {triage["alpha"][1]}]'
            return result

        beta = suite_angles[4]
        if not (triage['beta'][0] <= beta <= triage['beta'][1]):
            result['issue'] = f'beta={beta:.1f} outside [{triage["beta"][0]}, {triage["beta"][1]}]'
            return result

        zeta = suite_angles[2]
        if not (triage['zeta'][0] <= zeta <= triage['zeta'][1]):
            result['issue'] = f'zeta={zeta:.1f} outside [{triage["zeta"][0]}, {triage["zeta"][1]}]'
            return result

        # Sieve: classify sugar puckers and gamma
        puckerdm = self._classify_pucker(suite_angles[0])  # deltaMinus
        puckerd = self._classify_pucker(suite_angles[6])   # delta
        gamma_class = self._classify_gamma(suite_angles[5]) # gamma

        if puckerdm == 0 or puckerd == 0 or gamma_class == '':
            result['issue'] = f'sieve fail: puckerdm={puckerdm}, puckerd={puckerd}, gamma={gamma_class}'
            return result

        bin_key = f'{puckerdm}{puckerd}{gamma_class}'
        result['bin'] = bin_key

        conformers = self._conformers_by_bin.get(bin_key, [])
        if not conformers:
            result['issue'] = f'no conformers in bin {bin_key}'
            return result

        # 4D screening (epsilon, zeta, alpha, beta = indices 1-4)
        best_4d_dist = float('inf')
        best_4d_conf = None
        for conf in conformers:
            center = np.array(conf['angles'])
            d4 = self._l3_distance(suite_angles, center, self._normal_widths, range(1, 5))
            if d4 < best_4d_dist:
                best_4d_dist = d4
                best_4d_conf = conf

        if best_4d_conf is None or best_4d_dist >= 1.0:
            # No 4D match — outlier
            result['issue'] = f'no 4D match (best={best_4d_dist:.2f})'
            return result

        # 7D distance for suiteness
        center = np.array(best_4d_conf['angles'])
        d7 = self._l3_distance(suite_angles, center, self._normal_widths, range(7))

        if d7 > 1.0:
            # 4D matched but 7D didn't — outlier
            result['issue'] = f'7D reject {best_4d_conf["name"]} (d7={d7:.2f})'
            return result

        # Suiteness: raised cosine
        suiteness = (math.cos(math.pi * d7) + 1.0) / 2.0
        if suiteness < 0.01:
            suiteness = 0.01

        result['conformer'] = best_4d_conf['name']
        result['suiteness'] = round(suiteness, 3)
        result['is_outlier'] = False
        result['distance'] = round(d7, 3)
        return result

    def _score_backbone_suiteness(self, res_1: str, res_2: str,
                                   torsion_data: dict) -> Tuple[dict, float, list]:
        """Score backbone suiteness for both residues in a base pair.

        Uses Richardson suite classification (7 backbone torsion angles spanning
        two adjacent residues) to assess backbone quality.

        Returns:
            (geometry_issues dict, total penalty, backbone_details list of dicts)
        """
        geometry_issues = {}
        penalty = 0.0
        backbone_details = []

        if torsion_data is None or self.richardson_suites is None:
            return geometry_issues, penalty, backbone_details

        suiteness_scores = []

        for res_id in [res_1, res_2]:
            torsions = torsion_data.get(res_id, {})
            if not torsions:
                continue

            # Current residue needs: alpha, beta, gamma, delta
            alpha = torsions.get('alpha')
            beta = torsions.get('beta')
            gamma = torsions.get('gamma')
            delta = torsions.get('delta')

            if any(v is None for v in [alpha, beta, gamma, delta]):
                continue

            # Find predecessor to get deltaMinus, epsilon_prev, zeta_prev
            pred_id = self._find_predecessor(res_id, torsion_data)
            if pred_id is None:
                continue

            pred_torsions = torsion_data.get(pred_id, {})
            delta_prev = pred_torsions.get('delta')
            epsilon_prev = pred_torsions.get('epsilon')
            zeta_prev = pred_torsions.get('zeta')

            if any(v is None for v in [delta_prev, epsilon_prev, zeta_prev]):
                continue

            # Build 7-angle suite, converting from [-180,180] to [0,360]
            suite_angles = np.array([
                delta_prev % 360,
                epsilon_prev % 360,
                zeta_prev % 360,
                alpha % 360,
                beta % 360,
                gamma % 360,
                delta % 360,
            ])

            suite_result = self._score_single_suite(suite_angles)
            if suite_result is None:
                continue

            detail = {
                'residue': res_id,
                'conformer': suite_result['conformer'],
                'suiteness': suite_result['suiteness'],
                'bin': suite_result['bin'],
                'is_outlier': suite_result['is_outlier'],
            }
            if suite_result['distance'] is not None:
                detail['distance'] = suite_result['distance']
            if suite_result['issue'] is not None:
                detail['issue'] = suite_result['issue']

            backbone_details.append(detail)
            suiteness_scores.append(suite_result['suiteness'])

        # Calculate penalty from average suiteness
        if suiteness_scores:
            avg_suiteness = sum(suiteness_scores) / len(suiteness_scores)
            penalty = (1.0 - avg_suiteness) * self.config.PENALTY_WEIGHTS['backbone_suiteness']

            # Flag if any residue is an outlier
            if any(d['is_outlier'] for d in backbone_details):
                geometry_issues['backbone_outlier'] = True

        return geometry_issues, penalty, backbone_details

    def score_structure(self, basepair_data: list, hbond_data: pd.DataFrame,
                        torsion_data: dict = None) -> BaselineResult:
        """
        Score the entire structure based on base pair quality.

        Args:
            basepair_data: List of base pair dictionaries
            hbond_data: DataFrame of hydrogen bonds
            torsion_data: Optional dict of per-residue torsion angles (keyed by residue ID)

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
            'misaligned': 0, 'non_coplanar': 0, 'rotational_distortion': 0,
            'zero_hbond': 0, 'backbone_outlier': 0,
        }
        hbond_stats = {
            'poor_hbond_score': 0, 'bad_distance': 0, 'bad_angles': 0, 'bad_dihedral': 0,
            'incorrect_count': 0
        }
        
        all_suiteness = []
        for bp in basepair_data:
            bp_score = self._score_base_pair(bp, hbond_data, torsion_data)
            basepair_scores.append(bp_score)

            # Collect suiteness scores from backbone details
            for bd in bp_score.get('backbone', []):
                if bd.get('suiteness') is not None:
                    all_suiteness.append(bd['suiteness'])

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
        
        struct_avg_suiteness = round(sum(all_suiteness) / len(all_suiteness), 3) if all_suiteness else None

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
            detailed_issues=detailed_issues,
            avg_suiteness=struct_avg_suiteness,
        )
        
        self._print_results(result)
        
        return result
    
    def _score_base_pair(self, bp: Dict, hbond_data: pd.DataFrame,
                         torsion_data: dict = None) -> Dict:
        """
        Score a single base pair based on geometry, H-bonds, and backbone conformation.
        Uses edge-specific thresholds for geometry and H-bonds, global for dihedrals.

        Args:
            bp: Base pair dictionary with geometry info
            hbond_data: DataFrame of all hydrogen bonds
            torsion_data: Optional dict of per-residue torsion angles

        Returns:
            Dictionary with:
                - score: 0-100 quality score
                - geometry_penalty: penalty from geometry issues
                - hbond_penalty: penalty from H-bond issues
                - geometry_issues: dict of boolean flags
                - hbond_issues: dict of boolean flags
                - backbone: list of per-residue suiteness dicts
                - bp_info: identifier info
        """
        nt1_id = bp.get('res_1', '')
        nt2_id = bp.get('res_2', '')
        
        # Resolve edge and base-pair type (empty/unknown edge -> _OTHER)
        edge_type = bp.get('lw', '') or '_OTHER'
        bp_type = bp.get('bp_type', '')

        # Get base-pair + edge thresholds; fallback chain: (bp, edge) -> (bp, _OTHER) -> (_OTHERS, edge) -> (_OTHERS, _OTHER)
        geo_by_bp = getattr(self.config, 'GEOMETRY_THRESHOLDS_BY_BASE_PAIR_BY_EDGE', {})
        geo_thresh = (geo_by_bp.get(bp_type, {}).get(edge_type) or
                      geo_by_bp.get(bp_type, {}).get('_OTHER') or
                      geo_by_bp.get('_OTHERS', {}).get(edge_type) or
                      geo_by_bp.get('_OTHERS', {}).get('_OTHER') or {})

        hb_by_bp = getattr(self.config, 'HBOND_THRESHOLDS_BY_BASE_PAIR_BY_EDGE', {})
        hbond_thresh = (hb_by_bp.get(bp_type, {}).get(edge_type) or
                        hb_by_bp.get(bp_type, {}).get('_OTHER') or
                        hb_by_bp.get('_OTHERS', {}).get(edge_type) or
                        hb_by_bp.get('_OTHERS', {}).get('_OTHER') or {})
        
        # Get H-bonds for this base pair
        bp_hbonds = self._get_basepair_hbonds(nt1_id, nt2_id, hbond_data)
        
        # Analyze geometry using edge-specific thresholds
        geometry_issues = {}
        geometry_penalty = 0.0
        
        # Check alignment (shear and stretch) - use SHEAR_MIN/MAX range, no abs()
        shear = bp.get('shear', 0)
        stretch = bp.get('stretch', 0)
        shear_min = geo_thresh.get('SHEAR_MIN')
        shear_max = geo_thresh.get('SHEAR_MAX')
        stretch_min = geo_thresh.get('STRETCH_MIN')
        stretch_max = geo_thresh.get('STRETCH_MAX')

        misaligned = False
        if shear_min is not None and shear_max is not None and (shear < shear_min or shear > shear_max):
            misaligned = True
        if stretch_min is not None and stretch_max is not None and (stretch < stretch_min or stretch > stretch_max):
            misaligned = True
        if misaligned:
            geometry_issues['misaligned'] = True
            geometry_penalty += self.config.PENALTY_WEIGHTS['misaligned_pairs']

        # Check coplanarity (buckle and stagger) - use BUCKLE_MIN/MAX, STAGGER_MIN/MAX ranges
        buckle = bp.get('buckle', 0)
        stagger = bp.get('stagger', 0)
        buckle_min = geo_thresh.get('BUCKLE_MIN')
        buckle_max = geo_thresh.get('BUCKLE_MAX')
        stagger_min = geo_thresh.get('STAGGER_MIN')
        stagger_max = geo_thresh.get('STAGGER_MAX')

        non_coplanar = False
        if buckle_min is not None and buckle_max is not None and (buckle < buckle_min or buckle > buckle_max):
            non_coplanar = True
        if stagger_min is not None and stagger_max is not None and (stagger < stagger_min or stagger > stagger_max):
            non_coplanar = True
        if non_coplanar:
            geometry_issues['non_coplanar'] = True
            geometry_penalty += self.config.PENALTY_WEIGHTS['non_coplanar_pairs']

        # Check rotational distortion (propeller and opening)
        propeller = bp.get('propeller', 0)
        opening = bp.get('opening', 0)
        propeller_min = geo_thresh.get('PROPELLER_MIN')
        propeller_max = geo_thresh.get('PROPELLER_MAX')
        opening_min = geo_thresh.get('OPENING_MIN')
        opening_max = geo_thresh.get('OPENING_MAX')

        rot_distorted = False
        if (propeller_min is not None and propeller_max is not None and
                (propeller < propeller_min or propeller > propeller_max)):
            rot_distorted = True
        if (opening_min is not None and opening_max is not None and
                (opening < opening_min or opening > opening_max)):
            rot_distorted = True
        if rot_distorted:
            geometry_issues['rotational_distortion'] = True
            geometry_penalty += self.config.PENALTY_WEIGHTS['rotational_distortion_pairs']

        # Score backbone suiteness (Richardson suite classification)
        backbone_issues, backbone_penalty, backbone_details = self._score_backbone_suiteness(
            nt1_id, nt2_id, torsion_data
        )
        geometry_issues.update(backbone_issues)
        geometry_penalty += backbone_penalty

        # Check DSSR quality score (hbond_score in the JSON)
        dssr_quality = bp.get('hbond_score', 0.0)
        has_dssr_hbonds = dssr_quality > 0.0  # DSSR detected H-bonds if score > 0.0

        # Analyze H-bonds for this base pair
        hbond_issues = {}
        hbond_penalty = 0.0

        # Check DSSR hbond_score: penalize only when below HBOND_SCORE_MIN (not above max)
        # hbond_score is the average quality of all H-bonds in this base pair
        if dssr_quality > 0.0:
            hbond_score_min = geo_thresh.get('HBOND_SCORE_MIN')

            if hbond_score_min is not None and dssr_quality < hbond_score_min:
                hbond_issues['poor_hbond_score'] = True
                hbond_penalty += self.config.PENALTY_WEIGHTS['poor_hbond_score']
        
        # Reconcile DSSR hbond_score with CSV H-bond data
        # Only flag "zero_hbond" if BOTH sources agree there are no H-bonds
        if len(bp_hbonds) == 0:
            # If DSSR says no H-bonds (score == 0.0), then flag zero_hbond
            # Exception: when HBOND_SCORE_MIN and HBOND_SCORE_MAX are both 0, this edge type
            # typically does not form H-bonds (empirical data is all zeros) - do not penalize.
            geo_hb_min = geo_thresh.get('HBOND_SCORE_MIN')
            geo_hb_max = geo_thresh.get('HBOND_SCORE_MAX')
            expects_hbonds = not (geo_hb_min == 0 and geo_hb_max == 0)
            if not has_dssr_hbonds and expects_hbonds:
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
            
            # Check each H-bond using edge-specific thresholds
            # Collect flags first
            has_bad_distance = False
            has_bad_angles = False
            has_bad_dihedral = False

            # Get edge-specific H-bond thresholds (from _OTHER fallback; no hardcoded defaults)
            distance_min = hbond_thresh.get('DIST_MIN')
            distance_max = hbond_thresh.get('DIST_MAX')
            angle_min = hbond_thresh.get('ANGLE_MIN')

            # Dihedral uses config thresholds (structural cis/trans constraint)
            dihedral_cis_min = self.config.HBOND_DIHEDRAL_CIS_MIN
            dihedral_cis_max = self.config.HBOND_DIHEDRAL_CIS_MAX
            dihedral_trans_min = self.config.HBOND_DIHEDRAL_TRANS_MIN
            
            for _, hb in bp_hbonds.iterrows():
                # Distance check (only if thresholds from _OTHER are defined)
                if distance_min is not None and distance_max is not None:
                    distance = hb.get('distance', 0)
                    if distance < distance_min or distance > distance_max:
                        has_bad_distance = True

                # Angle checks (only if threshold from _OTHER is defined)
                if angle_min is not None:
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
        
        # Calculate base pair score
        total_penalty = geometry_penalty + hbond_penalty
        bp_score = max(0, min(100, self.config.BASE_SCORE - total_penalty))
        
        return {
            'score': bp_score,
            'geometry_penalty': geometry_penalty,
            'hbond_penalty': hbond_penalty,
            'geometry_issues': geometry_issues,
            'hbond_issues': hbond_issues,
            'backbone': backbone_details,
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
            'basepair_scores': result.basepair_scores,
            'avg_suiteness': result.avg_suiteness,
        }