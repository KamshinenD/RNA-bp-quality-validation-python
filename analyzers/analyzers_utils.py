"""Shared scoring logic for RNA structure quality assessment."""

import pandas as pd
from typing import Dict, List, Tuple


class BasePairScoring:
    """Core scoring logic for base pair geometry."""
    
    @staticmethod
    def check_misalignment(shear: float, stretch: float, stagger: float, config) -> bool:
        """Check if a base pair is misaligned based on translation parameters.
        
        Args:
            shear: Lateral displacement in Angstroms
            stretch: Extension/compression along long axis in Angstroms
            stagger: Displacement perpendicular to pair plane in Angstroms
            config: Configuration object with thresholds
            
        Returns:
            True if misaligned, False otherwise
        """
        return (abs(shear) > config.SHEAR_MAX or 
                not (config.STRETCH_MIN <= stretch <= config.STRETCH_MAX) or
                abs(stagger) > config.STAGGER_MAX)
    
    @staticmethod
    def check_twist(propeller: float, opening: float, config) -> bool:
        """Check if a base pair is twisted based on rotation parameters.
        Args:
            propeller: Rotation around short axis in degrees
            opening: Out-of-plane rotation in degrees
            config: Configuration object with thresholds
            
        Returns:
            True if twisted, False otherwise
        """
        return (not (config.PROPELLER_MIN <= propeller <= config.PROPELLER_MAX) or
                not (config.OPENING_MIN <= opening <= config.OPENING_MAX))
    
    
    @staticmethod
    def check_adjacent_pairing(res_1: str, res_2: str) -> bool:
        """Check if two adjacent residues are pairing (unusual)."""
        try:
            chain_1 = res_1.split('-')[0]
            chain_2 = res_2.split('-')[0]
            num_1 = int(res_1.split('-')[2])
            num_2 = int(res_2.split('-')[2])
            # Only adjacent if same chain AND consecutive numbers
            return chain_1 == chain_2 and abs(num_1 - num_2) == 1
        except (IndexError, ValueError):
            return False
    
    
    
    @staticmethod
    def check_self_pairing(res_1: str, res_2: str) -> bool:
        """Check if a residue is pairing with itself (serious error)."""
        try:
            chain_1 = res_1.split('-')[0]
            chain_2 = res_2.split('-')[0]
            num_1 = int(res_1.split('-')[2])
            num_2 = int(res_2.split('-')[2])
            # Self-pairing only if same chain AND same residue number
            return chain_1 == chain_2 and num_1 == num_2

        except (IndexError, ValueError):
            return False

    @staticmethod
    def check_coplanarity(buckle: float, config) -> bool:
        """Check if a base pair is non-coplanar.
        
        Args:
            buckle: Rotation around long axis in degrees
            config: Configuration object with thresholds
            
        Returns:
            True if non-coplanar, False otherwise
        """
        return abs(buckle) > config.BUCKLE_MAX
    
    @staticmethod
    def check_hbond_score(hbond_score: float, config) -> bool:
        """Check if a base pair has poor H-bonding.
        
        Args:
            hbond_score: H-bond score from DSSR/FR3D
            config: Configuration object with thresholds
            
        Returns:
            True if poor H-bonding, False otherwise or zero since we are checking for zero score independently already
        """
        return 0.0 < hbond_score < config.HBOND_SCORE_MIN
    
    @staticmethod
    def check_zero_hbond_score(hbond_score: float) -> bool:
        """Check if a base pair has zero H-bond score (complete absence of H-bonding)."""
        return hbond_score == 0.0

    
    @staticmethod
    def score_base_pair(bp: Dict, config) -> Dict:
        """Score a single base pair and identify its issues.
        
        Args:
            bp: Base pair dictionary with geometry parameters
            config: Configuration object
            
        Returns:
            Dictionary with boolean flags for each issue type
        """
        
        # Extract residues for adjacency check
        res_1 = bp.get('res_1', '')
        res_2 = bp.get('res_2', '')
        
        if BasePairScoring.check_adjacent_pairing(res_1, res_2):
            return None  # Skip adjacent pairs
        
        issues = {
            'misaligned': False,
            'twisted': False,
            'non_coplanar': False,
            'low_hbond_score': False,
            'zero_hbond': False,
            # 'adjacent_pairing': False,
            'self_pairing': False
        }
        
        # Extract parameters
        shear = bp.get('shear', 0)
        stretch = bp.get('stretch', 0)
        stagger = bp.get('stagger', 0)
        buckle = bp.get('buckle', 0)
        propeller = bp.get('propeller', 0)
        opening = bp.get('opening', 0)
        hbond_score = bp.get('hbond_score', 0)
        
        # Check each issue type
        issues['misaligned'] = BasePairScoring.check_misalignment(
            shear, stretch, stagger, config
        )
        issues['twisted'] = BasePairScoring.check_twist(
            propeller, opening, config
        )
        issues['non_coplanar'] = BasePairScoring.check_coplanarity(
            buckle, config
        )
        issues['low_hbond_score'] = BasePairScoring.check_hbond_score(
            hbond_score, config
        )
        issues['zero_hbond'] = BasePairScoring.check_zero_hbond_score(
            hbond_score
        )
        # issues['adjacent_pairing'] = BasePairScoring.check_adjacent_pairing(
        #     res_1, res_2
        # )
        issues['self_pairing'] = BasePairScoring.check_self_pairing(
            res_1, res_2
        )

        return issues
    
    @staticmethod
    def calculate_basepair_penalty(stats: Dict, total_pairs: int, config) -> float:
        """Calculate penalty score for base pairs.
        Returns:
            Total penalty score
        """
        if total_pairs == 0:
            return 0.0
        
        penalty = 0.0
        
        # Calculate fractions and apply weights
        if 'misaligned' in stats:
            misaligned_frac = stats['misaligned'] / total_pairs
            penalty += config.PENALTY_WEIGHTS['misaligned_pairs'] * misaligned_frac
        
        if 'twisted' in stats:
            twisted_frac = stats['twisted'] / total_pairs
            penalty += config.PENALTY_WEIGHTS['twisted_pairs'] * twisted_frac
        
        if 'non_coplanar' in stats:
            non_coplanar_frac = stats['non_coplanar'] / total_pairs
            penalty += config.PENALTY_WEIGHTS['non_coplanar_pairs'] * non_coplanar_frac
        
        if 'low_hbond_score' in stats:
            low_hbond_score_frac = stats['low_hbond_score'] / total_pairs
            penalty += config.PENALTY_WEIGHTS['low_hbond_score_pairs'] * low_hbond_score_frac
            
        if 'zero_hbond' in stats:
            zero_hbond_frac = stats['zero_hbond'] / total_pairs
            penalty += config.PENALTY_WEIGHTS['zero_hbond_pairs'] * zero_hbond_frac

            
        # if 'adjacent_pairing' in stats:
        #     adjacent_frac = stats['adjacent_pairing'] / total_pairs
        #     penalty += config.PENALTY_WEIGHTS['adjacent_pairing'] * adjacent_frac

        if 'self_pairing' in stats:
            selfpair_frac = stats['self_pairing'] / total_pairs
            penalty += config.PENALTY_WEIGHTS['self_pairing'] * selfpair_frac
    
        return penalty


class HBondScoring:
    """Core scoring logic for hydrogen bonds."""
    
    @staticmethod
    def check_distance(distance: float, config) -> bool:
        """Check if H-bond distance is outside acceptable range.
        
        Args:
            distance: H-bond distance in Angstroms
            config: Configuration object with thresholds
            
        Returns:
            True if bad distance, False otherwise
        """
        return not (config.HBOND_DISTANCE_MIN <= distance <= config.HBOND_DISTANCE_MAX)
    
    @staticmethod
    def check_angles(angle_1: float, angle_2: float, config) -> bool:
        """Check if H-bond angles are below minimum.
        
        Args:
            angle_1, angle_2: H-bond angles in degrees
            config: Configuration object with thresholds
            
        Returns:
            True if bad angles, False otherwise
        """
        return angle_1 < config.HBOND_ANGLE_MIN or angle_2 < config.HBOND_ANGLE_MIN
    
    @staticmethod
    def check_dihedral(dihedral: float, config) -> Tuple[bool, float]:
        """Check if dihedral angle is in forbidden zone.
        
        The "forbidden zone" (50° to 140° and -50° to -140°) represents 
        strained or unusual H-bond geometries.
        
        Args:
            dihedral: Dihedral angle in degrees
            config: Configuration object with thresholds
            
        Returns:
            Tuple of (is_bad, deviation_from_nearest_valid)
        """
        is_cis = config.HBOND_DIHEDRAL_CIS_MIN <= dihedral <= config.HBOND_DIHEDRAL_CIS_MAX
        is_trans = abs(dihedral) >= config.HBOND_DIHEDRAL_TRANS_MIN
        
        if is_cis or is_trans:
            return False, 0.0
        
        # Calculate distance to nearest acceptable range
        dist_to_cis_upper = abs(dihedral - config.HBOND_DIHEDRAL_CIS_MAX)
        dist_to_cis_lower = abs(dihedral - config.HBOND_DIHEDRAL_CIS_MIN)
        dist_to_trans = abs(abs(dihedral) - config.HBOND_DIHEDRAL_TRANS_MIN)
        
        deviation = min(dist_to_cis_upper, dist_to_cis_lower, dist_to_trans)
        
        return True, deviation
    
    @staticmethod
    def check_quality(quality_score: float, config) -> bool:
        """Check if H-bond quality score is weak.
        
        Args:
            quality_score: Overall H-bond quality score
            config: Configuration object with thresholds
            
        Returns:
            True if weak quality, False otherwise
        """
        return quality_score < config.HBOND_QUALITY_MIN
    
    @staticmethod
    def check_hbond_count(bp_type: str, actual_count: int, lw: str, config) -> Tuple[bool, str]:
        """Check if H-bond count matches expectations for base pair type.
        
        Args:
            bp_type: Base pair type (e.g., 'G-C', 'A-U')
            actual_count: Actual number of H-bonds
            lw: Leontis-Westhof notation
            config: Configuration object with expected counts
            
        Returns:
            Tuple of (is_incorrect, issue_description)
        """
        if bp_type not in config.EXPECTED_HBOND_COUNTS:
            return False, ""
        
        min_expected, max_expected, ideal = config.EXPECTED_HBOND_COUNTS[bp_type]
        
        # Check if count is outside expected range
        if actual_count < min_expected or actual_count > max_expected:
            return True, f"Wrong count: found {actual_count}, expected {min_expected}-{max_expected}"
        
        # For canonical Watson-Crick pairs, flag if not ideal
        if actual_count != ideal and lw == 'cWW':
            return True, f"Suboptimal cWW: found {actual_count}, expected {ideal}"
        
        return False, ""
    
    @staticmethod
    def score_hbond(hbond: pd.Series, config) -> Dict:
        """Score a single H-bond and identify its issues.
        
        Args:
            hbond: H-bond data as pandas Series
            config: Configuration object
            
        Returns:
            Dictionary with boolean flags for each issue type
        """
        issues = {
            'bad_distance': False,
            'bad_angles': False,
            'bad_dihedral': False,
            'weak_quality': False
        }
        
        # Check distance
        issues['bad_distance'] = HBondScoring.check_distance(
            hbond['distance'], config
        )
        
        # Check angles
        issues['bad_angles'] = HBondScoring.check_angles(
            hbond['angle_1'], hbond['angle_2'], config
        )
        
        # Check dihedral
        issues['bad_dihedral'], _ = HBondScoring.check_dihedral(
            hbond['dihedral_angle'], config
        )
        
        # Check quality
        issues['weak_quality'] = HBondScoring.check_quality(
            hbond['score'], config
        )
        
        return issues
    
    @staticmethod
    def calculate_hbond_penalty(stats: Dict, total_hbonds: int, 
                               total_pairs_checked: int, config) -> float:
        """Calculate penalty score for H-bonds.
        
        Args:
            stats: Dictionary with counts of each issue type
            total_hbonds: Total number of H-bonds
            total_pairs_checked: Total base pairs checked for H-bond count
            config: Configuration object with penalty weights
            
        Returns:
            Total penalty score
        """
        if total_hbonds == 0:
            return 0.0
        
        penalty = 0.0
        
        # Calculate fractions and apply weights
        if 'bad_distance' in stats:
            bad_distance_frac = stats['bad_distance'] / total_hbonds
            penalty += config.PENALTY_WEIGHTS['bad_hbond_distance'] * bad_distance_frac
        
        if 'bad_angles' in stats:
            bad_angles_frac = stats['bad_angles'] / total_hbonds
            penalty += config.PENALTY_WEIGHTS['bad_hbond_angles'] * bad_angles_frac
        
        if 'bad_dihedral' in stats:
            bad_dihedral_frac = stats['bad_dihedral'] / total_hbonds
            penalty += config.PENALTY_WEIGHTS['bad_hbond_dihedrals'] * bad_dihedral_frac
        
        if 'weak_quality' in stats:
            weak_quality_frac = stats['weak_quality'] / total_hbonds
            penalty += config.PENALTY_WEIGHTS['weak_hbond_quality'] * weak_quality_frac
        
        # H-bond count penalty uses different denominator
        if 'incorrect_hbond_count' in stats and total_pairs_checked > 0:
            incorrect_count_frac = stats['incorrect_hbond_count'] / total_pairs_checked
            penalty += config.PENALTY_WEIGHTS['incorrect_hbond_count'] * incorrect_count_frac
        
        return penalty


class RegionScoring:
    """Core scoring logic for regions/hotspots."""
    
    @staticmethod
    def calculate_region_score(bp_penalty: float, hb_penalty: float, config) -> float:
        """Calculate overall quality score for a region.
        
        Args:
            bp_penalty: Penalty from base pair issues
            hb_penalty: Penalty from H-bond issues
            config: Configuration object with base score
            
        Returns:
            Quality score for the region
        """
        return config.BASE_SCORE - bp_penalty - hb_penalty
    
    @staticmethod
    def calculate_residue_score(residue_issues: Dict, config) -> float:
        """Calculate quality score for a single residue.
        
        Args:
            residue_issues: Dictionary with issue fractions for the residue
            config: Configuration object
            
        Returns:
            Quality score for the residue
        """
        penalty = 0.0
        
        # Map issue types to penalty weights
        issue_weight_map = {
            'misaligned': 'misaligned_pairs',
            'twisted': 'twisted_pairs',
            'non_coplanar': 'non_coplanar_pairs',
            'low_hbond_score': 'low_hbond_score_pairs',
            # 'adjacent_pairing': 'adjacent_pairing',
            'self_pairing': 'self_pairing',
            'bad_distance': 'bad_hbond_distance',
            'bad_angles': 'bad_hbond_angles',
            'bad_dihedral': 'bad_hbond_dihedrals',
            'weak_quality': 'weak_hbond_quality',
            'incorrect_count': 'incorrect_hbond_count',
            'zero_hbond': 'zero_hbond_pairs',
            
        }
        
        for issue_type, weight_key in issue_weight_map.items():
            if issue_type in residue_issues and weight_key in config.PENALTY_WEIGHTS:
                penalty += config.PENALTY_WEIGHTS[weight_key] * residue_issues[issue_type]
        
        return config.BASE_SCORE - penalty


class ScoringUtils:
    """Additional utility functions for scoring."""
    
    
    @staticmethod
    def is_base_atom(atom_name: str) -> bool:
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

    @staticmethod
    def is_base_base_hbond(atom_1: str, atom_2: str) -> bool:
        """
        Check if H-bond is between two base atoms (not backbone/sugar).
        
        Args:
            atom_1: First atom name
            atom_2: Second atom name
            
        Returns:
            True if both atoms are base atoms
        """
        return (ScoringUtils.is_base_atom(atom_1) and 
                ScoringUtils.is_base_atom(atom_2))
        
    @staticmethod
    def build_base_pair_set(base_pairs: list) -> set:
        """
        Build a set of residue pairs that form base pairs.
        
        Args:
            base_pairs: List of base pair dictionaries
            
        Returns:
            Set of tuples representing base-paired residue pairs
        """
        base_pair_set = set()
        for bp in base_pairs:
            pair_key = tuple(sorted([bp['res_1'], bp['res_2']]))
            base_pair_set.add(pair_key)
        return base_pair_set    
    
    
    @staticmethod
    def extract_base_pair_type(res_1: str, res_2: str) -> str:
        """Extract base pair type from residue IDs.
        
        Format: CHAIN-BASE-NUMBER-
        Example: 'Q-G-22-', 'Q-C-45-' -> 'G-C'
        """
        parts1 = res_1.split('-')
        parts2 = res_2.split('-')
        
        base1 = parts1[1] if len(parts1) >= 2 else '?'
        base2 = parts2[1] if len(parts2) >= 2 else '?'
        
        return f"{base1}-{base2}"
    
    @staticmethod
    def count_hbonds_per_pair(hbond_df: pd.DataFrame) -> Dict[tuple, int]:
        """Count H-bonds between each pair of residues."""
        pair_counts = {}
        
        for _, row in hbond_df.iterrows():
            res_1 = row['res_1']
            res_2 = row['res_2']
            
            # Create canonical key (sorted to handle both directions)
            pair_key = tuple(sorted([res_1, res_2]))
            pair_counts[pair_key] = pair_counts.get(pair_key, 0) + 1
        
        return pair_counts