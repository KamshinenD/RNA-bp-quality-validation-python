"""Configuration file with all tunable thresholds and scoring parameters."""

class Config:
    """RNA quality scoring configuration."""
    
    BASE_SCORE = 100.0
    
#     PENALTY_WEIGHTS = {
#     'non_coplanar_pairs': 15.0,
#     'poor_hbond_pairs': 10.0,
#     'zero_hbond_pairs': 10.0,
#     'bad_hbond_distance': 15.0,
#     'bad_hbond_angles': 10.0,
    
#     'incorrect_hbond_count': 9.0,
    
#     'misaligned_pairs': 6.0,
#     'twisted_pairs': 7.0,
#     'bad_hbond_dihedrals': 9.0,
#     'weak_hbond_quality': 9.0,
    
#     'adjacent_pairing': 0.0,
#     'self_pairing': 0.0,
# }
    
    PENALTY_WEIGHTS = {
    # CRITICAL FAILURES (rare in good structures)
    'zero_hbond_pairs': 20.0,           # Almost never happens in good pairs
    'non_coplanar_pairs': 15.0,         # Severe stacking disruption
    
    # MAJOR DEFECTS (uncommon in good structures)
    'bad_hbond_distance': 15.0,         # Distance outside 2.3-3.7Å is serious
    'misaligned_pairs': 15.0,           # Shear/stretch problems
    'twisted_pairs': 10.0,              # Buckle/propeller issues
    
    # MODERATE ISSUES (somewhat common, especially in non-WC)
    'bad_hbond_dihedrals': 5.0,         # Dihedral problems
    'poor_hbond_pairs': 5.0,            # DSSR score < 2.2 (reduce - too strict!)
    'weak_hbond_quality': 5.0,          # Individual H-bond quality < 0.70
    
    # MINOR ISSUES (common even in good structures)
    'bad_hbond_angles': 5.0,            # One angle < 80° is very common
    'incorrect_hbond_count': 5.0,       # With tolerant ranges, less critical

    # MODERATE - Chi glycosidic bond angle in unexpected conformation
    'chi_conformation': 10.0,
}
    
    # ===== EXPECTED H-BOND COUNTS BY BASE PAIR TYPE =====
    # Format: (min_expected, max_expected, ideal_count)
    # - min_expected: Minimum acceptable H-bond count
    # - max_expected: Maximum acceptable H-bond count
    # - ideal_count: Ideal/optimal H-bond count
    # 
    # Usage (applies to both motifs and full RNA analysis):
    # 1. Range-based validation: Flags pairs if count < min_expected or count > max_expected
    # 2. Strict canonical checking: For cWW (canonical Watson-Crick) pairs, also flags if count != ideal_count
    # Used in scorer2.py (for scoring) and analyzers_utils.py (for analysis)
    EXPECTED_HBOND_COUNTS = {
        'G-C': (3, 3, 3),
        'C-G': (3, 3, 3),
        'A-U': (2, 2, 2),
        'U-A': (2, 2, 2),
        'G-U': (2, 2, 2),
        'U-G': (2, 2, 2),
        'A-A': (1, 2, 2),
        'G-G': (1, 2, 2),
        'U-U': (1, 2, 2),
        'C-C': (1, 2, 1),
        'A-C': (1, 2, 1),
        'C-A': (1, 2, 1),
        'A-G': (1, 2, 2),
        'G-A': (1, 2, 2),
        'U-C': (1, 2, 1),
        'C-U': (1, 2, 1),
    }
    
    # ===== GEOMETRY THRESHOLDS =====
    # OLD VALUES (commented out):
    # SHEAR_MAX = 1.2              # Very relaxed
    # STRETCH_MIN = -0.7           # Very tolerant
    # STRETCH_MAX = 0.3            # Very tolerant
    # STAGGER_MAX = 1.0            # Very relaxed
    # BUCKLE_MAX = 30.0            # Very relaxed
    # PROPELLER_MIN = -30.0        # Very relaxed
    # PROPELLER_MAX = 10.0         # Very relaxed
    # OPENING_MIN = -20.0          # Very relaxed
    # OPENING_MAX = 20.0           # Very relaxed
    # HBOND_SCORE_MIN = 2.2        # Higher threshold = more sensitive
    
    # NEW VALUES - DATA-DRIVEN (MIXED: 0.5 SD for translation, 1 SD for rotation)
    # Based on analysis of 3,794,390 base pairs
    # Translation properties (shear, stretch, stagger): 0.5 SD from ideal 0 (~38% coverage)
    # Rotation properties (buckle, propeller, opening): 1 SD from ideal 0 (~68% coverage)
    
    # OLD GLOBAL THRESHOLDS (COMMENTED OUT - NOW USING EDGE-BASED THRESHOLDS)
    # SHEAR_MAX = 1.5              # 0.5*std = 0.5*3.002 = 1.501 (translation property)
    # STRETCH_MIN = -1.3           # -0.5*std = -0.5*2.663 = -1.332 (translation property)
    # STRETCH_MAX = 1.3            # +0.5*std = 0.5*2.663 = 1.332 (translation property)
    # STAGGER_MAX = 0.4            # 0.5*std = 0.5*0.779 = 0.390 (translation property)
    # BUCKLE_MAX = 16.1            # 1*std = 16.058 (rotation property)
    # PROPELLER_MIN = -14.7        # -1*std = -14.680 (rotation property)
    # PROPELLER_MAX = 14.7         # +1*std = 14.680 (rotation property)
    # OPENING_MIN = -30.8          # -0.5*std = -0.5*61.517 = -30.759 (rotation property, stricter)
    # OPENING_MAX = 30.8           # +0.5*std = 0.5*61.517 = 30.759 (rotation property, stricter)
    
    HBOND_SCORE_MIN = 2.0        # Mean = 2.007 (flags pairs below average H-bond quality)
    
    # EDGE-BASED GEOMETRY THRESHOLDS 
    # Thresholds: ±1 SD from ideal 0 for translation/rotation parameters
    # Edges with < 1,000 pairs are grouped as 'uncommon_edges'
    GEOMETRY_THRESHOLDS_BY_EDGE = {
        '--': {
            'SHEAR_MAX': 2.75,
            'STRETCH_MIN': -3.89,
            'STRETCH_MAX': 3.89,
            'STAGGER_MAX': 0.99,
            'BUCKLE_MAX': 17.96,
            'PROPELLER_MIN': -20.78,
            'PROPELLER_MAX': 20.78,
            'OPENING_MIN': -110.11,
            'OPENING_MAX': 110.11,
        },
        '..H': {
            'SHEAR_MAX': 2.96,
            'STRETCH_MIN': -6.49,
            'STRETCH_MAX': 6.49,
            'STAGGER_MAX': 1.54,
            'BUCKLE_MAX': 26.22,
            'PROPELLER_MIN': -27.32,
            'PROPELLER_MAX': 27.32,
            'OPENING_MIN': -101.11,
            'OPENING_MAX': 101.11,
        },
        '..S': {
            'SHEAR_MAX': 5.12,
            'STRETCH_MIN': -4.29,
            'STRETCH_MAX': 4.29,
            'STAGGER_MAX': 0.50,
            'BUCKLE_MAX': 15.06,
            'PROPELLER_MIN': -11.78,
            'PROPELLER_MAX': 11.78,
            'OPENING_MIN': -37.81,
            'OPENING_MAX': 37.81,
        },
        '.H.': {
            'SHEAR_MAX': 4.28,
            'STRETCH_MIN': -3.86,
            'STRETCH_MAX': 3.86,
            'STAGGER_MAX': 0.76,
            'BUCKLE_MAX': 20.77,
            'PROPELLER_MIN': -23.49,
            'PROPELLER_MAX': 23.49,
            'OPENING_MIN': -155.15,
            'OPENING_MAX': 155.15,
        },
        '.HH': {
            'SHEAR_MAX': 4.78,
            'STRETCH_MIN': -7.14,
            'STRETCH_MAX': 7.14,
            'STAGGER_MAX': 1.46,
            'BUCKLE_MAX': 22.58,
            'PROPELLER_MIN': -21.06,
            'PROPELLER_MAX': 21.06,
            'OPENING_MIN': -107.70,
            'OPENING_MAX': 107.70,
        },
        '.HS': {
            'SHEAR_MAX': 7.31,
            'STRETCH_MIN': -3.99,
            'STRETCH_MAX': 3.99,
            'STAGGER_MAX': 0.59,
            'BUCKLE_MAX': 9.46,
            'PROPELLER_MIN': -6.68,
            'PROPELLER_MAX': 6.68,
            'OPENING_MIN': -55.87,
            'OPENING_MAX': 55.87,
        },
        '.HW': {
            'SHEAR_MAX': 1.42,
            'STRETCH_MIN': -3.72,
            'STRETCH_MAX': 3.72,
            'STAGGER_MAX': 0.97,
            'BUCKLE_MAX': 11.96,
            'PROPELLER_MIN': -32.38,
            'PROPELLER_MAX': 32.38,
            'OPENING_MIN': -83.13,
            'OPENING_MAX': 83.13,
        },
        '.SS': {
            'SHEAR_MAX': 3.27,
            'STRETCH_MIN': -7.60,
            'STRETCH_MAX': 7.60,
            'STAGGER_MAX': 1.25,
            'BUCKLE_MAX': 35.13,
            'PROPELLER_MIN': -16.16,
            'PROPELLER_MAX': 16.16,
            'OPENING_MIN': -161.17,
            'OPENING_MAX': 161.17,
        },
        '.SW': {
            'SHEAR_MAX': 2.32,
            'STRETCH_MIN': -3.37,
            'STRETCH_MAX': 3.37,
            'STAGGER_MAX': 0.54,
            'BUCKLE_MAX': 25.58,
            'PROPELLER_MIN': -10.73,
            'PROPELLER_MAX': 10.73,
            'OPENING_MIN': -73.74,
            'OPENING_MAX': 73.74,
        },
        '.W.': {
            'SHEAR_MAX': 1.56,
            'STRETCH_MIN': -2.14,
            'STRETCH_MAX': 2.14,
            'STAGGER_MAX': 0.31,
            'BUCKLE_MAX': 13.50,
            'PROPELLER_MIN': -42.11,
            'PROPELLER_MAX': 42.11,
            'OPENING_MIN': -77.57,
            'OPENING_MAX': 77.57,
        },
        '.WH': {
            'SHEAR_MAX': 4.36,
            'STRETCH_MIN': -2.76,
            'STRETCH_MAX': 2.76,
            'STAGGER_MAX': 1.03,
            'BUCKLE_MAX': 14.32,
            'PROPELLER_MIN': -16.90,
            'PROPELLER_MAX': 16.90,
            'OPENING_MIN': -116.66,
            'OPENING_MAX': 116.66,
        },
        '.WS': {
            'SHEAR_MAX': 2.85,
            'STRETCH_MIN': -3.94,
            'STRETCH_MAX': 3.94,
            'STAGGER_MAX': 1.02,
            'BUCKLE_MAX': 21.90,
            'PROPELLER_MIN': -11.75,
            'PROPELLER_MAX': 11.75,
            'OPENING_MIN': -101.47,
            'OPENING_MAX': 101.47,
        },
        '.WW': {
            'SHEAR_MAX': 0.73,
            'STRETCH_MIN': -1.26,
            'STRETCH_MAX': 1.26,
            'STAGGER_MAX': 0.46,
            'BUCKLE_MAX': 13.99,
            'PROPELLER_MIN': -14.62,
            'PROPELLER_MAX': 14.62,
            'OPENING_MIN': -74.98,
            'OPENING_MAX': 74.98,
        },
        'c.H': {
            'SHEAR_MAX': 5.09,
            'STRETCH_MIN': -3.13,
            'STRETCH_MAX': 3.13,
            'STAGGER_MAX': 1.14,
            'BUCKLE_MAX': 24.01,
            'PROPELLER_MIN': -22.50,
            'PROPELLER_MAX': 22.50,
            'OPENING_MIN': -90.99,
            'OPENING_MAX': 90.99,
        },
        'c.S': {
            'SHEAR_MAX': 3.68,
            'STRETCH_MIN': -5.28,
            'STRETCH_MAX': 5.28,
            'STAGGER_MAX': 1.82,
            'BUCKLE_MAX': 39.63,
            'PROPELLER_MIN': -22.55,
            'PROPELLER_MAX': 22.55,
            'OPENING_MIN': -122.56,
            'OPENING_MAX': 122.56,
        },
        'c.W': {
            'SHEAR_MAX': 4.50,
            'STRETCH_MIN': -1.97,
            'STRETCH_MAX': 1.97,
            'STAGGER_MAX': 1.29,
            'BUCKLE_MAX': 24.75,
            'PROPELLER_MIN': -20.46,
            'PROPELLER_MAX': 20.46,
            'OPENING_MIN': -58.45,
            'OPENING_MAX': 58.45,
        },
        'cH.': {
            'SHEAR_MAX': 3.72,
            'STRETCH_MIN': -4.75,
            'STRETCH_MAX': 4.75,
            'STAGGER_MAX': 1.03,
            'BUCKLE_MAX': 26.03,
            'PROPELLER_MIN': -18.07,
            'PROPELLER_MAX': 18.07,
            'OPENING_MIN': -130.71,
            'OPENING_MAX': 130.71,
        },
        'cHH': {
            'SHEAR_MAX': 3.06,
            'STRETCH_MIN': -5.89,
            'STRETCH_MAX': 5.89,
            'STAGGER_MAX': 1.13,
            'BUCKLE_MAX': 29.60,
            'PROPELLER_MIN': -22.16,
            'PROPELLER_MAX': 22.16,
            'OPENING_MIN': -154.02,
            'OPENING_MAX': 154.02,
        },
        'cHS': {
            'SHEAR_MAX': 6.46,
            'STRETCH_MIN': -1.88,
            'STRETCH_MAX': 1.88,
            'STAGGER_MAX': 1.44,
            'BUCKLE_MAX': 29.52,
            'PROPELLER_MIN': -27.76,
            'PROPELLER_MAX': 27.76,
            'OPENING_MIN': -28.99,
            'OPENING_MAX': 28.99,
        },
        'cHW': {
            'SHEAR_MAX': 2.13,
            'STRETCH_MIN': -4.04,
            'STRETCH_MAX': 4.04,
            'STAGGER_MAX': 0.93,
            'BUCKLE_MAX': 21.53,
            'PROPELLER_MIN': -19.13,
            'PROPELLER_MAX': 19.13,
            'OPENING_MIN': -85.56,
            'OPENING_MAX': 85.56,
        },
        'cS.': {
            'SHEAR_MAX': 4.76,
            'STRETCH_MIN': -4.60,
            'STRETCH_MAX': 4.60,
            'STAGGER_MAX': 1.77,
            'BUCKLE_MAX': 41.49,
            'PROPELLER_MIN': -23.80,
            'PROPELLER_MAX': 23.80,
            'OPENING_MIN': -103.05,
            'OPENING_MAX': 103.05,
        },
        'cSH': {
            'SHEAR_MAX': 7.01,
            'STRETCH_MIN': -1.79,
            'STRETCH_MAX': 1.79,
            'STAGGER_MAX': 1.14,
            'BUCKLE_MAX': 21.38,
            'PROPELLER_MIN': -21.48,
            'PROPELLER_MAX': 21.48,
            'OPENING_MIN': -22.55,
            'OPENING_MAX': 22.55,
        },
        'cSS': {
            'SHEAR_MAX': 2.35,
            'STRETCH_MIN': -7.23,
            'STRETCH_MAX': 7.23,
            'STAGGER_MAX': 1.28,
            'BUCKLE_MAX': 23.74,
            'PROPELLER_MIN': -20.86,
            'PROPELLER_MAX': 20.86,
            'OPENING_MIN': -152.05,
            'OPENING_MAX': 152.05,
        },
        'cSW': {
            'SHEAR_MAX': 5.01,
            'STRETCH_MIN': -3.18,
            'STRETCH_MAX': 3.18,
            'STAGGER_MAX': 1.18,
            'BUCKLE_MAX': 26.42,
            'PROPELLER_MIN': -26.30,
            'PROPELLER_MAX': 26.30,
            'OPENING_MIN': -88.84,
            'OPENING_MAX': 88.84,
        },
        'cW.': {
            'SHEAR_MAX': 4.27,
            'STRETCH_MIN': -2.30,
            'STRETCH_MAX': 2.30,
            'STAGGER_MAX': 1.13,
            'BUCKLE_MAX': 24.70,
            'PROPELLER_MIN': -19.05,
            'PROPELLER_MAX': 19.05,
            'OPENING_MIN': -69.46,
            'OPENING_MAX': 69.46,
        },
        'cWH': {
            'SHEAR_MAX': 2.36,
            'STRETCH_MIN': -3.39,
            'STRETCH_MAX': 3.39,
            'STAGGER_MAX': 1.16,
            'BUCKLE_MAX': 25.24,
            'PROPELLER_MIN': -18.47,
            'PROPELLER_MAX': 18.47,
            'OPENING_MIN': -75.41,
            'OPENING_MAX': 75.41,
        },
        'cWS': {
            'SHEAR_MAX': 5.16,
            'STRETCH_MIN': -3.10,
            'STRETCH_MAX': 3.10,
            'STAGGER_MAX': 1.29,
            'BUCKLE_MAX': 28.74,
            'PROPELLER_MIN': -23.40,
            'PROPELLER_MAX': 23.40,
            'OPENING_MIN': -88.94,
            'OPENING_MAX': 88.94,
        },
        'cWW': {
            'SHEAR_MAX': 1.08,
            'STRETCH_MIN': -0.51,
            'STRETCH_MAX': 0.51,
            'STAGGER_MAX': 0.56,
            'BUCKLE_MAX': 10.75,
            'PROPELLER_MIN': -11.97,
            'PROPELLER_MAX': 11.97,
            'OPENING_MIN': -10.19,
            'OPENING_MAX': 10.19,
        },
        't.H': {
            'SHEAR_MAX': 6.04,
            'STRETCH_MIN': -3.95,
            'STRETCH_MAX': 3.95,
            'STAGGER_MAX': 0.94,
            'BUCKLE_MAX': 19.73,
            'PROPELLER_MIN': -16.88,
            'PROPELLER_MAX': 16.88,
            'OPENING_MIN': -81.82,
            'OPENING_MAX': 81.82,
        },
        't.S': {
            'SHEAR_MAX': 3.31,
            'STRETCH_MIN': -6.52,
            'STRETCH_MAX': 6.52,
            'STAGGER_MAX': 1.88,
            'BUCKLE_MAX': 45.99,
            'PROPELLER_MIN': -12.54,
            'PROPELLER_MAX': 12.54,
            'OPENING_MIN': -137.84,
            'OPENING_MAX': 137.84,
        },
        't.W': {
            'SHEAR_MAX': 2.86,
            'STRETCH_MIN': -4.16,
            'STRETCH_MAX': 4.16,
            'STAGGER_MAX': 1.20,
            'BUCKLE_MAX': 24.90,
            'PROPELLER_MIN': -17.83,
            'PROPELLER_MAX': 17.83,
            'OPENING_MIN': -117.41,
            'OPENING_MAX': 117.41,
        },
        'tH.': {
            'SHEAR_MAX': 5.72,
            'STRETCH_MIN': -4.05,
            'STRETCH_MAX': 4.05,
            'STAGGER_MAX': 1.05,
            'BUCKLE_MAX': 16.13,
            'PROPELLER_MIN': -14.56,
            'PROPELLER_MAX': 14.56,
            'OPENING_MIN': -88.29,
            'OPENING_MAX': 88.29,
        },
        'tHH': {
            'SHEAR_MAX': 5.99,
            'STRETCH_MIN': -5.30,
            'STRETCH_MAX': 5.30,
            'STAGGER_MAX': 0.74,
            'BUCKLE_MAX': 19.91,
            'PROPELLER_MIN': -17.91,
            'PROPELLER_MAX': 17.91,
            'OPENING_MIN': -166.89,
            'OPENING_MAX': 166.89,
        },
        'tHS': {
            'SHEAR_MAX': 6.76,
            'STRETCH_MIN': -4.37,
            'STRETCH_MAX': 4.37,
            'STAGGER_MAX': 0.69,
            'BUCKLE_MAX': 14.71,
            'PROPELLER_MIN': -15.09,
            'PROPELLER_MAX': 15.09,
            'OPENING_MIN': -16.03,
            'OPENING_MAX': 16.03,
        },
        'tHW': {
            'SHEAR_MAX': 4.54,
            'STRETCH_MIN': -2.08,
            'STRETCH_MAX': 2.08,
            'STAGGER_MAX': 0.98,
            'BUCKLE_MAX': 16.18,
            'PROPELLER_MIN': -18.54,
            'PROPELLER_MAX': 18.54,
            'OPENING_MIN': -90.78,
            'OPENING_MAX': 90.78,
        },
        'tS.': {
            'SHEAR_MAX': 4.25,
            'STRETCH_MIN': -6.19,
            'STRETCH_MAX': 6.19,
            'STAGGER_MAX': 1.68,
            'BUCKLE_MAX': 41.81,
            'PROPELLER_MIN': -19.18,
            'PROPELLER_MAX': 19.18,
            'OPENING_MIN': -123.92,
            'OPENING_MAX': 123.92,
        },
        'tSH': {
            'SHEAR_MAX': 6.88,
            'STRETCH_MIN': -4.86,
            'STRETCH_MAX': 4.86,
            'STAGGER_MAX': 0.86,
            'BUCKLE_MAX': 17.61,
            'PROPELLER_MIN': -15.76,
            'PROPELLER_MAX': 15.76,
            'OPENING_MIN': -26.47,
            'OPENING_MAX': 26.47,
        },
        'tSS': {
            'SHEAR_MAX': 2.86,
            'STRETCH_MIN': -7.59,
            'STRETCH_MAX': 7.59,
            'STAGGER_MAX': 1.40,
            'BUCKLE_MAX': 29.58,
            'PROPELLER_MIN': -22.75,
            'PROPELLER_MAX': 22.75,
            'OPENING_MIN': -156.71,
            'OPENING_MAX': 156.71,
        },
        'tSW': {
            'SHEAR_MAX': 2.75,
            'STRETCH_MIN': -5.88,
            'STRETCH_MAX': 5.88,
            'STAGGER_MAX': 1.08,
            'BUCKLE_MAX': 24.32,
            'PROPELLER_MIN': -24.35,
            'PROPELLER_MAX': 24.35,
            'OPENING_MIN': -112.58,
            'OPENING_MAX': 112.58,
        },
        'tW.': {
            'SHEAR_MAX': 3.54,
            'STRETCH_MIN': -3.65,
            'STRETCH_MAX': 3.65,
            'STAGGER_MAX': 1.28,
            'BUCKLE_MAX': 25.78,
            'PROPELLER_MIN': -23.12,
            'PROPELLER_MAX': 23.12,
            'OPENING_MIN': -110.02,
            'OPENING_MAX': 110.02,
        },
        'tWH': {
            'SHEAR_MAX': 4.66,
            'STRETCH_MIN': -2.17,
            'STRETCH_MAX': 2.17,
            'STAGGER_MAX': 0.92,
            'BUCKLE_MAX': 16.21,
            'PROPELLER_MIN': -18.27,
            'PROPELLER_MAX': 18.27,
            'OPENING_MIN': -91.44,
            'OPENING_MAX': 91.44,
        },
        'tWS': {
            'SHEAR_MAX': 3.04,
            'STRETCH_MIN': -5.79,
            'STRETCH_MAX': 5.79,
            'STAGGER_MAX': 1.13,
            'BUCKLE_MAX': 23.60,
            'PROPELLER_MIN': -24.66,
            'PROPELLER_MAX': 24.66,
            'OPENING_MIN': -110.15,
            'OPENING_MAX': 110.15,
        },
        'tWW': {
            'SHEAR_MAX': 1.43,
            'STRETCH_MIN': -2.45,
            'STRETCH_MAX': 2.45,
            'STAGGER_MAX': 0.92,
            'BUCKLE_MAX': 25.47,
            'PROPELLER_MIN': -19.72,
            'PROPELLER_MAX': 19.72,
            'OPENING_MIN': -159.10,
            'OPENING_MAX': 159.10,
        },
        'uncommon_edges': {
            'SHEAR_MAX': 2.82,
            'STRETCH_MIN': -3.14,
            'STRETCH_MAX': 3.14,
            'STAGGER_MAX': 1.49,
            'BUCKLE_MAX': 30.15,
            'PROPELLER_MIN': -19.99,
            'PROPELLER_MAX': 19.99,
            'OPENING_MIN': -61.26,
            'OPENING_MAX': 61.26,
        },
    }
    
    # ===== EDGE-BASED H-BOND THRESHOLDS =====
    # Thresholds: ±1 SD from mean for H-bond parameters
    HBOND_THRESHOLDS_BY_EDGE = {
        '--': {
            'DIST_MIN': 2.73,
            'DIST_MAX': 3.36,
            'ANGLE_MIN': 119.14,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.58,
        },
        '.SW': {
            'DIST_MIN': 2.91,
            'DIST_MAX': 3.11,
            'ANGLE_MIN': 89.68,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.68,
        },
        '.WW': {
            'DIST_MIN': 2.74,
            'DIST_MAX': 2.89,
            'ANGLE_MIN': 110.26,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.82,
        },
        'c.H': {
            'DIST_MIN': 2.81,
            'DIST_MAX': 3.43,
            'ANGLE_MIN': 100.04,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.56,
        },
        'c.S': {
            'DIST_MIN': 2.75,
            'DIST_MAX': 3.43,
            'ANGLE_MIN': 93.08,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.60,
        },
        'c.W': {
            'DIST_MIN': 2.75,
            'DIST_MAX': 3.37,
            'ANGLE_MIN': 93.63,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.56,
        },
        'cH.': {
            'DIST_MIN': 2.69,
            'DIST_MAX': 3.36,
            'ANGLE_MIN': 112.03,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.52,
        },
        'cHH': {
            'DIST_MIN': 2.68,
            'DIST_MAX': 3.40,
            'ANGLE_MIN': 84.70,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.57,
        },
        'cHS': {
            'DIST_MIN': 2.78,
            'DIST_MAX': 3.44,
            'ANGLE_MIN': 78.59,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.57,
        },
        'cHW': {
            'DIST_MIN': 2.72,
            'DIST_MAX': 3.30,
            'ANGLE_MIN': 91.41,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.67,
        },
        'cS.': {
            'DIST_MIN': 2.66,
            'DIST_MAX': 3.27,
            'ANGLE_MIN': 91.13,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.67,
        },
        'cSH': {
            'DIST_MIN': 2.79,
            'DIST_MAX': 3.45,
            'ANGLE_MIN': 77.91,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.57,
        },
        'cSS': {
            'DIST_MIN': 2.69,
            'DIST_MAX': 3.30,
            'ANGLE_MIN': 93.72,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.65,
        },
        'cSW': {
            'DIST_MIN': 2.70,
            'DIST_MAX': 3.34,
            'ANGLE_MIN': 87.75,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.62,
        },
        'cW.': {
            'DIST_MIN': 2.77,
            'DIST_MAX': 3.37,
            'ANGLE_MIN': 93.24,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.56,
        },
        'cWH': {
            'DIST_MIN': 2.79,
            'DIST_MAX': 3.40,
            'ANGLE_MIN': 83.12,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.59,
        },
        'cWS': {
            'DIST_MIN': 2.73,
            'DIST_MAX': 3.35,
            'ANGLE_MIN': 90.54,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.63,
        },
        'cWW': {
            'DIST_MIN': 2.72,
            'DIST_MAX': 3.12,
            'ANGLE_MIN': 105.37,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.80,
        },
        't.H': {
            'DIST_MIN': 2.70,
            'DIST_MAX': 3.37,
            'ANGLE_MIN': 101.36,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.57,
        },
        't.S': {
            'DIST_MIN': 2.69,
            'DIST_MAX': 3.38,
            'ANGLE_MIN': 87.17,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.63,
        },
        't.W': {
            'DIST_MIN': 2.74,
            'DIST_MAX': 3.37,
            'ANGLE_MIN': 88.32,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.59,
        },
        'tH.': {
            'DIST_MIN': 2.68,
            'DIST_MAX': 3.35,
            'ANGLE_MIN': 102.09,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.56,
        },
        'tHH': {
            'DIST_MIN': 2.72,
            'DIST_MAX': 3.35,
            'ANGLE_MIN': 91.91,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.59,
        },
        'tHS': {
            'DIST_MIN': 2.78,
            'DIST_MAX': 3.30,
            'ANGLE_MIN': 100.15,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.64,
        },
        'tHW': {
            'DIST_MIN': 2.72,
            'DIST_MAX': 3.24,
            'ANGLE_MIN': 101.71,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.67,
        },
        'tS.': {
            'DIST_MIN': 2.65,
            'DIST_MAX': 3.37,
            'ANGLE_MIN': 86.87,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.62,
        },
        'tSH': {
            'DIST_MIN': 2.74,
            'DIST_MAX': 3.33,
            'ANGLE_MIN': 95.46,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.62,
        },
        'tSS': {
            'DIST_MIN': 2.73,
            'DIST_MAX': 3.33,
            'ANGLE_MIN': 89.46,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.69,
        },
        'tSW': {
            'DIST_MIN': 2.72,
            'DIST_MAX': 3.36,
            'ANGLE_MIN': 90.28,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.64,
        },
        'tW.': {
            'DIST_MIN': 2.80,
            'DIST_MAX': 3.45,
            'ANGLE_MIN': 89.24,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.55,
        },
        'tWH': {
            'DIST_MIN': 2.72,
            'DIST_MAX': 3.24,
            'ANGLE_MIN': 102.99,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.66,
        },
        'tWS': {
            'DIST_MIN': 2.77,
            'DIST_MAX': 3.38,
            'ANGLE_MIN': 91.76,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.64,
        },
        'tWW': {
            'DIST_MIN': 2.75,
            'DIST_MAX': 3.29,
            'ANGLE_MIN': 93.91,
            'DIHEDRAL_FORBID_MIN': 50.00,
            'DIHEDRAL_FORBID_MAX': 140.00,
            'QUALITY_MIN': 0.66,
        },
        'uncommon_edges': {
            'DISTANCE_MIN': 2.686610,
            'DISTANCE_MAX': 3.317934,
            'ANGLE_MIN': 91.81,
            'QUALITY_MIN': 0.65,
        },
    }
    
    # ===== HYDROGEN BOND DIHEDRAL ANGLE THRESHOLDS =====
    # Global thresholds for all edge types - only flag forbidden zone
    # 
    # Classification:
    #   - CIS range: -50° to +50° (acceptable)
    #   - TRANS range: |angle| >= 140° (acceptable)
    #   - Forbidden zone: 50° to 140° and -50° to -140° (penalize)
    #
    # Rationale:
    #   - Empirical analysis shows <5% of H-bonds fall in forbidden zone
    #   - CIS and TRANS ranges represent energetically favorable H-bond geometries
    #   - Intermediate angles create strained conformations
    #   - These ranges are consistent with observed distributions in high-quality RNA structures
    #
    # Note: We only penalize dihedral angles in the forbidden zone.
    # Values in CIS or TRANS ranges are acceptable regardless of edge type.
    HBOND_DIHEDRAL_CIS_MIN = -50.0    # Minimum acceptable CIS dihedral angle
    HBOND_DIHEDRAL_CIS_MAX = 50.0     # Maximum acceptable CIS dihedral angle
    HBOND_DIHEDRAL_TRANS_MIN = 140.0  # Minimum absolute value for TRANS dihedral angle
    
    # ===== OTHER HYDROGEN BOND THRESHOLDS =====
    HBOND_DISTANCE_MIN = 2.3
    HBOND_DISTANCE_MAX = 3.7
    HBOND_ANGLE_MIN = 80.0
    HBOND_QUALITY_MIN = 0.70
    
    
    
    
    
    
    
    
    
    
    
    # ===== DETAILED ISSUES THRESHOLD =====
    BASELINE = 75 
    # ===== DATA DIRECTORIES =====
    BASEPAIR_DIR = './data/basepairs'
    HBOND_DIR = './data/hbonds'
    
    # ===== OUTPUT CONFIGURATION =====
    MAX_ISSUES_DISPLAYED = 20
    MAX_SUMMARY_ISSUES = 20
    
    # ===== HOTSPOT DETECTION - MAXIMUM SENSITIVITY =====
    HOTSPOT_WINDOW_SIZE = 30
    HOTSPOT_SLIDE_STEP = 3
    HOTSPOT_MERGE_DISTANCE = 10

    # Severity levels - VERY LENIENT
    SEVERITY_CRITICAL = 40       # Easier to reach
    SEVERITY_SEVERE = 70         # Easier to reach
    SEVERITY_MODERATE = 80       # Easier to reach
    SEVERITY_MINOR = 90          # Easier to reach

    # Issue density thresholds - VERY SENSITIVE
    ISSUE_DENSITY_MEDIUM = 0.20  # Now 20% triggers reporting

    # Minimum residues - VERY SMALL
    MIN_HOTSPOT_RESIDUES = 4     # Only 2 residues needed!
    
    DAMAGE_THRESHOLD_MINOR = 90
    DAMAGE_THRESHOLD_MODERATE = 80
    DAMAGE_THRESHOLD_SEVERE = 70
    DAMAGE_THRESHOLD_CRITICAL = 40
    
    # MAXIMUM SENSITIVITY
    DAMAGE_THRESHOLD_FOR_HOTSPOTS = 80  # Was 75 - now catches even more
    
    # Larger connectivity
    SEQUENTIAL_NEIGHBOR_RANGE = 8  # Was 7 - even bigger regions
    
    MAX_SEVERITY_SPREAD = 3
    
    GAP_TOLERANCE_CRITICAL = 5
    GAP_TOLERANCE_SEVERE = 3
    GAP_TOLERANCE_MODERATE = 2
    
    # HOTSPOT VALIDATION - MAXIMUM SENSITIVITY
    MIN_HOTSPOT_PROBLEMATIC_COUNT = 1   # Only 1 bad pair needed!
    MIN_HOTSPOT_BP_FRACTION = 0.15      # Only 15% is enough
    MIN_HOTSPOT_ISSUE_DENSITY = 0.30    # Only 30% is enough
    MAX_HOTSPOT_SCORE = 70.0            # Higher = catches more
