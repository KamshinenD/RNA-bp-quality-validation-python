"""Configuration file with all tunable thresholds and scoring parameters."""

class Config:
    """RNA quality scoring configuration."""
    
    BASE_SCORE = 100.0
    
#     PENALTY_WEIGHTS = {
#     'non_coplanar_pairs': 15.0,
#     'low_hbond_score_pairs': 10.0,
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
    'low_hbond_score_pairs': 5.0,            # DSSR score < 2.2 (reduce - too strict!)
    'weak_hbond_quality': 5.0,          # Individual H-bond quality < 0.70
    
    # MINOR ISSUES (common even in good structures)
    'bad_hbond_angles': 5.0,            # One angle < 80° is very common
    'incorrect_hbond_count': 5.0,       # With tolerant ranges, less critical
}
# Total: 100.0 points
    
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

    # Base-pair + edge thresholds (mean ± SD, 2dp) generated from edge_type_analysis.json
    GEOMETRY_THRESHOLDS_BY_BASE_PAIR_BY_EDGE = {
        "A-C": {
            "cSS": {
                "SHEAR_MAX": 2.8,
                "STRETCH_MIN": 3.06,
                "STRETCH_MAX": 9.04,
                "STAGGER_MAX": 1.03,
                "BUCKLE_MAX": 20.67,
                "PROPELLER_MIN": -26.62,
                "PROPELLER_MAX": -2.24,
                "OPENING_MIN": 67.37,
                "OPENING_MAX": 206.25
            },
            "tHW": {
                "SHEAR_MAX": 5.24,
                "STRETCH_MIN": -2.16,
                "STRETCH_MAX": 0.78,
                "STAGGER_MAX": 1.28,
                "BUCKLE_MAX": 23.09,
                "PROPELLER_MIN": -22.76,
                "PROPELLER_MAX": 20.65,
                "OPENING_MIN": -107.79,
                "OPENING_MAX": -52.47
            },
            "cWS": {
                "SHEAR_MAX": 7.17,
                "STRETCH_MIN": -2.33,
                "STRETCH_MAX": 3.48,
                "STAGGER_MAX": 1.47,
                "BUCKLE_MAX": 36.87,
                "PROPELLER_MIN": -31.84,
                "PROPELLER_MAX": 3.69,
                "OPENING_MIN": 38.57,
                "OPENING_MAX": 111.64
            },
            "tSW": {
                "SHEAR_MAX": 3.4,
                "STRETCH_MIN": -8.46,
                "STRETCH_MAX": -4.7,
                "STAGGER_MAX": 1.1,
                "BUCKLE_MAX": 32.78,
                "PROPELLER_MIN": -19.11,
                "PROPELLER_MAX": 15.0,
                "OPENING_MIN": -165.65,
                "OPENING_MAX": -56.92
            },
            "cSH": {
                "SHEAR_MAX": 7.0,
                "STRETCH_MIN": -1.42,
                "STRETCH_MAX": 0.82,
                "STAGGER_MAX": 1.64,
                "BUCKLE_MAX": 33.75,
                "PROPELLER_MIN": -7.31,
                "PROPELLER_MAX": 28.73,
                "OPENING_MIN": -36.52,
                "OPENING_MAX": 4.93
            },
            "tWH": {
                "SHEAR_MAX": 4.58,
                "STRETCH_MIN": -1.67,
                "STRETCH_MAX": 1.1,
                "STAGGER_MAX": 1.18,
                "BUCKLE_MAX": 18.51,
                "PROPELLER_MIN": -21.65,
                "PROPELLER_MAX": 15.59,
                "OPENING_MIN": -114.42,
                "OPENING_MAX": -63.95
            },
            "cWH": {
                "SHEAR_MAX": 4.49,
                "STRETCH_MIN": -0.42,
                "STRETCH_MAX": 2.67,
                "STAGGER_MAX": 1.95,
                "BUCKLE_MAX": 55.75,
                "PROPELLER_MIN": -31.97,
                "PROPELLER_MAX": 10.19,
                "OPENING_MIN": -93.38,
                "OPENING_MAX": -43.74
            },
            "tH.": {
                "SHEAR_MAX": 7.67,
                "STRETCH_MIN": -5.74,
                "STRETCH_MAX": -0.97,
                "STAGGER_MAX": 0.9,
                "BUCKLE_MAX": 18.33,
                "PROPELLER_MIN": -14.81,
                "PROPELLER_MAX": 12.07,
                "OPENING_MIN": -71.31,
                "OPENING_MAX": 22.9
            },
            "cWW": {
                "SHEAR_MAX": 3.42,
                "STRETCH_MIN": -0.88,
                "STRETCH_MAX": 0.36,
                "STAGGER_MAX": 0.95,
                "BUCKLE_MAX": 20.02,
                "PROPELLER_MIN": -19.99,
                "PROPELLER_MAX": 9.23,
                "OPENING_MIN": -11.95,
                "OPENING_MAX": 20.47
            },
            "tWW": {
                "SHEAR_MAX": 1.19,
                "STRETCH_MIN": -2.22,
                "STRETCH_MAX": 1.14,
                "STAGGER_MAX": 0.98,
                "BUCKLE_MAX": 26.24,
                "PROPELLER_MIN": -17.03,
                "PROPELLER_MAX": 19.11,
                "OPENING_MIN": -145.16,
                "OPENING_MAX": 163.73
            },
            "cSW": {
                "SHEAR_MAX": 7.17,
                "STRETCH_MIN": -2.44,
                "STRETCH_MAX": 3.33,
                "STAGGER_MAX": 1.12,
                "BUCKLE_MAX": 25.63,
                "PROPELLER_MIN": -25.96,
                "PROPELLER_MAX": 1.65,
                "OPENING_MIN": 32.42,
                "OPENING_MAX": 108.25
            },
            "tSH": {
                "SHEAR_MAX": 7.45,
                "STRETCH_MIN": -5.81,
                "STRETCH_MAX": -2.89,
                "STAGGER_MAX": 1.25,
                "BUCKLE_MAX": 25.89,
                "PROPELLER_MIN": -23.15,
                "PROPELLER_MAX": 11.67,
                "OPENING_MIN": -43.19,
                "OPENING_MAX": 3.6
            },
            "tWS": {
                "SHEAR_MAX": 3.88,
                "STRETCH_MIN": 1.53,
                "STRETCH_MAX": 9.16,
                "STAGGER_MAX": 1.59,
                "BUCKLE_MAX": 21.71,
                "PROPELLER_MIN": -36.29,
                "PROPELLER_MAX": 23.13,
                "OPENING_MIN": 70.94,
                "OPENING_MAX": 147.74
            },
            "t.H": {
                "SHEAR_MAX": 8.49,
                "STRETCH_MIN": -5.12,
                "STRETCH_MAX": -0.89,
                "STAGGER_MAX": 0.9,
                "BUCKLE_MAX": 25.15,
                "PROPELLER_MIN": -21.81,
                "PROPELLER_MAX": 7.95,
                "OPENING_MIN": -53.29,
                "OPENING_MAX": 24.56
            },
            "tHS": {
                "SHEAR_MAX": 7.34,
                "STRETCH_MIN": -6.14,
                "STRETCH_MAX": -3.04,
                "STAGGER_MAX": 0.86,
                "BUCKLE_MAX": 20.75,
                "PROPELLER_MIN": -21.38,
                "PROPELLER_MAX": 15.77,
                "OPENING_MIN": -28.04,
                "OPENING_MAX": 11.55
            },
            "cHS": {
                "SHEAR_MAX": 7.26,
                "STRETCH_MIN": -0.83,
                "STRETCH_MAX": 2.08,
                "STAGGER_MAX": 2.64,
                "BUCKLE_MAX": 54.96,
                "PROPELLER_MIN": -35.59,
                "PROPELLER_MAX": 1.81,
                "OPENING_MIN": -20.41,
                "OPENING_MAX": 37.59
            },
            "c.W": {
                "SHEAR_MAX": 5.87,
                "STRETCH_MIN": -1.69,
                "STRETCH_MAX": 0.81,
                "STAGGER_MAX": 1.24,
                "BUCKLE_MAX": 35.41,
                "PROPELLER_MIN": -16.75,
                "PROPELLER_MAX": 11.67,
                "OPENING_MIN": -40.82,
                "OPENING_MAX": 30.02
            },
            "tHH": {
                "SHEAR_MAX": 7.54,
                "STRETCH_MIN": -4.28,
                "STRETCH_MAX": 6.73,
                "STAGGER_MAX": 1.27,
                "BUCKLE_MAX": 22.29,
                "PROPELLER_MIN": -45.93,
                "PROPELLER_MAX": 10.25,
                "OPENING_MIN": -191.58,
                "OPENING_MAX": 129.82
            },
            "cW.": {
                "SHEAR_MAX": 7.03,
                "STRETCH_MIN": -1.75,
                "STRETCH_MAX": 0.77,
                "STAGGER_MAX": 1.19,
                "BUCKLE_MAX": 24.23,
                "PROPELLER_MIN": -14.89,
                "PROPELLER_MAX": 12.86,
                "OPENING_MIN": -42.42,
                "OPENING_MAX": 29.89
            },
            "t.W": {
                "SHEAR_MAX": 5.97,
                "STRETCH_MIN": -3.33,
                "STRETCH_MAX": 0.89,
                "STAGGER_MAX": 2.13,
                "BUCKLE_MAX": 37.91,
                "PROPELLER_MIN": -12.38,
                "PROPELLER_MAX": 20.9,
                "OPENING_MIN": -106.4,
                "OPENING_MAX": 0.67
            },
            "cHW": {
                "SHEAR_MAX": 3.23,
                "STRETCH_MIN": -5.32,
                "STRETCH_MAX": -1.04,
                "STAGGER_MAX": 1.0,
                "BUCKLE_MAX": 29.12,
                "PROPELLER_MIN": -22.35,
                "PROPELLER_MAX": 19.48,
                "OPENING_MIN": 5.97,
                "OPENING_MAX": 96.58
            }
        },
        "A-G": {
            "tW.": {
                "SHEAR_MAX": 3.46,
                "STRETCH_MIN": 4.46,
                "STRETCH_MAX": 6.56,
                "STAGGER_MAX": 1.4,
                "BUCKLE_MAX": 43.05,
                "PROPELLER_MIN": -31.15,
                "PROPELLER_MAX": 8.04,
                "OPENING_MIN": 69.86,
                "OPENING_MAX": 183.89
            },
            "tSS": {
                "SHEAR_MAX": 3.06,
                "STRETCH_MIN": -8.57,
                "STRETCH_MAX": 6.21,
                "STAGGER_MAX": 1.55,
                "BUCKLE_MAX": 31.1,
                "PROPELLER_MIN": -24.04,
                "PROPELLER_MAX": 21.1,
                "OPENING_MIN": -175.29,
                "OPENING_MAX": 128.13
            },
            "tSH": {
                "SHEAR_MAX": 7.42,
                "STRETCH_MIN": -5.28,
                "STRETCH_MAX": -3.83,
                "STAGGER_MAX": 1.08,
                "BUCKLE_MAX": 21.77,
                "PROPELLER_MIN": -15.53,
                "PROPELLER_MAX": 13.08,
                "OPENING_MIN": -27.67,
                "OPENING_MAX": 4.44
            },
            "tS.": {
                "SHEAR_MAX": 2.85,
                "STRETCH_MIN": -7.55,
                "STRETCH_MAX": -6.52,
                "STAGGER_MAX": 2.28,
                "BUCKLE_MAX": 59.69,
                "PROPELLER_MIN": -12.35,
                "PROPELLER_MAX": 21.94,
                "OPENING_MIN": -158.53,
                "OPENING_MAX": -128.42
            },
            "tHS": {
                "SHEAR_MAX": 7.09,
                "STRETCH_MIN": -4.93,
                "STRETCH_MAX": -3.62,
                "STAGGER_MAX": 0.78,
                "BUCKLE_MAX": 16.42,
                "PROPELLER_MIN": -17.91,
                "PROPELLER_MAX": 10.03,
                "OPENING_MIN": -14.96,
                "OPENING_MAX": 5.61
            },
            "t.H": {
                "SHEAR_MAX": 6.46,
                "STRETCH_MIN": -1.51,
                "STRETCH_MAX": 5.14,
                "STAGGER_MAX": 1.19,
                "BUCKLE_MAX": 26.41,
                "PROPELLER_MIN": -15.97,
                "PROPELLER_MAX": 20.14,
                "OPENING_MIN": -130.33,
                "OPENING_MAX": -47.15
            },
            "cWW": {
                "SHEAR_MAX": 1.05,
                "STRETCH_MIN": 1.13,
                "STRETCH_MAX": 2.0,
                "STAGGER_MAX": 1.06,
                "BUCKLE_MAX": 17.68,
                "PROPELLER_MIN": -20.51,
                "PROPELLER_MAX": 12.22,
                "OPENING_MIN": -32.55,
                "OPENING_MAX": 1.93
            },
            "tH.": {
                "SHEAR_MAX": 4.72,
                "STRETCH_MIN": -0.52,
                "STRETCH_MAX": 5.55,
                "STAGGER_MAX": 1.48,
                "BUCKLE_MAX": 18.77,
                "PROPELLER_MIN": -17.56,
                "PROPELLER_MAX": 16.17,
                "OPENING_MIN": -164.96,
                "OPENING_MAX": -4.94
            },
            "cSW": {
                "SHEAR_MAX": 4.92,
                "STRETCH_MIN": 1.24,
                "STRETCH_MAX": 6.22,
                "STAGGER_MAX": 1.49,
                "BUCKLE_MAX": 25.63,
                "PROPELLER_MIN": -42.61,
                "PROPELLER_MAX": 14.07,
                "OPENING_MIN": 67.64,
                "OPENING_MAX": 139.77
            },
            "cSS": {
                "SHEAR_MAX": 2.31,
                "STRETCH_MIN": 5.86,
                "STRETCH_MAX": 8.57,
                "STAGGER_MAX": 1.52,
                "BUCKLE_MAX": 28.93,
                "PROPELLER_MIN": -28.35,
                "PROPELLER_MAX": 9.9,
                "OPENING_MIN": 113.96,
                "OPENING_MAX": 172.61
            },
            "cSH": {
                "SHEAR_MAX": 7.54,
                "STRETCH_MIN": -2.26,
                "STRETCH_MAX": 2.51,
                "STAGGER_MAX": 1.6,
                "BUCKLE_MAX": 30.29,
                "PROPELLER_MIN": -14.82,
                "PROPELLER_MAX": 31.46,
                "OPENING_MIN": -16.36,
                "OPENING_MAX": 26.83
            },
            "tSW": {
                "SHEAR_MAX": 3.7,
                "STRETCH_MIN": -7.26,
                "STRETCH_MAX": -3.83,
                "STAGGER_MAX": 1.1,
                "BUCKLE_MAX": 31.1,
                "PROPELLER_MIN": -7.86,
                "PROPELLER_MAX": 37.43,
                "OPENING_MIN": -141.51,
                "OPENING_MAX": -70.58
            },
            "tWS": {
                "SHEAR_MAX": 4.0,
                "STRETCH_MIN": 3.47,
                "STRETCH_MAX": 6.96,
                "STAGGER_MAX": 1.08,
                "BUCKLE_MAX": 25.74,
                "PROPELLER_MIN": -32.9,
                "PROPELLER_MAX": 3.33,
                "OPENING_MIN": 64.16,
                "OPENING_MAX": 133.83
            },
            "cWS": {
                "SHEAR_MAX": 5.32,
                "STRETCH_MIN": 0.86,
                "STRETCH_MAX": 6.02,
                "STAGGER_MAX": 1.72,
                "BUCKLE_MAX": 46.6,
                "PROPELLER_MIN": -32.93,
                "PROPELLER_MAX": 14.3,
                "OPENING_MIN": 68.67,
                "OPENING_MAX": 141.57
            },
            "cWH": {
                "SHEAR_MAX": 2.08,
                "STRETCH_MIN": 3.21,
                "STRETCH_MAX": 6.01,
                "STAGGER_MAX": 1.42,
                "BUCKLE_MAX": 35.52,
                "PROPELLER_MIN": -19.87,
                "PROPELLER_MAX": 19.26,
                "OPENING_MIN": -125.9,
                "OPENING_MAX": -73.65
            },
            "cHH": {
                "SHEAR_MAX": 3.47,
                "STRETCH_MIN": -6.81,
                "STRETCH_MAX": 3.77,
                "STAGGER_MAX": 1.7,
                "BUCKLE_MAX": 34.09,
                "PROPELLER_MIN": -11.63,
                "PROPELLER_MAX": 29.63,
                "OPENING_MIN": -127.18,
                "OPENING_MAX": 174.21
            },
            "tWW": {
                "SHEAR_MAX": 2.3,
                "STRETCH_MIN": -3.34,
                "STRETCH_MAX": 4.89,
                "STAGGER_MAX": 1.42,
                "BUCKLE_MAX": 47.14,
                "PROPELLER_MIN": -28.97,
                "PROPELLER_MAX": 16.02,
                "OPENING_MIN": -135.57,
                "OPENING_MAX": 176.53
            },
            "cHW": {
                "SHEAR_MAX": 2.07,
                "STRETCH_MIN": -6.25,
                "STRETCH_MAX": -2.55,
                "STAGGER_MAX": 0.85,
                "BUCKLE_MAX": 18.03,
                "PROPELLER_MIN": -30.88,
                "PROPELLER_MAX": 15.02,
                "OPENING_MIN": 33.1,
                "OPENING_MAX": 124.91
            },
            "t.S": {
                "SHEAR_MAX": 2.76,
                "STRETCH_MIN": 6.61,
                "STRETCH_MAX": 7.45,
                "STAGGER_MAX": 2.59,
                "BUCKLE_MAX": 65.24,
                "PROPELLER_MIN": -13.38,
                "PROPELLER_MAX": 6.71,
                "OPENING_MIN": 142.2,
                "OPENING_MAX": 155.72
            },
            "c.W": {
                "SHEAR_MAX": 4.02,
                "STRETCH_MIN": 0.45,
                "STRETCH_MAX": 3.54,
                "STAGGER_MAX": 1.79,
                "BUCKLE_MAX": 27.7,
                "PROPELLER_MIN": -29.82,
                "PROPELLER_MAX": 10.6,
                "OPENING_MIN": -77.66,
                "OPENING_MAX": 29.54
            },
            "tWH": {
                "SHEAR_MAX": 8.91,
                "STRETCH_MIN": -4.78,
                "STRETCH_MAX": -0.35,
                "STAGGER_MAX": 1.22,
                "BUCKLE_MAX": 24.78,
                "PROPELLER_MIN": -21.65,
                "PROPELLER_MAX": 25.97,
                "OPENING_MIN": -102.51,
                "OPENING_MAX": -49.17
            },
            "--": {
                "SHEAR_MAX": 3.32,
                "STRETCH_MIN": 1.91,
                "STRETCH_MAX": 5.11,
                "STAGGER_MAX": 1.59,
                "BUCKLE_MAX": 16.92,
                "PROPELLER_MIN": -15.94,
                "PROPELLER_MAX": 11.14,
                "OPENING_MIN": -125.15,
                "OPENING_MAX": -41.03
            },
            "tHW": {
                "SHEAR_MAX": 8.09,
                "STRETCH_MIN": -4.0,
                "STRETCH_MAX": 0.84,
                "STAGGER_MAX": 1.46,
                "BUCKLE_MAX": 20.27,
                "PROPELLER_MIN": -20.23,
                "PROPELLER_MAX": 19.49,
                "OPENING_MIN": -110.3,
                "OPENING_MAX": -5.64
            },
            "t.W": {
                "SHEAR_MAX": 3.59,
                "STRETCH_MIN": -6.33,
                "STRETCH_MAX": -4.11,
                "STAGGER_MAX": 1.65,
                "BUCKLE_MAX": 36.36,
                "PROPELLER_MIN": -2.55,
                "PROPELLER_MAX": 26.26,
                "OPENING_MIN": -155.67,
                "OPENING_MAX": -102.58
            },
            "cH.": {
                "SHEAR_MAX": 5.24,
                "STRETCH_MIN": -5.08,
                "STRETCH_MAX": 2.9,
                "STAGGER_MAX": 1.29,
                "BUCKLE_MAX": 29.2,
                "PROPELLER_MIN": -20.63,
                "PROPELLER_MAX": 7.77,
                "OPENING_MIN": -76.66,
                "OPENING_MAX": 154.0
            },
            "cW.": {
                "SHEAR_MAX": 5.94,
                "STRETCH_MIN": 0.97,
                "STRETCH_MAX": 3.25,
                "STAGGER_MAX": 2.09,
                "BUCKLE_MAX": 37.34,
                "PROPELLER_MIN": -11.38,
                "PROPELLER_MAX": 30.76,
                "OPENING_MIN": -31.86,
                "OPENING_MAX": 87.61
            },
            "cHS": {
                "SHEAR_MAX": 7.27,
                "STRETCH_MIN": -2.67,
                "STRETCH_MAX": 2.23,
                "STAGGER_MAX": 2.43,
                "BUCKLE_MAX": 50.07,
                "PROPELLER_MIN": -37.05,
                "PROPELLER_MAX": 7.69,
                "OPENING_MIN": -28.73,
                "OPENING_MAX": 34.84
            },
            "tHH": {
                "SHEAR_MAX": 7.41,
                "STRETCH_MIN": -1.63,
                "STRETCH_MAX": 7.05,
                "STAGGER_MAX": 1.35,
                "BUCKLE_MAX": 24.4,
                "PROPELLER_MIN": -18.51,
                "PROPELLER_MAX": 24.18,
                "OPENING_MIN": -185.24,
                "OPENING_MAX": 44.2
            },
            "c.H": {
                "SHEAR_MAX": 7.86,
                "STRETCH_MIN": -0.03,
                "STRETCH_MAX": 2.48,
                "STAGGER_MAX": 1.16,
                "BUCKLE_MAX": 26.23,
                "PROPELLER_MIN": -5.72,
                "PROPELLER_MAX": 24.96,
                "OPENING_MIN": -83.42,
                "OPENING_MAX": -32.63
            }
        },
        "C-G": {
            "cWW": {
                "SHEAR_MAX": 0.53,
                "STRETCH_MIN": -0.36,
                "STRETCH_MAX": 0.11,
                "STAGGER_MAX": 0.62,
                "BUCKLE_MAX": 10.41,
                "PROPELLER_MIN": -16.17,
                "PROPELLER_MAX": 2.73,
                "OPENING_MIN": -5.29,
                "OPENING_MAX": 6.94
            },
            "--": {
                "SHEAR_MAX": 2.06,
                "STRETCH_MIN": 1.71,
                "STRETCH_MAX": 4.83,
                "STAGGER_MAX": 0.89,
                "BUCKLE_MAX": 16.0,
                "PROPELLER_MIN": -17.55,
                "PROPELLER_MAX": 20.46,
                "OPENING_MIN": -140.97,
                "OPENING_MAX": -43.96
            },
            "tHH": {
                "SHEAR_MAX": 5.65,
                "STRETCH_MIN": -7.38,
                "STRETCH_MAX": 3.68,
                "STAGGER_MAX": 1.25,
                "BUCKLE_MAX": 25.16,
                "PROPELLER_MIN": -6.09,
                "PROPELLER_MAX": 20.68,
                "OPENING_MIN": -111.05,
                "OPENING_MAX": 202.56
            },
            "cSW": {
                "SHEAR_MAX": 6.01,
                "STRETCH_MIN": -1.96,
                "STRETCH_MAX": 3.65,
                "STAGGER_MAX": 1.63,
                "BUCKLE_MAX": 46.54,
                "PROPELLER_MIN": -23.8,
                "PROPELLER_MAX": 17.03,
                "OPENING_MIN": 30.2,
                "OPENING_MAX": 107.93
            },
            "cSH": {
                "SHEAR_MAX": 7.58,
                "STRETCH_MIN": -2.78,
                "STRETCH_MAX": 0.99,
                "STAGGER_MAX": 1.43,
                "BUCKLE_MAX": 26.35,
                "PROPELLER_MIN": -12.82,
                "PROPELLER_MAX": 20.55,
                "OPENING_MIN": -24.93,
                "OPENING_MAX": 21.39
            },
            "cSS": {
                "SHEAR_MAX": 2.94,
                "STRETCH_MIN": -1.41,
                "STRETCH_MAX": 10.56,
                "STAGGER_MAX": 1.57,
                "BUCKLE_MAX": 27.37,
                "PROPELLER_MIN": -30.71,
                "PROPELLER_MAX": 10.54,
                "OPENING_MIN": -30.21,
                "OPENING_MAX": 225.24
            },
            "cWH": {
                "SHEAR_MAX": 2.95,
                "STRETCH_MIN": 1.5,
                "STRETCH_MAX": 5.02,
                "STAGGER_MAX": 1.14,
                "BUCKLE_MAX": 37.38,
                "PROPELLER_MIN": -23.31,
                "PROPELLER_MAX": 19.04,
                "OPENING_MIN": -109.38,
                "OPENING_MAX": -48.07
            },
            "cW.": {
                "SHEAR_MAX": 3.21,
                "STRETCH_MIN": 1.36,
                "STRETCH_MAX": 4.43,
                "STAGGER_MAX": 0.94,
                "BUCKLE_MAX": 27.52,
                "PROPELLER_MIN": -22.77,
                "PROPELLER_MAX": 9.93,
                "OPENING_MIN": -134.16,
                "OPENING_MAX": -28.71
            },
            "t.W": {
                "SHEAR_MAX": 3.5,
                "STRETCH_MIN": -2.71,
                "STRETCH_MAX": 3.74,
                "STAGGER_MAX": 1.02,
                "BUCKLE_MAX": 20.19,
                "PROPELLER_MIN": -16.89,
                "PROPELLER_MAX": 19.6,
                "OPENING_MIN": -141.08,
                "OPENING_MAX": 18.13
            },
            "tWS": {
                "SHEAR_MAX": 2.67,
                "STRETCH_MIN": 1.76,
                "STRETCH_MAX": 9.56,
                "STAGGER_MAX": 1.17,
                "BUCKLE_MAX": 34.18,
                "PROPELLER_MIN": -21.54,
                "PROPELLER_MAX": 24.47,
                "OPENING_MIN": 41.04,
                "OPENING_MAX": 185.25
            },
            "tSW": {
                "SHEAR_MAX": 2.54,
                "STRETCH_MIN": -7.08,
                "STRETCH_MAX": -2.85,
                "STAGGER_MAX": 1.85,
                "BUCKLE_MAX": 43.25,
                "PROPELLER_MIN": -9.06,
                "PROPELLER_MAX": 35.85,
                "OPENING_MIN": -136.22,
                "OPENING_MAX": -64.29
            },
            "tWW": {
                "SHEAR_MAX": 0.56,
                "STRETCH_MIN": -3.73,
                "STRETCH_MAX": 3.47,
                "STAGGER_MAX": 1.04,
                "BUCKLE_MAX": 25.17,
                "PROPELLER_MIN": -22.74,
                "PROPELLER_MAX": 20.99,
                "OPENING_MIN": -164.08,
                "OPENING_MAX": 143.57
            },
            "tWH": {
                "SHEAR_MAX": 4.95,
                "STRETCH_MIN": -0.94,
                "STRETCH_MAX": 2.5,
                "STAGGER_MAX": 1.13,
                "BUCKLE_MAX": 20.76,
                "PROPELLER_MIN": -13.63,
                "PROPELLER_MAX": 26.47,
                "OPENING_MIN": -139.37,
                "OPENING_MAX": -76.97
            },
            "cHH": {
                "SHEAR_MAX": 2.9,
                "STRETCH_MIN": -8.73,
                "STRETCH_MAX": 0.37,
                "STAGGER_MAX": 1.24,
                "BUCKLE_MAX": 31.91,
                "PROPELLER_MIN": -27.04,
                "PROPELLER_MAX": 16.7,
                "OPENING_MIN": -30.21,
                "OPENING_MAX": 218.35
            },
            "cWS": {
                "SHEAR_MAX": 6.52,
                "STRETCH_MIN": -2.12,
                "STRETCH_MAX": 2.98,
                "STAGGER_MAX": 1.82,
                "BUCKLE_MAX": 46.08,
                "PROPELLER_MIN": -27.56,
                "PROPELLER_MAX": 15.55,
                "OPENING_MIN": 42.28,
                "OPENING_MAX": 98.61
            },
            "tW.": {
                "SHEAR_MAX": 4.54,
                "STRETCH_MIN": 0.92,
                "STRETCH_MAX": 4.83,
                "STAGGER_MAX": 1.77,
                "BUCKLE_MAX": 17.17,
                "PROPELLER_MIN": -19.86,
                "PROPELLER_MAX": 13.16,
                "OPENING_MIN": -155.73,
                "OPENING_MAX": 12.25
            },
            "tHW": {
                "SHEAR_MAX": 4.98,
                "STRETCH_MIN": -1.94,
                "STRETCH_MAX": 2.13,
                "STAGGER_MAX": 0.99,
                "BUCKLE_MAX": 25.9,
                "PROPELLER_MIN": -23.85,
                "PROPELLER_MAX": 15.82,
                "OPENING_MIN": -130.03,
                "OPENING_MAX": -19.15
            },
            "c.W": {
                "SHEAR_MAX": 3.47,
                "STRETCH_MIN": -1.6,
                "STRETCH_MAX": 3.27,
                "STAGGER_MAX": 1.28,
                "BUCKLE_MAX": 28.89,
                "PROPELLER_MIN": -24.76,
                "PROPELLER_MAX": 13.44,
                "OPENING_MIN": -81.26,
                "OPENING_MAX": 79.2
            },
            "cHW": {
                "SHEAR_MAX": 1.91,
                "STRETCH_MIN": -6.01,
                "STRETCH_MAX": -2.22,
                "STAGGER_MAX": 1.1,
                "BUCKLE_MAX": 28.62,
                "PROPELLER_MIN": -15.04,
                "PROPELLER_MAX": 19.07,
                "OPENING_MIN": 44.85,
                "OPENING_MAX": 132.81
            },
            "tSH": {
                "SHEAR_MAX": 7.44,
                "STRETCH_MIN": -6.12,
                "STRETCH_MAX": -3.28,
                "STAGGER_MAX": 1.98,
                "BUCKLE_MAX": 33.45,
                "PROPELLER_MIN": -22.89,
                "PROPELLER_MAX": 19.15,
                "OPENING_MIN": -50.14,
                "OPENING_MAX": 3.55
            }
        },
        "A-U": {
            "cWW": {
                "SHEAR_MAX": 0.73,
                "STRETCH_MIN": -0.33,
                "STRETCH_MAX": 0.17,
                "STAGGER_MAX": 0.55,
                "BUCKLE_MAX": 10.28,
                "PROPELLER_MIN": -15.95,
                "PROPELLER_MAX": 3.51,
                "OPENING_MIN": -6.58,
                "OPENING_MAX": 8.45
            },
            "cWH": {
                "SHEAR_MAX": 1.57,
                "STRETCH_MIN": 1.96,
                "STRETCH_MAX": 4.41,
                "STAGGER_MAX": 1.2,
                "BUCKLE_MAX": 25.04,
                "PROPELLER_MIN": -20.95,
                "PROPELLER_MAX": 9.85,
                "OPENING_MIN": -83.4,
                "OPENING_MAX": -47.33
            },
            "tWW": {
                "SHEAR_MAX": 0.89,
                "STRETCH_MIN": -1.51,
                "STRETCH_MAX": 1.72,
                "STAGGER_MAX": 0.76,
                "BUCKLE_MAX": 17.97,
                "PROPELLER_MIN": -14.23,
                "PROPELLER_MAX": 17.5,
                "OPENING_MIN": -134.45,
                "OPENING_MAX": 184.59
            },
            "cHS": {
                "SHEAR_MAX": 6.87,
                "STRETCH_MIN": -0.63,
                "STRETCH_MAX": 1.59,
                "STAGGER_MAX": 1.6,
                "BUCKLE_MAX": 32.48,
                "PROPELLER_MIN": -29.96,
                "PROPELLER_MAX": -5.77,
                "OPENING_MIN": -4.84,
                "OPENING_MAX": 25.97
            },
            "cSH": {
                "SHEAR_MAX": 6.98,
                "STRETCH_MIN": -1.67,
                "STRETCH_MAX": 0.87,
                "STAGGER_MAX": 1.69,
                "BUCKLE_MAX": 37.13,
                "PROPELLER_MIN": -20.6,
                "PROPELLER_MAX": 30.4,
                "OPENING_MIN": -28.68,
                "OPENING_MAX": 13.01
            },
            "cW.": {
                "SHEAR_MAX": 5.91,
                "STRETCH_MIN": -1.3,
                "STRETCH_MAX": 2.31,
                "STAGGER_MAX": 1.56,
                "BUCKLE_MAX": 24.79,
                "PROPELLER_MIN": -31.49,
                "PROPELLER_MAX": 7.75,
                "OPENING_MIN": -85.26,
                "OPENING_MAX": 32.52
            },
            "tHW": {
                "SHEAR_MAX": 4.9,
                "STRETCH_MIN": -2.72,
                "STRETCH_MAX": -1.17,
                "STAGGER_MAX": 1.09,
                "BUCKLE_MAX": 12.36,
                "PROPELLER_MIN": -19.47,
                "PROPELLER_MAX": 6.42,
                "OPENING_MIN": -111.59,
                "OPENING_MAX": -68.13
            },
            "tWH": {
                "SHEAR_MAX": 5.03,
                "STRETCH_MIN": -2.94,
                "STRETCH_MAX": -1.11,
                "STAGGER_MAX": 0.96,
                "BUCKLE_MAX": 13.25,
                "PROPELLER_MIN": -19.53,
                "PROPELLER_MAX": 10.62,
                "OPENING_MIN": -110.92,
                "OPENING_MAX": -67.02
            },
            "cSW": {
                "SHEAR_MAX": 6.74,
                "STRETCH_MIN": -1.68,
                "STRETCH_MAX": 2.62,
                "STAGGER_MAX": 1.64,
                "BUCKLE_MAX": 30.93,
                "PROPELLER_MIN": -22.62,
                "PROPELLER_MAX": 28.22,
                "OPENING_MIN": 39.63,
                "OPENING_MAX": 105.99
            },
            "cWS": {
                "SHEAR_MAX": 6.34,
                "STRETCH_MIN": -1.3,
                "STRETCH_MAX": 3.38,
                "STAGGER_MAX": 1.13,
                "BUCKLE_MAX": 27.07,
                "PROPELLER_MIN": -26.16,
                "PROPELLER_MAX": 23.52,
                "OPENING_MIN": 52.33,
                "OPENING_MAX": 114.91
            },
            "tSH": {
                "SHEAR_MAX": 6.87,
                "STRETCH_MIN": -9.03,
                "STRETCH_MAX": -5.02,
                "STAGGER_MAX": 1.46,
                "BUCKLE_MAX": 21.09,
                "PROPELLER_MIN": -20.47,
                "PROPELLER_MAX": 20.31,
                "OPENING_MIN": -77.67,
                "OPENING_MAX": -12.8
            },
            "cSS": {
                "SHEAR_MAX": 3.15,
                "STRETCH_MIN": 3.22,
                "STRETCH_MAX": 8.98,
                "STAGGER_MAX": 0.97,
                "BUCKLE_MAX": 21.89,
                "PROPELLER_MIN": -25.45,
                "PROPELLER_MAX": 13.49,
                "OPENING_MIN": 85.74,
                "OPENING_MAX": 198.2
            },
            "c.W": {
                "SHEAR_MAX": 6.94,
                "STRETCH_MIN": -2.11,
                "STRETCH_MAX": 1.41,
                "STAGGER_MAX": 1.53,
                "BUCKLE_MAX": 26.2,
                "PROPELLER_MIN": -13.62,
                "PROPELLER_MAX": 17.09,
                "OPENING_MIN": -45.62,
                "OPENING_MAX": 58.82
            },
            "t.W": {
                "SHEAR_MAX": 2.43,
                "STRETCH_MIN": -4.75,
                "STRETCH_MAX": -1.49,
                "STAGGER_MAX": 1.57,
                "BUCKLE_MAX": 27.75,
                "PROPELLER_MIN": -18.46,
                "PROPELLER_MAX": 13.71,
                "OPENING_MIN": -6.54,
                "OPENING_MAX": 200.89
            },
            "tSW": {
                "SHEAR_MAX": 2.27,
                "STRETCH_MIN": -7.76,
                "STRETCH_MAX": -3.98,
                "STAGGER_MAX": 0.95,
                "BUCKLE_MAX": 21.71,
                "PROPELLER_MIN": -10.54,
                "PROPELLER_MAX": 21.0,
                "OPENING_MIN": -147.22,
                "OPENING_MAX": -79.61
            },
            "cHW": {
                "SHEAR_MAX": 1.29,
                "STRETCH_MIN": -4.71,
                "STRETCH_MAX": -2.43,
                "STAGGER_MAX": 1.11,
                "BUCKLE_MAX": 22.99,
                "PROPELLER_MIN": -8.08,
                "PROPELLER_MAX": 16.09,
                "OPENING_MIN": 40.69,
                "OPENING_MAX": 95.7
            },
            "tWS": {
                "SHEAR_MAX": 2.33,
                "STRETCH_MIN": 2.82,
                "STRETCH_MAX": 7.41,
                "STAGGER_MAX": 1.04,
                "BUCKLE_MAX": 30.54,
                "PROPELLER_MIN": -31.19,
                "PROPELLER_MAX": 5.1,
                "OPENING_MIN": 83.62,
                "OPENING_MAX": 137.15
            },
            "t.H": {
                "SHEAR_MAX": 7.52,
                "STRETCH_MIN": -4.01,
                "STRETCH_MAX": 4.89,
                "STAGGER_MAX": 0.84,
                "BUCKLE_MAX": 16.59,
                "PROPELLER_MIN": -15.59,
                "PROPELLER_MAX": 16.28,
                "OPENING_MIN": -146.69,
                "OPENING_MAX": 0.43
            },
            "tHS": {
                "SHEAR_MAX": 7.3,
                "STRETCH_MIN": -9.61,
                "STRETCH_MAX": -3.24,
                "STAGGER_MAX": 1.34,
                "BUCKLE_MAX": 26.2,
                "PROPELLER_MIN": -25.74,
                "PROPELLER_MAX": 17.34,
                "OPENING_MIN": -66.42,
                "OPENING_MAX": -9.85
            },
            "tHH": {
                "SHEAR_MAX": 7.66,
                "STRETCH_MIN": -4.92,
                "STRETCH_MAX": 5.84,
                "STAGGER_MAX": 0.85,
                "BUCKLE_MAX": 16.21,
                "PROPELLER_MIN": -21.28,
                "PROPELLER_MAX": 15.39,
                "OPENING_MIN": -171.64,
                "OPENING_MAX": 157.66
            },
            "cHH": {
                "SHEAR_MAX": 3.25,
                "STRETCH_MIN": -8.61,
                "STRETCH_MAX": 2.2,
                "STAGGER_MAX": 1.09,
                "BUCKLE_MAX": 35.82,
                "PROPELLER_MIN": -25.63,
                "PROPELLER_MAX": 22.2,
                "OPENING_MIN": -75.18,
                "OPENING_MAX": 209.42
            },
            "tW.": {
                "SHEAR_MAX": 3.92,
                "STRETCH_MIN": 0.59,
                "STRETCH_MAX": 4.21,
                "STAGGER_MAX": 1.96,
                "BUCKLE_MAX": 23.22,
                "PROPELLER_MIN": -43.59,
                "PROPELLER_MAX": 13.16,
                "OPENING_MIN": -162.36,
                "OPENING_MAX": 13.01
            }
        },
        "G-U": {
            "cWW": {
                "SHEAR_MAX": 2.29,
                "STRETCH_MIN": -0.95,
                "STRETCH_MAX": 0.1,
                "STAGGER_MAX": 0.69,
                "BUCKLE_MAX": 11.15,
                "PROPELLER_MIN": -15.96,
                "PROPELLER_MAX": 3.47,
                "OPENING_MIN": -12.25,
                "OPENING_MAX": 12.39
            },
            "cSS": {
                "SHEAR_MAX": 2.69,
                "STRETCH_MIN": 5.11,
                "STRETCH_MAX": 8.97,
                "STAGGER_MAX": 1.11,
                "BUCKLE_MAX": 30.12,
                "PROPELLER_MIN": -28.82,
                "PROPELLER_MAX": -2.56,
                "OPENING_MIN": 107.32,
                "OPENING_MAX": 187.5
            },
            "cWS": {
                "SHEAR_MAX": 6.69,
                "STRETCH_MIN": -1.24,
                "STRETCH_MAX": 3.78,
                "STAGGER_MAX": 1.8,
                "BUCKLE_MAX": 29.27,
                "PROPELLER_MIN": -21.4,
                "PROPELLER_MAX": 30.4,
                "OPENING_MIN": 30.77,
                "OPENING_MAX": 109.12
            },
            "cHS": {
                "SHEAR_MAX": 8.01,
                "STRETCH_MIN": -3.26,
                "STRETCH_MAX": 0.21,
                "STAGGER_MAX": 1.5,
                "BUCKLE_MAX": 23.14,
                "PROPELLER_MIN": -51.33,
                "PROPELLER_MAX": 17.96,
                "OPENING_MIN": -31.38,
                "OPENING_MAX": 36.66
            },
            "tSW": {
                "SHEAR_MAX": 3.82,
                "STRETCH_MIN": -6.86,
                "STRETCH_MAX": -3.54,
                "STAGGER_MAX": 1.43,
                "BUCKLE_MAX": 35.74,
                "PROPELLER_MIN": -24.08,
                "PROPELLER_MAX": 16.89,
                "OPENING_MIN": -142.69,
                "OPENING_MAX": -46.83
            },
            "tWW": {
                "SHEAR_MAX": 1.83,
                "STRETCH_MIN": -4.03,
                "STRETCH_MAX": 2.36,
                "STAGGER_MAX": 1.4,
                "BUCKLE_MAX": 27.54,
                "PROPELLER_MIN": -31.4,
                "PROPELLER_MAX": 24.51,
                "OPENING_MIN": -183.49,
                "OPENING_MAX": 97.3
            },
            "c.W": {
                "SHEAR_MAX": 7.31,
                "STRETCH_MIN": -1.34,
                "STRETCH_MAX": 2.09,
                "STAGGER_MAX": 1.44,
                "BUCKLE_MAX": 20.9,
                "PROPELLER_MIN": -19.11,
                "PROPELLER_MAX": 15.74,
                "OPENING_MIN": -32.52,
                "OPENING_MAX": 53.3
            },
            "t.H": {
                "SHEAR_MAX": 8.87,
                "STRETCH_MIN": -3.81,
                "STRETCH_MAX": -1.65,
                "STAGGER_MAX": 1.28,
                "BUCKLE_MAX": 33.86,
                "PROPELLER_MIN": -13.78,
                "PROPELLER_MAX": 14.52,
                "OPENING_MIN": -45.01,
                "OPENING_MAX": -11.21
            },
            "cWH": {
                "SHEAR_MAX": 4.07,
                "STRETCH_MIN": 1.0,
                "STRETCH_MAX": 4.82,
                "STAGGER_MAX": 1.27,
                "BUCKLE_MAX": 38.63,
                "PROPELLER_MIN": -29.03,
                "PROPELLER_MAX": 5.11,
                "OPENING_MIN": -73.8,
                "OPENING_MAX": -26.82
            },
            "cSH": {
                "SHEAR_MAX": 8.42,
                "STRETCH_MIN": 0.9,
                "STRETCH_MAX": 2.43,
                "STAGGER_MAX": 1.06,
                "BUCKLE_MAX": 18.7,
                "PROPELLER_MIN": -9.25,
                "PROPELLER_MAX": 23.9,
                "OPENING_MIN": -2.54,
                "OPENING_MAX": 25.98
            },
            "tSH": {
                "SHEAR_MAX": 8.13,
                "STRETCH_MIN": -7.14,
                "STRETCH_MAX": -2.1,
                "STAGGER_MAX": 1.26,
                "BUCKLE_MAX": 24.05,
                "PROPELLER_MIN": -8.8,
                "PROPELLER_MAX": 21.53,
                "OPENING_MIN": -58.79,
                "OPENING_MAX": -2.8
            },
            "tHW": {
                "SHEAR_MAX": 6.85,
                "STRETCH_MIN": -1.85,
                "STRETCH_MAX": 1.39,
                "STAGGER_MAX": 1.58,
                "BUCKLE_MAX": 36.16,
                "PROPELLER_MIN": -31.98,
                "PROPELLER_MAX": 7.36,
                "OPENING_MIN": -128.11,
                "OPENING_MAX": -2.62
            },
            "cSW": {
                "SHEAR_MAX": 5.9,
                "STRETCH_MIN": -0.96,
                "STRETCH_MAX": 3.39,
                "STAGGER_MAX": 1.52,
                "BUCKLE_MAX": 30.55,
                "PROPELLER_MIN": -25.65,
                "PROPELLER_MAX": 31.94,
                "OPENING_MIN": 46.56,
                "OPENING_MAX": 115.7
            },
            "tWS": {
                "SHEAR_MAX": 2.62,
                "STRETCH_MIN": 3.16,
                "STRETCH_MAX": 8.09,
                "STAGGER_MAX": 1.18,
                "BUCKLE_MAX": 29.8,
                "PROPELLER_MIN": -24.7,
                "PROPELLER_MAX": 20.87,
                "OPENING_MIN": 72.58,
                "OPENING_MAX": 152.31
            },
            "cW.": {
                "SHEAR_MAX": 5.23,
                "STRETCH_MIN": -1.0,
                "STRETCH_MAX": 1.45,
                "STAGGER_MAX": 1.3,
                "BUCKLE_MAX": 20.73,
                "PROPELLER_MIN": -20.0,
                "PROPELLER_MAX": 16.8,
                "OPENING_MIN": -50.73,
                "OPENING_MAX": 28.81
            },
            "tW.": {
                "SHEAR_MAX": 4.47,
                "STRETCH_MIN": -0.3,
                "STRETCH_MAX": 2.67,
                "STAGGER_MAX": 1.35,
                "BUCKLE_MAX": 28.21,
                "PROPELLER_MIN": -32.23,
                "PROPELLER_MAX": 14.09,
                "OPENING_MIN": -65.62,
                "OPENING_MAX": 107.04
            },
            "t.W": {
                "SHEAR_MAX": 3.16,
                "STRETCH_MIN": -5.49,
                "STRETCH_MAX": -0.79,
                "STAGGER_MAX": 0.84,
                "BUCKLE_MAX": 16.29,
                "PROPELLER_MIN": -18.51,
                "PROPELLER_MAX": 3.86,
                "OPENING_MIN": -142.0,
                "OPENING_MAX": -11.48
            },
            "tSS": {
                "SHEAR_MAX": 4.21,
                "STRETCH_MIN": -8.04,
                "STRETCH_MAX": 9.12,
                "STAGGER_MAX": 1.6,
                "BUCKLE_MAX": 24.21,
                "PROPELLER_MIN": -30.03,
                "PROPELLER_MAX": 15.59,
                "OPENING_MIN": -138.91,
                "OPENING_MAX": 172.84
            },
            "tWH": {
                "SHEAR_MAX": 7.93,
                "STRETCH_MIN": -3.63,
                "STRETCH_MAX": 0.62,
                "STAGGER_MAX": 1.53,
                "BUCKLE_MAX": 27.24,
                "PROPELLER_MIN": -22.08,
                "PROPELLER_MAX": 23.21,
                "OPENING_MIN": -103.85,
                "OPENING_MAX": -28.86
            },
            "cHW": {
                "SHEAR_MAX": 3.86,
                "STRETCH_MIN": -4.65,
                "STRETCH_MAX": -0.34,
                "STAGGER_MAX": 1.52,
                "BUCKLE_MAX": 26.27,
                "PROPELLER_MIN": -27.0,
                "PROPELLER_MAX": 16.61,
                "OPENING_MIN": -9.71,
                "OPENING_MAX": 101.49
            },
            "tS.": {
                "SHEAR_MAX": 8.18,
                "STRETCH_MIN": -3.16,
                "STRETCH_MAX": -2.13,
                "STAGGER_MAX": 0.81,
                "BUCKLE_MAX": 34.28,
                "PROPELLER_MIN": -1.6,
                "PROPELLER_MAX": 28.55,
                "OPENING_MIN": -16.54,
                "OPENING_MAX": 17.02
            },
            "tHS": {
                "SHEAR_MAX": 8.24,
                "STRETCH_MIN": -5.18,
                "STRETCH_MAX": -1.61,
                "STAGGER_MAX": 1.04,
                "BUCKLE_MAX": 28.9,
                "PROPELLER_MIN": -19.22,
                "PROPELLER_MAX": 12.28,
                "OPENING_MIN": -41.41,
                "OPENING_MAX": 11.44
            }
        },
        "U-U": {
            "tHW": {
                "SHEAR_MAX": 4.55,
                "STRETCH_MIN": -1.94,
                "STRETCH_MAX": 0.78,
                "STAGGER_MAX": 0.96,
                "BUCKLE_MAX": 28.03,
                "PROPELLER_MIN": -21.11,
                "PROPELLER_MAX": 9.68,
                "OPENING_MIN": -107.73,
                "OPENING_MAX": -59.88
            },
            "cWW": {
                "SHEAR_MAX": 2.42,
                "STRETCH_MIN": -2.12,
                "STRETCH_MAX": -1.0,
                "STAGGER_MAX": 0.71,
                "BUCKLE_MAX": 14.65,
                "PROPELLER_MIN": -22.34,
                "PROPELLER_MAX": 3.17,
                "OPENING_MIN": -8.95,
                "OPENING_MAX": 25.68
            },
            "tWW": {
                "SHEAR_MAX": 2.54,
                "STRETCH_MIN": -3.28,
                "STRETCH_MAX": 2.17,
                "STAGGER_MAX": 0.73,
                "BUCKLE_MAX": 26.77,
                "PROPELLER_MIN": -18.26,
                "PROPELLER_MAX": 14.03,
                "OPENING_MIN": -228.97,
                "OPENING_MAX": 19.02
            },
            "cWS": {
                "SHEAR_MAX": 4.08,
                "STRETCH_MIN": 0.04,
                "STRETCH_MAX": 2.25,
                "STAGGER_MAX": 1.21,
                "BUCKLE_MAX": 20.51,
                "PROPELLER_MIN": 0.49,
                "PROPELLER_MAX": 45.79,
                "OPENING_MIN": 77.02,
                "OPENING_MAX": 112.9
            },
            "tWH": {
                "SHEAR_MAX": 5.73,
                "STRETCH_MIN": -4.02,
                "STRETCH_MAX": 0.46,
                "STAGGER_MAX": 1.06,
                "BUCKLE_MAX": 20.67,
                "PROPELLER_MIN": -14.64,
                "PROPELLER_MAX": 19.89,
                "OPENING_MIN": -97.53,
                "OPENING_MAX": -63.46
            },
            "cWH": {
                "SHEAR_MAX": 2.4,
                "STRETCH_MIN": 1.4,
                "STRETCH_MAX": 3.54,
                "STAGGER_MAX": 1.09,
                "BUCKLE_MAX": 21.45,
                "PROPELLER_MIN": -29.55,
                "PROPELLER_MAX": 7.46,
                "OPENING_MIN": -100.53,
                "OPENING_MAX": -48.63
            }
        },
        "A-A": {
            "cWW": {
                "SHEAR_MAX": 2.27,
                "STRETCH_MIN": 0.66,
                "STRETCH_MAX": 2.23,
                "STAGGER_MAX": 1.33,
                "BUCKLE_MAX": 39.81,
                "PROPELLER_MIN": -21.75,
                "PROPELLER_MAX": 5.34,
                "OPENING_MIN": -83.95,
                "OPENING_MAX": 45.38
            },
            "cWS": {
                "SHEAR_MAX": 7.15,
                "STRETCH_MIN": -2.06,
                "STRETCH_MAX": 3.59,
                "STAGGER_MAX": 1.31,
                "BUCKLE_MAX": 23.91,
                "PROPELLER_MIN": -17.14,
                "PROPELLER_MAX": 15.53,
                "OPENING_MIN": 33.14,
                "OPENING_MAX": 113.61
            },
            "tWS": {
                "SHEAR_MAX": 6.72,
                "STRETCH_MIN": -2.43,
                "STRETCH_MAX": 6.25,
                "STAGGER_MAX": 2.12,
                "BUCKLE_MAX": 34.47,
                "PROPELLER_MIN": -38.81,
                "PROPELLER_MAX": 2.74,
                "OPENING_MIN": 18.23,
                "OPENING_MAX": 126.66
            },
            "tWH": {
                "SHEAR_MAX": 5.23,
                "STRETCH_MIN": -0.11,
                "STRETCH_MAX": 2.12,
                "STAGGER_MAX": 1.26,
                "BUCKLE_MAX": 18.99,
                "PROPELLER_MIN": -20.91,
                "PROPELLER_MAX": 25.87,
                "OPENING_MIN": -118.56,
                "OPENING_MAX": -87.72
            },
            "cSW": {
                "SHEAR_MAX": 7.22,
                "STRETCH_MIN": -1.98,
                "STRETCH_MAX": 2.73,
                "STAGGER_MAX": 0.87,
                "BUCKLE_MAX": 31.37,
                "PROPELLER_MIN": -26.46,
                "PROPELLER_MAX": 2.6,
                "OPENING_MIN": 34.73,
                "OPENING_MAX": 98.15
            },
            "cSH": {
                "SHEAR_MAX": 7.14,
                "STRETCH_MIN": -1.09,
                "STRETCH_MAX": 0.99,
                "STAGGER_MAX": 1.39,
                "BUCKLE_MAX": 19.75,
                "PROPELLER_MIN": -9.31,
                "PROPELLER_MAX": 25.2,
                "OPENING_MIN": -10.79,
                "OPENING_MAX": 10.8
            },
            "cHS": {
                "SHEAR_MAX": 6.74,
                "STRETCH_MIN": -0.41,
                "STRETCH_MAX": 2.33,
                "STAGGER_MAX": 1.71,
                "BUCKLE_MAX": 35.09,
                "PROPELLER_MIN": -10.87,
                "PROPELLER_MAX": 55.41,
                "OPENING_MIN": -4.44,
                "OPENING_MAX": 29.01
            },
            "tHW": {
                "SHEAR_MAX": 5.47,
                "STRETCH_MIN": -0.67,
                "STRETCH_MAX": 2.33,
                "STAGGER_MAX": 1.06,
                "BUCKLE_MAX": 16.05,
                "PROPELLER_MIN": -23.38,
                "PROPELLER_MAX": 18.86,
                "OPENING_MIN": -130.58,
                "OPENING_MAX": -67.14
            },
            "cWH": {
                "SHEAR_MAX": 4.85,
                "STRETCH_MIN": 0.32,
                "STRETCH_MAX": 4.49,
                "STAGGER_MAX": 1.95,
                "BUCKLE_MAX": 53.06,
                "PROPELLER_MIN": -25.79,
                "PROPELLER_MAX": 11.39,
                "OPENING_MIN": -101.35,
                "OPENING_MAX": -33.24
            },
            "tWW": {
                "SHEAR_MAX": 1.73,
                "STRETCH_MIN": -1.46,
                "STRETCH_MAX": 1.38,
                "STAGGER_MAX": 1.09,
                "BUCKLE_MAX": 37.57,
                "PROPELLER_MIN": -15.4,
                "PROPELLER_MAX": 20.86,
                "OPENING_MIN": -171.1,
                "OPENING_MAX": 166.76
            },
            "tSH": {
                "SHEAR_MAX": 6.85,
                "STRETCH_MIN": -5.52,
                "STRETCH_MAX": -3.19,
                "STAGGER_MAX": 0.83,
                "BUCKLE_MAX": 20.88,
                "PROPELLER_MIN": -22.09,
                "PROPELLER_MAX": 12.16,
                "OPENING_MIN": -35.03,
                "OPENING_MAX": -5.01
            },
            "tHH": {
                "SHEAR_MAX": 6.22,
                "STRETCH_MIN": -4.94,
                "STRETCH_MAX": 5.41,
                "STAGGER_MAX": 0.65,
                "BUCKLE_MAX": 20.98,
                "PROPELLER_MIN": -15.29,
                "PROPELLER_MAX": 15.65,
                "OPENING_MIN": -181.53,
                "OPENING_MAX": 163.38
            },
            "tHS": {
                "SHEAR_MAX": 6.88,
                "STRETCH_MIN": -5.41,
                "STRETCH_MAX": -2.76,
                "STAGGER_MAX": 0.72,
                "BUCKLE_MAX": 11.23,
                "PROPELLER_MIN": -19.05,
                "PROPELLER_MAX": 10.89,
                "OPENING_MIN": -32.99,
                "OPENING_MAX": -0.68
            },
            "tSW": {
                "SHEAR_MAX": 3.2,
                "STRETCH_MIN": -7.72,
                "STRETCH_MAX": -3.35,
                "STAGGER_MAX": 0.93,
                "BUCKLE_MAX": 21.09,
                "PROPELLER_MIN": -7.09,
                "PROPELLER_MAX": 23.6,
                "OPENING_MIN": -148.62,
                "OPENING_MAX": -59.98
            },
            "cSS": {
                "SHEAR_MAX": 3.76,
                "STRETCH_MIN": 4.86,
                "STRETCH_MAX": 7.91,
                "STAGGER_MAX": 0.88,
                "BUCKLE_MAX": 25.19,
                "PROPELLER_MIN": -23.46,
                "PROPELLER_MAX": 1.89,
                "OPENING_MIN": 99.95,
                "OPENING_MAX": 176.05
            },
            "cHH": {
                "SHEAR_MAX": 6.47,
                "STRETCH_MIN": -8.46,
                "STRETCH_MAX": -2.2,
                "STAGGER_MAX": 0.95,
                "BUCKLE_MAX": 39.03,
                "PROPELLER_MIN": -6.94,
                "PROPELLER_MAX": 34.64,
                "OPENING_MIN": 46.89,
                "OPENING_MAX": 212.21
            }
        },
        "C-U": {
            "cHW": {
                "SHEAR_MAX": 4.08,
                "STRETCH_MIN": -6.43,
                "STRETCH_MAX": -1.16,
                "STAGGER_MAX": 1.12,
                "BUCKLE_MAX": 16.19,
                "PROPELLER_MIN": -12.95,
                "PROPELLER_MAX": 15.6,
                "OPENING_MIN": 26.98,
                "OPENING_MAX": 121.05
            },
            "cWW": {
                "SHEAR_MAX": 2.68,
                "STRETCH_MIN": -1.99,
                "STRETCH_MAX": -0.25,
                "STAGGER_MAX": 1.11,
                "BUCKLE_MAX": 21.06,
                "PROPELLER_MIN": -22.24,
                "PROPELLER_MAX": 8.99,
                "OPENING_MIN": -23.02,
                "OPENING_MAX": 29.4
            },
            "cSH": {
                "SHEAR_MAX": 6.97,
                "STRETCH_MIN": -0.69,
                "STRETCH_MAX": 1.36,
                "STAGGER_MAX": 1.58,
                "BUCKLE_MAX": 28.94,
                "PROPELLER_MIN": -17.93,
                "PROPELLER_MAX": 20.01,
                "OPENING_MIN": -39.3,
                "OPENING_MAX": 1.44
            },
            "cSW": {
                "SHEAR_MAX": 6.0,
                "STRETCH_MIN": -1.81,
                "STRETCH_MAX": 1.83,
                "STAGGER_MAX": 1.46,
                "BUCKLE_MAX": 40.54,
                "PROPELLER_MIN": -23.71,
                "PROPELLER_MAX": 42.69,
                "OPENING_MIN": 39.35,
                "OPENING_MAX": 118.92
            },
            "tWH": {
                "SHEAR_MAX": 6.91,
                "STRETCH_MIN": -5.08,
                "STRETCH_MAX": -0.52,
                "STAGGER_MAX": 1.4,
                "BUCKLE_MAX": 29.87,
                "PROPELLER_MIN": -19.12,
                "PROPELLER_MAX": 25.48,
                "OPENING_MIN": -91.68,
                "OPENING_MAX": -34.24
            },
            "tSW": {
                "SHEAR_MAX": 3.33,
                "STRETCH_MIN": -5.0,
                "STRETCH_MAX": -2.76,
                "STAGGER_MAX": 2.1,
                "BUCKLE_MAX": 36.96,
                "PROPELLER_MIN": -4.36,
                "PROPELLER_MAX": 44.7,
                "OPENING_MIN": -141.67,
                "OPENING_MAX": -48.02
            },
            "cWH": {
                "SHEAR_MAX": 5.21,
                "STRETCH_MIN": 0.3,
                "STRETCH_MAX": 3.1,
                "STAGGER_MAX": 1.69,
                "BUCKLE_MAX": 33.3,
                "PROPELLER_MIN": -26.18,
                "PROPELLER_MAX": 15.9,
                "OPENING_MIN": -82.31,
                "OPENING_MAX": -27.45
            },
            "cWS": {
                "SHEAR_MAX": 6.49,
                "STRETCH_MIN": -2.19,
                "STRETCH_MAX": 0.73,
                "STAGGER_MAX": 1.65,
                "BUCKLE_MAX": 28.53,
                "PROPELLER_MIN": -32.88,
                "PROPELLER_MAX": 18.72,
                "OPENING_MIN": 51.6,
                "OPENING_MAX": 89.63
            },
            "tW.": {
                "SHEAR_MAX": 5.77,
                "STRETCH_MIN": -2.48,
                "STRETCH_MAX": 1.72,
                "STAGGER_MAX": 1.07,
                "BUCKLE_MAX": 24.04,
                "PROPELLER_MIN": -10.33,
                "PROPELLER_MAX": 25.46,
                "OPENING_MIN": -124.48,
                "OPENING_MAX": -11.98
            },
            "tSH": {
                "SHEAR_MAX": 7.37,
                "STRETCH_MIN": -7.83,
                "STRETCH_MAX": -3.71,
                "STAGGER_MAX": 1.41,
                "BUCKLE_MAX": 22.77,
                "PROPELLER_MIN": -11.66,
                "PROPELLER_MAX": 28.63,
                "OPENING_MIN": -52.19,
                "OPENING_MAX": -6.69
            },
            "c.W": {
                "SHEAR_MAX": 3.6,
                "STRETCH_MIN": -1.91,
                "STRETCH_MAX": 1.59,
                "STAGGER_MAX": 2.37,
                "BUCKLE_MAX": 39.31,
                "PROPELLER_MIN": -9.29,
                "PROPELLER_MAX": 43.4,
                "OPENING_MIN": -103.93,
                "OPENING_MAX": 42.88
            },
            "t.W": {
                "SHEAR_MAX": 3.68,
                "STRETCH_MIN": -2.55,
                "STRETCH_MAX": 1.76,
                "STAGGER_MAX": 1.64,
                "BUCKLE_MAX": 27.65,
                "PROPELLER_MIN": -17.89,
                "PROPELLER_MAX": 20.43,
                "OPENING_MIN": -136.96,
                "OPENING_MAX": 17.13
            },
            "tHW": {
                "SHEAR_MAX": 7.16,
                "STRETCH_MIN": -5.24,
                "STRETCH_MAX": 0.28,
                "STAGGER_MAX": 1.29,
                "BUCKLE_MAX": 24.77,
                "PROPELLER_MIN": -20.36,
                "PROPELLER_MAX": 18.34,
                "OPENING_MIN": -94.1,
                "OPENING_MAX": -18.7
            },
            "tWW": {
                "SHEAR_MAX": 2.51,
                "STRETCH_MIN": -3.78,
                "STRETCH_MAX": 0.64,
                "STAGGER_MAX": 1.14,
                "BUCKLE_MAX": 30.28,
                "PROPELLER_MIN": -18.37,
                "PROPELLER_MAX": 32.9,
                "OPENING_MIN": -180.82,
                "OPENING_MAX": 62.76
            }
        },
        "G-G": {
            "cWH": {
                "SHEAR_MAX": 3.35,
                "STRETCH_MIN": 2.29,
                "STRETCH_MAX": 4.12,
                "STAGGER_MAX": 1.0,
                "BUCKLE_MAX": 28.27,
                "PROPELLER_MIN": -12.34,
                "PROPELLER_MAX": 16.27,
                "OPENING_MIN": -101.41,
                "OPENING_MAX": -62.17
            },
            "tSH": {
                "SHEAR_MAX": 7.87,
                "STRETCH_MIN": -5.43,
                "STRETCH_MAX": -2.27,
                "STAGGER_MAX": 1.04,
                "BUCKLE_MAX": 28.09,
                "PROPELLER_MIN": -18.07,
                "PROPELLER_MAX": 20.23,
                "OPENING_MIN": -51.94,
                "OPENING_MAX": -6.01
            },
            "tWH": {
                "SHEAR_MAX": 6.81,
                "STRETCH_MIN": -1.92,
                "STRETCH_MAX": 1.14,
                "STAGGER_MAX": 1.02,
                "BUCKLE_MAX": 34.19,
                "PROPELLER_MIN": -13.18,
                "PROPELLER_MAX": 24.08,
                "OPENING_MIN": -108.58,
                "OPENING_MAX": -52.51
            },
            "cSH": {
                "SHEAR_MAX": 7.62,
                "STRETCH_MIN": -1.63,
                "STRETCH_MAX": 3.57,
                "STAGGER_MAX": 1.81,
                "BUCKLE_MAX": 29.77,
                "PROPELLER_MIN": -28.79,
                "PROPELLER_MAX": 28.67,
                "OPENING_MIN": -30.49,
                "OPENING_MAX": 35.96
            },
            "cHW": {
                "SHEAR_MAX": 3.06,
                "STRETCH_MIN": -4.46,
                "STRETCH_MAX": -1.56,
                "STAGGER_MAX": 0.88,
                "BUCKLE_MAX": 23.57,
                "PROPELLER_MIN": -20.0,
                "PROPELLER_MAX": 22.34,
                "OPENING_MIN": 30.15,
                "OPENING_MAX": 120.93
            },
            "tHS": {
                "SHEAR_MAX": 7.86,
                "STRETCH_MIN": -4.84,
                "STRETCH_MAX": -1.78,
                "STAGGER_MAX": 0.95,
                "BUCKLE_MAX": 21.41,
                "PROPELLER_MIN": -17.01,
                "PROPELLER_MAX": 12.41,
                "OPENING_MIN": -37.71,
                "OPENING_MAX": 7.68
            },
            "tSS": {
                "SHEAR_MAX": 3.48,
                "STRETCH_MIN": -8.95,
                "STRETCH_MAX": 7.01,
                "STAGGER_MAX": 0.95,
                "BUCKLE_MAX": 28.82,
                "PROPELLER_MIN": -26.02,
                "PROPELLER_MAX": 19.62,
                "OPENING_MIN": -191.37,
                "OPENING_MAX": 151.16
            },
            "tWW": {
                "SHEAR_MAX": 1.95,
                "STRETCH_MIN": -3.25,
                "STRETCH_MAX": 2.43,
                "STAGGER_MAX": 1.04,
                "BUCKLE_MAX": 21.34,
                "PROPELLER_MIN": -16.04,
                "PROPELLER_MAX": 17.34,
                "OPENING_MIN": -224.21,
                "OPENING_MAX": 83.8
            },
            "cWW": {
                "SHEAR_MAX": 2.73,
                "STRETCH_MIN": 0.78,
                "STRETCH_MAX": 2.19,
                "STAGGER_MAX": 1.64,
                "BUCKLE_MAX": 23.33,
                "PROPELLER_MIN": -25.13,
                "PROPELLER_MAX": 19.88,
                "OPENING_MIN": -36.85,
                "OPENING_MAX": 10.24
            },
            "tSW": {
                "SHEAR_MAX": 5.24,
                "STRETCH_MIN": -8.46,
                "STRETCH_MAX": -0.12,
                "STAGGER_MAX": 1.58,
                "BUCKLE_MAX": 15.86,
                "PROPELLER_MIN": -7.88,
                "PROPELLER_MAX": 45.1,
                "OPENING_MIN": -152.99,
                "OPENING_MAX": 12.88
            },
            "tH.": {
                "SHEAR_MAX": 8.49,
                "STRETCH_MIN": -3.54,
                "STRETCH_MAX": 0.29,
                "STAGGER_MAX": 1.3,
                "BUCKLE_MAX": 11.01,
                "PROPELLER_MIN": -12.72,
                "PROPELLER_MAX": 6.33,
                "OPENING_MIN": -59.06,
                "OPENING_MAX": 1.18
            },
            "cWS": {
                "SHEAR_MAX": 6.9,
                "STRETCH_MIN": -0.5,
                "STRETCH_MAX": 2.45,
                "STAGGER_MAX": 1.68,
                "BUCKLE_MAX": 35.29,
                "PROPELLER_MIN": -16.62,
                "PROPELLER_MAX": 31.82,
                "OPENING_MIN": 19.72,
                "OPENING_MAX": 75.64
            },
            "cSS": {
                "SHEAR_MAX": 3.01,
                "STRETCH_MIN": -4.79,
                "STRETCH_MAX": 10.14,
                "STAGGER_MAX": 1.83,
                "BUCKLE_MAX": 33.8,
                "PROPELLER_MIN": -16.9,
                "PROPELLER_MAX": 21.21,
                "OPENING_MIN": -103.15,
                "OPENING_MAX": 212.17
            },
            "cHS": {
                "SHEAR_MAX": 7.27,
                "STRETCH_MIN": -2.67,
                "STRETCH_MAX": -0.41,
                "STAGGER_MAX": 2.05,
                "BUCKLE_MAX": 35.13,
                "PROPELLER_MIN": -19.53,
                "PROPELLER_MAX": 18.85,
                "OPENING_MIN": -11.85,
                "OPENING_MAX": 37.32
            },
            "tHW": {
                "SHEAR_MAX": 7.45,
                "STRETCH_MIN": -3.17,
                "STRETCH_MAX": 0.94,
                "STAGGER_MAX": 1.03,
                "BUCKLE_MAX": 25.02,
                "PROPELLER_MIN": -15.71,
                "PROPELLER_MAX": 32.8,
                "OPENING_MIN": -113.23,
                "OPENING_MAX": -41.97
            },
            "cW.": {
                "SHEAR_MAX": 5.91,
                "STRETCH_MIN": 0.76,
                "STRETCH_MAX": 2.46,
                "STAGGER_MAX": 1.04,
                "BUCKLE_MAX": 51.24,
                "PROPELLER_MIN": -26.17,
                "PROPELLER_MAX": 12.13,
                "OPENING_MIN": -59.1,
                "OPENING_MAX": 10.04
            },
            "t.H": {
                "SHEAR_MAX": 7.84,
                "STRETCH_MIN": -3.23,
                "STRETCH_MAX": -1.3,
                "STAGGER_MAX": 1.01,
                "BUCKLE_MAX": 32.32,
                "PROPELLER_MIN": -18.03,
                "PROPELLER_MAX": 7.63,
                "OPENING_MIN": -49.08,
                "OPENING_MAX": -5.21
            },
            "cSW": {
                "SHEAR_MAX": 6.89,
                "STRETCH_MIN": -1.31,
                "STRETCH_MAX": 3.93,
                "STAGGER_MAX": 1.75,
                "BUCKLE_MAX": 33.7,
                "PROPELLER_MIN": -28.13,
                "PROPELLER_MAX": 27.42,
                "OPENING_MIN": 11.63,
                "OPENING_MAX": 106.2
            }
        },
        "C-C": {
            "tHW": {
                "SHEAR_MAX": 7.27,
                "STRETCH_MIN": -5.1,
                "STRETCH_MAX": -1.94,
                "STAGGER_MAX": 1.84,
                "BUCKLE_MAX": 26.04,
                "PROPELLER_MIN": -24.96,
                "PROPELLER_MAX": 12.77,
                "OPENING_MIN": -88.76,
                "OPENING_MAX": 10.09
            },
            "cSW": {
                "SHEAR_MAX": 6.36,
                "STRETCH_MIN": -2.31,
                "STRETCH_MAX": -0.62,
                "STAGGER_MAX": 1.76,
                "BUCKLE_MAX": 55.54,
                "PROPELLER_MIN": -32.36,
                "PROPELLER_MAX": 3.3,
                "OPENING_MIN": 45.64,
                "OPENING_MAX": 79.81
            },
            "cSH": {
                "SHEAR_MAX": 7.01,
                "STRETCH_MIN": -1.66,
                "STRETCH_MAX": 0.94,
                "STAGGER_MAX": 2.55,
                "BUCKLE_MAX": 39.01,
                "PROPELLER_MIN": -3.55,
                "PROPELLER_MAX": 42.29,
                "OPENING_MIN": -41.53,
                "OPENING_MAX": 3.48
            },
            "cWW": {
                "SHEAR_MAX": 3.58,
                "STRETCH_MIN": -1.95,
                "STRETCH_MAX": -0.72,
                "STAGGER_MAX": 1.21,
                "BUCKLE_MAX": 17.0,
                "PROPELLER_MIN": -18.36,
                "PROPELLER_MAX": 11.86,
                "OPENING_MIN": -16.1,
                "OPENING_MAX": 30.01
            },
            "cWS": {
                "SHEAR_MAX": 6.53,
                "STRETCH_MIN": -1.96,
                "STRETCH_MAX": -0.47,
                "STAGGER_MAX": 0.95,
                "BUCKLE_MAX": 52.44,
                "PROPELLER_MIN": -23.46,
                "PROPELLER_MAX": 7.82,
                "OPENING_MIN": 55.18,
                "OPENING_MAX": 78.42
            },
            "tSH": {
                "SHEAR_MAX": 7.52,
                "STRETCH_MIN": -6.0,
                "STRETCH_MAX": -3.43,
                "STAGGER_MAX": 1.69,
                "BUCKLE_MAX": 31.73,
                "PROPELLER_MIN": -19.02,
                "PROPELLER_MAX": 31.47,
                "OPENING_MIN": -35.38,
                "OPENING_MAX": 2.94
            },
            "cWH": {
                "SHEAR_MAX": 4.08,
                "STRETCH_MIN": 0.54,
                "STRETCH_MAX": 2.77,
                "STAGGER_MAX": 1.57,
                "BUCKLE_MAX": 26.34,
                "PROPELLER_MIN": -27.87,
                "PROPELLER_MAX": 13.2,
                "OPENING_MIN": -79.24,
                "OPENING_MAX": -35.52
            },
            "tHS": {
                "SHEAR_MAX": 7.19,
                "STRETCH_MIN": -6.34,
                "STRETCH_MAX": -4.58,
                "STAGGER_MAX": 1.01,
                "BUCKLE_MAX": 24.7,
                "PROPELLER_MIN": -23.83,
                "PROPELLER_MAX": 2.28,
                "OPENING_MIN": -30.76,
                "OPENING_MAX": -8.04
            },
            "tWW": {
                "SHEAR_MAX": 2.18,
                "STRETCH_MIN": -0.43,
                "STRETCH_MAX": 2.96,
                "STAGGER_MAX": 1.01,
                "BUCKLE_MAX": 40.73,
                "PROPELLER_MIN": -6.65,
                "PROPELLER_MAX": 28.39,
                "OPENING_MIN": -97.29,
                "OPENING_MAX": 152.77
            }
        },
        "A-PSU": {
            "cWW": {
                "SHEAR_MAX": 1.66,
                "STRETCH_MIN": -0.32,
                "STRETCH_MAX": 0.36,
                "STAGGER_MAX": 0.6,
                "BUCKLE_MAX": 14.74,
                "PROPELLER_MIN": -15.49,
                "PROPELLER_MAX": 2.87,
                "OPENING_MIN": -9.05,
                "OPENING_MAX": 10.82
            }
        },
        "C-OMG": {
            "cWW": {
                "SHEAR_MAX": 0.53,
                "STRETCH_MIN": -0.29,
                "STRETCH_MAX": 0.16,
                "STAGGER_MAX": 0.49,
                "BUCKLE_MAX": 10.31,
                "PROPELLER_MIN": -14.5,
                "PROPELLER_MAX": 3.7,
                "OPENING_MIN": -4.34,
                "OPENING_MAX": 7.12
            }
        },
        "G-PSU": {
            "cWW": {
                "SHEAR_MAX": 2.26,
                "STRETCH_MIN": -0.61,
                "STRETCH_MAX": 0.1,
                "STAGGER_MAX": 0.49,
                "BUCKLE_MAX": 9.07,
                "PROPELLER_MIN": -12.46,
                "PROPELLER_MAX": 5.42,
                "OPENING_MIN": -10.52,
                "OPENING_MAX": 7.87
            }
        },
        "G-OMC": {
            "cWW": {
                "SHEAR_MAX": 0.53,
                "STRETCH_MIN": -0.29,
                "STRETCH_MAX": 0.14,
                "STAGGER_MAX": 0.53,
                "BUCKLE_MAX": 8.45,
                "PROPELLER_MIN": -13.52,
                "PROPELLER_MAX": 7.36,
                "OPENING_MIN": -3.73,
                "OPENING_MAX": 6.98
            }
        },
        "A-DT": {
            "cWW": {
                "SHEAR_MAX": 0.51,
                "STRETCH_MIN": -0.34,
                "STRETCH_MAX": 0.15,
                "STAGGER_MAX": 0.63,
                "BUCKLE_MAX": 12.37,
                "PROPELLER_MIN": -15.32,
                "PROPELLER_MAX": 4.47,
                "OPENING_MIN": -7.09,
                "OPENING_MAX": 8.42
            }
        },
        "C-DG": {
            "cWW": {
                "SHEAR_MAX": 0.55,
                "STRETCH_MIN": -0.36,
                "STRETCH_MAX": 0.11,
                "STAGGER_MAX": 0.5,
                "BUCKLE_MAX": 8.92,
                "PROPELLER_MIN": -13.36,
                "PROPELLER_MAX": 5.82,
                "OPENING_MIN": -4.42,
                "OPENING_MAX": 6.12
            }
        },
        "DA-U": {
            "cWW": {
                "SHEAR_MAX": 0.57,
                "STRETCH_MIN": -0.33,
                "STRETCH_MAX": 0.21,
                "STAGGER_MAX": 0.56,
                "BUCKLE_MAX": 11.34,
                "PROPELLER_MIN": -14.89,
                "PROPELLER_MAX": 7.37,
                "OPENING_MIN": -5.56,
                "OPENING_MAX": 9.02
            }
        },
        "DC-G": {
            "cWW": {
                "SHEAR_MAX": 0.63,
                "STRETCH_MIN": -0.38,
                "STRETCH_MAX": 0.16,
                "STAGGER_MAX": 0.64,
                "BUCKLE_MAX": 12.4,
                "PROPELLER_MIN": -13.25,
                "PROPELLER_MAX": 5.14,
                "OPENING_MIN": -5.11,
                "OPENING_MAX": 8.19
            }
        },
        "DC-DG": {
            "cWW": {
                "SHEAR_MAX": 0.58,
                "STRETCH_MIN": -0.45,
                "STRETCH_MAX": 0.14,
                "STAGGER_MAX": 0.65,
                "BUCKLE_MAX": 9.57,
                "PROPELLER_MIN": -13.16,
                "PROPELLER_MAX": 3.64,
                "OPENING_MIN": -5.79,
                "OPENING_MAX": 7.7
            }
        },
        "DA-DT": {
            "cWW": {
                "SHEAR_MAX": 0.53,
                "STRETCH_MIN": -0.4,
                "STRETCH_MAX": 0.14,
                "STAGGER_MAX": 0.71,
                "BUCKLE_MAX": 9.81,
                "PROPELLER_MIN": -16.73,
                "PROPELLER_MAX": 1.8,
                "OPENING_MIN": -6.98,
                "OPENING_MAX": 9.06
            }
        }
    }
    
    # Base-pair + edge H-bond thresholds (mean ± SD, 2dp) generated from edge_type_analysis.json
    HBOND_THRESHOLDS_BY_BASE_PAIR_BY_EDGE = {
        "A-C": {
            "c.W": {
                "DIST_MIN": 2.72,
                "DIST_MAX": 3.46,
                "ANGLE_MIN": 85.07,
                "QUALITY_MIN": 0.52,
                "DIHEDRAL_FORBID_MIN": -115.15,
                "DIHEDRAL_FORBID_MAX": 151.02
            },
            "cHH": {
                "DIST_MIN": 2.82,
                "DIST_MAX": 3.48,
                "ANGLE_MIN": 77.45,
                "QUALITY_MIN": 0.55,
                "DIHEDRAL_FORBID_MIN": -68.16,
                "DIHEDRAL_FORBID_MAX": 94.78
            },
            "cHS": {
                "DIST_MIN": 2.76,
                "DIST_MAX": 3.43,
                "ANGLE_MIN": 76.9,
                "QUALITY_MIN": 0.57,
                "DIHEDRAL_FORBID_MIN": -88.46,
                "DIHEDRAL_FORBID_MAX": 120.83
            },
            "cHW": {
                "DIST_MIN": 2.74,
                "DIST_MAX": 3.39,
                "ANGLE_MIN": 81.82,
                "QUALITY_MIN": 0.6,
                "DIHEDRAL_FORBID_MIN": -94.49,
                "DIHEDRAL_FORBID_MAX": 89.3
            },
            "cSH": {
                "DIST_MIN": 2.89,
                "DIST_MAX": 3.47,
                "ANGLE_MIN": 73.84,
                "QUALITY_MIN": 0.54,
                "DIHEDRAL_FORBID_MIN": -116.83,
                "DIHEDRAL_FORBID_MAX": 111.72
            },
            "cSS": {
                "DIST_MIN": 2.63,
                "DIST_MAX": 3.22,
                "ANGLE_MIN": 94.04,
                "QUALITY_MIN": 0.7,
                "DIHEDRAL_FORBID_MIN": -21.02,
                "DIHEDRAL_FORBID_MAX": 137.08
            },
            "cSW": {
                "DIST_MIN": 2.62,
                "DIST_MAX": 3.22,
                "ANGLE_MIN": 92.66,
                "QUALITY_MIN": 0.67,
                "DIHEDRAL_FORBID_MIN": -49.37,
                "DIHEDRAL_FORBID_MAX": 160.75
            },
            "cW.": {
                "DIST_MIN": 2.77,
                "DIST_MAX": 3.33,
                "ANGLE_MIN": 85.66,
                "QUALITY_MIN": 0.59,
                "DIHEDRAL_FORBID_MIN": -128.59,
                "DIHEDRAL_FORBID_MAX": 154.54
            },
            "cWH": {
                "DIST_MIN": 2.82,
                "DIST_MAX": 3.45,
                "ANGLE_MIN": 75.72,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -94.11,
                "DIHEDRAL_FORBID_MAX": 106.29
            },
            "cWS": {
                "DIST_MIN": 2.69,
                "DIST_MAX": 3.27,
                "ANGLE_MIN": 92.26,
                "QUALITY_MIN": 0.66,
                "DIHEDRAL_FORBID_MIN": -67.14,
                "DIHEDRAL_FORBID_MAX": 147.28
            },
            "cWW": {
                "DIST_MIN": 2.82,
                "DIST_MAX": 3.38,
                "ANGLE_MIN": 91.72,
                "QUALITY_MIN": 0.62,
                "DIHEDRAL_FORBID_MIN": -103.35,
                "DIHEDRAL_FORBID_MAX": 148.89
            },
            "t.H": {
                "DIST_MIN": 2.66,
                "DIST_MAX": 3.39,
                "ANGLE_MIN": 97.43,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -90.51,
                "DIHEDRAL_FORBID_MAX": 188.1
            },
            "t.W": {
                "DIST_MIN": 2.85,
                "DIST_MAX": 3.41,
                "ANGLE_MIN": 93.01,
                "QUALITY_MIN": 0.54,
                "DIHEDRAL_FORBID_MIN": -111.22,
                "DIHEDRAL_FORBID_MAX": 148.08
            },
            "tH.": {
                "DIST_MIN": 2.69,
                "DIST_MAX": 3.32,
                "ANGLE_MIN": 101.93,
                "QUALITY_MIN": 0.6,
                "DIHEDRAL_FORBID_MIN": -68.93,
                "DIHEDRAL_FORBID_MAX": 206.28
            },
            "tHH": {
                "DIST_MIN": 2.86,
                "DIST_MAX": 3.38,
                "ANGLE_MIN": 96.44,
                "QUALITY_MIN": 0.6,
                "DIHEDRAL_FORBID_MIN": -56.73,
                "DIHEDRAL_FORBID_MAX": 175.93
            },
            "tHS": {
                "DIST_MIN": 2.77,
                "DIST_MAX": 3.39,
                "ANGLE_MIN": 92.96,
                "QUALITY_MIN": 0.6,
                "DIHEDRAL_FORBID_MIN": -75.77,
                "DIHEDRAL_FORBID_MAX": 157.42
            },
            "tHW": {
                "DIST_MIN": 2.78,
                "DIST_MAX": 3.3,
                "ANGLE_MIN": 101.65,
                "QUALITY_MIN": 0.66,
                "DIHEDRAL_FORBID_MIN": -87.5,
                "DIHEDRAL_FORBID_MAX": 153.42
            },
            "tSH": {
                "DIST_MIN": 2.78,
                "DIST_MAX": 3.4,
                "ANGLE_MIN": 96.51,
                "QUALITY_MIN": 0.62,
                "DIHEDRAL_FORBID_MIN": -60.87,
                "DIHEDRAL_FORBID_MAX": 146.79
            },
            "tSW": {
                "DIST_MIN": 2.67,
                "DIST_MAX": 3.32,
                "ANGLE_MIN": 88.49,
                "QUALITY_MIN": 0.65,
                "DIHEDRAL_FORBID_MIN": -71.8,
                "DIHEDRAL_FORBID_MAX": 45.71
            },
            "tWH": {
                "DIST_MIN": 2.77,
                "DIST_MAX": 3.29,
                "ANGLE_MIN": 103.85,
                "QUALITY_MIN": 0.64,
                "DIHEDRAL_FORBID_MIN": -111.06,
                "DIHEDRAL_FORBID_MAX": 134.74
            },
            "tWS": {
                "DIST_MIN": 2.69,
                "DIST_MAX": 3.37,
                "ANGLE_MIN": 90.06,
                "QUALITY_MIN": 0.63,
                "DIHEDRAL_FORBID_MIN": -53.41,
                "DIHEDRAL_FORBID_MAX": 65.07
            },
            "tWW": {
                "DIST_MIN": 2.85,
                "DIST_MAX": 3.36,
                "ANGLE_MIN": 97.96,
                "QUALITY_MIN": 0.65,
                "DIHEDRAL_FORBID_MIN": -79.58,
                "DIHEDRAL_FORBID_MAX": 141.94
            },
            "_OTHER": {
                "DIST_MIN": 2.76,
                "DIST_MAX": 3.38,
                "ANGLE_MIN": 90.56,
                "QUALITY_MIN": 0.58,
                "DIHEDRAL_FORBID_MIN": -108.78,
                "DIHEDRAL_FORBID_MAX": 117.29
            }
        },
        "A-G": {
            "--": {
                "DIST_MIN": 2.97,
                "DIST_MAX": 3.49,
                "ANGLE_MIN": 128.63,
                "QUALITY_MIN": 0.49,
                "DIHEDRAL_FORBID_MIN": -141.81,
                "DIHEDRAL_FORBID_MAX": 75.27
            },
            "c.H": {
                "DIST_MIN": 2.88,
                "DIST_MAX": 3.38,
                "ANGLE_MIN": 114.37,
                "QUALITY_MIN": 0.62,
                "DIHEDRAL_FORBID_MIN": -31.62,
                "DIHEDRAL_FORBID_MAX": 177.27
            },
            "c.W": {
                "DIST_MIN": 2.8,
                "DIST_MAX": 3.36,
                "ANGLE_MIN": 89.91,
                "QUALITY_MIN": 0.63,
                "DIHEDRAL_FORBID_MIN": -50.77,
                "DIHEDRAL_FORBID_MAX": 124.28
            },
            "cH.": {
                "DIST_MIN": 2.67,
                "DIST_MAX": 3.33,
                "ANGLE_MIN": 120.6,
                "QUALITY_MIN": 0.54,
                "DIHEDRAL_FORBID_MIN": -181.34,
                "DIHEDRAL_FORBID_MAX": 62.02
            },
            "cHH": {
                "DIST_MIN": 2.72,
                "DIST_MAX": 3.42,
                "ANGLE_MIN": 83.13,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -116.63,
                "DIHEDRAL_FORBID_MAX": 101.55
            },
            "cHS": {
                "DIST_MIN": 2.88,
                "DIST_MAX": 3.47,
                "ANGLE_MIN": 75.95,
                "QUALITY_MIN": 0.55,
                "DIHEDRAL_FORBID_MIN": -87.06,
                "DIHEDRAL_FORBID_MAX": 121.12
            },
            "cHW": {
                "DIST_MIN": 2.75,
                "DIST_MAX": 3.34,
                "ANGLE_MIN": 90.18,
                "QUALITY_MIN": 0.65,
                "DIHEDRAL_FORBID_MIN": -63.0,
                "DIHEDRAL_FORBID_MAX": 79.32
            },
            "cSH": {
                "DIST_MIN": 2.87,
                "DIST_MAX": 3.48,
                "ANGLE_MIN": 78.07,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -80.58,
                "DIHEDRAL_FORBID_MAX": 111.72
            },
            "cSS": {
                "DIST_MIN": 2.76,
                "DIST_MAX": 3.38,
                "ANGLE_MIN": 91.24,
                "QUALITY_MIN": 0.6,
                "DIHEDRAL_FORBID_MIN": -55.25,
                "DIHEDRAL_FORBID_MAX": 157.37
            },
            "cSW": {
                "DIST_MIN": 2.74,
                "DIST_MAX": 3.37,
                "ANGLE_MIN": 87.62,
                "QUALITY_MIN": 0.63,
                "DIHEDRAL_FORBID_MIN": -72.77,
                "DIHEDRAL_FORBID_MAX": 107.37
            },
            "cW.": {
                "DIST_MIN": 2.88,
                "DIST_MAX": 3.45,
                "ANGLE_MIN": 89.78,
                "QUALITY_MIN": 0.52,
                "DIHEDRAL_FORBID_MIN": -95.16,
                "DIHEDRAL_FORBID_MAX": 156.31
            },
            "cWH": {
                "DIST_MIN": 2.77,
                "DIST_MAX": 3.39,
                "ANGLE_MIN": 94.32,
                "QUALITY_MIN": 0.6,
                "DIHEDRAL_FORBID_MIN": -98.54,
                "DIHEDRAL_FORBID_MAX": 71.43
            },
            "cWS": {
                "DIST_MIN": 2.73,
                "DIST_MAX": 3.37,
                "ANGLE_MIN": 90.22,
                "QUALITY_MIN": 0.64,
                "DIHEDRAL_FORBID_MIN": -75.36,
                "DIHEDRAL_FORBID_MAX": 94.39
            },
            "cWW": {
                "DIST_MIN": 2.77,
                "DIST_MAX": 3.23,
                "ANGLE_MIN": 96.88,
                "QUALITY_MIN": 0.76,
                "DIHEDRAL_FORBID_MIN": -41.8,
                "DIHEDRAL_FORBID_MAX": 50.36
            },
            "t.S": {
                "DIST_MIN": 2.69,
                "DIST_MAX": 3.39,
                "ANGLE_MIN": 86.22,
                "QUALITY_MIN": 0.63,
                "DIHEDRAL_FORBID_MIN": -64.94,
                "DIHEDRAL_FORBID_MAX": 14.32
            },
            "t.W": {
                "DIST_MIN": 2.92,
                "DIST_MAX": 3.47,
                "ANGLE_MIN": 77.65,
                "QUALITY_MIN": 0.53,
                "DIHEDRAL_FORBID_MIN": -91.12,
                "DIHEDRAL_FORBID_MAX": 192.83
            },
            "tHH": {
                "DIST_MIN": 2.81,
                "DIST_MAX": 3.39,
                "ANGLE_MIN": 84.7,
                "QUALITY_MIN": 0.6,
                "DIHEDRAL_FORBID_MIN": -57.65,
                "DIHEDRAL_FORBID_MAX": 175.18
            },
            "tHS": {
                "DIST_MIN": 2.78,
                "DIST_MAX": 3.3,
                "ANGLE_MIN": 100.74,
                "QUALITY_MIN": 0.64,
                "DIHEDRAL_FORBID_MIN": -101.19,
                "DIHEDRAL_FORBID_MAX": 155.07
            },
            "tHW": {
                "DIST_MIN": 2.76,
                "DIST_MAX": 3.38,
                "ANGLE_MIN": 87.31,
                "QUALITY_MIN": 0.58,
                "DIHEDRAL_FORBID_MIN": -99.1,
                "DIHEDRAL_FORBID_MAX": 152.17
            },
            "tS.": {
                "DIST_MIN": 2.66,
                "DIST_MAX": 3.39,
                "ANGLE_MIN": 86.08,
                "QUALITY_MIN": 0.63,
                "DIHEDRAL_FORBID_MIN": -61.29,
                "DIHEDRAL_FORBID_MAX": 32.76
            },
            "tSH": {
                "DIST_MIN": 2.76,
                "DIST_MAX": 3.32,
                "ANGLE_MIN": 97.14,
                "QUALITY_MIN": 0.62,
                "DIHEDRAL_FORBID_MIN": -99.97,
                "DIHEDRAL_FORBID_MAX": 161.55
            },
            "tSS": {
                "DIST_MIN": 2.73,
                "DIST_MAX": 3.32,
                "ANGLE_MIN": 89.69,
                "QUALITY_MIN": 0.69,
                "DIHEDRAL_FORBID_MIN": -56.5,
                "DIHEDRAL_FORBID_MAX": 30.31
            },
            "tSW": {
                "DIST_MIN": 2.75,
                "DIST_MAX": 3.37,
                "ANGLE_MIN": 89.8,
                "QUALITY_MIN": 0.64,
                "DIHEDRAL_FORBID_MIN": -75.66,
                "DIHEDRAL_FORBID_MAX": 114.24
            },
            "tW.": {
                "DIST_MIN": 3.09,
                "DIST_MAX": 3.64,
                "ANGLE_MIN": 69.71,
                "QUALITY_MIN": 0.42,
                "DIHEDRAL_FORBID_MIN": -95.61,
                "DIHEDRAL_FORBID_MAX": 187.29
            },
            "tWH": {
                "DIST_MIN": 2.77,
                "DIST_MAX": 3.38,
                "ANGLE_MIN": 89.63,
                "QUALITY_MIN": 0.59,
                "DIHEDRAL_FORBID_MIN": -92.64,
                "DIHEDRAL_FORBID_MAX": 144.97
            },
            "tWS": {
                "DIST_MIN": 2.8,
                "DIST_MAX": 3.38,
                "ANGLE_MIN": 92.36,
                "QUALITY_MIN": 0.66,
                "DIHEDRAL_FORBID_MIN": -52.51,
                "DIHEDRAL_FORBID_MAX": 129.88
            },
            "tWW": {
                "DIST_MIN": 2.89,
                "DIST_MAX": 3.52,
                "ANGLE_MIN": 77.24,
                "QUALITY_MIN": 0.51,
                "DIHEDRAL_FORBID_MIN": -90.91,
                "DIHEDRAL_FORBID_MAX": 127.93
            },
            "_OTHER": {
                "DIST_MIN": 2.77,
                "DIST_MAX": 3.41,
                "ANGLE_MIN": 99.31,
                "QUALITY_MIN": 0.59,
                "DIHEDRAL_FORBID_MIN": -69.74,
                "DIHEDRAL_FORBID_MAX": 178.86
            }
        },
        "A-U": {
            "c.W": {
                "DIST_MIN": 2.79,
                "DIST_MAX": 3.37,
                "ANGLE_MIN": 105.05,
                "QUALITY_MIN": 0.53,
                "DIHEDRAL_FORBID_MIN": -181.21,
                "DIHEDRAL_FORBID_MAX": 64.87
            },
            "cHH": {
                "DIST_MIN": 2.59,
                "DIST_MAX": 3.31,
                "ANGLE_MIN": 96.76,
                "QUALITY_MIN": 0.59,
                "DIHEDRAL_FORBID_MIN": -99.1,
                "DIHEDRAL_FORBID_MAX": 71.33
            },
            "cHS": {
                "DIST_MIN": 2.68,
                "DIST_MAX": 3.28,
                "ANGLE_MIN": 96.0,
                "QUALITY_MIN": 0.62,
                "DIHEDRAL_FORBID_MIN": -126.39,
                "DIHEDRAL_FORBID_MAX": 121.0
            },
            "cHW": {
                "DIST_MIN": 2.74,
                "DIST_MAX": 3.24,
                "ANGLE_MIN": 98.59,
                "QUALITY_MIN": 0.74,
                "DIHEDRAL_FORBID_MIN": -54.28,
                "DIHEDRAL_FORBID_MAX": 56.88
            },
            "cSH": {
                "DIST_MIN": 2.75,
                "DIST_MAX": 3.41,
                "ANGLE_MIN": 81.9,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -114.21,
                "DIHEDRAL_FORBID_MAX": 103.66
            },
            "cSS": {
                "DIST_MIN": 2.7,
                "DIST_MAX": 3.27,
                "ANGLE_MIN": 96.28,
                "QUALITY_MIN": 0.69,
                "DIHEDRAL_FORBID_MIN": -17.29,
                "DIHEDRAL_FORBID_MAX": 126.7
            },
            "cSW": {
                "DIST_MIN": 2.65,
                "DIST_MAX": 3.28,
                "ANGLE_MIN": 88.73,
                "QUALITY_MIN": 0.62,
                "DIHEDRAL_FORBID_MIN": -99.13,
                "DIHEDRAL_FORBID_MAX": 128.64
            },
            "cW.": {
                "DIST_MIN": 2.67,
                "DIST_MAX": 3.27,
                "ANGLE_MIN": 106.94,
                "QUALITY_MIN": 0.63,
                "DIHEDRAL_FORBID_MIN": -56.56,
                "DIHEDRAL_FORBID_MAX": 182.52
            },
            "cWH": {
                "DIST_MIN": 2.76,
                "DIST_MAX": 3.31,
                "ANGLE_MIN": 92.51,
                "QUALITY_MIN": 0.67,
                "DIHEDRAL_FORBID_MIN": -69.88,
                "DIHEDRAL_FORBID_MAX": 49.33
            },
            "cWS": {
                "DIST_MIN": 2.74,
                "DIST_MAX": 3.34,
                "ANGLE_MIN": 91.86,
                "QUALITY_MIN": 0.62,
                "DIHEDRAL_FORBID_MIN": -82.12,
                "DIHEDRAL_FORBID_MAX": 130.08
            },
            "cWW": {
                "DIST_MIN": 2.76,
                "DIST_MAX": 3.13,
                "ANGLE_MIN": 106.33,
                "QUALITY_MIN": 0.84,
                "DIHEDRAL_FORBID_MIN": -29.05,
                "DIHEDRAL_FORBID_MAX": 15.9
            },
            "t.H": {
                "DIST_MIN": 2.72,
                "DIST_MAX": 3.36,
                "ANGLE_MIN": 104.85,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -150.15,
                "DIHEDRAL_FORBID_MAX": 149.28
            },
            "t.W": {
                "DIST_MIN": 2.86,
                "DIST_MAX": 3.43,
                "ANGLE_MIN": 99.62,
                "QUALITY_MIN": 0.57,
                "DIHEDRAL_FORBID_MIN": -82.85,
                "DIHEDRAL_FORBID_MAX": 184.97
            },
            "tHH": {
                "DIST_MIN": 2.61,
                "DIST_MAX": 3.31,
                "ANGLE_MIN": 97.22,
                "QUALITY_MIN": 0.59,
                "DIHEDRAL_FORBID_MIN": -68.28,
                "DIHEDRAL_FORBID_MAX": 179.41
            },
            "tHS": {
                "DIST_MIN": 2.72,
                "DIST_MAX": 3.41,
                "ANGLE_MIN": 90.18,
                "QUALITY_MIN": 0.58,
                "DIHEDRAL_FORBID_MIN": -102.81,
                "DIHEDRAL_FORBID_MAX": 63.64
            },
            "tHW": {
                "DIST_MIN": 2.7,
                "DIST_MAX": 3.18,
                "ANGLE_MIN": 105.69,
                "QUALITY_MIN": 0.7,
                "DIHEDRAL_FORBID_MIN": -104.84,
                "DIHEDRAL_FORBID_MAX": 134.08
            },
            "tSH": {
                "DIST_MIN": 2.68,
                "DIST_MAX": 3.37,
                "ANGLE_MIN": 84.25,
                "QUALITY_MIN": 0.63,
                "DIHEDRAL_FORBID_MIN": -33.43,
                "DIHEDRAL_FORBID_MAX": 120.75
            },
            "tSW": {
                "DIST_MIN": 2.71,
                "DIST_MAX": 3.36,
                "ANGLE_MIN": 87.56,
                "QUALITY_MIN": 0.64,
                "DIHEDRAL_FORBID_MIN": -52.07,
                "DIHEDRAL_FORBID_MAX": 54.0
            },
            "tW.": {
                "DIST_MIN": 2.67,
                "DIST_MAX": 3.36,
                "ANGLE_MIN": 104.15,
                "QUALITY_MIN": 0.61,
                "DIHEDRAL_FORBID_MIN": -55.12,
                "DIHEDRAL_FORBID_MAX": 190.77
            },
            "tWH": {
                "DIST_MIN": 2.7,
                "DIST_MAX": 3.19,
                "ANGLE_MIN": 106.06,
                "QUALITY_MIN": 0.69,
                "DIHEDRAL_FORBID_MIN": -109.4,
                "DIHEDRAL_FORBID_MAX": 131.08
            },
            "tWS": {
                "DIST_MIN": 2.71,
                "DIST_MAX": 3.41,
                "ANGLE_MIN": 87.33,
                "QUALITY_MIN": 0.59,
                "DIHEDRAL_FORBID_MIN": -74.65,
                "DIHEDRAL_FORBID_MAX": 49.45
            },
            "tWW": {
                "DIST_MIN": 2.72,
                "DIST_MAX": 3.2,
                "ANGLE_MIN": 97.98,
                "QUALITY_MIN": 0.69,
                "DIHEDRAL_FORBID_MIN": -108.32,
                "DIHEDRAL_FORBID_MAX": 125.51
            },
            "_OTHER": {
                "DIST_MIN": 2.67,
                "DIST_MAX": 3.34,
                "ANGLE_MIN": 96.81,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -143.11,
                "DIHEDRAL_FORBID_MAX": 113.56
            }
        },
        "C-G": {
            "--": {
                "DIST_MIN": 2.74,
                "DIST_MAX": 3.29,
                "ANGLE_MIN": 123.47,
                "QUALITY_MIN": 0.62,
                "DIHEDRAL_FORBID_MIN": -91.59,
                "DIHEDRAL_FORBID_MAX": 175.69
            },
            "c.H": {
                "DIST_MIN": 2.85,
                "DIST_MAX": 3.51,
                "ANGLE_MIN": 101.72,
                "QUALITY_MIN": 0.55,
                "DIHEDRAL_FORBID_MIN": -28.29,
                "DIHEDRAL_FORBID_MAX": 165.94
            },
            "c.W": {
                "DIST_MIN": 2.82,
                "DIST_MAX": 3.43,
                "ANGLE_MIN": 87.89,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -80.35,
                "DIHEDRAL_FORBID_MAX": 148.95
            },
            "cHH": {
                "DIST_MIN": 2.7,
                "DIST_MAX": 3.39,
                "ANGLE_MIN": 80.95,
                "QUALITY_MIN": 0.59,
                "DIHEDRAL_FORBID_MIN": -80.28,
                "DIHEDRAL_FORBID_MAX": 79.16
            },
            "cHS": {
                "DIST_MIN": 2.94,
                "DIST_MAX": 3.55,
                "ANGLE_MIN": 69.13,
                "QUALITY_MIN": 0.51,
                "DIHEDRAL_FORBID_MIN": -90.01,
                "DIHEDRAL_FORBID_MAX": 82.12
            },
            "cHW": {
                "DIST_MIN": 2.69,
                "DIST_MAX": 3.34,
                "ANGLE_MIN": 88.54,
                "QUALITY_MIN": 0.63,
                "DIHEDRAL_FORBID_MIN": -73.56,
                "DIHEDRAL_FORBID_MAX": 100.19
            },
            "cSH": {
                "DIST_MIN": 2.89,
                "DIST_MAX": 3.52,
                "ANGLE_MIN": 75.33,
                "QUALITY_MIN": 0.55,
                "DIHEDRAL_FORBID_MIN": -72.52,
                "DIHEDRAL_FORBID_MAX": 110.23
            },
            "cSS": {
                "DIST_MIN": 2.71,
                "DIST_MAX": 3.31,
                "ANGLE_MIN": 97.34,
                "QUALITY_MIN": 0.64,
                "DIHEDRAL_FORBID_MIN": -69.75,
                "DIHEDRAL_FORBID_MAX": 107.15
            },
            "cSW": {
                "DIST_MIN": 2.83,
                "DIST_MAX": 3.45,
                "ANGLE_MIN": 78.52,
                "QUALITY_MIN": 0.57,
                "DIHEDRAL_FORBID_MIN": -77.29,
                "DIHEDRAL_FORBID_MAX": 112.13
            },
            "cW.": {
                "DIST_MIN": 2.75,
                "DIST_MAX": 3.35,
                "ANGLE_MIN": 95.31,
                "QUALITY_MIN": 0.61,
                "DIHEDRAL_FORBID_MIN": -94.52,
                "DIHEDRAL_FORBID_MAX": 130.64
            },
            "cWH": {
                "DIST_MIN": 2.86,
                "DIST_MAX": 3.49,
                "ANGLE_MIN": 73.44,
                "QUALITY_MIN": 0.54,
                "DIHEDRAL_FORBID_MIN": -88.32,
                "DIHEDRAL_FORBID_MAX": 82.67
            },
            "cWS": {
                "DIST_MIN": 2.83,
                "DIST_MAX": 3.38,
                "ANGLE_MIN": 83.94,
                "QUALITY_MIN": 0.6,
                "DIHEDRAL_FORBID_MIN": -92.86,
                "DIHEDRAL_FORBID_MAX": 117.08
            },
            "cWW": {
                "DIST_MIN": 2.72,
                "DIST_MAX": 3.1,
                "ANGLE_MIN": 107.03,
                "QUALITY_MIN": 0.81,
                "DIHEDRAL_FORBID_MIN": -64.74,
                "DIHEDRAL_FORBID_MAX": 121.83
            },
            "t.W": {
                "DIST_MIN": 2.71,
                "DIST_MAX": 3.29,
                "ANGLE_MIN": 94.42,
                "QUALITY_MIN": 0.62,
                "DIHEDRAL_FORBID_MIN": -109.9,
                "DIHEDRAL_FORBID_MAX": 142.62
            },
            "tHH": {
                "DIST_MIN": 2.62,
                "DIST_MAX": 3.29,
                "ANGLE_MIN": 93.4,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -98.46,
                "DIHEDRAL_FORBID_MAX": 173.26
            },
            "tHS": {
                "DIST_MIN": 2.77,
                "DIST_MAX": 3.38,
                "ANGLE_MIN": 86.14,
                "QUALITY_MIN": 0.61,
                "DIHEDRAL_FORBID_MIN": -42.96,
                "DIHEDRAL_FORBID_MAX": 124.95
            },
            "tHW": {
                "DIST_MIN": 2.73,
                "DIST_MAX": 3.35,
                "ANGLE_MIN": 89.3,
                "QUALITY_MIN": 0.62,
                "DIHEDRAL_FORBID_MIN": -51.95,
                "DIHEDRAL_FORBID_MAX": 183.68
            },
            "tSH": {
                "DIST_MIN": 2.8,
                "DIST_MAX": 3.39,
                "ANGLE_MIN": 86.81,
                "QUALITY_MIN": 0.62,
                "DIHEDRAL_FORBID_MIN": -52.81,
                "DIHEDRAL_FORBID_MAX": 137.35
            },
            "tSS": {
                "DIST_MIN": 2.72,
                "DIST_MAX": 3.39,
                "ANGLE_MIN": 92.32,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -101.93,
                "DIHEDRAL_FORBID_MAX": 138.2
            },
            "tSW": {
                "DIST_MIN": 2.78,
                "DIST_MAX": 3.37,
                "ANGLE_MIN": 97.0,
                "QUALITY_MIN": 0.61,
                "DIHEDRAL_FORBID_MIN": -109.82,
                "DIHEDRAL_FORBID_MAX": 58.73
            },
            "tW.": {
                "DIST_MIN": 2.8,
                "DIST_MAX": 3.36,
                "ANGLE_MIN": 102.15,
                "QUALITY_MIN": 0.57,
                "DIHEDRAL_FORBID_MIN": -121.53,
                "DIHEDRAL_FORBID_MAX": 140.36
            },
            "tWH": {
                "DIST_MIN": 2.76,
                "DIST_MAX": 3.31,
                "ANGLE_MIN": 96.67,
                "QUALITY_MIN": 0.59,
                "DIHEDRAL_FORBID_MIN": -133.2,
                "DIHEDRAL_FORBID_MAX": 172.06
            },
            "tWS": {
                "DIST_MIN": 2.73,
                "DIST_MAX": 3.38,
                "ANGLE_MIN": 98.78,
                "QUALITY_MIN": 0.62,
                "DIHEDRAL_FORBID_MIN": -79.99,
                "DIHEDRAL_FORBID_MAX": 129.55
            },
            "tWW": {
                "DIST_MIN": 2.74,
                "DIST_MAX": 3.24,
                "ANGLE_MIN": 97.2,
                "QUALITY_MIN": 0.71,
                "DIHEDRAL_FORBID_MIN": -68.64,
                "DIHEDRAL_FORBID_MAX": 144.18
            },
            "_OTHER": {
                "DIST_MIN": 2.73,
                "DIST_MAX": 3.41,
                "ANGLE_MIN": 98.52,
                "QUALITY_MIN": 0.54,
                "DIHEDRAL_FORBID_MIN": -113.27,
                "DIHEDRAL_FORBID_MAX": 159.03
            }
        },
        "C-U": {
            "c.H": {
                "DIST_MIN": 2.7,
                "DIST_MAX": 3.31,
                "ANGLE_MIN": 101.76,
                "QUALITY_MIN": 0.58,
                "DIHEDRAL_FORBID_MIN": -123.33,
                "DIHEDRAL_FORBID_MAX": 159.93
            },
            "c.W": {
                "DIST_MIN": 2.69,
                "DIST_MAX": 3.25,
                "ANGLE_MIN": 101.46,
                "QUALITY_MIN": 0.57,
                "DIHEDRAL_FORBID_MIN": -193.08,
                "DIHEDRAL_FORBID_MAX": 92.95
            },
            "cHH": {
                "DIST_MIN": 2.57,
                "DIST_MAX": 3.26,
                "ANGLE_MIN": 97.52,
                "QUALITY_MIN": 0.58,
                "DIHEDRAL_FORBID_MIN": -126.38,
                "DIHEDRAL_FORBID_MAX": 30.25
            },
            "cHS": {
                "DIST_MIN": 2.79,
                "DIST_MAX": 3.5,
                "ANGLE_MIN": 69.86,
                "QUALITY_MIN": 0.54,
                "DIHEDRAL_FORBID_MIN": -86.51,
                "DIHEDRAL_FORBID_MAX": 76.18
            },
            "cHW": {
                "DIST_MIN": 2.65,
                "DIST_MAX": 3.26,
                "ANGLE_MIN": 98.16,
                "QUALITY_MIN": 0.66,
                "DIHEDRAL_FORBID_MIN": -34.34,
                "DIHEDRAL_FORBID_MAX": 162.75
            },
            "cSH": {
                "DIST_MIN": 2.77,
                "DIST_MAX": 3.4,
                "ANGLE_MIN": 83.26,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -127.59,
                "DIHEDRAL_FORBID_MAX": 95.49
            },
            "cSW": {
                "DIST_MIN": 2.71,
                "DIST_MAX": 3.36,
                "ANGLE_MIN": 90.43,
                "QUALITY_MIN": 0.6,
                "DIHEDRAL_FORBID_MIN": -92.22,
                "DIHEDRAL_FORBID_MAX": 99.62
            },
            "cW.": {
                "DIST_MIN": 2.78,
                "DIST_MAX": 3.4,
                "ANGLE_MIN": 98.48,
                "QUALITY_MIN": 0.57,
                "DIHEDRAL_FORBID_MIN": -87.35,
                "DIHEDRAL_FORBID_MAX": 174.54
            },
            "cWH": {
                "DIST_MIN": 2.85,
                "DIST_MAX": 3.46,
                "ANGLE_MIN": 77.38,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -86.58,
                "DIHEDRAL_FORBID_MAX": 126.39
            },
            "cWS": {
                "DIST_MIN": 2.71,
                "DIST_MAX": 3.28,
                "ANGLE_MIN": 100.79,
                "QUALITY_MIN": 0.65,
                "DIHEDRAL_FORBID_MIN": -89.92,
                "DIHEDRAL_FORBID_MAX": 65.99
            },
            "cWW": {
                "DIST_MIN": 2.79,
                "DIST_MAX": 3.36,
                "ANGLE_MIN": 98.08,
                "QUALITY_MIN": 0.64,
                "DIHEDRAL_FORBID_MIN": -82.73,
                "DIHEDRAL_FORBID_MAX": 131.41
            },
            "t.W": {
                "DIST_MIN": 2.71,
                "DIST_MAX": 3.31,
                "ANGLE_MIN": 96.64,
                "QUALITY_MIN": 0.61,
                "DIHEDRAL_FORBID_MIN": -81.53,
                "DIHEDRAL_FORBID_MAX": 190.39
            },
            "tHS": {
                "DIST_MIN": 2.7,
                "DIST_MAX": 3.3,
                "ANGLE_MIN": 95.73,
                "QUALITY_MIN": 0.65,
                "DIHEDRAL_FORBID_MIN": -78.59,
                "DIHEDRAL_FORBID_MAX": 79.41
            },
            "tHW": {
                "DIST_MIN": 2.77,
                "DIST_MAX": 3.39,
                "ANGLE_MIN": 86.19,
                "QUALITY_MIN": 0.57,
                "DIHEDRAL_FORBID_MIN": -114.55,
                "DIHEDRAL_FORBID_MAX": 135.7
            },
            "tSH": {
                "DIST_MIN": 2.66,
                "DIST_MAX": 3.33,
                "ANGLE_MIN": 93.27,
                "QUALITY_MIN": 0.58,
                "DIHEDRAL_FORBID_MIN": -115.8,
                "DIHEDRAL_FORBID_MAX": 142.6
            },
            "tSW": {
                "DIST_MIN": 2.63,
                "DIST_MAX": 3.35,
                "ANGLE_MIN": 95.59,
                "QUALITY_MIN": 0.57,
                "DIHEDRAL_FORBID_MIN": -97.75,
                "DIHEDRAL_FORBID_MAX": 104.56
            },
            "tW.": {
                "DIST_MIN": 2.8,
                "DIST_MAX": 3.39,
                "ANGLE_MIN": 103.56,
                "QUALITY_MIN": 0.58,
                "DIHEDRAL_FORBID_MIN": -77.23,
                "DIHEDRAL_FORBID_MAX": 170.31
            },
            "tWH": {
                "DIST_MIN": 2.75,
                "DIST_MAX": 3.35,
                "ANGLE_MIN": 95.59,
                "QUALITY_MIN": 0.59,
                "DIHEDRAL_FORBID_MIN": -132.5,
                "DIHEDRAL_FORBID_MAX": 109.45
            },
            "tWW": {
                "DIST_MIN": 2.73,
                "DIST_MAX": 3.35,
                "ANGLE_MIN": 87.94,
                "QUALITY_MIN": 0.61,
                "DIHEDRAL_FORBID_MIN": -81.16,
                "DIHEDRAL_FORBID_MAX": 75.66
            },
            "_OTHER": {
                "DIST_MIN": 2.65,
                "DIST_MAX": 3.32,
                "ANGLE_MIN": 95.48,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -135.62,
                "DIHEDRAL_FORBID_MAX": 132.15
            }
        },
        "G-U": {
            "c.W": {
                "DIST_MIN": 2.71,
                "DIST_MAX": 3.3,
                "ANGLE_MIN": 98.2,
                "QUALITY_MIN": 0.59,
                "DIHEDRAL_FORBID_MIN": -148.29,
                "DIHEDRAL_FORBID_MAX": 124.11
            },
            "cHS": {
                "DIST_MIN": 2.77,
                "DIST_MAX": 3.48,
                "ANGLE_MIN": 84.16,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -75.1,
                "DIHEDRAL_FORBID_MAX": 98.79
            },
            "cHW": {
                "DIST_MIN": 2.69,
                "DIST_MAX": 3.34,
                "ANGLE_MIN": 80.67,
                "QUALITY_MIN": 0.62,
                "DIHEDRAL_FORBID_MIN": -77.75,
                "DIHEDRAL_FORBID_MAX": 114.12
            },
            "cSH": {
                "DIST_MIN": 2.73,
                "DIST_MAX": 3.43,
                "ANGLE_MIN": 78.47,
                "QUALITY_MIN": 0.6,
                "DIHEDRAL_FORBID_MIN": -62.37,
                "DIHEDRAL_FORBID_MAX": 54.7
            },
            "cSS": {
                "DIST_MIN": 2.7,
                "DIST_MAX": 3.23,
                "ANGLE_MIN": 100.82,
                "QUALITY_MIN": 0.69,
                "DIHEDRAL_FORBID_MIN": -53.39,
                "DIHEDRAL_FORBID_MAX": 136.77
            },
            "cSW": {
                "DIST_MIN": 2.76,
                "DIST_MAX": 3.43,
                "ANGLE_MIN": 87.83,
                "QUALITY_MIN": 0.57,
                "DIHEDRAL_FORBID_MIN": -95.96,
                "DIHEDRAL_FORBID_MAX": 100.76
            },
            "cW.": {
                "DIST_MIN": 2.73,
                "DIST_MAX": 3.35,
                "ANGLE_MIN": 94.49,
                "QUALITY_MIN": 0.54,
                "DIHEDRAL_FORBID_MIN": -142.66,
                "DIHEDRAL_FORBID_MAX": 113.2
            },
            "cWH": {
                "DIST_MIN": 2.83,
                "DIST_MAX": 3.48,
                "ANGLE_MIN": 75.82,
                "QUALITY_MIN": 0.52,
                "DIHEDRAL_FORBID_MIN": -128.85,
                "DIHEDRAL_FORBID_MAX": 41.12
            },
            "cWS": {
                "DIST_MIN": 2.81,
                "DIST_MAX": 3.42,
                "ANGLE_MIN": 93.57,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -112.41,
                "DIHEDRAL_FORBID_MAX": 85.02
            },
            "cWW": {
                "DIST_MIN": 2.71,
                "DIST_MAX": 3.16,
                "ANGLE_MIN": 97.29,
                "QUALITY_MIN": 0.77,
                "DIHEDRAL_FORBID_MIN": -51.35,
                "DIHEDRAL_FORBID_MAX": 149.51
            },
            "t.H": {
                "DIST_MIN": 2.59,
                "DIST_MAX": 3.28,
                "ANGLE_MIN": 102.39,
                "QUALITY_MIN": 0.59,
                "DIHEDRAL_FORBID_MIN": -86.23,
                "DIHEDRAL_FORBID_MAX": 203.12
            },
            "t.S": {
                "DIST_MIN": 2.65,
                "DIST_MAX": 3.21,
                "ANGLE_MIN": 102.97,
                "QUALITY_MIN": 0.69,
                "DIHEDRAL_FORBID_MIN": -7.54,
                "DIHEDRAL_FORBID_MAX": 228.04
            },
            "t.W": {
                "DIST_MIN": 2.62,
                "DIST_MAX": 3.22,
                "ANGLE_MIN": 96.41,
                "QUALITY_MIN": 0.68,
                "DIHEDRAL_FORBID_MIN": -65.61,
                "DIHEDRAL_FORBID_MAX": 125.05
            },
            "tHH": {
                "DIST_MIN": 2.57,
                "DIST_MAX": 3.19,
                "ANGLE_MIN": 94.2,
                "QUALITY_MIN": 0.64,
                "DIHEDRAL_FORBID_MIN": -43.07,
                "DIHEDRAL_FORBID_MAX": 197.75
            },
            "tHS": {
                "DIST_MIN": 2.63,
                "DIST_MAX": 3.35,
                "ANGLE_MIN": 89.54,
                "QUALITY_MIN": 0.57,
                "DIHEDRAL_FORBID_MIN": -91.4,
                "DIHEDRAL_FORBID_MAX": 162.55
            },
            "tHW": {
                "DIST_MIN": 2.73,
                "DIST_MAX": 3.32,
                "ANGLE_MIN": 88.47,
                "QUALITY_MIN": 0.6,
                "DIHEDRAL_FORBID_MIN": -99.07,
                "DIHEDRAL_FORBID_MAX": 150.63
            },
            "tS.": {
                "DIST_MIN": 2.62,
                "DIST_MAX": 3.25,
                "ANGLE_MIN": 102.61,
                "QUALITY_MIN": 0.6,
                "DIHEDRAL_FORBID_MIN": -129.13,
                "DIHEDRAL_FORBID_MAX": 187.98
            },
            "tSH": {
                "DIST_MIN": 2.58,
                "DIST_MAX": 3.23,
                "ANGLE_MIN": 91.57,
                "QUALITY_MIN": 0.65,
                "DIHEDRAL_FORBID_MIN": -69.61,
                "DIHEDRAL_FORBID_MAX": 136.23
            },
            "tSS": {
                "DIST_MIN": 2.77,
                "DIST_MAX": 3.39,
                "ANGLE_MIN": 92.35,
                "QUALITY_MIN": 0.59,
                "DIHEDRAL_FORBID_MIN": -111.3,
                "DIHEDRAL_FORBID_MAX": 121.55
            },
            "tSW": {
                "DIST_MIN": 2.66,
                "DIST_MAX": 3.35,
                "ANGLE_MIN": 91.78,
                "QUALITY_MIN": 0.61,
                "DIHEDRAL_FORBID_MIN": -88.38,
                "DIHEDRAL_FORBID_MAX": 108.31
            },
            "tW.": {
                "DIST_MIN": 2.73,
                "DIST_MAX": 3.3,
                "ANGLE_MIN": 94.92,
                "QUALITY_MIN": 0.63,
                "DIHEDRAL_FORBID_MIN": -87.01,
                "DIHEDRAL_FORBID_MAX": 175.95
            },
            "tWH": {
                "DIST_MIN": 2.67,
                "DIST_MAX": 3.37,
                "ANGLE_MIN": 88.78,
                "QUALITY_MIN": 0.56,
                "DIHEDRAL_FORBID_MIN": -105.8,
                "DIHEDRAL_FORBID_MAX": 147.78
            },
            "tWS": {
                "DIST_MIN": 2.75,
                "DIST_MAX": 3.36,
                "ANGLE_MIN": 94.39,
                "QUALITY_MIN": 0.6,
                "DIHEDRAL_FORBID_MIN": -82.63,
                "DIHEDRAL_FORBID_MAX": 132.11
            },
            "tWW": {
                "DIST_MIN": 2.78,
                "DIST_MAX": 3.35,
                "ANGLE_MIN": 89.21,
                "QUALITY_MIN": 0.64,
                "DIHEDRAL_FORBID_MIN": -70.26,
                "DIHEDRAL_FORBID_MAX": 151.84
            },
            "_OTHER": {
                "DIST_MIN": 2.75,
                "DIST_MAX": 3.44,
                "ANGLE_MIN": 84.96,
                "QUALITY_MIN": 0.54,
                "DIHEDRAL_FORBID_MIN": -111.02,
                "DIHEDRAL_FORBID_MAX": 139.85
            }
        }
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
    
    
    
    
    
    
    
    
    
    
    
    # ===== QUALITY GRADE THRESHOLDS =====
    GRADE_EXCELLENT = 85
    GRADE_GOOD = 70
    GRADE_FAIR = 50
    GRADE_POOR = 0
    
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
