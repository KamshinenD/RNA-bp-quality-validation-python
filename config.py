"""Configuration file with all tunable thresholds and scoring parameters."""

class Config:
    """RNA quality scoring configuration."""
    
    BASE_SCORE = 100.0
# }
    
    PENALTY_WEIGHTS = {
        # CRITICAL - No hydrogen bonds detected (DSSR + CSV both report 0)
        'zero_hbond_pairs': 20.0,

        # SEVERE - H-bond distance outside edge-specific range (DIST_MIN to DIST_MAX)
        'bad_hbond_distance': 18.0,

        # MAJOR - In-plane translational displacement (shear or stretch outside thresholds)
        'misaligned_pairs': 14.0,

        # MAJOR - Out-of-plane issues (buckle rotation or stagger displacement)
        'non_coplanar_pairs': 12.0,

        # MODERATE - H-bond quality score below edge-specific QUALITY_MIN threshold
        'weak_hbond_quality': 10.0,

        # MODERATE - H-bond dihedral in forbidden zone (not cis: -50 to 50, not trans: |angle| >= 140)
        'bad_hbond_dihedrals': 8.0,

        # MODERATE - H-bond count outside expected range for base pair type (or != ideal for cWW)
        'incorrect_hbond_count': 8.0,

        # MINOR - H-bond donor/acceptor angle below edge-specific ANGLE_MIN threshold
        'bad_hbond_angles': 5.0,

        # MINOR - Rotational distortion (propeller or opening outside thresholds)
        'rotational_distortion_pairs': 5.0,
    }
    # Total: 100.0
    
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
    
    HBOND_SCORE_MIN = 2.0        # Mean = 2.007 (flags pairs below average H-bond quality)
    
    # Base-pair + edge thresholds (mean ± SD, 2dp) generated from edge_type_analysis.json
    GEOMETRY_THRESHOLDS_BY_BASE_PAIR_BY_EDGE = {
        "A-A": {
            "cSH": {"SHEAR_MAX": 7.15, "STRETCH_MIN": -1.13, "STRETCH_MAX": 0.99, "STAGGER_MAX": 1.45, "BUCKLE_MAX": 20.85, "PROPELLER_MIN": -10.0, "PROPELLER_MAX": 25.2, "OPENING_MIN": -11.69, "OPENING_MAX": 11.95},
            "cSW": {"SHEAR_MAX": 7.24, "STRETCH_MIN": -2.02, "STRETCH_MAX": 2.96, "STAGGER_MAX": 0.86, "BUCKLE_MAX": 30.43, "PROPELLER_MIN": -25.29, "PROPELLER_MAX": 3.43, "OPENING_MIN": 33.86, "OPENING_MAX": 101.05},
            "cWS": {"SHEAR_MAX": 7.16, "STRETCH_MIN": -2.06, "STRETCH_MAX": 3.72, "STAGGER_MAX": 1.26, "BUCKLE_MAX": 23.35, "PROPELLER_MIN": -18.26, "PROPELLER_MAX": 14.26, "OPENING_MIN": 32.57, "OPENING_MAX": 115.01},
            "cWW": {"SHEAR_MAX": 2.25, "STRETCH_MIN": 0.75, "STRETCH_MAX": 2.21, "STAGGER_MAX": 1.25, "BUCKLE_MAX": 37.29, "PROPELLER_MIN": -21.26, "PROPELLER_MAX": 5.39, "OPENING_MIN": -77.53, "OPENING_MAX": 42.63},
            "tHH": {"SHEAR_MAX": 6.31, "STRETCH_MIN": -4.85, "STRETCH_MAX": 5.5, "STAGGER_MAX": 0.66, "BUCKLE_MAX": 20.93, "PROPELLER_MIN": -15.57, "PROPELLER_MAX": 15.29, "OPENING_MIN": -184.34, "OPENING_MAX": 159.88},
            "tHS": {"SHEAR_MAX": 6.9, "STRETCH_MIN": -5.4, "STRETCH_MAX": -2.71, "STAGGER_MAX": 0.73, "BUCKLE_MAX": 11.74, "PROPELLER_MIN": -18.55, "PROPELLER_MAX": 11.5, "OPENING_MIN": -33.41, "OPENING_MAX": -0.99},
            "tHW": {"SHEAR_MAX": 5.51, "STRETCH_MIN": -0.76, "STRETCH_MAX": 2.36, "STAGGER_MAX": 1.1, "BUCKLE_MAX": 15.92, "PROPELLER_MIN": -22.64, "PROPELLER_MAX": 19.24, "OPENING_MIN": -130.35, "OPENING_MAX": -65.64},
            "tSH": {"SHEAR_MAX": 6.89, "STRETCH_MIN": -5.37, "STRETCH_MAX": -3.16, "STAGGER_MAX": 0.87, "BUCKLE_MAX": 21.6, "PROPELLER_MIN": -20.96, "PROPELLER_MAX": 10.97, "OPENING_MIN": -34.76, "OPENING_MAX": -4.62},
            "tSW": {"SHEAR_MAX": 3.29, "STRETCH_MIN": -7.76, "STRETCH_MAX": -3.39, "STAGGER_MAX": 0.99, "BUCKLE_MAX": 23.54, "PROPELLER_MIN": -8.39, "PROPELLER_MAX": 23.12, "OPENING_MIN": -150.0, "OPENING_MAX": -58.0},
            "tWH": {"SHEAR_MAX": 5.19, "STRETCH_MIN": -0.04, "STRETCH_MAX": 2.15, "STAGGER_MAX": 1.2, "BUCKLE_MAX": 18.37, "PROPELLER_MIN": -21.06, "PROPELLER_MAX": 24.82, "OPENING_MIN": -118.3, "OPENING_MAX": -88.86},
            "tWW": {"SHEAR_MAX": 1.76, "STRETCH_MIN": -1.48, "STRETCH_MAX": 1.39, "STAGGER_MAX": 1.07, "BUCKLE_MAX": 37.46, "PROPELLER_MIN": -15.91, "PROPELLER_MAX": 20.7, "OPENING_MIN": -170.78, "OPENING_MAX": 167.13}
        },
        "A-C": {
            "_OTHER": {"SHEAR_MAX": 5.63, "STRETCH_MIN": -4.03, "STRETCH_MAX": 1.19, "STAGGER_MAX": 1.32, "BUCKLE_MAX": 25.58, "PROPELLER_MIN": -21.45, "PROPELLER_MAX": 15.58, "OPENING_MIN": -75.27, "OPENING_MAX": 50.43},
            "cSH": {"SHEAR_MAX": 7.03, "STRETCH_MIN": -1.44, "STRETCH_MAX": 0.83, "STAGGER_MAX": 1.68, "BUCKLE_MAX": 34.6, "PROPELLER_MIN": -7.69, "PROPELLER_MAX": 28.85, "OPENING_MIN": -35.95, "OPENING_MAX": 5.6},
            "cSS": {"SHEAR_MAX": 2.76, "STRETCH_MIN": 3.23, "STRETCH_MAX": 8.96, "STAGGER_MAX": 1.03, "BUCKLE_MAX": 21.32, "PROPELLER_MIN": -26.75, "PROPELLER_MAX": -1.68, "OPENING_MIN": 71.23, "OPENING_MAX": 204.54},
            "cSW": {"SHEAR_MAX": 7.18, "STRETCH_MIN": -2.48, "STRETCH_MAX": 3.15, "STAGGER_MAX": 1.14, "BUCKLE_MAX": 25.38, "PROPELLER_MIN": -26.02, "PROPELLER_MAX": 2.58, "OPENING_MIN": 32.25, "OPENING_MAX": 106.26},
            "cWH": {"SHEAR_MAX": 4.57, "STRETCH_MIN": -0.41, "STRETCH_MAX": 2.82, "STAGGER_MAX": 1.89, "BUCKLE_MAX": 55.67, "PROPELLER_MIN": -31.86, "PROPELLER_MAX": 10.77, "OPENING_MIN": -94.83, "OPENING_MAX": -45.14},
            "cWS": {"SHEAR_MAX": 7.14, "STRETCH_MIN": -2.29, "STRETCH_MAX": 3.73, "STAGGER_MAX": 1.47, "BUCKLE_MAX": 35.82, "PROPELLER_MIN": -31.65, "PROPELLER_MAX": 4.95, "OPENING_MIN": 38.4, "OPENING_MAX": 114.21},
            "cWW": {"SHEAR_MAX": 3.36, "STRETCH_MIN": -0.83, "STRETCH_MAX": 0.41, "STAGGER_MAX": 0.93, "BUCKLE_MAX": 18.73, "PROPELLER_MIN": -20.97, "PROPELLER_MAX": 9.16, "OPENING_MIN": -12.99, "OPENING_MAX": 20.29},
            "tHW": {"SHEAR_MAX": 5.24, "STRETCH_MIN": -2.22, "STRETCH_MAX": 0.8, "STAGGER_MAX": 1.32, "BUCKLE_MAX": 21.86, "PROPELLER_MIN": -22.38, "PROPELLER_MAX": 20.13, "OPENING_MIN": -107.66, "OPENING_MAX": -51.96},
            "tSW": {"SHEAR_MAX": 3.49, "STRETCH_MIN": -8.51, "STRETCH_MAX": -4.87, "STAGGER_MAX": 1.16, "BUCKLE_MAX": 32.95, "PROPELLER_MIN": -20.14, "PROPELLER_MAX": 15.15, "OPENING_MIN": -167.77, "OPENING_MAX": -61.31},
            "tWH": {"SHEAR_MAX": 4.62, "STRETCH_MIN": -1.69, "STRETCH_MAX": 1.14, "STAGGER_MAX": 1.16, "BUCKLE_MAX": 17.79, "PROPELLER_MIN": -20.53, "PROPELLER_MAX": 15.88, "OPENING_MIN": -114.57, "OPENING_MAX": -63.97},
            "tWS": {"SHEAR_MAX": 3.9, "STRETCH_MIN": 1.54, "STRETCH_MAX": 9.16, "STAGGER_MAX": 1.55, "BUCKLE_MAX": 22.79, "PROPELLER_MIN": -35.71, "PROPELLER_MAX": 23.76, "OPENING_MIN": 70.28, "OPENING_MAX": 147.86},
            "tWW": {"SHEAR_MAX": 1.27, "STRETCH_MIN": -2.35, "STRETCH_MAX": 1.33, "STAGGER_MAX": 0.96, "BUCKLE_MAX": 26.25, "PROPELLER_MIN": -17.66, "PROPELLER_MAX": 20.32, "OPENING_MIN": -149.89, "OPENING_MAX": 156.84}
        },
        "A-G": {
            "_OTHER": {"SHEAR_MAX": 4.08, "STRETCH_MIN": -5.39, "STRETCH_MAX": 5.43, "STAGGER_MAX": 1.78, "BUCKLE_MAX": 35.73, "PROPELLER_MIN": -15.35, "PROPELLER_MAX": 20.17, "OPENING_MIN": -144.99, "OPENING_MAX": 102.07},
            "cHW": {"SHEAR_MAX": 1.98, "STRETCH_MIN": -6.18, "STRETCH_MAX": -2.89, "STAGGER_MAX": 0.84, "BUCKLE_MAX": 16.86, "PROPELLER_MIN": -29.79, "PROPELLER_MAX": 14.2, "OPENING_MIN": 39.41, "OPENING_MAX": 123.4},
            "cSH": {"SHEAR_MAX": 7.5, "STRETCH_MIN": -2.2, "STRETCH_MAX": 2.48, "STAGGER_MAX": 1.61, "BUCKLE_MAX": 29.75, "PROPELLER_MIN": -15.25, "PROPELLER_MAX": 30.89, "OPENING_MIN": -15.75, "OPENING_MAX": 26.37},
            "cSS": {"SHEAR_MAX": 2.3, "STRETCH_MIN": 5.89, "STRETCH_MAX": 8.55, "STAGGER_MAX": 1.54, "BUCKLE_MAX": 29.33, "PROPELLER_MIN": -28.18, "PROPELLER_MAX": 10.11, "OPENING_MIN": 113.53, "OPENING_MAX": 172.51},
            "cSW": {"SHEAR_MAX": 4.89, "STRETCH_MIN": 1.38, "STRETCH_MAX": 6.16, "STAGGER_MAX": 1.52, "BUCKLE_MAX": 24.87, "PROPELLER_MIN": -41.57, "PROPELLER_MAX": 13.47, "OPENING_MIN": 68.83, "OPENING_MAX": 138.5},
            "cWH": {"SHEAR_MAX": 2.08, "STRETCH_MIN": 3.23, "STRETCH_MAX": 5.9, "STAGGER_MAX": 1.38, "BUCKLE_MAX": 34.55, "PROPELLER_MIN": -19.3, "PROPELLER_MAX": 19.15, "OPENING_MIN": -124.84, "OPENING_MAX": -73.3},
            "cWS": {"SHEAR_MAX": 5.37, "STRETCH_MIN": 0.81, "STRETCH_MAX": 5.99, "STAGGER_MAX": 1.74, "BUCKLE_MAX": 46.2, "PROPELLER_MIN": -33.9, "PROPELLER_MAX": 14.26, "OPENING_MIN": 67.39, "OPENING_MAX": 140.52},
            "cWW": {"SHEAR_MAX": 1.05, "STRETCH_MIN": 1.11, "STRETCH_MAX": 2.01, "STAGGER_MAX": 1.08, "BUCKLE_MAX": 16.89, "PROPELLER_MIN": -20.23, "PROPELLER_MAX": 12.31, "OPENING_MIN": -32.84, "OPENING_MAX": 2.32},
            "tHS": {"SHEAR_MAX": 7.08, "STRETCH_MIN": -4.92, "STRETCH_MAX": -3.56, "STAGGER_MAX": 0.79, "BUCKLE_MAX": 16.56, "PROPELLER_MIN": -16.77, "PROPELLER_MAX": 10.44, "OPENING_MIN": -14.99, "OPENING_MAX": 5.44},
            "tSH": {"SHEAR_MAX": 7.42, "STRETCH_MIN": -5.3, "STRETCH_MAX": -3.81, "STAGGER_MAX": 1.08, "BUCKLE_MAX": 21.57, "PROPELLER_MIN": -14.95, "PROPELLER_MAX": 12.82, "OPENING_MIN": -28.13, "OPENING_MAX": 4.02},
            "tSS": {"SHEAR_MAX": 3.04, "STRETCH_MIN": -8.58, "STRETCH_MAX": 6.22, "STAGGER_MAX": 1.57, "BUCKLE_MAX": 31.53, "PROPELLER_MIN": -24.33, "PROPELLER_MAX": 21.22, "OPENING_MIN": -175.39, "OPENING_MAX": 128.34},
            "tSW": {"SHEAR_MAX": 3.69, "STRETCH_MIN": -7.3, "STRETCH_MAX": -3.8, "STAGGER_MAX": 1.1, "BUCKLE_MAX": 31.15, "PROPELLER_MIN": -8.55, "PROPELLER_MAX": 36.3, "OPENING_MIN": -142.13, "OPENING_MAX": -69.61},
            "tWH": {"SHEAR_MAX": 8.9, "STRETCH_MIN": -4.84, "STRETCH_MAX": -0.28, "STAGGER_MAX": 1.17, "BUCKLE_MAX": 25.81, "PROPELLER_MIN": -20.52, "PROPELLER_MAX": 25.09, "OPENING_MIN": -102.5, "OPENING_MAX": -49.44},
            "tWS": {"SHEAR_MAX": 4.03, "STRETCH_MIN": 3.42, "STRETCH_MAX": 6.96, "STAGGER_MAX": 1.04, "BUCKLE_MAX": 25.26, "PROPELLER_MIN": -32.8, "PROPELLER_MAX": 3.46, "OPENING_MIN": 64.09, "OPENING_MAX": 133.69}
        },
        "A-U": {
            "_OTHER": {"SHEAR_MAX": 5.72, "STRETCH_MIN": -4.08, "STRETCH_MAX": 3.15, "STAGGER_MAX": 1.26, "BUCKLE_MAX": 20.37, "PROPELLER_MIN": -21.94, "PROPELLER_MAX": 16.55, "OPENING_MIN": -113.22, "OPENING_MAX": 100.64},
            "cHW": {"SHEAR_MAX": 1.37, "STRETCH_MIN": -4.82, "STRETCH_MAX": -2.47, "STAGGER_MAX": 1.03, "BUCKLE_MAX": 20.52, "PROPELLER_MIN": -9.55, "PROPELLER_MAX": 16.46, "OPENING_MIN": 41.83, "OPENING_MAX": 97.65},
            "cSH": {"SHEAR_MAX": 6.95, "STRETCH_MIN": -1.53, "STRETCH_MAX": 0.9, "STAGGER_MAX": 1.65, "BUCKLE_MAX": 34.0, "PROPELLER_MIN": -17.94, "PROPELLER_MAX": 30.62, "OPENING_MIN": -27.71, "OPENING_MAX": 11.98},
            "cSW": {"SHEAR_MAX": 6.72, "STRETCH_MIN": -1.71, "STRETCH_MAX": 2.8, "STAGGER_MAX": 1.64, "BUCKLE_MAX": 30.09, "PROPELLER_MIN": -22.45, "PROPELLER_MAX": 27.97, "OPENING_MIN": 37.1, "OPENING_MAX": 108.37},
            "cWH": {"SHEAR_MAX": 1.51, "STRETCH_MIN": 2.0, "STRETCH_MAX": 4.4, "STAGGER_MAX": 1.18, "BUCKLE_MAX": 24.19, "PROPELLER_MIN": -20.1, "PROPELLER_MAX": 9.58, "OPENING_MIN": -83.65, "OPENING_MAX": -47.71},
            "cWS": {"SHEAR_MAX": 6.24, "STRETCH_MIN": -1.2, "STRETCH_MAX": 3.61, "STAGGER_MAX": 1.16, "BUCKLE_MAX": 27.4, "PROPELLER_MIN": -24.92, "PROPELLER_MAX": 22.56, "OPENING_MIN": 52.48, "OPENING_MAX": 116.89},
            "cWW": {"SHEAR_MAX": 0.71, "STRETCH_MIN": -0.33, "STRETCH_MAX": 0.19, "STAGGER_MAX": 0.55, "BUCKLE_MAX": 9.79, "PROPELLER_MIN": -15.79, "PROPELLER_MAX": 3.6, "OPENING_MIN": -6.88, "OPENING_MAX": 8.16},
            "tHW": {"SHEAR_MAX": 4.91, "STRETCH_MIN": -2.71, "STRETCH_MAX": -1.13, "STAGGER_MAX": 1.1, "BUCKLE_MAX": 12.25, "PROPELLER_MIN": -19.13, "PROPELLER_MAX": 6.94, "OPENING_MIN": -111.53, "OPENING_MAX": -67.58},
            "tSH": {"SHEAR_MAX": 6.87, "STRETCH_MIN": -9.02, "STRETCH_MAX": -5.04, "STAGGER_MAX": 1.47, "BUCKLE_MAX": 22.02, "PROPELLER_MIN": -20.46, "PROPELLER_MAX": 21.4, "OPENING_MIN": -77.79, "OPENING_MAX": -12.85},
            "tSW": {"SHEAR_MAX": 2.33, "STRETCH_MIN": -8.0, "STRETCH_MAX": -3.95, "STAGGER_MAX": 0.9, "BUCKLE_MAX": 21.61, "PROPELLER_MIN": -10.9, "PROPELLER_MAX": 21.98, "OPENING_MIN": -153.68, "OPENING_MAX": -74.04},
            "tWH": {"SHEAR_MAX": 5.04, "STRETCH_MIN": -2.9, "STRETCH_MAX": -1.1, "STAGGER_MAX": 0.92, "BUCKLE_MAX": 13.13, "PROPELLER_MIN": -18.52, "PROPELLER_MAX": 11.1, "OPENING_MIN": -110.92, "OPENING_MAX": -66.99},
            "tWW": {"SHEAR_MAX": 0.79, "STRETCH_MIN": -1.5, "STRETCH_MAX": 1.74, "STAGGER_MAX": 0.78, "BUCKLE_MAX": 17.55, "PROPELLER_MIN": -13.85, "PROPELLER_MAX": 17.42, "OPENING_MIN": -136.28, "OPENING_MAX": 183.84}
        },
        "C-C": {
            "cWW": {"SHEAR_MAX": 3.51, "STRETCH_MIN": -1.94, "STRETCH_MAX": -0.71, "STAGGER_MAX": 1.19, "BUCKLE_MAX": 16.47, "PROPELLER_MIN": -19.28, "PROPELLER_MAX": 10.29, "OPENING_MIN": -15.21, "OPENING_MAX": 29.74}
        },
        "C-G": {
            "_OTHER": {"SHEAR_MAX": 3.49, "STRETCH_MIN": -1.53, "STRETCH_MAX": 5.08, "STAGGER_MAX": 1.24, "BUCKLE_MAX": 24.37, "PROPELLER_MIN": -19.83, "PROPELLER_MAX": 15.3, "OPENING_MIN": -144.29, "OPENING_MAX": 50.5},
            "cSS": {"SHEAR_MAX": 2.74, "STRETCH_MIN": -1.47, "STRETCH_MAX": 10.57, "STAGGER_MAX": 1.59, "BUCKLE_MAX": 30.31, "PROPELLER_MIN": -30.75, "PROPELLER_MAX": 10.93, "OPENING_MIN": -30.92, "OPENING_MAX": 224.84},
            "cWH": {"SHEAR_MAX": 2.66, "STRETCH_MIN": 1.63, "STRETCH_MAX": 5.08, "STAGGER_MAX": 1.15, "BUCKLE_MAX": 36.28, "PROPELLER_MIN": -21.66, "PROPELLER_MAX": 18.58, "OPENING_MIN": -109.72, "OPENING_MAX": -46.42},
            "cWW": {"SHEAR_MAX": 0.51, "STRETCH_MIN": -0.37, "STRETCH_MAX": 0.11, "STAGGER_MAX": 0.61, "BUCKLE_MAX": 10.12, "PROPELLER_MIN": -15.72, "PROPELLER_MAX": 3.07, "OPENING_MIN": -5.4, "OPENING_MAX": 7.0},
            "tSW": {"SHEAR_MAX": 2.53, "STRETCH_MIN": -7.17, "STRETCH_MAX": -2.92, "STAGGER_MAX": 1.79, "BUCKLE_MAX": 42.79, "PROPELLER_MIN": -10.88, "PROPELLER_MAX": 36.36, "OPENING_MIN": -134.57, "OPENING_MAX": -65.92},
            "tWW": {"SHEAR_MAX": 0.62, "STRETCH_MIN": -3.67, "STRETCH_MAX": 3.49, "STAGGER_MAX": 1.06, "BUCKLE_MAX": 24.12, "PROPELLER_MIN": -21.92, "PROPELLER_MAX": 20.86, "OPENING_MIN": -163.3, "OPENING_MAX": 144.12}
        },
        "C-U": {
            "_OTHER": {"SHEAR_MAX": 4.7, "STRETCH_MIN": -3.4, "STRETCH_MAX": 2.25, "STAGGER_MAX": 1.61, "BUCKLE_MAX": 25.8, "PROPELLER_MIN": -15.91, "PROPELLER_MAX": 25.84, "OPENING_MIN": -125.51, "OPENING_MAX": 36.96},
            "cSH": {"SHEAR_MAX": 7.01, "STRETCH_MIN": -0.68, "STRETCH_MAX": 1.35, "STAGGER_MAX": 1.63, "BUCKLE_MAX": 29.43, "PROPELLER_MIN": -18.13, "PROPELLER_MAX": 20.75, "OPENING_MIN": -38.98, "OPENING_MAX": 1.74},
            "cSW": {"SHEAR_MAX": 6.02, "STRETCH_MIN": -1.8, "STRETCH_MAX": 1.95, "STAGGER_MAX": 1.54, "BUCKLE_MAX": 41.19, "PROPELLER_MIN": -25.48, "PROPELLER_MAX": 38.03, "OPENING_MIN": 38.86, "OPENING_MAX": 119.01},
            "cWW": {"SHEAR_MAX": 2.62, "STRETCH_MIN": -2.0, "STRETCH_MAX": -0.2, "STAGGER_MAX": 1.12, "BUCKLE_MAX": 20.46, "PROPELLER_MIN": -21.85, "PROPELLER_MAX": 9.03, "OPENING_MIN": -23.77, "OPENING_MAX": 28.8}
        },
        "G-G": {
            "_OTHER": {"SHEAR_MAX": 5.66, "STRETCH_MIN": -3.8, "STRETCH_MAX": 2.97, "STAGGER_MAX": 1.2, "BUCKLE_MAX": 26.62, "PROPELLER_MIN": -21.39, "PROPELLER_MAX": 15.96, "OPENING_MIN": -94.1, "OPENING_MAX": 48.57},
            "cHW": {"SHEAR_MAX": 3.13, "STRETCH_MIN": -4.43, "STRETCH_MAX": -1.72, "STAGGER_MAX": 0.91, "BUCKLE_MAX": 23.15, "PROPELLER_MIN": -18.3, "PROPELLER_MAX": 19.43, "OPENING_MIN": 34.84, "OPENING_MAX": 120.35},
            "cSH": {"SHEAR_MAX": 7.67, "STRETCH_MIN": -1.6, "STRETCH_MAX": 3.62, "STAGGER_MAX": 1.8, "BUCKLE_MAX": 29.49, "PROPELLER_MIN": -27.29, "PROPELLER_MAX": 28.95, "OPENING_MIN": -29.26, "OPENING_MAX": 36.69},
            "cWH": {"SHEAR_MAX": 3.23, "STRETCH_MIN": 2.33, "STRETCH_MAX": 4.17, "STAGGER_MAX": 0.97, "BUCKLE_MAX": 26.98, "PROPELLER_MIN": -11.2, "PROPELLER_MAX": 16.56, "OPENING_MIN": -102.78, "OPENING_MAX": -63.59},
            "cWW": {"SHEAR_MAX": 2.47, "STRETCH_MIN": 0.77, "STRETCH_MAX": 2.25, "STAGGER_MAX": 1.61, "BUCKLE_MAX": 23.45, "PROPELLER_MIN": -24.21, "PROPELLER_MAX": 19.22, "OPENING_MIN": -37.95, "OPENING_MAX": 10.32},
            "tHS": {"SHEAR_MAX": 7.85, "STRETCH_MIN": -4.82, "STRETCH_MAX": -1.79, "STAGGER_MAX": 1.01, "BUCKLE_MAX": 21.06, "PROPELLER_MIN": -17.85, "PROPELLER_MAX": 12.19, "OPENING_MIN": -37.74, "OPENING_MAX": 7.35},
            "tHW": {"SHEAR_MAX": 7.38, "STRETCH_MIN": -3.1, "STRETCH_MAX": 0.95, "STAGGER_MAX": 0.95, "BUCKLE_MAX": 23.63, "PROPELLER_MIN": -13.68, "PROPELLER_MAX": 30.55, "OPENING_MIN": -116.52, "OPENING_MAX": -47.89},
            "tSH": {"SHEAR_MAX": 7.88, "STRETCH_MIN": -5.48, "STRETCH_MAX": -2.3, "STAGGER_MAX": 1.05, "BUCKLE_MAX": 28.24, "PROPELLER_MIN": -17.86, "PROPELLER_MAX": 21.45, "OPENING_MIN": -53.07, "OPENING_MAX": -5.42},
            "tSS": {"SHEAR_MAX": 3.53, "STRETCH_MIN": -9.16, "STRETCH_MAX": 6.74, "STAGGER_MAX": 1.02, "BUCKLE_MAX": 28.2, "PROPELLER_MIN": -26.51, "PROPELLER_MAX": 19.1, "OPENING_MIN": -195.65, "OPENING_MAX": 144.94},
            "tWH": {"SHEAR_MAX": 6.88, "STRETCH_MIN": -1.99, "STRETCH_MAX": 1.08, "STAGGER_MAX": 0.99, "BUCKLE_MAX": 32.57, "PROPELLER_MIN": -11.67, "PROPELLER_MAX": 23.86, "OPENING_MIN": -109.36, "OPENING_MAX": -52.92}
        },
        "G-U": {
            "_OTHER": {"SHEAR_MAX": 6.47, "STRETCH_MIN": -3.58, "STRETCH_MAX": 1.54, "STAGGER_MAX": 1.17, "BUCKLE_MAX": 27.81, "PROPELLER_MIN": -21.27, "PROPELLER_MAX": 14.46, "OPENING_MIN": -81.84, "OPENING_MAX": 51.11},
            "cSH": {"SHEAR_MAX": 8.44, "STRETCH_MIN": 0.85, "STRETCH_MAX": 2.47, "STAGGER_MAX": 1.05, "BUCKLE_MAX": 18.8, "PROPELLER_MIN": -9.62, "PROPELLER_MAX": 23.65, "OPENING_MIN": -3.75, "OPENING_MAX": 26.26},
            "cSS": {"SHEAR_MAX": 2.63, "STRETCH_MIN": 4.99, "STRETCH_MAX": 8.99, "STAGGER_MAX": 1.11, "BUCKLE_MAX": 29.77, "PROPELLER_MIN": -28.35, "PROPELLER_MAX": -1.53, "OPENING_MIN": 104.32, "OPENING_MAX": 188.11},
            "cSW": {"SHEAR_MAX": 5.91, "STRETCH_MIN": -1.02, "STRETCH_MAX": 3.6, "STAGGER_MAX": 1.53, "BUCKLE_MAX": 31.58, "PROPELLER_MIN": -23.12, "PROPELLER_MAX": 32.86, "OPENING_MIN": 47.35, "OPENING_MAX": 116.45},
            "cWH": {"SHEAR_MAX": 4.32, "STRETCH_MIN": 0.75, "STRETCH_MAX": 4.83, "STAGGER_MAX": 1.35, "BUCKLE_MAX": 40.63, "PROPELLER_MIN": -28.65, "PROPELLER_MAX": 4.58, "OPENING_MIN": -76.33, "OPENING_MAX": -27.83},
            "cWS": {"SHEAR_MAX": 6.53, "STRETCH_MIN": -1.13, "STRETCH_MAX": 3.96, "STAGGER_MAX": 1.7, "BUCKLE_MAX": 29.54, "PROPELLER_MIN": -19.23, "PROPELLER_MAX": 31.88, "OPENING_MIN": 31.46, "OPENING_MAX": 111.3},
            "cWW": {"SHEAR_MAX": 2.27, "STRETCH_MIN": -0.94, "STRETCH_MAX": 0.12, "STAGGER_MAX": 0.69, "BUCKLE_MAX": 10.82, "PROPELLER_MIN": -15.57, "PROPELLER_MAX": 3.8, "OPENING_MIN": -12.63, "OPENING_MAX": 12.24},
            "tSW": {"SHEAR_MAX": 3.51, "STRETCH_MIN": -6.92, "STRETCH_MAX": -3.59, "STAGGER_MAX": 1.43, "BUCKLE_MAX": 36.03, "PROPELLER_MIN": -24.69, "PROPELLER_MAX": 16.72, "OPENING_MIN": -144.98, "OPENING_MAX": -51.63},
            "tWW": {"SHEAR_MAX": 1.69, "STRETCH_MIN": -3.86, "STRETCH_MAX": 2.43, "STAGGER_MAX": 1.36, "BUCKLE_MAX": 27.21, "PROPELLER_MIN": -30.36, "PROPELLER_MAX": 23.05, "OPENING_MIN": -178.1, "OPENING_MAX": 105.7}
        },
        "U-U": {
            "cWW": {"SHEAR_MAX": 2.39, "STRETCH_MIN": -2.12, "STRETCH_MAX": -0.96, "STAGGER_MAX": 0.78, "BUCKLE_MAX": 14.77, "PROPELLER_MIN": -22.22, "PROPELLER_MAX": 3.94, "OPENING_MIN": -10.22, "OPENING_MAX": 25.76},
            "tWW": {"SHEAR_MAX": 2.5, "STRETCH_MIN": -3.17, "STRETCH_MAX": 2.35, "STAGGER_MAX": 0.79, "BUCKLE_MAX": 27.29, "PROPELLER_MIN": -19.1, "PROPELLER_MAX": 13.63, "OPENING_MIN": -226.16, "OPENING_MAX": 33.38}
        }
    }
    
    # Base-pair + edge H-bond thresholds (mean ± SD, 2dp) generated from edge_type_analysis.json
    HBOND_THRESHOLDS_BY_BASE_PAIR_BY_EDGE = {
        "A-C": {
            "_OTHER": {"DIST_MIN": 2.75, "DIST_MAX": 3.38, "ANGLE_MIN": 91.37, "QUALITY_MIN": 0.6},
            "cSH": {"DIST_MIN": 2.87, "DIST_MAX": 3.46, "ANGLE_MIN": 75.46, "QUALITY_MIN": 0.54},
            "cSS": {"DIST_MIN": 2.63, "DIST_MAX": 3.23, "ANGLE_MIN": 94.03, "QUALITY_MIN": 0.69},
            "cSW": {"DIST_MIN": 2.62, "DIST_MAX": 3.22, "ANGLE_MIN": 93.1, "QUALITY_MIN": 0.67},
            "cWH": {"DIST_MIN": 2.78, "DIST_MAX": 3.38, "ANGLE_MIN": 80.34, "QUALITY_MIN": 0.59},
            "cWS": {"DIST_MIN": 2.67, "DIST_MAX": 3.26, "ANGLE_MIN": 93.6, "QUALITY_MIN": 0.66},
            "cWW": {"DIST_MIN": 2.83, "DIST_MAX": 3.38, "ANGLE_MIN": 93.74, "QUALITY_MIN": 0.63},
            "tHW": {"DIST_MIN": 2.78, "DIST_MAX": 3.3, "ANGLE_MIN": 102.07, "QUALITY_MIN": 0.65},
            "tSH": {"DIST_MIN": 2.77, "DIST_MAX": 3.37, "ANGLE_MIN": 97.65, "QUALITY_MIN": 0.63},
            "tSW": {"DIST_MIN": 2.66, "DIST_MAX": 3.29, "ANGLE_MIN": 88.99, "QUALITY_MIN": 0.67},
            "tWH": {"DIST_MIN": 2.78, "DIST_MAX": 3.27, "ANGLE_MIN": 105.37, "QUALITY_MIN": 0.65},
            "tWS": {"DIST_MIN": 2.7, "DIST_MAX": 3.37, "ANGLE_MIN": 90.7, "QUALITY_MIN": 0.63},
            "tWW": {"DIST_MIN": 2.84, "DIST_MAX": 3.34, "ANGLE_MIN": 98.67, "QUALITY_MIN": 0.67}
        },
        "A-G": {
            "_OTHER": {"DIST_MIN": 2.8, "DIST_MAX": 3.43, "ANGLE_MIN": 84.84, "QUALITY_MIN": 0.56},
            "cHW": {"DIST_MIN": 2.73, "DIST_MAX": 3.29, "ANGLE_MIN": 95.59, "QUALITY_MIN": 0.69},
            "cSH": {"DIST_MIN": 2.86, "DIST_MAX": 3.47, "ANGLE_MIN": 79.32, "QUALITY_MIN": 0.57},
            "cSS": {"DIST_MIN": 2.77, "DIST_MAX": 3.39, "ANGLE_MIN": 91.27, "QUALITY_MIN": 0.6},
            "cSW": {"DIST_MIN": 2.75, "DIST_MAX": 3.36, "ANGLE_MIN": 89.18, "QUALITY_MIN": 0.65},
            "cWH": {"DIST_MIN": 2.73, "DIST_MAX": 3.33, "ANGLE_MIN": 98.87, "QUALITY_MIN": 0.63},
            "cWS": {"DIST_MIN": 2.73, "DIST_MAX": 3.36, "ANGLE_MIN": 91.0, "QUALITY_MIN": 0.65},
            "cWW": {"DIST_MIN": 2.76, "DIST_MAX": 3.21, "ANGLE_MIN": 97.68, "QUALITY_MIN": 0.78},
            "tHS": {"DIST_MIN": 2.78, "DIST_MAX": 3.29, "ANGLE_MIN": 101.3, "QUALITY_MIN": 0.65},
            "tSH": {"DIST_MIN": 2.75, "DIST_MAX": 3.31, "ANGLE_MIN": 98.12, "QUALITY_MIN": 0.62},
            "tSS": {"DIST_MIN": 2.74, "DIST_MAX": 3.32, "ANGLE_MIN": 90.16, "QUALITY_MIN": 0.7},
            "tSW": {"DIST_MIN": 2.76, "DIST_MAX": 3.36, "ANGLE_MIN": 89.82, "QUALITY_MIN": 0.65},
            "tWH": {"DIST_MIN": 2.75, "DIST_MAX": 3.36, "ANGLE_MIN": 91.35, "QUALITY_MIN": 0.61},
            "tWS": {"DIST_MIN": 2.81, "DIST_MAX": 3.38, "ANGLE_MIN": 92.22, "QUALITY_MIN": 0.66}
        },
        "A-U": {
            "_OTHER": {"DIST_MIN": 2.63, "DIST_MAX": 3.34, "ANGLE_MIN": 94.59, "QUALITY_MIN": 0.6},
            "cHS": {"DIST_MIN": 2.66, "DIST_MAX": 3.22, "ANGLE_MIN": 102.47, "QUALITY_MIN": 0.64},
            "cHW": {"DIST_MIN": 2.75, "DIST_MAX": 3.23, "ANGLE_MIN": 100.22, "QUALITY_MIN": 0.75},
            "cSH": {"DIST_MIN": 2.74, "DIST_MAX": 3.4, "ANGLE_MIN": 85.5, "QUALITY_MIN": 0.57},
            "cSS": {"DIST_MIN": 2.7, "DIST_MAX": 3.28, "ANGLE_MIN": 95.25, "QUALITY_MIN": 0.68},
            "cSW": {"DIST_MIN": 2.62, "DIST_MAX": 3.23, "ANGLE_MIN": 91.61, "QUALITY_MIN": 0.62},
            "cWH": {"DIST_MIN": 2.75, "DIST_MAX": 3.26, "ANGLE_MIN": 98.03, "QUALITY_MIN": 0.72},
            "cWS": {"DIST_MIN": 2.73, "DIST_MAX": 3.32, "ANGLE_MIN": 93.49, "QUALITY_MIN": 0.62},
            "cWW": {"DIST_MIN": 2.77, "DIST_MAX": 3.13, "ANGLE_MIN": 106.94, "QUALITY_MIN": 0.85},
            "tHW": {"DIST_MIN": 2.71, "DIST_MAX": 3.16, "ANGLE_MIN": 106.97, "QUALITY_MIN": 0.71},
            "tSH": {"DIST_MIN": 2.69, "DIST_MAX": 3.38, "ANGLE_MIN": 83.91, "QUALITY_MIN": 0.63},
            "tSW": {"DIST_MIN": 2.73, "DIST_MAX": 3.36, "ANGLE_MIN": 88.23, "QUALITY_MIN": 0.65},
            "tWH": {"DIST_MIN": 2.7, "DIST_MAX": 3.17, "ANGLE_MIN": 107.65, "QUALITY_MIN": 0.7},
            "tWS": {"DIST_MIN": 2.7, "DIST_MAX": 3.4, "ANGLE_MIN": 88.11, "QUALITY_MIN": 0.6},
            "tWW": {"DIST_MIN": 2.72, "DIST_MAX": 3.19, "ANGLE_MIN": 98.53, "QUALITY_MIN": 0.7}
        },
        "C-G": {
            "_OTHER": {"DIST_MIN": 2.71, "DIST_MAX": 3.35, "ANGLE_MIN": 91.77, "QUALITY_MIN": 0.6},
            "cSH": {"DIST_MIN": 2.92, "DIST_MAX": 3.52, "ANGLE_MIN": 74.9, "QUALITY_MIN": 0.54},
            "cSS": {"DIST_MIN": 2.72, "DIST_MAX": 3.31, "ANGLE_MIN": 97.63, "QUALITY_MIN": 0.64},
            "cSW": {"DIST_MIN": 2.79, "DIST_MAX": 3.41, "ANGLE_MIN": 81.27, "QUALITY_MIN": 0.59},
            "cWH": {"DIST_MIN": 2.77, "DIST_MAX": 3.44, "ANGLE_MIN": 80.47, "QUALITY_MIN": 0.56},
            "cWS": {"DIST_MIN": 2.8, "DIST_MAX": 3.36, "ANGLE_MIN": 86.72, "QUALITY_MIN": 0.61},
            "cWW": {"DIST_MIN": 2.72, "DIST_MAX": 3.1, "ANGLE_MIN": 108.1, "QUALITY_MIN": 0.81},
            "tSW": {"DIST_MIN": 2.78, "DIST_MAX": 3.37, "ANGLE_MIN": 98.32, "QUALITY_MIN": 0.61},
            "tWH": {"DIST_MIN": 2.77, "DIST_MAX": 3.32, "ANGLE_MIN": 96.94, "QUALITY_MIN": 0.59},
            "tWS": {"DIST_MIN": 2.72, "DIST_MAX": 3.35, "ANGLE_MIN": 100.54, "QUALITY_MIN": 0.63},
            "tWW": {"DIST_MIN": 2.74, "DIST_MAX": 3.24, "ANGLE_MIN": 98.08, "QUALITY_MIN": 0.71}
        },
        "C-U": {
            "_OTHER": {"DIST_MIN": 2.68, "DIST_MAX": 3.32, "ANGLE_MIN": 96.42, "QUALITY_MIN": 0.6},
            "cSH": {"DIST_MIN": 2.72, "DIST_MAX": 3.37, "ANGLE_MIN": 86.84, "QUALITY_MIN": 0.57},
            "cSW": {"DIST_MIN": 2.7, "DIST_MAX": 3.33, "ANGLE_MIN": 95.39, "QUALITY_MIN": 0.62},
            "cWW": {"DIST_MIN": 2.77, "DIST_MAX": 3.36, "ANGLE_MIN": 100.78, "QUALITY_MIN": 0.65}
        },
        "G-U": {
            "_OTHER": {"DIST_MIN": 2.67, "DIST_MAX": 3.33, "ANGLE_MIN": 91.25, "QUALITY_MIN": 0.59},
            "cSH": {"DIST_MIN": 2.72, "DIST_MAX": 3.42, "ANGLE_MIN": 79.47, "QUALITY_MIN": 0.61},
            "cSS": {"DIST_MIN": 2.69, "DIST_MAX": 3.2, "ANGLE_MIN": 104.47, "QUALITY_MIN": 0.7},
            "cSW": {"DIST_MIN": 2.72, "DIST_MAX": 3.42, "ANGLE_MIN": 90.77, "QUALITY_MIN": 0.57},
            "cWH": {"DIST_MIN": 2.74, "DIST_MAX": 3.42, "ANGLE_MIN": 81.98, "QUALITY_MIN": 0.53},
            "cWS": {"DIST_MIN": 2.79, "DIST_MAX": 3.41, "ANGLE_MIN": 99.37, "QUALITY_MIN": 0.56},
            "cWW": {"DIST_MIN": 2.71, "DIST_MAX": 3.15, "ANGLE_MIN": 98.66, "QUALITY_MIN": 0.77},
            "tSH": {"DIST_MIN": 2.59, "DIST_MAX": 3.25, "ANGLE_MIN": 96.72, "QUALITY_MIN": 0.64},
            "tSW": {"DIST_MIN": 2.64, "DIST_MAX": 3.33, "ANGLE_MIN": 96.18, "QUALITY_MIN": 0.62},
            "tWS": {"DIST_MIN": 2.72, "DIST_MAX": 3.33, "ANGLE_MIN": 94.54, "QUALITY_MIN": 0.63},
            "tWW": {"DIST_MIN": 2.77, "DIST_MAX": 3.35, "ANGLE_MIN": 90.63, "QUALITY_MIN": 0.64}
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
    MAX_HOTSPOT_SCORE = 70.0




    # ===== BACKBONE TORSION THRESHOLDS (mean ± 1SD, outside = outlier) =====
    # Structure: edge -> bp_type -> angle -> (min, max)
    # Only edge/bp_type with >=1000 data; edges with '.', '--', or <1k go to _OTHER
    # Generated from torsion_thresholds_analysis.json
    TORSION_THRESHOLDS_BY_EDGE_BY_BASE_PAIR = {
        "cHW": {
            "A-G": {"alpha": (-97.43, 91.09), "beta": (-118.08, 182.01), "chi": (-171.1, 52.46), "delta": (72.38, 129.56), "epsilon": (-161.58, -85.2), "eta": (-75.75, 195.67), "gamma": (-98.12, 110.48), "v1": (-34.76, 21.49), "v2": (-17.39, 47.01), "v3": (-44.19, 8.35), "v4": (3.77, 24.68), "zeta": (-115.64, 36.98)},
            "A-U": {"alpha": (-89.14, 125.58), "beta": (-147.35, 163.39), "chi": (-181.69, 13.85), "delta": (78.79, 143.39), "epsilon": (-162.46, -73.45), "eta": (-71.96, 179.78), "gamma": (-58.82, 105.45), "v1": (-25.18, 35.24), "v2": (-33.1, 37.43), "v3": (-37.68, 21.02), "v4": (-0.5, 24.2), "zeta": (-99.45, 67.79)},
            "G-G": {"alpha": (-112.32, 79.53), "beta": (-125.86, 179.88), "chi": (-191.75, 27.67), "delta": (72.22, 133.73), "epsilon": (-171.97, -66.87), "eta": (-69.47, 183.15), "gamma": (-47.38, 132.34), "v1": (-33.2, 23.55), "v2": (-21.18, 46.06), "v3": (-44.01, 12.49), "v4": (1.25, 25.06), "zeta": (-103.15, 57.91)},
            "_OTHER": {"alpha": (-160.14, 168.52), "beta": (-211.04, 243.28), "chi": (-220.35, 190.03), "delta": (34.66, 159.32), "epsilon": (-214.35, 41.48), "eta": (-178.33, 233.04), "gamma": (-169.47, 170.89), "v1": (-64.03, 43.98), "v2": (-48.88, 53.38), "v3": (-50.35, 39.05), "v4": (-34.34, 40.12), "zeta": (-177.99, 84.56)},
        },
        "cSH": {
            "A-A": {"alpha": (-135.07, 94.93), "beta": (-126.82, 172.77), "chi": (-203.37, 37.7), "delta": (80.27, 142.21), "epsilon": (-154.86, -78.09), "eta": (-119.34, 166.97), "gamma": (-6.84, 110.37), "v1": (-29.02, 30.99), "v2": (-30.52, 39.23), "v3": (-36.81, 20.68), "v4": (-2.59, 20.5), "zeta": (-134.38, 92.31)},
            "A-C": {"alpha": (-121.34, 94.18), "beta": (-132.8, 176.74), "chi": (-194.01, -54.33), "delta": (73.5, 133.94), "epsilon": (-169.22, -58.16), "eta": (-99.3, 181.81), "gamma": (8.86, 96.44), "v1": (-32.77, 23.14), "v2": (-20.87, 45.08), "v3": (-42.83, 12.39), "v4": (1.07, 24.11), "zeta": (-121.11, 85.89)},
            "A-G": {"alpha": (-120.88, 76.66), "beta": (-134.99, 170.23), "chi": (-206.32, 28.99), "delta": (72.61, 130.16), "epsilon": (-172.19, -66.97), "eta": (-102.19, 181.8), "gamma": (-50.75, 119.67), "v1": (-33.51, 19.88), "v2": (-16.98, 45.97), "v3": (-43.78, 9.26), "v4": (1.71, 25.16), "zeta": (-117.8, 68.27)},
            "A-U": {"alpha": (-119.98, 94.8), "beta": (-137.58, 167.97), "chi": (-202.12, 22.93), "delta": (75.4, 134.83), "epsilon": (-167.35, -60.9), "eta": (-87.82, 190.25), "gamma": (-6.48, 115.33), "v1": (-30.82, 26.6), "v2": (-24.0, 42.8), "v3": (-40.97, 14.16), "v4": (1.49, 23.51), "zeta": (-120.4, 86.98)},
            "C-U": {"alpha": (-115.99, 51.55), "beta": (-111.59, 201.97), "chi": (-184.74, -91.93), "delta": (72.16, 133.91), "epsilon": (-192.81, -23.56), "eta": (-95.26, 181.83), "gamma": (0.68, 97.56), "v1": (-31.32, 25.69), "v2": (-22.44, 45.08), "v3": (-44.08, 12.38), "v4": (3.17, 25.94), "zeta": (-118.67, 97.6)},
            "G-G": {"alpha": (-120.71, 67.29), "beta": (-132.7, 162.35), "chi": (-196.43, 8.94), "delta": (78.54, 142.94), "epsilon": (-175.76, -45.8), "eta": (-128.1, 162.2), "gamma": (-44.66, 125.97), "v1": (-26.9, 32.63), "v2": (-31.83, 38.47), "v3": (-37.68, 21.32), "v4": (-2.27, 22.86), "zeta": (-128.54, 62.26)},
            "G-U": {"alpha": (-121.3, 37.13), "beta": (-77.22, 185.69), "chi": (-201.42, -8.72), "delta": (80.44, 144.44), "epsilon": (-189.46, -59.33), "eta": (-87.65, 175.86), "gamma": (-58.08, 126.67), "v1": (-25.12, 36.33), "v2": (-34.35, 36.72), "v3": (-36.46, 21.83), "v4": (-0.09, 22.43), "zeta": (-98.62, 118.62)},
            "_OTHER": {"alpha": (-176.39, 173.01), "beta": (-201.63, 201.35), "chi": (-212.63, 173.01), "delta": (-38.37, 164.27), "epsilon": (-173.47, 172.89), "eta": (-173.27, 191.61), "gamma": (-175.89, 177.27), "v1": (-40.11, 42.4), "v2": (-43.28, 45.36), "v3": (-47.78, 40.29), "v4": (-33.0, 37.21), "zeta": (-173.79, 178.95)},
        },
        "cSS": {
            "A-C": {"alpha": (-115.7, 31.69), "beta": (-88.1, 211.22), "chi": (-189.43, -111.48), "delta": (73.4, 94.71), "epsilon": (-172.97, -102.54), "eta": (-13.33, 219.49), "gamma": (-3.82, 100.84), "v1": (-30.81, -12.41), "v2": (22.13, 44.01), "v3": (-43.28, -23.84), "v4": (15.0, 26.46), "zeta": (-110.79, -32.36)},
            "A-G": {"alpha": (-107.72, 59.5), "beta": (-100.66, 197.47), "chi": (-213.49, -38.05), "delta": (68.58, 115.13), "epsilon": (-171.94, -87.98), "eta": (-62.38, 206.81), "gamma": (-31.39, 115.34), "v1": (-36.49, 7.6), "v2": (-1.11, 50.15), "v3": (-47.73, -5.04), "v4": (8.65, 27.25), "zeta": (-123.3, -0.97)},
            "C-G": {"alpha": (-112.71, 12.08), "beta": (-82.1, 217.65), "chi": (-196.93, -99.13), "delta": (70.25, 100.74), "epsilon": (-168.96, -100.55), "eta": (4.02, 227.5), "gamma": (5.28, 97.22), "v1": (-34.67, -6.47), "v2": (14.96, 48.11), "v3": (-46.04, -17.95), "v4": (13.11, 26.45), "zeta": (-131.73, -21.53)},
            "G-U": {"alpha": (-113.35, 16.19), "beta": (-114.49, 205.8), "chi": (-193.51, -94.3), "delta": (70.11, 105.35), "epsilon": (-175.29, -100.24), "eta": (4.15, 231.55), "gamma": (11.95, 99.15), "v1": (-34.08, -1.92), "v2": (9.98, 47.87), "v3": (-46.28, -14.15), "v4": (12.12, 27.11), "zeta": (-117.63, -32.69)},
            "_OTHER": {"alpha": (-129.82, 214.69), "beta": (-238.6, 240.88), "chi": (-228.81, 160.73), "delta": (73.34, 149.08), "epsilon": (-202.6, 170.53), "eta": (-185.07, 220.36), "gamma": (-177.9, 148.06), "v1": (-31.88, 41.11), "v2": (-36.29, 42.75), "v3": (-43.04, 25.82), "v4": (-13.67, 33.13), "zeta": (-144.76, 129.25)},
        },
        "cSW": {
            "A-A": {"alpha": (-117.76, 53.59), "beta": (-104.14, 200.3), "chi": (-220.98, -16.02), "delta": (69.5, 103.83), "epsilon": (-165.46, -100.92), "eta": (-15.09, 214.55), "gamma": (-27.86, 115.85), "v1": (-35.96, -3.45), "v2": (11.47, 49.42), "v3": (-46.97, -15.16), "v4": (12.2, 26.56), "zeta": (-121.41, 16.45)},
            "A-C": {"alpha": (-111.95, 82.73), "beta": (-135.42, 175.23), "chi": (-222.41, -15.38), "delta": (69.14, 110.79), "epsilon": (-171.46, -95.76), "eta": (-20.62, 202.43), "gamma": (-22.91, 108.05), "v1": (-36.89, 2.26), "v2": (4.31, 50.16), "v3": (-47.23, -8.92), "v4": (9.58, 26.06), "zeta": (-107.81, -33.16)},
            "A-G": {"alpha": (-101.12, 60.61), "beta": (-100.68, 199.79), "chi": (-211.55, 0.22), "delta": (69.42, 126.65), "epsilon": (-166.97, -85.87), "eta": (-79.79, 192.95), "gamma": (-53.1, 121.8), "v1": (-35.75, 16.82), "v2": (-12.89, 49.22), "v3": (-46.76, 5.32), "v4": (4.38, 26.21), "zeta": (-133.6, 39.0)},
            "A-U": {"alpha": (-120.33, 79.64), "beta": (-104.56, 203.31), "chi": (-206.49, -54.88), "delta": (68.96, 122.82), "epsilon": (-167.23, -94.04), "eta": (-23.3, 189.29), "gamma": (-40.84, 111.77), "v1": (-35.67, 14.65), "v2": (-9.55, 49.49), "v3": (-47.32, 1.93), "v4": (6.39, 26.95), "zeta": (-116.26, 24.86)},
            "C-U": {"alpha": (-124.43, 76.57), "beta": (-99.21, 197.4), "chi": (-193.14, -93.82), "delta": (70.51, 127.19), "epsilon": (-158.5, -98.82), "eta": (25.04, 197.62), "gamma": (-24.29, 122.9), "v1": (-31.97, 17.57), "v2": (-13.65, 46.64), "v3": (-46.21, 5.83), "v4": (4.42, 28.04), "zeta": (-128.69, 56.85)},
            "G-U": {"alpha": (-116.83, 70.21), "beta": (-104.41, 194.52), "chi": (-204.14, -24.78), "delta": (70.5, 130.02), "epsilon": (-172.73, -71.62), "eta": (-18.3, 203.51), "gamma": (-35.28, 120.58), "v1": (-34.15, 21.41), "v2": (-17.63, 47.38), "v3": (-45.23, 8.68), "v4": (3.91, 25.62), "zeta": (-118.9, 26.23)},
            "_OTHER": {"alpha": (-199.59, 177.03), "beta": (-206.66, 233.96), "chi": (-224.26, 178.77), "delta": (70.29, 161.1), "epsilon": (-174.9, 115.35), "eta": (-179.82, 228.54), "gamma": (-179.58, 189.65), "v1": (-36.28, 40.65), "v2": (-45.73, 44.44), "v3": (-45.44, 39.23), "v4": (-21.16, 36.12), "zeta": (-180.85, 172.41)},
        },
        "cWH": {
            "A-C": {"alpha": (-127.71, 85.95), "beta": (-144.93, 156.78), "chi": (-222.61, 1.22), "delta": (67.95, 117.07), "epsilon": (-174.35, -61.91), "eta": (-45.68, 201.28), "gamma": (-105.2, 136.07), "v1": (-37.48, 8.37), "v2": (-2.44, 50.93), "v3": (-47.96, -3.65), "v4": (7.96, 26.55), "zeta": (-129.08, 30.59)},
            "A-G": {"alpha": (-92.63, 83.74), "beta": (-160.06, 144.58), "chi": (-183.21, 31.55), "delta": (77.93, 141.95), "epsilon": (-160.58, -64.12), "eta": (-141.43, 149.01), "gamma": (-62.34, 112.8), "v1": (-27.66, 33.63), "v2": (-31.18, 39.41), "v3": (-38.57, 19.41), "v4": (0.02, 23.58), "zeta": (-117.97, 40.74)},
            "A-U": {"alpha": (-112.28, 71.24), "beta": (-126.62, 179.06), "chi": (-198.91, -29.05), "delta": (75.98, 140.9), "epsilon": (-166.25, -69.02), "eta": (-63.56, 188.51), "gamma": (-56.33, 117.58), "v1": (-27.8, 31.51), "v2": (-29.56, 40.82), "v3": (-40.61, 18.61), "v4": (0.15, 24.97), "zeta": (-108.24, 72.32)},
            "C-G": {"alpha": (-109.92, 76.64), "beta": (-129.12, 172.75), "chi": (-194.13, 47.12), "delta": (68.86, 124.37), "epsilon": (-182.66, -61.28), "eta": (-66.31, 200.88), "gamma": (-56.25, 126.76), "v1": (-35.78, 15.62), "v2": (-11.08, 49.13), "v3": (-46.76, 3.77), "v4": (4.37, 27.04), "zeta": (-125.21, 39.02)},
            "G-G": {"alpha": (-110.0, 95.24), "beta": (-105.56, 193.42), "chi": (-195.71, -13.41), "delta": (72.67, 135.79), "epsilon": (-167.75, -73.28), "eta": (-88.49, 182.28), "gamma": (-44.37, 123.35), "v1": (-32.72, 25.12), "v2": (-22.89, 45.57), "v3": (-43.69, 13.79), "v4": (0.95, 25.02), "zeta": (-104.81, 68.79)},
            "G-U": {"alpha": (-100.47, 93.57), "beta": (-108.99, 197.46), "chi": (-218.27, -22.61), "delta": (68.95, 105.71), "epsilon": (-179.81, -71.92), "eta": (-36.53, 220.01), "gamma": (-69.18, 122.08), "v1": (-36.04, -1.04), "v2": (9.08, 49.43), "v3": (-47.11, -13.3), "v4": (11.11, 27.39), "zeta": (-119.39, 64.66)},
            "_OTHER": {"alpha": (-188.64, 199.98), "beta": (-234.88, 214.55), "chi": (-234.89, 177.04), "delta": (63.26, 152.1), "epsilon": (-200.31, 179.28), "eta": (-179.5, 240.52), "gamma": (-169.54, 195.22), "v1": (-36.81, 39.51), "v2": (-38.61, 51.43), "v3": (-51.15, 28.99), "v4": (-13.13, 36.03), "zeta": (-171.04, 154.66)},
        },
        "cWS": {
            "A-A": {"alpha": (-110.24, 97.41), "beta": (-98.43, 198.38), "chi": (-215.8, 8.52), "delta": (70.75, 121.46), "epsilon": (-168.1, -86.04), "eta": (-19.94, 199.12), "gamma": (-48.93, 128.99), "v1": (-35.62, 15.3), "v2": (-9.24, 48.48), "v3": (-45.8, 0.91), "v4": (7.42, 25.86), "zeta": (-123.7, 5.83)},
            "A-C": {"alpha": (-109.1, 71.81), "beta": (-113.8, 193.64), "chi": (-225.32, -20.1), "delta": (69.79, 102.17), "epsilon": (-176.04, -95.17), "eta": (-13.47, 208.14), "gamma": (-47.58, 121.99), "v1": (-35.14, -5.38), "v2": (13.64, 48.67), "v3": (-46.52, -16.84), "v4": (12.67, 26.61), "zeta": (-104.75, -35.99)},
            "A-G": {"alpha": (-107.28, 91.73), "beta": (-107.46, 195.69), "chi": (-224.28, -10.79), "delta": (68.71, 109.98), "epsilon": (-172.33, -94.32), "eta": (-41.57, 206.48), "gamma": (-28.0, 116.21), "v1": (-36.62, 2.26), "v2": (4.87, 50.17), "v3": (-47.55, -9.8), "v4": (10.34, 26.74), "zeta": (-119.15, -17.89)},
            "A-U": {"alpha": (-105.81, 52.9), "beta": (-85.03, 200.79), "chi": (-213.16, -55.69), "delta": (69.24, 109.66), "epsilon": (-188.26, -69.27), "eta": (-45.41, 188.29), "gamma": (-33.46, 117.77), "v1": (-35.73, 2.54), "v2": (4.84, 49.35), "v3": (-46.99, -10.14), "v4": (11.15, 26.47), "zeta": (-130.72, 32.84)},
            "G-U": {"alpha": (-111.51, 67.23), "beta": (-121.5, 185.64), "chi": (-205.23, -53.84), "delta": (69.09, 125.05), "epsilon": (-169.81, -81.77), "eta": (-19.64, 200.16), "gamma": (-21.32, 120.35), "v1": (-34.25, 17.88), "v2": (-12.88, 48.5), "v3": (-47.05, 4.31), "v4": (6.01, 27.58), "zeta": (-130.39, 20.97)},
            "_OTHER": {"alpha": (-159.52, 176.58), "beta": (-176.43, 238.7), "chi": (-239.55, 178.4), "delta": (68.0, 165.02), "epsilon": (-179.35, 169.35), "eta": (-162.6, 230.9), "gamma": (-172.72, 178.69), "v1": (-34.32, 37.48), "v2": (-39.93, 46.53), "v3": (-46.53, 40.75), "v4": (-29.28, 43.05), "zeta": (-165.16, 101.06)},
        },
        "cWW": {
            "A-A": {"alpha": (-114.38, 35.91), "beta": (-55.49, 216.2), "chi": (-227.04, -2.45), "delta": (69.12, 105.9), "epsilon": (-174.06, -96.94), "eta": (-7.7, 223.9), "gamma": (-7.97, 117.77), "v1": (-38.12, -2.56), "v2": (10.01, 50.66), "v3": (-46.91, -13.55), "v4": (10.92, 25.25), "zeta": (-112.43, -7.96)},
            "A-C": {"alpha": (-112.57, 35.49), "beta": (-66.95, 216.51), "chi": (-227.26, -30.93), "delta": (70.41, 99.42), "epsilon": (-171.19, -96.69), "eta": (9.36, 231.04), "gamma": (-10.68, 116.45), "v1": (-35.18, -7.59), "v2": (16.35, 48.3), "v3": (-45.82, -19.18), "v4": (13.69, 25.82), "zeta": (-110.44, -25.4)},
            "A-DT": {"alpha": (-109.26, 29.9), "beta": (-64.1, 219.28), "chi": (-194.69, -79.22), "delta": (71.3, 115.62), "epsilon": (-203.69, -59.09), "eta": (16.32, 232.44), "gamma": (-25.26, 112.31), "v1": (-34.08, 13.55), "v2": (-4.63, 47.92), "v3": (-47.01, -4.38), "v4": (9.47, 30.58), "zeta": (-116.23, -39.12)},
            "A-G": {"alpha": (-111.85, 47.51), "beta": (-81.89, 201.91), "chi": (-223.67, 24.71), "delta": (68.06, 112.36), "epsilon": (-174.68, -88.26), "eta": (-35.65, 224.48), "gamma": (-33.27, 122.66), "v1": (-38.35, 4.66), "v2": (2.07, 51.4), "v3": (-47.84, -7.5), "v4": (9.52, 25.88), "zeta": (-121.46, 0.63)},
            "A-PSU": {"alpha": (-105.6, 28.28), "beta": (-41.78, 229.44), "chi": (-226.04, -67.89), "delta": (73.5, 88.52), "epsilon": (-165.41, -115.57), "eta": (75.42, 223.31), "gamma": (-26.35, 111.67), "v1": (-31.21, -16.92), "v2": (27.99, 43.82), "v3": (-42.58, -28.98), "v4": (17.08, 25.71), "zeta": (-97.76, -60.65)},
            "A-U": {"alpha": (-111.45, 29.19), "beta": (-80.03, 217.64), "chi": (-218.8, -67.39), "delta": (70.15, 96.77), "epsilon": (-175.27, -103.87), "eta": (38.0, 227.46), "gamma": (-11.55, 111.97), "v1": (-34.77, -9.27), "v2": (18.82, 48.08), "v3": (-45.88, -21.6), "v4": (15.02, 26.25), "zeta": (-109.3, -33.52)},
            "C-C": {"alpha": (-115.19, 40.01), "beta": (-90.02, 208.42), "chi": (-203.94, -93.57), "delta": (70.71, 97.17), "epsilon": (-170.39, -94.13), "eta": (14.9, 225.35), "gamma": (-10.94, 118.05), "v1": (-33.48, -8.62), "v2": (18.41, 46.98), "v3": (-45.37, -21.51), "v4": (15.22, 26.68), "zeta": (-111.74, -28.87)},
            "C-DG": {"alpha": (-110.37, 46.59), "beta": (-75.25, 210.04), "chi": (-201.49, -68.16), "delta": (71.12, 117.07), "epsilon": (-201.38, -59.66), "eta": (18.03, 229.34), "gamma": (-40.57, 120.14), "v1": (-34.54, 13.04), "v2": (-5.04, 48.15), "v3": (-47.14, -2.99), "v4": (7.37, 30.7), "zeta": (-119.92, -34.03)},
            "C-G": {"alpha": (-112.22, 40.13), "beta": (-89.98, 211.53), "chi": (-225.94, -50.4), "delta": (70.73, 94.99), "epsilon": (-174.3, -105.43), "eta": (27.65, 228.58), "gamma": (-19.96, 117.9), "v1": (-34.2, -11.24), "v2": (21.07, 47.37), "v3": (-45.35, -23.29), "v4": (15.25, 26.27), "zeta": (-108.44, -31.15)},
            "C-OMG": {"alpha": (-109.23, 20.21), "beta": (-64.1, 218.67), "chi": (-226.13, -44.46), "delta": (68.41, 96.67), "epsilon": (-162.13, -104.66), "eta": (-19.43, 222.83), "gamma": (-3.57, 103.43), "v1": (-33.88, -9.86), "v2": (19.97, 47.17), "v3": (-45.51, -22.65), "v4": (15.07, 27.08), "zeta": (-110.68, -28.14)},
            "C-U": {"alpha": (-111.87, 42.73), "beta": (-101.58, 204.86), "chi": (-204.72, -89.02), "delta": (67.95, 103.55), "epsilon": (-171.86, -89.91), "eta": (2.66, 226.36), "gamma": (-18.08, 115.61), "v1": (-35.47, -2.3), "v2": (10.86, 49.67), "v3": (-47.79, -15.26), "v4": (13.11, 27.62), "zeta": (-116.04, -11.35)},
            "DA-DT": {"alpha": (-96.51, 25.18), "beta": (-59.06, 213.61), "chi": (-136.47, -88.94), "delta": (111.61, 152.99), "epsilon": (-221.08, -2.35), "eta": (-27.84, 233.27), "gamma": (-29.85, 101.13), "v1": (14.96, 45.72), "v2": (-44.22, -6.87), "v3": (-6.08, 31.63), "v4": (-8.58, 21.5), "zeta": (-163.32, -23.41)},
            "DA-U": {"alpha": (-109.15, 36.13), "beta": (-66.12, 215.73), "chi": (-189.17, -87.22), "delta": (71.67, 120.19), "epsilon": (-196.49, -70.58), "eta": (14.5, 229.48), "gamma": (-37.45, 113.44), "v1": (-34.04, 15.98), "v2": (-8.84, 47.43), "v3": (-46.24, 0.23), "v4": (6.48, 29.57), "zeta": (-124.69, -32.42)},
            "DC-G": {"alpha": (-109.2, 48.76), "beta": (-76.41, 211.3), "chi": (-204.42, -58.68), "delta": (71.28, 116.43), "epsilon": (-197.88, -66.84), "eta": (18.56, 228.95), "gamma": (-43.67, 118.96), "v1": (-34.47, 14.12), "v2": (-5.26, 48.16), "v3": (-47.06, -3.83), "v4": (9.12, 30.52), "zeta": (-124.01, -21.57)},
            "G-G": {"alpha": (-114.59, 52.89), "beta": (-79.49, 201.01), "chi": (-221.67, 43.07), "delta": (68.72, 100.62), "epsilon": (-178.11, -79.22), "eta": (-16.29, 226.53), "gamma": (-11.37, 131.16), "v1": (-38.1, -7.41), "v2": (15.69, 50.39), "v3": (-46.8, -17.79), "v4": (10.98, 26.35), "zeta": (-123.38, -18.85)},
            "G-OMC": {"alpha": (-106.59, 55.75), "beta": (-43.05, 221.43), "chi": (-234.86, -33.38), "delta": (69.53, 90.71), "epsilon": (-163.78, -117.9), "eta": (17.58, 228.35), "gamma": (-29.8, 116.16), "v1": (-31.41, -16.45), "v2": (27.34, 44.51), "v3": (-43.55, -28.33), "v4": (16.83, 26.46), "zeta": (-102.65, -35.68)},
            "G-PSU": {"alpha": (-107.75, 25.32), "beta": (-65.54, 220.77), "chi": (-227.07, -65.97), "delta": (75.41, 86.4), "epsilon": (-160.69, -123.95), "eta": (27.5, 234.44), "gamma": (-20.48, 109.57), "v1": (-29.35, -17.17), "v2": (29.63, 41.55), "v3": (-40.84, -31.25), "v4": (18.5, 25.96), "zeta": (-103.78, -46.81)},
            "G-U": {"alpha": (-113.34, 40.82), "beta": (-69.22, 218.8), "chi": (-229.0, -38.36), "delta": (71.14, 94.77), "epsilon": (-172.23, -105.53), "eta": (29.53, 230.6), "gamma": (-18.81, 119.89), "v1": (-33.65, -11.1), "v2": (21.11, 46.85), "v3": (-45.01, -23.53), "v4": (15.65, 26.25), "zeta": (-111.92, -33.42)},
            "U-U": {"alpha": (-114.3, 29.95), "beta": (-82.3, 217.09), "chi": (-206.05, -96.26), "delta": (71.23, 95.89), "epsilon": (-169.37, -101.39), "eta": (-6.76, 229.98), "gamma": (-3.52, 109.08), "v1": (-33.16, -9.62), "v2": (19.59, 46.61), "v3": (-44.99, -22.56), "v4": (15.9, 26.26), "zeta": (-116.79, -29.05)},
            "_OTHER": {"alpha": (-179.27, 211.9), "beta": (-238.83, 248.94), "chi": (-243.11, 214.22), "delta": (-91.5, 167.13), "epsilon": (-231.27, 184.08), "eta": (-234.1, 243.77), "gamma": (-221.72, 186.95), "v1": (-79.77, 48.98), "v2": (-51.94, 69.49), "v3": (-66.62, 51.18), "v4": (-37.57, 45.57), "zeta": (-175.75, 156.1)},
        },
        "tHH": {
            "A-A": {"alpha": (-96.25, 117.4), "beta": (-46.97, 212.96), "chi": (-198.6, -62.67), "delta": (80.68, 144.42), "epsilon": (-157.96, -79.06), "eta": (-91.97, 164.86), "gamma": (-6.0, 104.51), "v1": (-27.04, 35.14), "v2": (-33.45, 38.01), "v3": (-36.63, 21.41), "v4": (-0.21, 21.14), "zeta": (-125.24, 91.43)},
            "_OTHER": {"alpha": (-178.18, 237.34), "beta": (-191.46, 230.77), "chi": (-208.92, 179.23), "delta": (71.62, 162.54), "epsilon": (-161.93, 18.17), "eta": (-182.86, 184.26), "gamma": (-168.93, 173.83), "v1": (-35.86, 39.18), "v2": (-43.5, 48.73), "v3": (-44.02, 33.15), "v4": (-16.05, 34.54), "zeta": (-171.71, 177.38)},
        },
        "tHS": {
            "A-A": {"alpha": (-108.82, 9.31), "beta": (-25.56, 225.89), "chi": (-232.04, 13.93), "delta": (70.52, 100.38), "epsilon": (-171.08, -88.54), "eta": (-38.01, 227.78), "gamma": (3.45, 106.25), "v1": (-35.6, -7.11), "v2": (15.59, 48.63), "v3": (-45.94, -18.44), "v4": (13.39, 25.53), "zeta": (-115.96, 28.1)},
            "A-G": {"alpha": (-110.38, 13.19), "beta": (-51.78, 220.92), "chi": (-212.21, -52.6), "delta": (68.63, 114.19), "epsilon": (-172.67, -77.6), "eta": (-43.73, 227.79), "gamma": (-0.04, 108.69), "v1": (-37.48, 5.85), "v2": (0.35, 50.72), "v3": (-47.53, -5.91), "v4": (8.92, 25.85), "zeta": (-134.39, 26.97)},
            "G-G": {"alpha": (-109.43, 45.81), "beta": (-96.03, 208.3), "chi": (-218.41, -24.69), "delta": (68.87, 119.3), "epsilon": (-170.84, -83.64), "eta": (-39.75, 227.24), "gamma": (-27.85, 122.23), "v1": (-36.9, 10.08), "v2": (-4.93, 50.23), "v3": (-47.22, -1.38), "v4": (7.2, 25.69), "zeta": (-119.82, 71.51)},
            "_OTHER": {"alpha": (-155.75, 162.81), "beta": (-177.95, 234.16), "chi": (-205.13, 166.32), "delta": (73.99, 153.56), "epsilon": (-166.92, -94.23), "eta": (-226.09, 229.21), "gamma": (-146.44, 164.17), "v1": (-34.2, 40.04), "v2": (-39.36, 47.73), "v3": (-45.74, 29.19), "v4": (-6.42, 48.8), "zeta": (-145.11, 169.31)},
        },
        "tHW": {
            "A-A": {"alpha": (-101.2, 79.47), "beta": (-61.38, 211.67), "chi": (-212.14, 0.23), "delta": (68.53, 115.64), "epsilon": (-166.09, -89.29), "eta": (-85.54, 177.84), "gamma": (-37.79, 119.83), "v1": (-36.69, 8.6), "v2": (-1.91, 50.21), "v3": (-47.56, -4.74), "v4": (9.17, 26.72), "zeta": (-111.67, 31.82)},
            "A-C": {"alpha": (-115.68, 46.33), "beta": (-120.36, 183.74), "chi": (-208.41, -26.21), "delta": (69.68, 124.45), "epsilon": (-166.78, -86.13), "eta": (-125.49, 174.02), "gamma": (-27.68, 107.6), "v1": (-33.78, 17.65), "v2": (-12.34, 47.82), "v3": (-46.52, 3.77), "v4": (6.0, 27.8), "zeta": (-119.49, 64.28)},
            "A-U": {"alpha": (-118.54, 41.77), "beta": (-114.89, 177.66), "chi": (-217.29, -5.01), "delta": (69.83, 123.69), "epsilon": (-167.25, -81.96), "eta": (-134.91, 179.2), "gamma": (-44.04, 138.46), "v1": (-35.24, 16.9), "v2": (-11.49, 48.6), "v3": (-46.33, 3.06), "v4": (6.33, 26.53), "zeta": (-130.58, 23.19)},
            "G-G": {"alpha": (-114.57, 67.21), "beta": (-104.62, 188.65), "chi": (-199.1, -20.28), "delta": (72.78, 132.1), "epsilon": (-173.78, -63.24), "eta": (-120.52, 179.46), "gamma": (-52.19, 136.43), "v1": (-32.55, 24.16), "v2": (-20.55, 45.41), "v3": (-43.67, 10.94), "v4": (2.92, 25.5), "zeta": (-129.19, 54.27)},
            "_OTHER": {"alpha": (-175.12, 167.85), "beta": (-173.74, 220.12), "chi": (-228.8, 179.38), "delta": (71.75, 168.32), "epsilon": (-200.9, 27.16), "eta": (-239.11, 246.21), "gamma": (-150.06, 173.26), "v1": (-33.53, 49.54), "v2": (-54.38, 46.67), "v3": (-46.52, 42.33), "v4": (-13.3, 36.95), "zeta": (-168.5, 161.46)},
        },
        "tSH": {
            "A-A": {"alpha": (-110.77, 45.43), "beta": (-51.83, 216.71), "chi": (-217.93, -7.12), "delta": (69.97, 116.88), "epsilon": (-188.33, -37.41), "eta": (-77.42, 216.08), "gamma": (-9.97, 121.65), "v1": (-36.96, 9.77), "v2": (-3.62, 49.74), "v3": (-46.44, -3.14), "v4": (8.48, 25.2), "zeta": (-138.83, 20.12)},
            "A-G": {"alpha": (-109.52, 9.9), "beta": (-49.71, 223.16), "chi": (-210.97, -59.06), "delta": (68.59, 108.17), "epsilon": (-172.82, -81.03), "eta": (1.91, 230.49), "gamma": (1.49, 105.29), "v1": (-37.22, 0.69), "v2": (6.79, 50.42), "v3": (-47.36, -11.46), "v4": (11.06, 26.11), "zeta": (-124.44, 29.34)},
            "A-U": {"alpha": (-111.82, 63.52), "beta": (-72.76, 210.07), "chi": (-214.32, -33.87), "delta": (69.09, 104.84), "epsilon": (-165.28, -81.31), "eta": (-36.08, 215.34), "gamma": (-53.07, 121.29), "v1": (-35.56, -1.76), "v2": (10.03, 49.26), "v3": (-47.15, -14.38), "v4": (12.39, 27.11), "zeta": (-128.62, -14.64)},
            "G-G": {"alpha": (-109.85, 31.36), "beta": (-39.82, 225.24), "chi": (-220.26, -38.72), "delta": (69.28, 104.16), "epsilon": (-173.32, -81.33), "eta": (-14.56, 224.11), "gamma": (-1.34, 110.89), "v1": (-36.91, -3.28), "v2": (11.37, 49.86), "v3": (-46.72, -15.17), "v4": (12.39, 25.6), "zeta": (-116.8, 0.75)},
            "_OTHER": {"alpha": (-99.35, 138.46), "beta": (-179.18, 242.27), "chi": (-224.25, 179.31), "delta": (-177.75, 159.25), "epsilon": (-186.22, 170.52), "eta": (-178.16, 240.18), "gamma": (-190.12, 169.9), "v1": (-35.54, 41.13), "v2": (-45.62, 48.38), "v3": (-44.34, 44.54), "v4": (-32.26, 28.62), "zeta": (-172.68, 173.04)},
        },
        "tSS": {
            "A-G": {"alpha": (-111.57, 61.95), "beta": (-86.66, 202.31), "chi": (-223.35, -6.91), "delta": (70.91, 99.51), "epsilon": (-167.44, -98.75), "eta": (-14.92, 220.65), "gamma": (-39.27, 126.01), "v1": (-35.26, -8.55), "v2": (17.19, 48.12), "v3": (-45.61, -19.54), "v4": (13.14, 25.88), "zeta": (-120.9, 12.9)},
            "G-G": {"alpha": (-118.56, 88.29), "beta": (-172.46, 125.05), "chi": (-189.67, -47.47), "delta": (77.67, 140.7), "epsilon": (-165.76, -79.6), "eta": (-65.39, 188.82), "gamma": (-21.31, 108.24), "v1": (-28.81, 32.8), "v2": (-30.18, 40.7), "v3": (-39.36, 18.38), "v4": (1.23, 22.99), "zeta": (-114.69, 97.49)},
            "_OTHER": {"alpha": (-142.05, 213.87), "beta": (-109.78, 221.19), "chi": (-223.1, 107.44), "delta": (67.44, 156.97), "epsilon": (-178.26, 121.85), "eta": (-176.15, 224.98), "gamma": (-170.06, 144.97), "v1": (-35.53, 34.69), "v2": (-37.94, 49.82), "v3": (-48.07, 28.93), "v4": (-7.36, 32.69), "zeta": (-110.17, 140.99)},
        },
        "tSW": {
            "A-A": {"alpha": (-120.1, 32.99), "beta": (-9.2, 219.88), "chi": (-211.32, -64.48), "delta": (68.29, 104.83), "epsilon": (-168.55, -95.43), "eta": (-12.06, 212.2), "gamma": (-60.19, 117.35), "v1": (-37.67, -4.08), "v2": (11.51, 50.92), "v3": (-47.68, -14.61), "v4": (11.28, 26.1), "zeta": (-102.38, -28.02)},
            "A-C": {"alpha": (-110.36, 59.88), "beta": (-99.88, 196.05), "chi": (-208.19, -58.2), "delta": (69.76, 122.76), "epsilon": (-172.74, -82.41), "eta": (-21.77, 212.13), "gamma": (-29.0, 122.0), "v1": (-35.22, 14.46), "v2": (-9.47, 48.88), "v3": (-46.7, 1.95), "v4": (6.36, 26.45), "zeta": (-125.39, 3.0)},
            "A-G": {"alpha": (-111.5, 72.96), "beta": (-84.88, 203.97), "chi": (-214.53, -8.17), "delta": (69.16, 120.75), "epsilon": (-165.64, -84.41), "eta": (-40.64, 206.97), "gamma": (-29.46, 130.84), "v1": (-36.5, 10.81), "v2": (-6.3, 49.88), "v3": (-47.19, 0.32), "v4": (5.49, 26.34), "zeta": (-128.57, 20.59)},
            "A-U": {"alpha": (-108.4, 49.68), "beta": (-135.94, 170.13), "chi": (-193.46, -62.68), "delta": (70.76, 125.97), "epsilon": (-170.91, -80.63), "eta": (-100.36, 181.96), "gamma": (-30.72, 110.72), "v1": (-35.25, 19.84), "v2": (-14.73, 48.39), "v3": (-45.85, 5.45), "v4": (6.03, 25.75), "zeta": (-138.98, 32.76)},
            "C-G": {"alpha": (-105.08, 54.11), "beta": (-110.16, 195.59), "chi": (-197.7, -50.55), "delta": (69.72, 126.01), "epsilon": (-168.33, -91.02), "eta": (-21.89, 201.47), "gamma": (-51.75, 123.11), "v1": (-34.18, 16.94), "v2": (-12.8, 48.15), "v3": (-46.52, 5.06), "v4": (4.72, 26.98), "zeta": (-110.8, 56.23)},
            "G-U": {"alpha": (-97.41, 93.32), "beta": (-143.01, 167.66), "chi": (-192.38, 21.03), "delta": (68.41, 123.31), "epsilon": (-172.6, -74.93), "eta": (-51.97, 198.84), "gamma": (-46.57, 143.14), "v1": (-34.31, 14.36), "v2": (-9.53, 48.62), "v3": (-47.25, 2.2), "v4": (5.88, 27.82), "zeta": (-138.12, 20.07)},
            "_OTHER": {"alpha": (-173.72, 158.91), "beta": (-223.15, 235.94), "chi": (-228.98, 178.4), "delta": (70.8, 150.95), "epsilon": (-174.72, -9.36), "eta": (-169.73, 216.87), "gamma": (-207.28, 178.82), "v1": (-35.05, 45.11), "v2": (-41.85, 48.61), "v3": (-48.42, 27.53), "v4": (-11.64, 38.01), "zeta": (-177.82, 152.34)},
        },
        "tWH": {
            "A-A": {"alpha": (-106.95, 62.73), "beta": (-99.33, 194.13), "chi": (-200.86, 30.72), "delta": (72.86, 132.87), "epsilon": (-161.39, -76.34), "eta": (-114.68, 178.03), "gamma": (-26.66, 120.18), "v1": (-32.94, 24.86), "v2": (-21.28, 45.55), "v3": (-43.47, 11.41), "v4": (3.1, 24.81), "zeta": (-130.06, 11.61)},
            "A-C": {"alpha": (-115.73, 50.75), "beta": (-112.64, 187.3), "chi": (-211.05, -4.29), "delta": (73.36, 134.17), "epsilon": (-167.57, -76.38), "eta": (-122.69, 173.0), "gamma": (-50.66, 127.09), "v1": (-31.28, 25.52), "v2": (-22.27, 44.4), "v3": (-43.22, 12.43), "v4": (2.48, 25.61), "zeta": (-126.36, 58.19)},
            "A-G": {"alpha": (-110.33, 53.38), "beta": (-75.82, 197.47), "chi": (-207.99, -10.53), "delta": (68.66, 115.32), "epsilon": (-167.41, -75.23), "eta": (-65.79, 218.03), "gamma": (-21.89, 130.95), "v1": (-36.93, 7.24), "v2": (-0.9, 50.19), "v3": (-47.33, -5.05), "v4": (8.51, 26.44), "zeta": (-136.61, 23.49)},
            "A-U": {"alpha": (-119.76, 37.99), "beta": (-103.93, 184.07), "chi": (-212.6, 9.78), "delta": (70.87, 127.16), "epsilon": (-171.13, -72.07), "eta": (-130.47, 176.82), "gamma": (-44.31, 136.7), "v1": (-33.88, 19.87), "v2": (-15.06, 47.13), "v3": (-45.28, 6.09), "v4": (5.03, 26.43), "zeta": (-133.98, 26.08)},
            "G-G": {"alpha": (-106.9, 56.06), "beta": (-116.82, 188.76), "chi": (-204.19, -34.93), "delta": (69.82, 126.28), "epsilon": (-166.16, -81.9), "eta": (-68.52, 212.08), "gamma": (-44.77, 127.37), "v1": (-35.04, 19.45), "v2": (-14.13, 48.55), "v3": (-46.37, 4.91), "v4": (6.26, 26.52), "zeta": (-127.21, 72.95)},
            "_OTHER": {"alpha": (-178.98, 181.87), "beta": (-206.26, 241.2), "chi": (-230.05, 171.53), "delta": (12.16, 160.45), "epsilon": (-209.9, 77.29), "eta": (-212.79, 241.27), "gamma": (-142.42, 176.24), "v1": (-43.4, 49.36), "v2": (-50.46, 67.36), "v3": (-69.42, 35.67), "v4": (-11.32, 46.85), "zeta": (-187.81, 169.75)},
        },
        "tWS": {
            "A-C": {"alpha": (-108.17, 51.9), "beta": (-91.6, 204.65), "chi": (-208.81, -56.99), "delta": (69.13, 108.34), "epsilon": (-158.6, -111.41), "eta": (-29.61, 213.67), "gamma": (-9.13, 111.32), "v1": (-36.82, 1.1), "v2": (6.42, 50.23), "v3": (-47.44, -11.22), "v4": (11.04, 26.5), "zeta": (-117.28, 8.24)},
            "A-G": {"alpha": (-115.04, 60.77), "beta": (-99.43, 196.62), "chi": (-208.8, -27.95), "delta": (69.65, 120.86), "epsilon": (-167.32, -80.82), "eta": (-75.86, 195.69), "gamma": (-35.66, 120.75), "v1": (-36.1, 13.48), "v2": (-7.93, 49.27), "v3": (-46.59, 0.47), "v4": (6.88, 26.17), "zeta": (-121.32, 22.93)},
            "_OTHER": {"alpha": (-177.42, 173.3), "beta": (-178.65, 226.53), "chi": (-224.66, 177.01), "delta": (49.39, 167.24), "epsilon": (-166.57, 178.84), "eta": (-178.5, 209.62), "gamma": (-165.85, 174.41), "v1": (-40.85, 41.96), "v2": (-41.86, 58.95), "v3": (-61.91, 39.03), "v4": (-22.16, 45.3), "zeta": (-174.04, 142.81)},
        },
        "tWW": {
            "A-A": {"alpha": (-117.42, 79.67), "beta": (-126.87, 181.8), "chi": (-204.87, -6.51), "delta": (73.63, 133.04), "epsilon": (-168.44, -71.38), "eta": (-106.62, 162.29), "gamma": (-46.57, 130.78), "v1": (-32.72, 24.55), "v2": (-20.91, 44.84), "v3": (-42.64, 11.25), "v4": (2.6, 24.55), "zeta": (-130.39, 38.61)},
            "A-C": {"alpha": (-120.1, 79.54), "beta": (-91.35, 195.63), "chi": (-193.37, 36.76), "delta": (79.53, 141.34), "epsilon": (-157.92, -69.83), "eta": (-114.1, 142.8), "gamma": (-58.49, 118.51), "v1": (-27.9, 31.55), "v2": (-29.8, 38.94), "v3": (-37.58, 19.06), "v4": (-0.82, 22.3), "zeta": (-114.2, 88.4)},
            "A-U": {"alpha": (-115.05, 59.53), "beta": (-127.14, 174.51), "chi": (-180.54, 43.06), "delta": (73.5, 133.65), "epsilon": (-167.66, -69.73), "eta": (-137.33, 148.67), "gamma": (-46.23, 131.13), "v1": (-29.06, 27.12), "v2": (-22.86, 42.88), "v3": (-42.98, 11.99), "v4": (3.56, 27.26), "zeta": (-126.23, 81.39)},
            "C-G": {"alpha": (-122.5, 65.21), "beta": (-143.79, 156.97), "chi": (-198.11, 1.68), "delta": (75.54, 135.97), "epsilon": (-172.66, -67.14), "eta": (-78.98, 190.64), "gamma": (-48.08, 130.82), "v1": (-30.5, 26.92), "v2": (-24.22, 42.72), "v3": (-41.2, 14.31), "v4": (1.33, 24.16), "zeta": (-122.65, 87.53)},
            "G-U": {"alpha": (-107.49, 89.33), "beta": (-134.35, 170.26), "chi": (-196.51, -14.64), "delta": (71.94, 133.04), "epsilon": (-167.75, -78.02), "eta": (-51.44, 194.83), "gamma": (-74.23, 129.78), "v1": (-31.38, 24.56), "v2": (-21.03, 45.25), "v3": (-44.53, 11.32), "v4": (2.97, 26.95), "zeta": (-123.99, 74.08)},
            "U-U": {"alpha": (-109.59, 30.61), "beta": (-140.89, 160.2), "chi": (-191.54, -38.75), "delta": (80.59, 145.64), "epsilon": (-158.36, -74.49), "eta": (-159.36, 123.25), "gamma": (-29.59, 107.08), "v1": (-24.74, 36.63), "v2": (-35.04, 36.58), "v3": (-36.51, 22.66), "v4": (-0.57, 22.53), "zeta": (-101.62, 68.43)},
            "_OTHER": {"alpha": (-170.33, 155.98), "beta": (-227.46, 234.13), "chi": (-214.72, 182.6), "delta": (67.26, 172.17), "epsilon": (-178.6, 81.52), "eta": (-208.2, 245.28), "gamma": (-175.84, 213.58), "v1": (-37.08, 56.07), "v2": (-59.37, 50.71), "v3": (-48.07, 44.1), "v4": (-25.98, 41.41), "zeta": (-190.39, 179.58)},
        },
        "_OTHER": {
            "alpha": (-178.45, 179.71), "beta": (-232.47, 230.47), "chi": (-230.41, 179.43), "delta": (8.26, 163.54), "epsilon": (-182.96, 169.53), "eta": (-193.33, 232.77), "gamma": (-182.6, 182.0), "v1": (-65.87, 46.73), "v2": (-43.89, 50.73), "v3": (-51.56, 39.51), "v4": (-21.51, 47.16), "zeta": (-181.73, 185.21),
        },
    }

            # Higher = catches more