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
    MAX_HOTSPOT_SCORE = 70.0            # Higher = catches more