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

        # MODERATE - Chi glycosidic bond angle in unexpected conformation
        # Full penalty if both residues wrong; half if only one
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







    # ===== BACKBONE TORSION THRESHOLDS (2.5-97.5 pct; unimodal or bimodal GMM) =====
    # angle -> (min, max) or [(m1,M1),(m2,M2)] for bimodal
    # Only edge/bp_type with >=1000 data; edges with '.', '--', or <1k go to _OTHER
    # Generated from torsion_thresholds_analysis.json
    TORSION_THRESHOLDS_BY_EDGE_BY_BASE_PAIR = {
        "cHW": {
            "A-G": {"alpha": [(-170.65, -29.36), (28.63, 172.01)], "beta": (-177.77, 178.31), "chi": [(-179.25, -88.74), (-85.23, 179.41)], "delta": [(74.86, 96.12), (135.35, 154.03)], "epsilon": (-168.12, -65.41), "eta": [(14.79, 178.99), (-179.6, -118.27)], "gamma": [(5.46, 177.57), (-179.18, -48.68)], "v1": [(-33.29, -13.16), (23.62, 42.66)], "v2": [(26.47, 41.91), (-39.87, -28.35)], "v3": (-40.5, 27.18), "v4": (-6.54, 29.76), "zeta": [(51.12, 160.05), (-133.88, -37.62)]},
            "A-U": {"alpha": [(-175.79, -16.05), (29.18, 177.03)], "beta": (-177.89, 178.07), "chi": [(-177.51, -79.82), (-49.26, 178.72)], "delta": [(73.78, 90.29), (115.94, 157.24)], "epsilon": (-166.71, 8.5), "eta": [(-179.17, -62.1), (19.59, 177.84)], "gamma": [(-177.94, -36.79), (-14.36, 174.8)], "v1": (-29.87, 42.46), "v2": [(27.47, 41.41), (-40.61, -18.69)], "v3": (-40.87, 29.75), "v4": (-9.18, 29.98), "zeta": [(-143.61, -28.51), (40.47, 156.0)]},
            "G-G": {"alpha": [(37.23, 176.54), (-168.55, -30.08)], "beta": (-178.18, 178.54), "chi": [(-177.8, -89.1), (-79.81, 179.04)], "delta": [(74.95, 91.2), (131.95, 156.13)], "epsilon": (-169.29, 64.71), "eta": [(147.74, 179.03), (-178.77, 130.93)], "gamma": [(-179.74, -130.73), (-58.35, 177.18)], "v1": (-32.34, 39.91), "v2": [(27.92, 42.59), (-40.68, -23.85)], "v3": (-41.1, 28.42), "v4": (-8.68, 27.83), "zeta": [(25.63, 153.0), (-132.5, -15.28)]},
            "_OTHER": {"alpha": (-176.4, 171.19), "beta": (-179.4, 179.41), "chi": (-177.83, 178.9), "delta": (-10.8, 154.72), "epsilon": (-178.47, 99.23), "eta": (-178.55, 179.45), "gamma": (-177.37, 174.45), "v1": (-63.01, 48.12), "v2": (-40.43, 45.79), "v3": (-44.79, 32.51), "v4": (-41.63, 39.68), "zeta": (-168.44, 161.13)},
        },
        "cSH": {
            "A-A": {"alpha": [(-172.76, 75.5), (124.42, 178.39)], "beta": [(100.37, 177.84), (-178.27, -105.19)], "chi": [(-178.57, -77.05), (46.06, 179.67)], "delta": [(126.22, 153.23), (76.32, 92.34)], "epsilon": (-163.72, -69.22), "eta": [(53.74, 177.56), (-178.7, -60.49)], "gamma": [(-3.66, 168.95), (-178.77, -44.43)], "v1": (-34.83, 37.35), "v2": [(-37.94, -21.26), (26.73, 43.65)], "v3": (-38.7, 28.61), "v4": (-9.17, 25.63), "zeta": [(-175.71, -39.79), (92.08, 178.49)]},
            "A-C": {"alpha": [(-173.77, 63.13), (123.85, 178.36)], "beta": (-178.28, 178.23), "chi": [(-176.1, -84.47), (4.89, 179.55)], "delta": [(74.15, 92.92), (129.87, 155.65)], "epsilon": (-166.1, 83.87), "eta": [(-178.06, -102.39), (41.98, 177.64)], "gamma": [(20.57, 157.98), (-175.52, -17.31)], "v1": (-31.95, 37.67), "v2": [(28.3, 42.06), (-38.74, -22.56)], "v3": (-39.93, 28.1), "v4": (-7.87, 26.84), "zeta": [(-166.19, -38.99), (68.52, 174.98)]},
            "A-G": {"alpha": [(-171.34, -25.5), (35.85, 176.71)], "beta": (-177.88, 178.06), "chi": [(-178.42, -70.26), (29.76, 179.69)], "delta": [(75.36, 93.22), (121.34, 155.78)], "epsilon": (-171.44, 65.32), "eta": [(52.36, 177.23), (-178.66, -21.77)], "gamma": [(-59.85, 173.12), (-179.26, -130.64)], "v1": (-32.85, 39.29), "v2": [(22.66, 42.09), (-39.4, -4.39)], "v3": (-40.21, 28.32), "v4": (-9.99, 29.4), "zeta": [(-169.59, -28.03), (34.54, 175.68)]},
            "A-U": {"alpha": [(48.02, 177.38), (-174.64, -3.24)], "beta": (-177.78, 178.03), "chi": [(-178.93, -83.4), (33.5, 179.68)], "delta": [(76.43, 90.96), (111.2, 152.86)], "epsilon": (-166.39, 75.2), "eta": [(-179.34, -57.72), (51.04, 178.58)], "gamma": [(-156.98, 123.67), (161.14, 179.23)], "v1": (-31.85, 38.41), "v2": [(25.4, 40.49), (-38.88, -19.0)], "v3": (-39.04, 26.88), "v4": (-6.17, 27.9), "zeta": [(-172.89, -43.52), (45.24, 176.84)]},
            "C-U": {"alpha": [(-169.2, -28.86), (30.5, 177.8)], "beta": (-178.73, 178.98), "chi": (-174.28, -82.68), "delta": [(129.98, 153.26), (74.54, 87.08)], "epsilon": [(-176.98, -79.73), (41.33, 179.75)], "eta": [(-177.91, 58.9), (103.15, 177.88)], "gamma": [(21.83, 159.88), (-177.48, -19.81)], "v1": (-29.24, 38.98), "v2": [(-39.02, -28.12), (29.3, 41.32)], "v3": (-40.58, 27.11), "v4": (-5.35, 28.03), "zeta": [(16.52, 177.91), (-175.43, -15.65)]},
            "G-G": {"alpha": [(27.77, 177.17), (-174.29, -27.85)], "beta": [(-178.5, -79.57), (97.38, 178.27)], "chi": [(-177.78, -70.06), (-44.46, 179.35)], "delta": [(118.04, 158.73), (75.06, 90.68)], "epsilon": [(-171.82, -71.91), (12.48, 178.77)], "eta": [(24.71, 177.38), (-178.31, -66.65)], "gamma": [(7.29, 177.1), (-179.03, -39.21)], "v1": (-31.05, 40.96), "v2": [(-41.05, -16.23), (27.6, 41.93)], "v3": (-39.75, 31.46), "v4": (-12.05, 27.93), "zeta": [(-172.25, -36.85), (44.74, 177.2)]},
            "G-U": {"alpha": [(-166.18, -34.76), (25.83, 177.03)], "beta": [(71.32, 176.77), (-178.46, -81.33)], "chi": (-178.62, 178.4), "delta": [(128.56, 155.63), (75.8, 92.81)], "epsilon": [(-173.59, -81.55), (7.88, 179.44)], "eta": [(128.13, 178.23), (-178.45, 107.29)], "gamma": [(-179.47, -64.9), (-40.54, 177.4)], "v1": [(26.29, 45.69), (-30.91, -12.66)], "v2": [(-43.28, -23.46), (27.43, 40.15)], "v3": (-38.99, 29.05), "v4": (-6.06, 27.85), "zeta": [(-152.13, -42.33), (69.05, 172.21)]},
            "_OTHER": {"alpha": (-172.87, 177.74), "beta": (-178.8, 178.89), "chi": (-179.6, 179.72), "delta": (-33.74, 164.18), "epsilon": (-173.61, 164.62), "eta": (-175.92, 176.13), "gamma": (-171.16, 177.21), "v1": (-38.62, 42.89), "v2": (-41.96, 54.65), "v3": (-49.58, 40.23), "v4": (-32.78, 38.81), "zeta": (-174.07, 172.8)},
        },
        "cSS": {
            "A-C": {"alpha": [(-161.92, -35.4), (31.14, 177.37)], "beta": (-178.75, 178.77), "chi": (-173.01, -109.61), "delta": (75.27, 102.49), "epsilon": (-164.59, -91.34), "eta": [(-179.35, 142.33), (151.65, 178.78)], "gamma": [(30.68, 169.09), (-178.59, -57.05)], "v1": (-30.3, 0.27), "v2": (3.28, 41.18), "v3": (-41.06, -13.96), "v4": (2.67, 28.27), "zeta": (-114.45, 87.44)},
            "A-G": {"alpha": [(25.45, 174.29), (-163.05, -34.84)], "beta": (-178.25, 178.43), "chi": (-178.29, 177.55), "delta": [(119.91, 157.63), (75.08, 90.74)], "epsilon": (-166.63, -65.4), "eta": [(147.2, 178.82), (-179.06, 130.86)], "gamma": [(20.99, 173.9), (-178.43, -46.74)], "v1": (-31.92, 38.71), "v2": [(-40.15, -25.12), (25.16, 41.93)], "v3": (-41.1, 26.26), "v4": (-4.96, 30.18), "zeta": [(24.59, 175.65), (-155.37, -42.88)]},
            "C-G": {"alpha": [(-139.64, -38.35), (13.12, 177.32)], "beta": (-178.89, 179.11), "chi": (-174.86, -72.8), "delta": [(75.87, 89.93), (134.73, 156.2)], "epsilon": (-167.77, -87.07), "eta": [(-179.58, 142.62), (150.87, 178.95)], "gamma": (-148.45, 156.33), "v1": (-31.22, 35.27), "v2": [(29.33, 41.43), (-40.04, -22.68)], "v3": (-40.3, 24.33), "v4": (-2.52, 27.8), "zeta": (-166.27, 139.72)},
            "G-U": {"alpha": [(-149.32, -41.48), (52.89, 178.59)], "beta": (-179.07, 179.06), "chi": (-176.28, -59.76), "delta": [(74.62, 90.2), (114.31, 155.39)], "epsilon": [(-169.05, -114.31), (-106.03, 175.96)], "eta": [(156.3, 179.28), (-179.85, 143.81)], "gamma": (-93.86, 159.71), "v1": (-29.47, 36.54), "v2": [(28.54, 40.75), (-39.3, -24.4)], "v3": (-41.6, 24.99), "v4": (-3.2, 29.49), "zeta": (-123.06, 93.65)},
            "_OTHER": {"alpha": (-177.19, 177.94), "beta": (-178.81, 179.42), "chi": (-178.3, 178.95), "delta": (67.42, 151.28), "epsilon": (-176.17, 162.26), "eta": (-179.5, 179.87), "gamma": (-177.65, 176.42), "v1": (-34.27, 48.65), "v2": (-38.69, 41.38), "v3": (-44.44, 27.6), "v4": (-19.14, 49.85), "zeta": (-147.71, 130.14)},
        },
        "cSW": {
            "A-A": {"alpha": [(-172.92, -39.13), (38.66, 177.6)], "beta": (-178.11, 178.28), "chi": [(-178.96, -145.05), (-136.03, 179.5)], "delta": [(75.84, 88.66), (126.06, 152.97)], "epsilon": (-168.64, -92.16), "eta": [(-178.27, 144.91), (150.4, 178.06)], "gamma": [(31.22, 169.55), (-178.79, -62.14)], "v1": (-31.23, 36.58), "v2": [(28.98, 41.22), (-38.87, -28.21)], "v3": (-40.9, 24.62), "v4": (-3.44, 28.44), "zeta": [(-140.71, -27.63), (5.66, 175.05)]},
            "A-C": {"alpha": [(-167.56, -37.38), (29.78, 176.45)], "beta": (-178.24, 178.38), "chi": (-178.55, 178.43), "delta": [(75.78, 90.77), (128.68, 155.48)], "epsilon": (-166.84, -86.91), "eta": [(128.69, 177.87), (-178.18, 122.38)], "gamma": [(16.54, 170.49), (-179.44, -34.25)], "v1": (-32.16, 37.38), "v2": [(28.81, 41.35), (-40.23, -20.85)], "v3": (-40.22, 26.16), "v4": (-4.55, 27.26), "zeta": (-125.42, 48.42)},
            "A-G": {"alpha": [(-158.68, -33.52), (21.17, 174.37)], "beta": (-178.42, 178.58), "chi": [(-178.76, -132.28), (-124.4, 179.4)], "delta": [(74.52, 90.72), (132.21, 157.25)], "epsilon": (-168.05, -68.18), "eta": [(26.73, 177.91), (-178.86, -16.49)], "gamma": [(-179.12, -50.85), (32.15, 177.12)], "v1": [(26.22, 43.13), (-32.61, -14.75)], "v2": [(-40.9, -23.56), (28.56, 42.05)], "v3": (-41.15, 28.99), "v4": (-7.53, 28.16), "zeta": [(-169.75, -38.1), (43.72, 175.92)]},
            "A-U": {"alpha": [(33.4, 178.39), (-175.36, -36.37)], "beta": (-178.78, 178.77), "chi": [(-177.42, -94.42), (-57.12, 179.72)], "delta": [(73.43, 90.4), (131.55, 156.22)], "epsilon": (-168.02, -76.78), "eta": [(148.26, 177.48), (-176.79, 135.24)], "gamma": [(25.61, 175.18), (-178.85, -39.37)], "v1": [(-33.09, -11.48), (26.61, 42.46)], "v2": [(26.19, 42.39), (-41.0, -25.81)], "v3": (-41.69, 27.71), "v4": (-5.53, 28.93), "zeta": [(-141.24, -38.39), (25.19, 176.55)]},
            "C-U": {"alpha": [(-178.56, -30.69), (25.26, 179.21)], "beta": (-178.39, 178.33), "chi": [(-174.69, -119.84), (-107.48, 179.18)], "delta": [(76.78, 89.48), (136.22, 162.79)], "epsilon": (-166.46, -77.53), "eta": [(-179.73, -25.93), (17.72, 178.08)], "gamma": [(-179.79, -129.87), (-43.79, 176.69)], "v1": (-28.13, 38.17), "v2": [(25.33, 39.88), (-39.64, -30.06)], "v3": (-39.63, 28.91), "v4": (-8.29, 30.8), "zeta": [(-160.28, -42.79), (86.22, 176.0)]},
            "G-U": {"alpha": [(32.61, 174.47), (-169.37, -36.72)], "beta": (-177.19, 177.64), "chi": [(-178.05, -92.32), (-50.35, 179.61)], "delta": [(129.43, 155.4), (74.12, 90.43)], "epsilon": (-170.42, 55.58), "eta": [(30.31, 177.13), (-179.52, -15.19)], "gamma": [(-179.13, -53.64), (7.36, 173.7)], "v1": (-31.47, 39.46), "v2": [(-40.69, -23.82), (29.97, 41.72)], "v3": (-41.04, 27.34), "v4": (-5.82, 28.16), "zeta": [(-155.49, -32.91), (24.41, 172.65)]},
            "_OTHER": {"alpha": (-179.28, 178.7), "beta": (-178.94, 178.77), "chi": (-177.81, 177.77), "delta": (75.27, 159.68), "epsilon": (-176.25, 164.36), "eta": (-175.13, 179.7), "gamma": (-178.88, 177.24), "v1": (-35.39, 40.05), "v2": (-42.48, 46.8), "v3": (-45.1, 39.12), "v4": (-25.43, 34.97), "zeta": (-169.85, 166.34)},
        },
        "cWH": {
            "A-C": {"alpha": [(66.71, 173.26), (-162.66, -35.69)], "beta": [(88.12, 178.98), (-178.81, -89.03)], "chi": [(-178.5, -121.08), (-98.78, 179.54)], "delta": [(72.49, 90.74), (137.49, 156.53)], "epsilon": (-170.25, 156.9), "eta": [(52.94, 179.17), (-179.57, 13.49)], "gamma": [(-77.65, 178.42), (-179.44, -136.57)], "v1": [(-34.23, -16.57), (29.68, 42.17)], "v2": [(29.95, 41.93), (-40.84, -31.09)], "v3": (-41.86, 26.62), "v4": (-4.47, 29.1), "zeta": [(35.78, 176.3), (-172.72, -36.97)]},
            "A-G": {"alpha": [(-174.93, -31.52), (21.56, 174.61)], "beta": (-178.33, 178.47), "chi": [(-27.27, 179.38), (-178.4, -77.44)], "delta": [(126.83, 158.54), (75.39, 91.69)], "epsilon": (-165.82, 60.52), "eta": [(-178.91, 85.78), (125.87, 179.29)], "gamma": [(8.36, 173.79), (-178.16, -35.41)], "v1": [(26.42, 45.71), (-33.78, -7.43)], "v2": [(-43.42, -18.85), (26.11, 41.24)], "v3": (-39.52, 29.06), "v4": (-9.57, 28.38), "zeta": [(-152.95, -37.25), (34.12, 161.22)]},
            "A-U": {"alpha": [(27.15, 176.36), (-171.49, -35.35)], "beta": (-177.68, 178.17), "chi": [(-177.53, -81.41), (-66.47, 179.47)], "delta": [(131.54, 154.95), (74.05, 89.81)], "epsilon": (-169.58, -31.09), "eta": [(-177.67, 134.57), (149.38, 178.07)], "gamma": [(13.39, 175.7), (-178.74, -45.81)], "v1": (-29.62, 41.23), "v2": [(-40.51, -26.59), (28.24, 41.43)], "v3": (-41.02, 28.65), "v4": (-7.92, 29.08), "zeta": [(30.03, 172.89), (-140.5, -39.07)]},
            "C-G": {"alpha": [(21.83, 176.42), (-170.98, -27.75)], "beta": (-177.27, 178.07), "chi": [(-178.86, -135.01), (-127.46, 179.02)], "delta": [(70.66, 94.41), (124.87, 154.79)], "epsilon": (-171.74, 153.39), "eta": [(-178.65, 133.87), (143.84, 178.61)], "gamma": [(-178.68, -56.36), (-4.32, 175.0)], "v1": (-36.99, 38.75), "v2": [(22.76, 45.46), (-39.76, -20.22)], "v3": (-42.27, 27.73), "v4": (-7.47, 30.69), "zeta": [(-163.38, -29.38), (29.6, 174.36)]},
            "G-G": {"alpha": [(44.67, 177.35), (-173.98, -34.32)], "beta": (-177.79, 178.27), "chi": [(-178.01, -71.17), (-6.81, 179.55)], "delta": [(73.87, 91.36), (134.41, 159.07)], "epsilon": (-168.5, 1.37), "eta": [(146.64, 178.75), (-178.1, 102.69)], "gamma": [(11.4, 175.52), (-178.72, -47.72)], "v1": [(-33.02, -15.55), (26.01, 42.66)], "v2": [(28.56, 42.84), (-41.88, -27.92)], "v3": [(-41.87, -28.12), (14.38, 34.46)], "v4": (-9.5, 28.35), "zeta": [(-131.59, -45.14), (38.28, 168.04)]},
            "G-U": {"alpha": [(28.1, 174.63), (-164.5, -35.82)], "beta": (-177.76, 178.37), "chi": [(-128.53, 179.61), (-178.4, -137.89)], "delta": [(74.0, 89.86), (106.89, 156.17)], "epsilon": (-168.72, 65.26), "eta": [(150.88, 178.69), (-179.37, 138.4)], "gamma": [(22.79, 177.3), (-178.76, -44.4)], "v1": (-34.11, 37.01), "v2": [(29.29, 42.72), (-39.96, 10.58)], "v3": (-40.7, 25.48), "v4": (-4.04, 29.53), "zeta": [(-159.18, -39.72), (38.54, 172.3)]},
            "_OTHER": {"alpha": (-169.4, 173.8), "beta": (-178.07, 179.56), "chi": (-178.44, 176.14), "delta": (71.13, 172.6), "epsilon": (-175.97, 174.22), "eta": (-179.58, 179.8), "gamma": (-177.9, 177.72), "v1": (-33.89, 45.34), "v2": (-47.96, 46.83), "v3": (-46.55, 45.34), "v4": (-31.2, 45.8), "zeta": (-171.59, 161.69)},
        },
        "cWS": {
            "A-A": {"alpha": [(-173.45, -38.21), (36.32, 177.47)], "beta": (-178.29, 178.5), "chi": [(65.41, 179.8), (-178.88, -75.48)], "delta": (76.06, 150.7), "epsilon": (-167.69, -66.66), "eta": [(145.74, 177.95), (-177.92, 108.47)], "gamma": [(-13.79, 176.86), (-179.45, -80.04)], "v1": (-31.36, 39.38), "v2": [(26.6, 41.04), (-39.73, -19.03)], "v3": (-40.32, 26.8), "v4": (-5.09, 28.55), "zeta": [(-163.93, -45.9), (46.11, 169.63)]},
            "A-C": {"alpha": [(33.1, 175.05), (-165.61, -38.98)], "beta": (-178.32, 178.45), "chi": (-178.3, 178.42), "delta": [(75.09, 91.1), (131.19, 156.09)], "epsilon": (-168.67, -89.45), "eta": [(149.72, 178.41), (-179.0, 139.82)], "gamma": [(165.28, 179.75), (-176.29, 146.06)], "v1": (-31.35, 36.11), "v2": [(27.2, 41.3), (-40.28, -23.16)], "v3": (-40.92, 25.05), "v4": (-3.18, 28.48), "zeta": (-118.01, 64.23)},
            "A-G": {"alpha": [(42.04, 178.08), (-169.76, -30.98)], "beta": (-178.14, 178.08), "chi": (-178.2, 177.86), "delta": [(73.98, 90.59), (127.42, 157.52)], "epsilon": (-169.04, -82.14), "eta": [(151.65, 179.05), (-179.09, 137.79)], "gamma": [(-179.41, -66.22), (26.4, 173.54)], "v1": (-31.71, 37.7), "v2": [(28.29, 41.83), (-39.23, -26.9)], "v3": (-42.03, 24.87), "v4": (-2.84, 29.37), "zeta": [(-157.51, -44.89), (20.7, 179.32)]},
            "A-U": {"alpha": [(21.01, 173.65), (-156.59, -39.24)], "beta": (-177.5, 177.79), "chi": (-177.43, 177.42), "delta": [(76.09, 90.95), (135.58, 154.98)], "epsilon": (-175.47, 170.22), "eta": [(148.69, 178.07), (-178.13, 122.78)], "gamma": [(32.46, 175.19), (-179.45, -58.66)], "v1": (-29.98, 37.87), "v2": [(26.65, 40.93), (-39.5, -30.26)], "v3": (-40.17, 24.49), "v4": (-2.41, 27.7), "zeta": [(-158.79, -47.71), (37.47, 177.68)]},
            "G-U": {"alpha": [(-160.08, -35.86), (26.99, 176.25)], "beta": (-178.19, 177.94), "chi": [(-177.82, -105.92), (-68.01, 179.84)], "delta": [(74.45, 89.19), (130.95, 153.97)], "epsilon": (-168.43, -67.4), "eta": [(-177.09, 127.37), (138.81, 177.7)], "gamma": [(24.68, 171.96), (-179.31, -50.34)], "v1": (-31.08, 38.86), "v2": [(25.79, 41.91), (-40.14, -22.79)], "v3": (-41.07, 26.24), "v4": (-3.68, 30.38), "zeta": [(-165.04, -40.39), (24.95, 171.56)]},
            "_OTHER": {"alpha": (-176.91, 178.34), "beta": (-178.63, 179.81), "chi": (-179.67, 179.19), "delta": (67.61, 164.52), "epsilon": (-172.49, 160.63), "eta": (-179.64, 179.52), "gamma": (-178.94, 174.89), "v1": (-33.67, 37.46), "v2": (-38.09, 46.98), "v3": (-55.73, 40.3), "v4": (-28.57, 43.0), "zeta": (-162.95, 133.71)},
        },
        "cWW": {
            "A-A": {"alpha": [(-143.12, -37.83), (24.78, 176.69)], "beta": (-178.12, 177.6), "chi": (-178.32, 178.39), "delta": [(132.02, 156.11), (75.64, 90.08)], "epsilon": (-171.05, -78.1), "eta": [(144.4, 177.29), (-179.02, 125.17)], "gamma": [(32.23, 172.54), (-179.5, -63.53)], "v1": (-32.98, 38.15), "v2": [(-40.46, -26.97), (29.09, 42.13)], "v3": (-40.5, 25.22), "v4": (-2.56, 27.64), "zeta": [(-137.68, -43.72), (15.28, 173.15)]},
            "A-C": {"alpha": [(-149.24, -38.02), (30.82, 176.18)], "beta": (-177.71, 178.39), "chi": (-178.4, 178.09), "delta": [(103.21, 154.47), (75.7, 89.84)], "epsilon": (-168.02, -86.5), "eta": [(144.21, 178.56), (-179.58, 134.23)], "gamma": [(28.46, 172.55), (-179.2, -82.18)], "v1": (-32.0, 35.6), "v2": (-35.38, 41.48), "v3": (-40.53, 23.76), "v4": (-1.71, 27.35), "zeta": (-130.99, 91.7)},
            "A-DT": {"alpha": [(-152.0, -30.82), (24.88, 176.72)], "beta": (-178.38, 178.55), "chi": (-175.81, 154.24), "delta": (75.19, 149.69), "epsilon": (-175.65, 176.24), "eta": (-177.74, 178.32), "gamma": [(15.62, 171.58), (-178.67, -52.2)], "v1": (-32.53, 39.99), "v2": (-36.91, 43.04), "v3": (-43.1, 27.74), "v4": (-7.68, 36.95), "zeta": (-145.44, 55.18)},
            "A-G": {"alpha": [(35.41, 175.07), (-156.49, -29.42)], "beta": (-177.1, 177.92), "chi": (-178.61, 178.52), "delta": [(116.73, 155.21), (75.25, 90.14)], "epsilon": (-171.18, -69.21), "eta": [(146.2, 178.7), (-179.5, 130.02)], "gamma": [(20.86, 173.34), (-179.23, -62.91)], "v1": [(25.63, 43.03), (-33.57, -16.76)], "v2": [(-40.77, -22.92), (29.81, 42.27)], "v3": (-40.75, 25.68), "v4": (-3.49, 27.48), "zeta": [(-148.99, -40.05), (21.48, 175.77)]},
            "A-PSU": {"alpha": [(-97.49, -44.78), (50.43, 170.43)], "beta": (-177.38, 178.72), "chi": (-178.75, 176.44), "delta": (75.45, 90.4), "epsilon": (-165.42, -113.5), "eta": [(153.35, 178.17), (-179.81, 133.49)], "gamma": [(39.24, 171.13), (-178.67, -118.48)], "v1": (-31.31, -15.05), "v2": (29.1, 41.4), "v3": (-40.82, -28.01), "v4": (11.68, 27.49), "zeta": (-131.35, -57.99)},
            "A-U": {"alpha": [(-141.84, -39.62), (28.39, 176.74)], "beta": (-178.74, 178.86), "chi": [(-177.6, -147.17), (-135.89, 179.69)], "delta": (74.94, 144.67), "epsilon": [(-101.42, 178.88), (-169.91, -111.43)], "eta": [(-179.52, 142.76), (152.24, 177.94)], "gamma": [(-179.34, -75.42), (32.34, 170.79)], "v1": (-31.56, 35.84), "v2": (-34.93, 41.77), "v3": (-41.35, 22.42), "v4": (0.21, 28.0), "zeta": (-127.91, 69.06)},
            "C-C": {"alpha": [(-162.54, -34.72), (30.73, 176.66)], "beta": (-178.63, 178.66), "chi": (-176.26, 162.08), "delta": (75.2, 144.64), "epsilon": (-167.76, -82.74), "eta": [(-179.58, 129.47), (146.28, 178.8)], "gamma": [(25.48, 172.55), (-179.19, -71.84)], "v1": (-30.57, 35.29), "v2": (-34.66, 40.82), "v3": (-41.36, 22.45), "v4": (0.42, 29.68), "zeta": (-143.34, 51.69)},
            "C-DG": {"alpha": [(-148.65, -31.45), (22.28, 176.36)], "beta": (-178.13, 177.92), "chi": (-176.46, 174.11), "delta": (74.81, 151.12), "epsilon": [(-175.35, -100.77), (61.23, 179.7)], "eta": [(142.99, 178.43), (-179.71, 127.44)], "gamma": [(-15.41, 170.32), (-179.18, -66.25)], "v1": (-33.1, 40.49), "v2": (-37.38, 43.04), "v3": (-43.69, 29.02), "v4": (-10.74, 36.32), "zeta": (-140.08, 76.11)},
            "C-G": {"alpha": [(-145.29, -38.37), (32.85, 175.92)], "beta": (-178.65, 178.79), "chi": (-178.12, 177.52), "delta": (74.49, 142.8), "epsilon": [(-169.57, -112.47), (-102.13, 178.86)], "eta": [(150.89, 178.22), (-179.59, 140.28)], "gamma": (-171.52, 172.77), "v1": (-32.17, 34.2), "v2": (-33.57, 42.23), "v3": (-41.79, 21.14), "v4": (1.22, 28.22), "zeta": (-122.08, 75.0)},
            "C-OMG": {"alpha": [(-158.26, -38.41), (18.76, 172.71)], "beta": (-178.28, 178.52), "chi": (-177.75, 177.04), "delta": (74.33, 145.21), "epsilon": (-161.71, -84.05), "eta": [(145.19, 178.1), (-179.45, 125.68)], "gamma": (-163.87, 145.68), "v1": (-31.89, 35.38), "v2": (-35.05, 41.54), "v3": (-41.21, 22.6), "v4": (1.04, 30.89), "zeta": (-118.58, 55.56)},
            "C-U": {"alpha": [(25.54, 175.43), (-156.87, -36.98)], "beta": (-178.68, 178.73), "chi": [(-176.81, -136.61), (-128.95, 179.64)], "delta": [(73.36, 87.88), (118.52, 156.26)], "epsilon": (-168.07, -76.08), "eta": [(151.5, 178.78), (-179.62, 141.21)], "gamma": [(25.76, 171.7), (-178.96, -61.9)], "v1": (-30.39, 37.25), "v2": [(30.01, 41.52), (-39.72, 5.06)], "v3": (-42.05, 24.57), "v4": (-2.43, 29.44), "zeta": [(9.09, 171.92), (-146.17, -41.42)]},
            "DA-DT": {"alpha": [(-142.32, -14.78), (13.39, 176.45)], "beta": (-178.21, 178.17), "chi": (-160.41, -77.14), "delta": (85.22, 161.14), "epsilon": (-178.55, 178.19), "eta": (-178.26, 178.69), "gamma": [(-36.12, 169.1), (-179.08, -59.49)], "v1": (-20.75, 45.32), "v2": (-42.93, 32.0), "v3": (-35.25, 38.17), "v4": (-25.39, 35.47), "zeta": [(-167.79, -69.25), (28.81, 179.15)]},
            "DA-U": {"alpha": [(-157.96, -28.12), (31.83, 175.83)], "beta": (-177.98, 178.33), "chi": (-176.52, -58.13), "delta": (75.35, 150.82), "epsilon": (-175.58, 166.72), "eta": [(142.01, 178.41), (-179.11, 125.72)], "gamma": [(15.54, 173.03), (-177.77, -53.85)], "v1": (-32.44, 41.64), "v2": (-37.69, 41.76), "v3": (-42.11, 29.88), "v4": (-12.02, 33.79), "zeta": (-150.97, 77.18)},
            "DC-G": {"alpha": [(-149.25, -26.6), (23.88, 175.25)], "beta": (-177.77, 178.2), "chi": [(-176.79, -70.79), (56.52, 179.9)], "delta": (74.84, 149.65), "epsilon": (-175.75, 169.71), "eta": [(149.92, 178.09), (-179.66, 140.09)], "gamma": [(-178.7, -63.13), (-8.35, 173.6)], "v1": (-33.76, 40.87), "v2": (-37.72, 43.05), "v3": (-43.53, 26.83), "v4": (-6.1, 37.27), "zeta": [(-147.83, -52.72), (21.32, 170.99)]},
            "G-G": {"alpha": [(-165.52, -28.86), (45.29, 176.21)], "beta": [(-179.09, -90.09), (86.14, 178.71)], "chi": [(-122.91, 179.43), (-179.07, -137.7)], "delta": [(70.34, 91.6), (122.63, 153.47)], "epsilon": (-171.12, 26.79), "eta": [(136.82, 179.29), (-179.65, 117.21)], "gamma": [(19.63, 174.55), (-178.8, -92.87)], "v1": (-37.93, 35.28), "v2": [(29.38, 44.88), (-38.64, 3.43)], "v3": (-43.43, 24.2), "v4": (-4.02, 29.31), "zeta": [(12.21, 172.47), (-144.7, -42.82)]},
            "G-OMC": {"alpha": [(-150.61, -41.0), (41.05, 174.24)], "beta": (-177.09, 178.13), "chi": (-178.8, 178.8), "delta": (74.5, 88.43), "epsilon": (-165.66, -114.76), "eta": [(152.24, 177.89), (-179.34, 142.0)], "gamma": (-169.33, 174.03), "v1": (-30.63, -16.16), "v2": (29.54, 41.05), "v3": (-41.41, -28.86), "v4": (11.06, 28.7), "zeta": [(-100.82, -57.75), (26.91, 74.97)]},
            "G-PSU": {"alpha": [(-102.25, -44.87), (69.04, 174.9)], "beta": (-178.49, 178.35), "chi": (-177.87, 176.46), "delta": (75.27, 89.31), "epsilon": (-164.68, -114.85), "eta": (-174.22, 178.46), "gamma": [(40.72, 162.24), (-178.45, -96.88)], "v1": (-30.58, -11.44), "v2": (25.59, 41.25), "v3": (-41.0, -29.79), "v4": (15.99, 28.77), "zeta": (-120.53, -43.66)},
            "G-U": {"alpha": [(-147.52, -39.6), (39.42, 176.43)], "beta": (-178.36, 178.54), "chi": (-178.34, 177.89), "delta": (74.72, 142.26), "epsilon": (-169.08, -93.67), "eta": [(151.24, 178.46), (-179.67, 138.55)], "gamma": [(35.89, 173.44), (-179.43, -87.62)], "v1": (-31.63, 34.01), "v2": (-33.51, 41.74), "v3": (-41.49, 20.66), "v4": (1.74, 28.36), "zeta": (-134.76, 62.91)},
            "U-U": {"alpha": [(-156.0, -37.51), (26.2, 177.72)], "beta": (-178.74, 178.92), "chi": [(-176.73, -144.37), (-135.26, 179.58)], "delta": (74.66, 143.36), "epsilon": (-168.44, -90.19), "eta": [(151.67, 178.8), (-179.66, 138.19)], "gamma": [(31.84, 169.0), (-179.25, -58.23)], "v1": (-29.95, 35.69), "v2": (-34.67, 40.95), "v3": (-41.08, 21.5), "v4": (1.78, 28.69), "zeta": (-150.47, 76.22)},
            "_OTHER": {"alpha": (-179.75, 178.7), "beta": (-179.97, 179.95), "chi": (-179.63, 179.84), "delta": (-103.83, 169.37), "epsilon": (-179.42, 179.24), "eta": (-179.83, 179.93), "gamma": (-179.69, 179.58), "v1": (-110.56, 49.72), "v2": (-51.45, 60.97), "v3": (-63.01, 55.41), "v4": (-42.48, 45.07), "zeta": (-171.86, 175.56)},
        },
        "tHH": {
            "A-A": {"alpha": [(-175.97, -30.07), (44.99, 178.29)], "beta": (-177.16, 177.46), "chi": [(-175.34, -117.2), (-88.52, 179.47)], "delta": [(133.37, 155.19), (77.12, 91.28)], "epsilon": (-163.06, -64.59), "eta": [(-179.28, -121.57), (-44.85, 177.91)], "gamma": [(29.92, 166.77), (-179.4, -33.28)], "v1": [(29.34, 43.17), (-31.51, -16.11)], "v2": [(-40.44, -28.61), (28.24, 41.3)], "v3": (-38.75, 28.05), "v4": (-5.94, 24.68), "zeta": [(50.51, 176.19), (-171.23, -48.37)]},
            "_OTHER": {"alpha": (-177.71, 179.04), "beta": (-174.56, 178.3), "chi": (-179.01, 177.85), "delta": (74.75, 160.53), "epsilon": (-165.26, 53.04), "eta": (-178.56, 179.66), "gamma": (-175.64, 172.57), "v1": (-32.06, 42.55), "v2": (-41.47, 48.43), "v3": (-43.73, 32.54), "v4": (-18.78, 37.85), "zeta": (-173.9, 172.3)},
        },
        "tHS": {
            "A-A": {"alpha": [(-94.66, -43.48), (-170.16, 173.52)], "beta": (-178.4, 178.26), "chi": (-179.13, 179.03), "delta": [(102.71, 154.75), (76.05, 88.41)], "epsilon": (-170.21, -87.66), "eta": (-178.84, 178.96), "gamma": (-154.76, 164.61), "v1": (-31.04, 36.05), "v2": [(-39.46, -14.35), (31.2, 41.16)], "v3": (-40.05, 24.11), "v4": (-2.46, 26.35), "zeta": [(-152.22, -36.86), (34.82, 172.54)]},
            "A-G": {"alpha": [(-127.88, -33.92), (20.08, 176.51)], "beta": (-178.27, 178.57), "chi": (-177.57, 176.46), "delta": [(75.9, 89.43), (126.6, 155.46)], "epsilon": (-170.08, -30.27), "eta": (-178.2, 178.57), "gamma": (-161.39, 169.36), "v1": (-31.18, 38.6), "v2": [(29.48, 41.34), (-40.38, -21.24)], "v3": (-40.22, 26.51), "v4": (-4.72, 27.15), "zeta": [(-167.12, -40.24), (30.42, 177.92)]},
            "G-G": {"alpha": [(-122.7, -37.18), (40.15, 176.69)], "beta": (-178.56, 178.93), "chi": (-178.25, 178.01), "delta": [(76.68, 90.42), (133.71, 154.13)], "epsilon": (-169.42, -78.55), "eta": [(145.75, 179.24), (-179.48, 114.22)], "gamma": [(37.95, 173.66), (-179.53, -99.54)], "v1": (-31.26, 38.43), "v2": [(28.52, 40.92), (-39.81, -28.8)], "v3": (-39.61, 26.73), "v4": (-4.66, 26.39), "zeta": [(-174.34, -31.71), (26.11, 177.34)]},
            "_OTHER": {"alpha": (-152.96, 157.03), "beta": (-177.7, 178.48), "chi": (-179.51, 158.15), "delta": (74.36, 151.82), "epsilon": (-170.74, -95.22), "eta": (-178.33, 178.92), "gamma": (-165.7, 168.36), "v1": (-33.75, 41.66), "v2": (-38.28, 44.67), "v3": (-51.42, 27.64), "v4": (-5.79, 51.73), "zeta": (-156.6, 176.91)},
        },
        "tHW": {
            "A-A": {"alpha": [(-174.43, -30.85), (23.34, 176.22)], "beta": (-176.85, 177.67), "chi": [(-38.88, 179.72), (-178.68, -79.7)], "delta": [(124.7, 157.67), (74.44, 91.39)], "epsilon": (-166.08, -67.35), "eta": [(143.99, 178.22), (-178.44, 116.69)], "gamma": [(-47.84, 175.41), (-179.21, -117.12)], "v1": [(9.8, 47.25), (-31.54, -14.06)], "v2": [(-41.79, -17.77), (26.91, 42.04)], "v3": (-40.99, 26.6), "v4": (-4.44, 29.21), "zeta": [(28.82, 171.23), (-158.14, -43.81)]},
            "A-C": {"alpha": [(-162.38, -36.3), (25.14, 176.83)], "beta": (-177.84, 177.94), "chi": [(-178.2, -82.71), (2.76, 179.68)], "delta": [(75.25, 92.09), (130.05, 154.17)], "epsilon": (-165.69, -68.95), "eta": [(58.82, 178.98), (-179.19, -32.83)], "gamma": [(17.03, 169.81), (-178.39, -31.89)], "v1": (-30.73, 39.4), "v2": [(15.24, 41.25), (-39.28, -23.54)], "v3": (-40.75, 26.58), "v4": (-5.51, 32.96), "zeta": [(-168.14, -32.44), (38.31, 172.71)]},
            "A-U": {"alpha": [(22.51, 175.83), (-163.28, -34.57)], "beta": [(-179.19, -93.56), (73.92, 178.64)], "chi": (-178.41, 178.03), "delta": [(74.66, 91.23), (124.57, 154.3)], "epsilon": (-167.28, -60.1), "eta": (-178.27, 178.41), "gamma": [(-3.53, 176.9), (-179.2, -73.09)], "v1": (-31.83, 39.9), "v2": [(26.07, 41.19), (-40.09, -20.86)], "v3": (-41.25, 26.44), "v4": (-4.55, 30.44), "zeta": [(-164.77, -37.2), (36.06, 174.35)]},
            "G-G": {"alpha": [(19.44, 176.9), (-170.38, -26.32)], "beta": (-177.93, 178.61), "chi": [(-177.69, -67.25), (-20.43, 179.69)], "delta": [(123.13, 153.2), (75.6, 91.96)], "epsilon": (-170.34, 71.51), "eta": [(-178.87, -32.81), (73.31, 178.72)], "gamma": [(157.77, 179.35), (-175.43, 147.12)], "v1": (-31.7, 40.4), "v2": [(-39.23, -24.33), (23.95, 41.72)], "v3": (-39.81, 27.26), "v4": (-6.74, 29.68), "zeta": [(-166.89, -26.81), (25.32, 173.34)]},
            "_OTHER": {"alpha": (-171.96, 172.45), "beta": (-179.58, 178.4), "chi": (-179.03, 179.99), "delta": (71.14, 163.98), "epsilon": (-177.59, 174.39), "eta": (-179.6, 179.52), "gamma": (-178.51, 176.61), "v1": (-33.72, 43.0), "v2": (-44.21, 50.35), "v3": (-50.84, 36.7), "v4": (-21.21, 40.38), "zeta": (-165.93, 178.78)},
        },
        "tSH": {
            "A-A": {"alpha": [(-160.04, -32.86), (35.07, 176.94)], "beta": (-177.54, 178.37), "chi": (-178.6, 178.4), "delta": (77.02, 147.86), "epsilon": [(-175.53, -87.8), (46.35, 179.61)], "eta": [(-179.63, 114.79), (148.0, 179.21)], "gamma": [(30.39, 172.12), (-179.33, -67.11)], "v1": (-30.82, 38.41), "v2": [(29.31, 40.26), (-39.03, -17.98)], "v3": (-39.14, 25.19), "v4": (-3.75, 26.48), "zeta": [(-173.11, -41.23), (37.41, 177.12)]},
            "A-G": {"alpha": [(-129.07, -36.57), (14.13, 176.38)], "beta": (-178.39, 178.53), "chi": (-177.38, 175.48), "delta": [(74.54, 89.59), (120.05, 154.95)], "epsilon": (-169.02, -67.76), "eta": [(-179.61, 128.4), (150.29, 178.57)], "gamma": (-148.46, 167.66), "v1": [(-31.64, -15.45), (23.9, 42.39)], "v2": [(29.24, 41.46), (-40.02, -20.5)], "v3": (-41.13, 24.9), "v4": (-2.81, 28.88), "zeta": [(-164.81, -35.36), (25.56, 177.59)]},
            "A-U": {"alpha": [(40.73, 178.02), (-153.22, -32.93)], "beta": (-177.31, 178.17), "chi": (-178.32, 176.91), "delta": [(75.54, 89.16), (128.52, 151.92)], "epsilon": (-165.72, -81.5), "eta": [(84.85, 178.38), (-179.17, 36.98)], "gamma": [(-179.32, -56.09), (30.54, 173.83)], "v1": (-30.91, 37.34), "v2": [(29.27, 40.66), (-38.87, -27.85)], "v3": (-41.45, 24.61), "v4": (-3.18, 29.36), "zeta": [(-169.42, -28.79), (32.89, 178.05)]},
            "G-G": {"alpha": [(-138.11, -33.09), (12.34, 176.61)], "beta": (-177.5, 178.46), "chi": (-177.95, 177.49), "delta": [(74.72, 89.07), (106.26, 153.63)], "epsilon": (-172.65, -62.27), "eta": [(149.76, 178.73), (-179.53, 110.61)], "gamma": (-163.44, 170.29), "v1": [(-31.69, -17.76), (21.25, 44.07)], "v2": [(30.02, 41.89), (-39.97, -18.05)], "v3": (-41.31, 24.39), "v4": (-1.3, 27.39), "zeta": [(-160.08, -25.04), (44.33, 178.06)]},
            "_OTHER": {"alpha": (-139.48, 157.61), "beta": (-179.14, 178.1), "chi": (-177.89, 179.22), "delta": (-171.24, 157.23), "epsilon": (-173.34, 170.26), "eta": (-178.94, 178.88), "gamma": (-179.22, 178.08), "v1": (-36.07, 44.1), "v2": (-43.59, 52.1), "v3": (-53.24, 42.59), "v4": (-31.05, 33.43), "zeta": (-170.1, 176.75)},
        },
        "tSS": {
            "A-G": {"alpha": [(-155.48, -36.1), (35.02, 176.08)], "beta": (-177.55, 177.78), "chi": (-178.53, 178.43), "delta": [(75.18, 90.52), (99.86, 157.95)], "epsilon": (-167.69, -76.28), "eta": [(148.44, 178.54), (-179.31, 136.54)], "gamma": [(31.9, 175.67), (-179.57, -102.94)], "v1": (-32.89, 35.49), "v2": (-35.8, 42.14), "v3": (-41.32, 24.28), "v4": (-2.19, 27.87), "zeta": [(-161.39, -41.14), (37.01, 177.54)]},
            "G-G": {"alpha": [(-174.14, -12.4), (46.83, 178.45)], "beta": [(107.42, 178.93), (-178.68, -86.78)], "chi": (-176.85, 173.08), "delta": [(75.45, 91.01), (132.15, 154.08)], "epsilon": (-168.16, -24.49), "eta": [(-178.77, 84.84), (129.78, 179.01)], "gamma": [(23.42, 169.17), (-178.83, -33.2)], "v1": [(-31.39, -13.33), (27.84, 43.56)], "v2": [(27.02, 41.5), (-40.61, -27.84)], "v3": (-40.07, 27.82), "v4": (-5.75, 26.94), "zeta": [(70.84, 174.72), (-172.09, -48.09)]},
            "_OTHER": {"alpha": (-160.56, 172.21), "beta": (-175.98, 177.94), "chi": (-178.87, 179.58), "delta": (70.83, 154.99), "epsilon": (-177.67, 114.99), "eta": (-179.32, 179.63), "gamma": (-178.37, 177.99), "v1": (-33.63, 37.14), "v2": (-35.96, 42.51), "v3": (-44.12, 27.12), "v4": (-15.72, 40.35), "zeta": (-97.81, 135.96)},
        },
        "tSW": {
            "A-A": {"alpha": [(-148.0, -35.99), (22.79, 174.27)], "beta": (-174.59, 177.27), "chi": [(-178.29, -112.85), (-100.38, 179.62)], "delta": [(72.3, 88.62), (133.92, 159.33)], "epsilon": [(-167.22, -97.78), (-83.2, 177.68)], "eta": [(152.38, 178.37), (-179.12, 131.15)], "gamma": [(32.08, 171.61), (-178.13, -68.22)], "v1": (-32.97, 36.63), "v2": [(30.77, 43.66), (-39.85, -24.15)], "v3": [(-42.33, -30.8), (14.02, 33.21)], "v4": (-3.83, 27.96), "zeta": (-138.09, 29.32)},
            "A-C": {"alpha": [(-164.48, -34.45), (24.24, 176.52)], "beta": (-178.38, 178.01), "chi": (-177.0, 175.76), "delta": [(127.48, 154.09), (75.57, 90.95)], "epsilon": (-166.08, -65.24), "eta": [(150.51, 178.51), (-178.68, 139.69)], "gamma": [(27.48, 174.78), (-179.6, -49.95)], "v1": (-30.48, 39.0), "v2": [(-39.17, -23.41), (27.87, 40.77)], "v3": (-40.53, 26.58), "v4": (-5.31, 27.47), "zeta": [(-160.54, -38.88), (33.65, 177.11)]},
            "A-G": {"alpha": [(-170.82, -34.88), (35.82, 177.35)], "beta": (-177.75, 178.1), "chi": [(-178.23, -103.33), (-45.3, 179.58)], "delta": [(75.31, 90.28), (130.64, 155.01)], "epsilon": (-167.17, -69.96), "eta": [(-178.93, 139.12), (150.01, 178.7)], "gamma": [(-175.86, 153.57), (164.79, 179.75)], "v1": (-32.4, 38.29), "v2": [(28.74, 41.92), (-40.09, -23.12)], "v3": (-40.53, 27.23), "v4": (-8.34, 28.32), "zeta": [(-153.8, -42.68), (38.21, 175.95)]},
            "A-U": {"alpha": [(-162.44, -36.08), (20.12, 176.45)], "beta": (-176.79, 177.62), "chi": [(-175.73, -68.3), (-10.95, 179.11)], "delta": [(74.72, 90.86), (122.74, 152.26)], "epsilon": (-167.11, -54.84), "eta": [(7.37, 177.65), (-178.42, -55.04)], "gamma": [(5.81, 172.36), (-179.36, -39.2)], "v1": [(-32.1, -16.37), (30.16, 41.3)], "v2": [(29.49, 42.36), (-39.25, -20.87)], "v3": (-41.13, 26.29), "v4": (-5.13, 27.42), "zeta": [(-174.11, -50.04), (47.95, 176.96)]},
            "C-G": {"alpha": [(19.86, 174.62), (-152.8, -35.28)], "beta": (-178.43, 178.71), "chi": [(-176.25, -60.97), (18.95, 179.45)], "delta": [(73.56, 91.32), (133.77, 157.53)], "epsilon": (-166.84, -67.93), "eta": [(144.44, 177.46), (-177.66, 132.15)], "gamma": [(32.92, 174.23), (-178.9, -44.71)], "v1": (-30.0, 39.52), "v2": [(25.22, 41.32), (-40.22, -26.64)], "v3": (-41.29, 27.65), "v4": (-7.66, 29.29), "zeta": [(-142.64, -43.71), (25.65, 166.0)]},
            "G-U": {"alpha": [(-174.98, -28.61), (29.6, 176.27)], "beta": (-178.24, 178.62), "chi": [(-178.03, -106.41), (-66.11, 179.3)], "delta": [(134.02, 158.3), (74.55, 91.08)], "epsilon": (-171.4, -56.22), "eta": [(-179.29, -22.96), (54.19, 178.1)], "gamma": [(-179.36, -128.88), (-55.05, 176.74)], "v1": (-30.86, 38.54), "v2": [(-40.69, -27.7), (25.65, 41.13)], "v3": (-41.12, 28.41), "v4": (-6.88, 30.35), "zeta": [(-160.88, -36.66), (30.62, 178.83)]},
            "_OTHER": {"alpha": (-174.13, 172.04), "beta": (-179.9, 179.68), "chi": (-179.78, 178.07), "delta": (35.12, 157.43), "epsilon": (-170.4, 144.85), "eta": (-177.07, 179.03), "gamma": (-177.02, 178.78), "v1": (-38.59, 44.49), "v2": (-40.62, 45.63), "v3": (-44.97, 31.2), "v4": (-10.87, 43.23), "zeta": (-174.7, 170.84)},
        },
        "tWH": {
            "A-A": {"alpha": [(-168.36, -28.57), (23.33, 175.36)], "beta": (-177.46, 177.6), "chi": [(-178.87, -67.37), (25.81, 179.7)], "delta": [(74.96, 90.62), (128.63, 155.36)], "epsilon": (-165.66, -15.58), "eta": [(-179.14, 115.49), (147.85, 179.26)], "gamma": [(8.24, 175.0), (-179.34, -48.23)], "v1": [(-31.97, -14.41), (26.35, 45.1)], "v2": [(28.56, 41.12), (-40.88, -23.79)], "v3": (-39.97, 27.77), "v4": (-5.96, 28.48), "zeta": [(-168.28, -34.65), (37.86, 177.52)]},
            "A-C": {"alpha": [(-166.03, -28.26), (31.77, 176.3)], "beta": (-178.15, 178.1), "chi": (-177.36, 176.78), "delta": [(130.43, 157.61), (75.29, 90.94)], "epsilon": (-170.13, -57.33), "eta": [(-179.01, -86.55), (50.76, 178.94)], "gamma": [(-64.97, 176.57), (-179.43, -126.85)], "v1": (-30.33, 40.35), "v2": [(-40.44, -28.4), (27.13, 40.55)], "v3": (-40.25, 28.5), "v4": (-7.3, 29.39), "zeta": [(33.91, 174.36), (-162.56, -34.74)]},
            "A-G": {"alpha": [(-157.67, -33.25), (36.31, 175.55)], "beta": [(-179.04, -98.02), (84.83, 178.03)], "chi": [(-178.22, -119.11), (-110.51, 179.25)], "delta": [(74.32, 91.78), (130.33, 163.61)], "epsilon": (-167.85, 15.35), "eta": [(-179.29, 96.29), (145.69, 178.91)], "gamma": [(26.22, 175.59), (-178.41, -53.48)], "v1": [(-32.5, -14.94), (19.53, 42.89)], "v2": [(26.81, 42.42), (-43.0, -23.93)], "v3": (-41.21, 27.16), "v4": (-4.68, 28.38), "zeta": [(-173.86, -7.45), (29.41, 179.07)]},
            "A-U": {"alpha": [(-164.02, -31.73), (23.57, 175.71)], "beta": [(73.09, 178.45), (-179.15, -93.64)], "chi": (-178.48, 178.48), "delta": [(75.62, 91.28), (119.99, 155.68)], "epsilon": (-167.96, 47.86), "eta": [(126.87, 179.01), (-179.01, 60.11)], "gamma": [(-13.01, 176.1), (-179.35, -84.1)], "v1": (-31.83, 40.77), "v2": [(24.51, 41.07), (-40.84, -21.22)], "v3": (-40.23, 27.36), "v4": (-5.71, 30.79), "zeta": [(-165.5, -45.61), (42.56, 177.2)]},
            "G-G": {"alpha": [(-152.23, -29.79), (24.47, 174.18)], "beta": (-178.0, 177.92), "chi": [(-176.88, -65.04), (0.26, 179.57)], "delta": [(74.3, 89.09), (126.71, 154.24)], "epsilon": (-167.45, -68.74), "eta": (-174.46, 176.85), "gamma": [(-174.85, 149.52), (162.31, 179.47)], "v1": [(-31.88, -7.88), (29.24, 43.98)], "v2": [(25.43, 41.82), (-40.69, -23.91)], "v3": (-41.08, 26.42), "v4": (-3.84, 31.3), "zeta": [(-169.19, -45.64), (67.48, 174.52)]},
            "_OTHER": {"alpha": (-175.67, 179.75), "beta": (-179.39, 179.69), "chi": (-179.43, 179.33), "delta": (-26.01, 155.67), "epsilon": (-178.2, 178.44), "eta": (-178.65, 179.93), "gamma": (-179.43, 178.13), "v1": (-36.45, 42.98), "v2": (-40.42, 56.13), "v3": (-63.95, 31.95), "v4": (-16.08, 49.29), "zeta": (-179.08, 179.81)},
        },
        "tWS": {
            "A-C": {"alpha": [(-162.65, -38.33), (31.67, 174.73)], "beta": (-177.64, 178.28), "chi": (-177.0, 176.0), "delta": [(75.16, 90.86), (122.75, 154.34)], "epsilon": (-165.53, -93.09), "eta": [(154.18, 178.81), (-179.51, 139.26)], "gamma": [(33.56, 173.41), (-179.14, -55.11)], "v1": [(-31.68, -15.65), (28.55, 42.71)], "v2": [(28.04, 41.59), (-39.03, -25.94)], "v3": (-41.18, 25.05), "v4": (-3.16, 27.5), "zeta": [(20.86, 175.72), (-137.7, -41.53)]},
            "A-G": {"alpha": [(-172.51, -30.95), (30.38, 177.53)], "beta": (-178.03, 178.26), "chi": [(-178.02, -80.83), (-35.05, 179.72)], "delta": [(74.74, 91.12), (118.06, 155.58)], "epsilon": (-167.61, -55.64), "eta": [(149.97, 178.94), (-179.12, 134.8)], "gamma": [(25.47, 175.52), (-179.43, -45.37)], "v1": (-31.34, 39.67), "v2": [(27.59, 41.79), (-40.05, -19.61)], "v3": (-41.03, 27.05), "v4": (-6.11, 28.23), "zeta": [(34.4, 175.37), (-162.27, -32.97)]},
            "_OTHER": {"alpha": (-173.81, 178.84), "beta": (-179.51, 178.15), "chi": (-179.17, 177.96), "delta": (4.5, 171.69), "epsilon": (-169.74, 171.01), "eta": (-178.21, 177.62), "gamma": (-177.4, 176.13), "v1": (-40.39, 44.23), "v2": (-46.6, 56.79), "v3": (-60.78, 44.03), "v4": (-27.31, 44.39), "zeta": (-175.69, 176.16)},
        },
        "tWW": {
            "A-A": {"alpha": [(-173.09, -29.52), (28.71, 176.63)], "beta": (-177.93, 177.54), "chi": [(6.0, 179.73), (-178.41, -77.1)], "delta": [(75.65, 91.28), (114.03, 156.53)], "epsilon": (-171.2, 41.74), "eta": [(-175.89, 3.38), (55.01, 177.28)], "gamma": [(-36.04, 177.1), (-179.51, -68.72)], "v1": (-31.38, 41.46), "v2": [(24.86, 40.23), (-41.54, -17.09)], "v3": (-39.76, 28.32), "v4": (-6.84, 31.01), "zeta": [(-167.99, -45.52), (50.05, 174.22)]},
            "A-C": {"alpha": [(45.32, 179.22), (-177.42, -5.67)], "beta": (-177.95, 177.61), "chi": [(-177.84, -76.45), (-0.74, 179.14)], "delta": [(134.51, 154.32), (75.01, 95.53)], "epsilon": (-167.88, -58.15), "eta": [(-101.13, 176.66), (-177.89, -139.96)], "gamma": [(-17.77, 176.62), (-179.25, -63.39)], "v1": (-31.77, 41.27), "v2": [(-40.07, -28.41), (26.91, 40.3)], "v3": (-39.81, 28.2), "v4": (-7.38, 29.74), "zeta": [(45.62, 173.72), (-165.69, -23.88)]},
            "A-U": {"alpha": [(-171.29, -29.96), (26.28, 176.12)], "beta": (-178.07, 178.06), "chi": [(-178.44, -98.31), (-58.78, 178.81)], "delta": [(73.76, 91.93), (125.26, 157.12)], "epsilon": (-169.51, 52.74), "eta": [(-178.08, 50.64), (104.83, 178.85)], "gamma": [(-8.81, 176.67), (-179.43, -59.53)], "v1": (-29.82, 42.05), "v2": [(18.49, 40.77), (-41.26, -22.67)], "v3": (-41.08, 28.02), "v4": (-5.86, 33.51), "zeta": [(-173.07, -28.18), (50.2, 176.83)]},
            "C-G": {"alpha": [(-169.93, -31.48), (33.1, 174.86)], "beta": [(92.55, 178.38), (-178.42, -100.24)], "chi": [(-15.39, 179.12), (-175.95, -98.98)], "delta": [(74.48, 95.3), (127.74, 155.51)], "epsilon": (-169.58, 63.83), "eta": [(143.85, 178.15), (-178.53, 126.92)], "gamma": [(-52.22, 174.09), (-179.28, -132.96)], "v1": (-31.14, 40.0), "v2": [(23.44, 41.49), (-39.94, -24.92)], "v3": (-40.49, 27.92), "v4": (-7.63, 28.32), "zeta": [(-170.68, -37.01), (43.55, 177.09)]},
            "G-U": {"alpha": [(-177.01, -29.76), (22.05, 177.73)], "beta": (-177.46, 178.29), "chi": [(-177.01, -83.76), (-70.0, 179.14)], "delta": [(74.56, 89.37), (130.2, 158.22)], "epsilon": (-167.63, -60.66), "eta": [(-178.96, 34.96), (63.03, 177.22)], "gamma": [(-2.49, 177.81), (-179.22, -53.67)], "v1": (-30.5, 41.02), "v2": [(21.47, 41.34), (-40.3, -27.05)], "v3": (-41.29, 28.1), "v4": (-7.23, 31.83), "zeta": [(46.22, 176.37), (-166.74, -29.28)]},
            "U-U": {"alpha": [(14.3, 176.52), (-155.13, -24.25)], "beta": [(81.48, 178.06), (-178.39, -87.01)], "chi": [(-174.69, -72.59), (-12.13, 178.62)], "delta": [(132.32, 156.09), (75.9, 89.62)], "epsilon": (-167.18, -29.79), "eta": [(25.17, 176.18), (-176.24, -79.58)], "gamma": [(0.02, 172.52), (-178.93, -51.41)], "v1": [(27.06, 43.47), (-30.89, -13.28)], "v2": [(-41.62, -27.05), (28.07, 41.39)], "v3": (-39.59, 28.71), "v4": (-6.2, 27.26), "zeta": [(48.24, 152.88), (-144.06, -40.63)]},
            "_OTHER": {"alpha": (-168.34, 178.05), "beta": (-177.91, 179.54), "chi": (-179.83, 179.49), "delta": (46.38, 169.17), "epsilon": (-177.17, 88.55), "eta": (-179.03, 179.56), "gamma": (-175.76, 177.74), "v1": (-37.06, 49.79), "v2": (-47.91, 50.13), "v3": (-50.66, 43.6), "v4": (-28.62, 53.49), "zeta": (-177.62, 174.54)},
        },
        "_OTHER": {
            "alpha": (-175.87, 178.85), "beta": (-179.9, 178.97), "chi": (-178.6, 179.7), "delta": (-18.62, 161.52), "epsilon": (-173.78, 170.44), "eta": (-179.26, 179.63), "gamma": (-179.79, 177.83), "v1": (-70.08, 47.48), "v2": (-44.96, 52.19), "v3": (-56.76, 37.55), "v4": (-30.2, 46.7), "zeta": (-176.84, 177.46),
        },
    }




            # Higher = catches more