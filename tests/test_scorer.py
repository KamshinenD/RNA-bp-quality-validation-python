"""Tests for scorer2.py - RNA structure quality scoring."""

import pytest
import pandas as pd
from scorer2 import Scorer, BaselineResult


class TestScorer:
    """Tests for the Scorer class."""

    def test_scorer_initialization(self, config):
        """Test that Scorer initializes correctly with config."""
        scorer = Scorer(config)
        assert scorer.config == config

    def test_empty_structure_returns_empty_result(self, config):
        """Test scoring an empty structure returns appropriate empty result."""
        scorer = Scorer(config)
        result = scorer.score_structure([], pd.DataFrame())

        assert result.overall_score == "N/A"
        assert result.total_base_pairs == 0
        assert result.summary == "No base pairs found for analysis"

    def test_score_perfect_base_pair(self, config, sample_base_pair, sample_hbond_data):
        """Test scoring a perfect base pair returns high score."""
        scorer = Scorer(config)
        bp_score = scorer._score_base_pair(sample_base_pair, sample_hbond_data)

        # Perfect geometry and 3 good H-bonds should score high
        assert bp_score['score'] >= 85
        assert bp_score['geometry_penalty'] == 0
        assert bp_score['bp_info']['bp_type'] == 'G-C'
        assert bp_score['bp_info']['num_hbonds'] == 3

    def test_score_misaligned_base_pair(self, config, sample_misaligned_base_pair, sample_hbond_data):
        """Test that misaligned base pair is penalized."""
        scorer = Scorer(config)
        bp_score = scorer._score_base_pair(sample_misaligned_base_pair, sample_hbond_data)

        # Should have misalignment penalty
        assert bp_score['geometry_issues']['misaligned'] is True
        assert bp_score['geometry_penalty'] > 0
        assert bp_score['score'] < 100

    def test_score_rotational_distortion_base_pair(self, config, sample_twisted_base_pair, sample_hbond_data):
        """Test that rotationally distorted base pair is penalized."""
        scorer = Scorer(config)
        bp_score = scorer._score_base_pair(sample_twisted_base_pair, sample_hbond_data)

        # Should have rotational distortion penalty
        assert bp_score['geometry_issues']['rotational_distortion'] is True
        assert bp_score['geometry_penalty'] > 0
        assert bp_score['score'] < 100

    def test_score_non_coplanar_base_pair(self, config, sample_non_coplanar_base_pair, sample_hbond_data):
        """Test that non-coplanar base pair is penalized."""
        scorer = Scorer(config)
        bp_score = scorer._score_base_pair(sample_non_coplanar_base_pair, sample_hbond_data)

        # Should have coplanarity penalty
        assert bp_score['geometry_issues']['non_coplanar'] is True
        assert bp_score['geometry_penalty'] > 0
        assert bp_score['score'] < 100

    def test_score_base_pair_with_bad_hbonds(self, config, sample_base_pair, sample_bad_hbond_data):
        """Test that base pair with poor H-bonds is penalized."""
        scorer = Scorer(config)
        bp_score = scorer._score_base_pair(sample_base_pair, sample_bad_hbond_data)

        # Should have H-bond penalties
        assert bp_score['hbond_penalty'] > 0
        assert bp_score['score'] < 100
        # Check for specific H-bond issues
        assert (bp_score['hbond_issues'].get('bad_distance') or
                bp_score['hbond_issues'].get('bad_angles') or
                bp_score['hbond_issues'].get('bad_dihedral') or
                bp_score['hbond_issues'].get('weak_quality'))

    def test_score_base_pair_with_no_hbonds(self, config, sample_base_pair, empty_hbond_data):
        """Test that base pair with no H-bonds is heavily penalized."""
        # Set hbond_score to 0 to indicate no H-bonds detected by DSSR
        sample_base_pair_no_hbonds = sample_base_pair.copy()
        sample_base_pair_no_hbonds['hbond_score'] = 0.0

        scorer = Scorer(config)
        bp_score = scorer._score_base_pair(sample_base_pair_no_hbonds, empty_hbond_data)

        # Should flag zero_hbond issue
        assert bp_score['geometry_issues']['zero_hbond'] is True
        assert bp_score['geometry_penalty'] >= config.PENALTY_WEIGHTS['zero_hbond_pairs']
        assert bp_score['score'] < 85

    def test_score_structure_with_multiple_base_pairs(self, config, sample_basepair_list, sample_hbond_data):
        """Test scoring a structure with multiple base pairs."""
        scorer = Scorer(config)
        result = scorer.score_structure(sample_basepair_list, sample_hbond_data)

        assert isinstance(result, BaselineResult)
        assert result.total_base_pairs == 3
        assert result.overall_score > 0
        assert len(result.basepair_scores) == 3

    def test_is_base_atom(self, config):
        """Test base atom detection."""
        scorer = Scorer(config)

        # Base atoms
        assert scorer._is_base_atom('N1') is True
        assert scorer._is_base_atom('N2') is True
        assert scorer._is_base_atom('O6') is True
        assert scorer._is_base_atom('N4') is True

        # Backbone/sugar atoms
        assert scorer._is_base_atom('P') is False
        assert scorer._is_base_atom('OP1') is False
        assert scorer._is_base_atom('OP2') is False
        assert scorer._is_base_atom("C5'") is False
        assert scorer._is_base_atom("O2'") is False
        assert scorer._is_base_atom("C1'") is False
        assert scorer._is_base_atom('C5*') is False
        assert scorer._is_base_atom('O4*') is False

    def test_get_basepair_hbonds(self, config, sample_hbond_data):
        """Test extraction of H-bonds for a specific base pair."""
        scorer = Scorer(config)

        bp_hbonds = scorer._get_basepair_hbonds('A-G-10-', 'A-C-20-', sample_hbond_data)

        # Should get all 3 H-bonds for this base pair
        assert len(bp_hbonds) == 3
        assert all(bp_hbonds['res_1'].isin(['A-G-10-', 'A-C-20-']))
        assert all(bp_hbonds['res_2'].isin(['A-G-10-', 'A-C-20-']))

    def test_get_basepair_hbonds_bidirectional(self, config):
        """Test that H-bond extraction works bidirectionally."""
        scorer = Scorer(config)

        # Create H-bonds in both directions
        hbond_data = pd.DataFrame([
            {
                'res_1': 'A-G-10-',
                'res_2': 'A-C-20-',
                'atom_1': 'N1',
                'atom_2': 'N3',
                'distance': 2.9,
                'angle_1': 150.0,
                'angle_2': 145.0,
                'dihedral_angle': 10.0,
                'score': 0.85
            },
            {
                'res_1': 'A-C-20-',
                'res_2': 'A-G-10-',
                'atom_1': 'N4',
                'atom_2': 'O6',
                'distance': 2.85,
                'angle_1': 160.0,
                'angle_2': 155.0,
                'dihedral_angle': 5.0,
                'score': 0.92
            }
        ])

        bp_hbonds = scorer._get_basepair_hbonds('A-G-10-', 'A-C-20-', hbond_data)

        # Should get both H-bonds regardless of direction
        assert len(bp_hbonds) == 2

    def test_export_to_dict(self, config, sample_basepair_list, sample_hbond_data):
        """Test exporting results to dictionary format."""
        scorer = Scorer(config)
        result = scorer.score_structure(sample_basepair_list, sample_hbond_data)

        export_dict = scorer.export_to_dict(result)

        assert 'overall_score' in export_dict
        assert 'total_base_pairs' in export_dict
        assert 'geometry_issues' in export_dict
        assert 'hbond_issues' in export_dict
        assert export_dict['total_base_pairs'] == 3

    def test_incorrect_hbond_count_penalty(self, config):
        """Test that incorrect H-bond count is penalized."""
        scorer = Scorer(config)

        # G-C pair should have 3 H-bonds
        gc_pair = {
            'res_1': 'A-G-10-',
            'res_2': 'A-C-20-',
            'bp_type': 'G-C',
            'lw': 'cWW',
            'shear': 0.1,
            'stretch': 0.05,
            'stagger': 0.15,
            'buckle': 5.0,
            'propeller': -3.0,
            'opening': 2.0,
            'hbond_score': 3.5
        }

        # Create H-bond data with only 2 bonds (incorrect for G-C)
        hbond_data = pd.DataFrame([
            {
                'res_1': 'A-G-10-',
                'res_2': 'A-C-20-',
                'atom_1': 'N1',
                'atom_2': 'N3',
                'distance': 2.9,
                'angle_1': 150.0,
                'angle_2': 145.0,
                'dihedral_angle': 10.0,
                'score': 0.85
            },
            {
                'res_1': 'A-G-10-',
                'res_2': 'A-C-20-',
                'atom_1': 'N2',
                'atom_2': 'O2',
                'distance': 3.0,
                'angle_1': 155.0,
                'angle_2': 150.0,
                'dihedral_angle': -15.0,
                'score': 0.90
            }
        ])

        bp_score = scorer._score_base_pair(gc_pair, hbond_data)

        # Should flag incorrect H-bond count
        assert bp_score['hbond_issues']['incorrect_count'] is True
        assert bp_score['hbond_penalty'] >= config.PENALTY_WEIGHTS['incorrect_hbond_count']

    def test_score_never_exceeds_100(self, config, sample_base_pair, sample_hbond_data):
        """Test that score is capped at 100."""
        scorer = Scorer(config)
        bp_score = scorer._score_base_pair(sample_base_pair, sample_hbond_data)

        assert bp_score['score'] <= 100

    def test_score_never_below_0(self, config):
        """Test that score is floored at 0 even with multiple issues."""
        scorer = Scorer(config)

        # Create terrible base pair with all issues
        terrible_bp = {
            'res_1': 'A-G-10-',
            'res_2': 'A-C-20-',
            'bp_type': 'G-C',
            'lw': 'cWW',
            'shear': 5.0,  # Way too high
            'stretch': 3.0,  # Way too high
            'stagger': 2.0,  # Way too high
            'buckle': 50.0,  # Way too high
            'propeller': -40.0,  # Way too low
            'opening': 60.0,  # Way too high
            'hbond_score': 0.0  # No H-bonds
        }

        bp_score = scorer._score_base_pair(terrible_bp, pd.DataFrame())

        assert bp_score['score'] >= 0
