"""Pytest configuration and shared fixtures."""

import pytest
import pandas as pd
from pathlib import Path
import sys

# Add parent directory to path so we can import modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from config import Config


@pytest.fixture
def config():
    """Provide a standard configuration object for tests."""
    return Config()


@pytest.fixture
def sample_base_pair():
    """Sample base pair data with good geometry."""
    return {
        'res_1': 'A-G-10-',
        'res_2': 'A-C-20-',
        'bp_type': 'G-C',
        'lw': 'cWW',  # canonical Watson-Crick
        'shear': 0.1,
        'stretch': 0.05,
        'stagger': 0.15,
        'buckle': 5.0,
        'propeller': -3.0,
        'opening': 2.0,
        'hbond_score': 3.5  # DSSR quality score
    }


@pytest.fixture
def sample_misaligned_base_pair():
    """Sample base pair with misalignment issues."""
    return {
        'res_1': 'A-G-10-',
        'res_2': 'A-C-20-',
        'bp_type': 'G-C',
        'lw': 'cWW',
        'shear': 2.5,  # Exceeds threshold (1.5)
        'stretch': 0.05,
        'stagger': 0.15,
        'buckle': 5.0,
        'propeller': -3.0,
        'opening': 2.0,
        'hbond_score': 3.5
    }


@pytest.fixture
def sample_twisted_base_pair():
    """Sample base pair with twist issues."""
    return {
        'res_1': 'A-G-10-',
        'res_2': 'A-C-20-',
        'bp_type': 'G-C',
        'lw': 'cWW',
        'shear': 0.1,
        'stretch': 0.05,
        'stagger': 0.15,
        'buckle': 25.0,  # Exceeds threshold (16.1)
        'propeller': -3.0,
        'opening': 2.0,
        'hbond_score': 3.5
    }


@pytest.fixture
def sample_non_coplanar_base_pair():
    """Sample base pair with coplanarity issues."""
    return {
        'res_1': 'A-G-10-',
        'res_2': 'A-C-20-',
        'bp_type': 'G-C',
        'lw': 'cWW',
        'shear': 0.1,
        'stretch': 0.05,
        'stagger': 0.8,  # Exceeds threshold (0.4)
        'buckle': 5.0,
        'propeller': -3.0,
        'opening': 2.0,
        'hbond_score': 3.5
    }


@pytest.fixture
def sample_hbond_data():
    """Sample hydrogen bond DataFrame with good bonds."""
    return pd.DataFrame([
        {
            'res_1': 'A-G-10-',
            'res_2': 'A-C-20-',
            'atom_1': 'N1',  # Base atom
            'atom_2': 'N3',  # Base atom
            'distance': 2.9,  # Good distance
            'angle_1': 150.0,  # Good angle
            'angle_2': 145.0,  # Good angle
            'dihedral_angle': 10.0,  # CIS orientation
            'score': 0.85,  # Good quality
            'res_type_1': 'RNA',
            'res_type_2': 'RNA'
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
            'score': 0.90,
            'res_type_1': 'RNA',
            'res_type_2': 'RNA'
        },
        {
            'res_1': 'A-G-10-',
            'res_2': 'A-C-20-',
            'atom_1': 'O6',
            'atom_2': 'N4',
            'distance': 2.85,
            'angle_1': 160.0,
            'angle_2': 155.0,
            'dihedral_angle': 5.0,
            'score': 0.92,
            'res_type_1': 'RNA',
            'res_type_2': 'RNA'
        }
    ])


@pytest.fixture
def sample_bad_hbond_data():
    """Sample hydrogen bond DataFrame with problematic bonds."""
    return pd.DataFrame([
        {
            'res_1': 'A-G-10-',
            'res_2': 'A-C-20-',
            'atom_1': 'N1',
            'atom_2': 'N3',
            'distance': 4.2,  # Too long (max is 3.7)
            'angle_1': 60.0,  # Too small (min is 80)
            'angle_2': 145.0,
            'dihedral_angle': 90.0,  # Forbidden zone (not CIS or TRANS)
            'score': 0.50,  # Weak quality (min is 0.70)
            'res_type_1': 'RNA',
            'res_type_2': 'RNA'
        }
    ])


@pytest.fixture
def sample_basepair_list():
    """Sample list of base pairs for structure scoring."""
    return [
        {
            'res_1': 'A-G-1-',
            'res_2': 'A-C-24-',
            'bp_type': 'G-C',
            'lw': 'cWW',
            'shear': 0.1,
            'stretch': 0.05,
            'stagger': 0.15,
            'buckle': 5.0,
            'propeller': -3.0,
            'opening': 2.0,
            'hbond_score': 3.5
        },
        {
            'res_1': 'A-A-2-',
            'res_2': 'A-U-23-',
            'bp_type': 'A-U',
            'lw': 'cWW',
            'shear': -0.2,
            'stretch': 0.1,
            'stagger': 0.2,
            'buckle': 4.0,
            'propeller': -2.0,
            'opening': 1.5,
            'hbond_score': 3.2
        },
        {
            'res_1': 'A-C-3-',
            'res_2': 'A-G-22-',
            'bp_type': 'C-G',
            'lw': 'cWW',
            'shear': 0.15,
            'stretch': -0.1,
            'stagger': 0.1,
            'buckle': 6.0,
            'propeller': -4.0,
            'opening': 2.5,
            'hbond_score': 3.8
        }
    ]


@pytest.fixture
def empty_hbond_data():
    """Empty hydrogen bond DataFrame."""
    return pd.DataFrame()
