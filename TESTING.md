# Testing Guide for RNA Structure Quality Scoring

## Overview

A test suite has been created for your RNA structure quality scoring system. The tests validate the correctness of scoring logic in `scorer2.py` and data loading in `utils/data_loader.py`.

## Test Suite Structure

```
score_structure_motifs/
├── tests/
│   ├── __init__.py                    # Test package
│   ├── conftest.py                    # Shared fixtures and test data
│   ├── test_scorer.py                 # Tests for scorer2.py (16 test scenarios)
│   ├── test_data_loader.py            # Tests for data loading (15 test scenarios)
│   └── README.md                      # Test documentation
├── pytest.ini                         # Pytest configuration
├── requirements-test.txt              # Testing dependencies
├── run_tests.sh                       # Helper script to run tests
└── TESTING.md                         # This file
```

## Quick Start

### 1. Install Testing Dependencies

```bash
pip install -r requirements-test.txt
```

This installs:
- `pytest` - Testing framework
- `pytest-cov` - Coverage reporting
- `pytest-mock` - Mocking utilities

### 2. Run All Tests

```bash
# Using the helper script (recommended)
./run_tests.sh

# Or directly with pytest
pytest

# With verbose output
pytest -v
```

### 3. Run Specific Tests

```bash
# Test only the scorer
pytest tests/test_scorer.py

# Test only data loader
pytest tests/test_data_loader.py

# Test a specific function
pytest tests/test_scorer.py::TestScorer::test_score_perfect_base_pair
```

## What's Being Tested

### Scorer Tests ([test_scorer.py](tests/test_scorer.py))
Tests for [scorer2.py](scorer2.py) - the main scoring engine:

- ✓ Empty structure handling
- ✓ Perfect base pair scoring (should score ≥85)
- ✓ Misalignment detection and penalty
- ✓ Twist detection and penalty
- ✓ Non-coplanarity detection and penalty
- ✓ H-bond quality assessment
- ✓ Zero H-bond penalty (heavily penalized)
- ✓ Incorrect H-bond count detection
- ✓ Grade assignment (EXCELLENT/GOOD/FAIR/POOR)
- ✓ Score bounds (0-100)
- ✓ Base vs backbone atom detection
- ✓ Bidirectional H-bond extraction
- ✓ Export to dictionary format

### Data Loader Tests ([test_data_loader.py](tests/test_data_loader.py))
Tests for [utils/data_loader.py](utils/data_loader.py):

- ✓ Base pair loading from JSON
- ✓ H-bond loading from CSV
- ✓ Stacking interaction filtering
- ✓ RNA-RNA vs RNA-protein filtering
- ✓ File not found handling
- ✓ CIF file download (mocked)
- ✓ Nucleotide counting (mocked)
- ✓ Validation metrics fetching (mocked)
- ✓ Backward compatibility (dict format)
- ✓ Error handling

## Test Fixtures

Reusable test data is defined in [tests/conftest.py](tests/conftest.py):

- `config` - Standard configuration object from [config.py](config.py)
- `sample_base_pair` - Perfect base pair (G-C, cWW)
- `sample_misaligned_base_pair` - High shear value
- `sample_twisted_base_pair` - High buckle
- `sample_non_coplanar_base_pair` - High stagger
- `sample_hbond_data` - 3 good quality H-bonds
- `sample_bad_hbond_data` - H-bond with multiple issues
- `sample_basepair_list` - 3 base pairs for structure testing
- `empty_hbond_data` - Empty DataFrame

## Coverage Report

Generate HTML coverage report:

```bash
pytest --cov=. --cov-report=html
open htmlcov/index.html  # View coverage in browser
```

## Common Test Patterns

### Testing Score Ranges
```python
def test_score_perfect_base_pair(self, config, sample_base_pair, sample_hbond_data):
    scorer = Scorer(config)
    bp_score = scorer._score_base_pair(sample_base_pair, sample_hbond_data)

    assert bp_score['score'] >= 85  # Perfect geometry + 3 good H-bonds
    assert bp_score['geometry_penalty'] == 0
```

### Testing Issue Detection
```python
def test_score_misaligned_base_pair(self, config, sample_misaligned_base_pair, sample_hbond_data):
    scorer = Scorer(config)
    bp_score = scorer._score_base_pair(sample_misaligned_base_pair, sample_hbond_data)

    assert bp_score['geometry_issues']['misaligned'] is True
    assert bp_score['geometry_penalty'] > 0
```

### Testing Data Loading
```python
def test_load_basepairs_success(self, config, tmp_path):
    # Create test JSON file
    test_json = tmp_path / "basepairs" / "TEST.json"
    test_json.parent.mkdir(parents=True)
    test_json.write_text(json.dumps([{"bp_type": "G-C"}]))

    loader = DataLoader(config)
    result = loader.load_basepairs("TEST")
    assert len(result) == 1
```

## Why These Tests Matter

### 1. **Prevent Regressions**
When you modify scoring logic in [scorer2.py](scorer2.py) or thresholds in [config.py](config.py), tests ensure you don't break existing functionality:

```bash
# After changing code
pytest tests/test_scorer.py -v
```

### 2. **Validate Scoring Logic**
Tests verify that penalties are applied correctly:
- Misalignment → 15 points penalty
- Zero H-bonds → 20 points penalty
- Bad H-bond distance → 15 points penalty

### 3. **Scientific Reproducibility**
Tests document expected behavior:
- G-C pairs should have 3 H-bonds
- H-bond distances should be 2.3-3.7Å
- Dihedral angles should be CIS (-50° to +50°) or TRANS (|angle| ≥ 140°)

### 4. **Catch Edge Cases**
Tests ensure proper handling of:
- Empty structures
- Missing data files
- Base vs backbone atoms
- Bidirectional H-bond lookup

## Continuous Integration

Add to `.github/workflows/test.yml`:

```yaml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-python@v4
        with:
          python-version: '3.9'
      - run: pip install -r requirements-test.txt
      - run: pytest --cov=. --cov-report=xml
```

## Next Steps

### Running Your First Test
```bash
# Run a simple test to verify everything works
pytest tests/test_scorer.py::TestScorer::test_scorer_initialization -v
```

### Adding New Tests
When you add new features to scorer2.py, create tests following this pattern:

```python
def test_new_feature(self, config):
    """Test that new feature works correctly."""
    # Setup
    scorer = Scorer(config)

    # Execute
    result = scorer.new_method(test_data)

    # Assert
    assert result.expected == actual
```

### Test-Driven Development
1. Write a failing test for new functionality
2. Implement the feature
3. Run tests to verify it works
4. Refactor with confidence

## Troubleshooting

### "ModuleNotFoundError"
```bash
# Ensure you're in the project root
cd /Users/kdewan2/Desktop/Projects/score_structure_motifs

# Install dependencies
pip install -r requirements-test.txt
```

### "No tests ran"
```bash
# Verify test files exist
ls tests/test_*.py

# Run with discovery
pytest --collect-only
```

### Tests Failing
```bash
# Run with detailed output
pytest -vv --tb=long

# Run single test to isolate issue
pytest tests/test_scorer.py::TestScorer::test_specific_test -vv
```

## Summary

You now have **31 automated tests** covering:
- Core scoring engine (scorer2.py) - 16 tests
- Data loading (utils/data_loader.py) - 15 tests
- Edge cases and error handling

Run `./run_tests.sh` anytime to verify your code works correctly!

## Key Files Tested

- [scorer2.py](scorer2.py) - Main scoring engine (Scorer class)
- [utils/data_loader.py](utils/data_loader.py) - Data loading utilities
- [config.py](config.py) - Configuration with thresholds and weights

## What's NOT Tested

The following files are NOT tested because they're not used in the main execution path:
- `analyzers/base_pair_analyzer_bc.py` - Not used by scorer2.py
- `analyzers/hbond_analyzer_bc.py` - Not used by scorer2.py
- `config_bc_thresholds.py` - Alternative config, not used by app.py

The main entry points ([app.py](app.py) and [scorer2.py](scorer2.py)) use the files listed in "Key Files Tested" above.
