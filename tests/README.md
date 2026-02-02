# RNA Structure Quality Scoring - Test Suite

This directory contains automated tests for the RNA structure quality scoring system.

## Test Structure

```
tests/
├── __init__.py              # Test package initialization
├── conftest.py              # Pytest fixtures and shared test data
├── test_scorer.py           # Tests for scorer2.py (16 tests)
└── test_data_loader.py      # Tests for data loading (15 tests)
```

## Running Tests

### Install Dependencies

First, install pytest and required testing packages:

```bash
pip install -r requirements-test.txt
```

### Run All Tests

```bash
# From project root directory
pytest

# With verbose output
pytest -v

# With coverage report
pytest --cov=. --cov-report=html
```

### Run Specific Test Files

```bash
# Test only the scorer
pytest tests/test_scorer.py

# Test only data loader
pytest tests/test_data_loader.py
```

### Run Specific Tests

```bash
# Run a specific test class
pytest tests/test_scorer.py::TestScorer

# Run a specific test function
pytest tests/test_scorer.py::TestScorer::test_score_perfect_base_pair

# Run tests matching a pattern
pytest -k "hbond"
```

## Test Coverage

The test suite covers:

### Scorer ([scorer2.py](../scorer2.py))
- Empty structure handling
- Perfect base pair scoring
- Geometry issue detection (misalignment, twist, coplanarity)
- H-bond quality assessment
- Zero H-bond penalty
- Incorrect H-bond count detection
- Grade assignment
- Score bounds (0-100)
- Base atom detection
- Bidirectional H-bond lookup
- Export functionality

### Data Loader ([utils/data_loader.py](../utils/data_loader.py))
- Base pair loading from JSON
- H-bond loading from CSV
- Stacking interaction filtering
- RNA-RNA vs RNA-protein filtering
- CIF file download (mocked)
- Nucleotide counting (mocked)
- Validation metrics fetching (mocked)
- Error handling

## Test Fixtures

Key fixtures defined in [conftest.py](conftest.py):

- `config`: Standard configuration object from [config.py](../config.py)
- `sample_base_pair`: Good quality base pair (G-C, cWW)
- `sample_misaligned_base_pair`: Base pair with alignment issues
- `sample_twisted_base_pair`: Base pair with twist issues
- `sample_non_coplanar_base_pair`: Base pair with coplanarity issues
- `sample_hbond_data`: Good quality H-bonds (3 bonds)
- `sample_bad_hbond_data`: Problematic H-bonds
- `sample_basepair_list`: Multiple base pairs for structure testing
- `empty_hbond_data`: Empty H-bond DataFrame

## Writing New Tests

When adding new functionality, follow this pattern:

```python
def test_new_feature(self, config):
    """Test description explaining what is being tested."""
    # Setup
    scorer = Scorer(config)

    # Execute
    result = scorer.new_method(test_data)

    # Assert
    assert result.expected_value == actual_value
```

## Files Being Tested

These tests focus on the actual execution path through:
- [scorer2.py](../scorer2.py) - Main scoring engine (Scorer class)
- [utils/data_loader.py](../utils/data_loader.py) - Data loading utilities
- [config.py](../config.py) - Configuration with thresholds and weights

**Note**: The `_bc` analyzer files (base_pair_analyzer_bc.py, hbond_analyzer_bc.py) are NOT tested because they're not used in the main execution path through scorer2.py.

## Test Results

All 32 tests pass successfully:
- 16 tests for scorer2.py
- 15 tests for data_loader.py

Run `pytest -v` to see detailed test results.
