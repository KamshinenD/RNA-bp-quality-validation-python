#!/bin/bash
# Helper script to run tests for RNA structure quality scoring

echo "================================================"
echo "RNA Structure Quality Scoring - Test Suite"
echo "================================================"
echo ""

# Check if pytest is installed
if ! command -v pytest &> /dev/null; then
    echo "Error: pytest is not installed"
    echo "Install it with: pip install -r requirements-test.txt"
    exit 1
fi

echo "Running all tests..."
echo ""

# Run pytest with coverage
pytest -v --tb=short

echo ""
echo "================================================"
echo "Test run complete!"
echo ""
echo "For more options:"
echo "  - Run specific test file: pytest tests/test_scorer.py"
echo "  - Run with coverage: pytest --cov=. --cov-report=html"
echo "  - Run specific test: pytest tests/test_scorer.py::TestScorer::test_score_perfect_base_pair"
echo "================================================"
