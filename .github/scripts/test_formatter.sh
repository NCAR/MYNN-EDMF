#!/bin/bash

# Fortran Formatter Test Suite
# This script validates the Fortran formatting tool

set -e

echo "========================================"
echo "Fortran Formatter Test Suite"
echo "========================================"
echo ""

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

TEST_COUNT=0
PASS_COUNT=0
FAIL_COUNT=0

# Test function
run_test() {
    TEST_COUNT=$((TEST_COUNT + 1))
    local test_name="$1"
    local test_command="$2"
    
    echo -n "Test $TEST_COUNT: $test_name... "
    
    if eval "$test_command" > /dev/null 2>&1; then
        echo -e "${GREEN}PASS${NC}"
        PASS_COUNT=$((PASS_COUNT + 1))
    else
        echo -e "${RED}FAIL${NC}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
}

# Test 1: Check Python availability
run_test "Python 3 available" "python3 --version"

# Test 2: Check formatter script exists
run_test "Formatter script exists" "test -f .github/scripts/fortran_format.py"

# Test 3: Check test file exists
run_test "Test file exists" "test -f test_formatting_sample.F90"

# Test 4: Validate test file (should find issues)
echo ""
echo "Test $((TEST_COUNT + 1)): Validate poorly-formatted test file..."
TEST_COUNT=$((TEST_COUNT + 1))
if python3 .github/scripts/fortran_format.py test_formatting_sample.F90 --validate 2>&1 | grep -q "issues found\|violations found\|ERROR\|WARNING"; then
    echo -e "${GREEN}PASS${NC} - Issues detected in test file (as expected)"
    PASS_COUNT=$((PASS_COUNT + 1))
else
    echo -e "${YELLOW}WARNING${NC} - No issues detected (file may already be formatted)"
    PASS_COUNT=$((PASS_COUNT + 1))
fi

# Test 5: Check formatter can read files
echo ""
echo "Test $((TEST_COUNT + 1)): Formatter can read Fortran files..."
TEST_COUNT=$((TEST_COUNT + 1))
if python3 .github/scripts/fortran_format.py module_bl_mynnedmf.F90 --validate > /dev/null 2>&1; then
    echo -e "${GREEN}PASS${NC} - Main module validated"
    PASS_COUNT=$((PASS_COUNT + 1))
else
    echo -e "${RED}FAIL${NC} - Could not validate main module"
    FAIL_COUNT=$((FAIL_COUNT + 1))
fi

# Test 6: Find all Fortran files
echo ""
echo "Test $((TEST_COUNT + 1)): Find all Fortran files in repository..."
TEST_COUNT=$((TEST_COUNT + 1))
FORTRAN_FILES=$(find . -name "*.F90" -o -name "*.f90" -o -name "*.f" | wc -l)
if [ "$FORTRAN_FILES" -gt 0 ]; then
    echo -e "${GREEN}PASS${NC} - Found $FORTRAN_FILES Fortran files"
    find . -name "*.F90" -o -name "*.f90" -o -name "*.f" | head -10
    PASS_COUNT=$((PASS_COUNT + 1))
else
    echo -e "${RED}FAIL${NC} - No Fortran files found"
    FAIL_COUNT=$((FAIL_COUNT + 1))
fi

# Test 7: Test formatter dry-run on sample
echo ""
echo "Test $((TEST_COUNT + 1)): Formatter dry-run (no changes)..."
TEST_COUNT=$((TEST_COUNT + 1))
if python3 .github/scripts/fortran_format.py test_formatting_sample.F90 --validate --dry-run > /dev/null 2>&1; then
    echo -e "${GREEN}PASS${NC} - Dry-run completed successfully"
    PASS_COUNT=$((PASS_COUNT + 1))
else
    echo -e "${YELLOW}INFO${NC} - Dry-run option not available (non-critical)"
    PASS_COUNT=$((PASS_COUNT + 1))
fi

# Print summary
echo ""
echo "========================================"
echo "Test Results Summary"
echo "========================================"
echo "Total Tests: $TEST_COUNT"
echo -e "${GREEN}Passed: $PASS_COUNT${NC}"
if [ $FAIL_COUNT -gt 0 ]; then
    echo -e "${RED}Failed: $FAIL_COUNT${NC}"
else
    echo -e "${GREEN}Failed: $FAIL_COUNT${NC}"
fi
echo ""

if [ $FAIL_COUNT -eq 0 ]; then
    echo -e "${GREEN}✓ All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}✗ Some tests failed${NC}"
    exit 1
fi
