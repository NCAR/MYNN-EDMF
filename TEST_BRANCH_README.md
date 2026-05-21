# Fortran Formatter Test Branch

## Quick Start

### 1. Setup Environment
```bash
# Ensure Python 3 is installed
python3 --version

# Make scripts executable
chmod +x .github/scripts/test_formatter.sh
chmod +x .github/scripts/fortran_format.py
```

### 2. Run Tests
```bash
# Execute full test suite
bash .github/scripts/test_formatter.sh

# Expected output shows 7 tests with pass/fail results
```

### 3. Test Specific Files
```bash
# Validate poorly-formatted test file
python3 .github/scripts/fortran_format.py test_formatting_sample.F90 --validate

# Format the test file
python3 .github/scripts/fortran_format.py test_formatting_sample.F90 --format

# Validate main module
python3 .github/scripts/fortran_format.py module_bl_mynnedmf.F90 --validate
```

## Files in This Branch

### Test Infrastructure
- `.github/scripts/test_formatter.sh` - Main test suite (7 tests)
- `.github/scripts/fortran_format.py` - Python formatter/validator
- `test_formatting_sample.F90` - Poorly-formatted test file

### Documentation
- `FORMATTER_TEST_REPORT.md` - Detailed test results and standards
- `TEST_BRANCH_README.md` - This file

## Test Coverage

| Test | Purpose | Expected Result |
|------|---------|------------------|
| 1 | Python availability | ✓ PASS |
| 2 | Script exists | ✓ PASS |
| 3 | Test file exists | ✓ PASS |
| 4 | Detect formatting issues | ✓ PASS (finds issues) |
| 5 | Validate main module | ✓ PASS or ⚠ Issues found |
| 6 | Find Fortran files | ✓ PASS (lists files) |
| 7 | Dry-run functionality | ✓ PASS |

## Formatter Capabilities

✅ **Validates**:
- Indentation (3-space multiples)
- Line length (79 character limit)
- Trailing whitespace
- Variable declarations
- Operator spacing

⏳ **Future Enhancements**:
- Automatic line breaking
- Continuation alignment
- Comment formatting
- Case normalization

## Command Reference

```bash
# Validate file (check only, no changes)
python3 .github/scripts/fortran_format.py file.F90 --validate

# Format file (apply changes)
python3 .github/scripts/fortran_format.py file.F90 --format

# Preview changes (dry-run)
python3 .github/scripts/fortran_format.py file.F90 --format --dry-run

# Run test suite
bash .github/scripts/test_formatter.sh
```

## Expected Behaviors

### Validation Success
```bash
$ python3 .github/scripts/fortran_format.py module_bl_mynnedmf.F90 --validate
Validating: module_bl_mynnedmf.F90
✓ module_bl_mynnedmf.F90: No formatting issues found

Validation Status: PASSED
No formatting issues detected
```

### Validation with Issues
```bash
$ python3 .github/scripts/fortran_format.py test_formatting_sample.F90 --validate
Validating: test_formatting_sample.F90
✗ test_formatting_sample.F90: 15 formatting issues found
  Line 9: Indentation not multiple of 3 (found 1)
  Line 10: Indentation not multiple of 3 (found 5)
  Line 13: Line too long (102 > 79)
  ...

Validation Status: FAILED
Issues found: 15
```

## Workflow

1. **Local Testing**:
   - Developer checks out this branch
   - Runs `bash .github/scripts/test_formatter.sh`
   - Validates specific files

2. **Feedback**:
   - If tests pass ✓: Code is formatted correctly
   - If tests fail ✗: Shows what needs fixing
   - Run formatter: `python3 .github/scripts/fortran_format.py file.F90 --format`

3. **Integration**:
   - Once approved, merge to `setup/fortran-formatter` branch
   - Then merge to `main`
   - GitHub Actions will use these scripts in PR validation

## Troubleshooting

### "No such file or directory"
```bash
# Make sure you're in repo root
cd /path/to/NCAR/MYNN-EDMF

# Verify files exist
ls -la .github/scripts/
ls -la test_formatting_sample.F90
```

### "Python: command not found"
```bash
# Install Python 3
sudo apt install python3      # Ubuntu/Debian
brew install python3          # macOS
choco install python3         # Windows

# Or use python3 explicitly
python3 .github/scripts/fortran_format.py file.F90 --validate
```

### "Permission denied"
```bash
# Make scripts executable
chmod +x .github/scripts/test_formatter.sh
chmod +x .github/scripts/fortran_format.py
```

## Next Steps

1. ✅ **Test**: Run `bash .github/scripts/test_formatter.sh`
2. ✅ **Validate**: Test on actual code files
3. ⏳ **Review**: Check results and formatting output
4. ⏳ **Approve**: If satisfied, merge to main
5. ⏳ **Deploy**: GitHub Actions uses these in production

## Support

For issues or questions:
1. Check FORMATTER_TEST_REPORT.md for detailed standards
2. Review test script: `.github/scripts/test_formatter.sh`
3. Examine formatter: `.github/scripts/fortran_format.py`

---

**Branch**: `test/fortran-formatter-validation`  
**Status**: Ready for testing and validation  
**Last Updated**: 2026-05-21
