#!/usr/bin/env python3
"""
Fortran Code Formatter and Validator
Enforces consistent Fortran 90 free-form code formatting standards.

Usage:
    python3 fortran_format.py <file.f90> [--validate] [--format] [--dry-run]

Options:
    --validate   : Check for formatting issues without modifying
    --format     : Apply formatting changes to file
    --dry-run    : Show what would be formatted without applying changes
"""

import sys
import re
from pathlib import Path
from typing import List, Tuple

class FortranFormatter:
    """Validates and formats Fortran 90 code."""
    
    # Configuration
    INDENT_SIZE = 3
    MAX_LINE_LENGTH = 79
    
    def __init__(self, filename: str):
        self.filename = filename
        self.lines = []
        self.issues = []
        self.load_file()
    
    def load_file(self):
        """Load the Fortran file."""
        try:
            with open(self.filename, 'r', encoding='utf-8') as f:
                self.lines = f.readlines()
        except FileNotFoundError:
            print(f"ERROR: File not found: {self.filename}")
            sys.exit(1)
        except Exception as e:
            print(f"ERROR: Could not read file: {e}")
            sys.exit(1)
    
    def save_file(self, lines: List[str]):
        """Save formatted lines to file."""
        try:
            with open(self.filename, 'w', encoding='utf-8') as f:
                f.writelines(lines)
        except Exception as e:
            print(f"ERROR: Could not write file: {e}")
            sys.exit(1)
    
    def validate(self) -> bool:
        """Validate file formatting. Returns True if valid."""
        self.issues = []
        
        for line_num, line in enumerate(self.lines, 1):
            self._check_line(line_num, line)
        
        return len(self.issues) == 0
    
    def _check_line(self, line_num: int, line: str):
        """Check a single line for formatting issues."""
        # Skip empty lines and comments
        if not line.strip() or line.strip().startswith('!'):
            return
        
        # Check for trailing whitespace
        if line.rstrip() != line.rstrip('\n'):
            self.issues.append((line_num, "Trailing whitespace"))
        
        # Check line length (excluding newline)
        line_content = line.rstrip('\n')
        if len(line_content) > self.MAX_LINE_LENGTH:
            self.issues.append((line_num, f"Line too long ({len(line_content)} > {self.MAX_LINE_LENGTH})"))
        
        # Check indentation (lines that aren't continuations)
        if not line_content.startswith('&'):
            indent = len(line) - len(line.lstrip())
            if indent > 0 and indent % self.INDENT_SIZE != 0:
                self.issues.append((line_num, f"Indentation not multiple of {self.INDENT_SIZE} (found {indent})"))
    
    def format_file(self) -> List[str]:
        """Format the file and return formatted lines."""
        formatted_lines = []
        
        for line in self.lines:
            formatted_line = self._format_line(line)
            formatted_lines.append(formatted_line)
        
        return formatted_lines
    
    def _format_line(self, line: str) -> str:
        """Format a single line."""
        # Remove trailing whitespace but preserve newline
        has_newline = line.endswith('\n')
        line_content = line.rstrip('\n').rstrip()
        
        # Fix indentation if not a continuation line
        if line_content and not line_content.startswith('&'):
            # Count leading spaces
            indent = len(line_content) - len(line_content.lstrip())
            if indent > 0:
                # Round to nearest multiple of INDENT_SIZE
                new_indent = (indent // self.INDENT_SIZE) * self.INDENT_SIZE
                if indent != new_indent:
                    line_content = ' ' * new_indent + line_content.lstrip()
        
        # Add newline back
        if has_newline:
            return line_content + '\n'
        else:
            return line_content
    
    def print_issues(self):
        """Print all formatting issues found."""
        if not self.issues:
            print(f"✓ {self.filename}: No formatting issues found")
            return
        
        print(f"✗ {self.filename}: {len(self.issues)} formatting issues found\n")
        for line_num, issue in self.issues:
            print(f"  Line {line_num}: {issue}")
    
    def print_summary(self):
        """Print validation summary."""
        if self.issues:
            print(f"\nValidation Status: FAILED")
            print(f"Issues found: {len(self.issues)}")
            return False
        else:
            print(f"\nValidation Status: PASSED")
            print(f"No formatting issues detected")
            return True


def main():
    """Main entry point."""
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    
    filename = sys.argv[1]
    validate_mode = '--validate' in sys.argv
    format_mode = '--format' in sys.argv
    dry_run = '--dry-run' in sys.argv
    
    # Default to validate if no mode specified
    if not validate_mode and not format_mode:
        validate_mode = True
    
    formatter = FortranFormatter(filename)
    
    if validate_mode:
        print(f"Validating: {filename}")
        is_valid = formatter.validate()
        formatter.print_issues()
        formatter.print_summary()
        sys.exit(0 if is_valid else 1)
    
    if format_mode:
        print(f"Formatting: {filename}")
        formatted_lines = formatter.format_file()
        
        if dry_run:
            print("DRY RUN - No changes applied")
            print("Preview of changes:")
            for i, (orig, formatted) in enumerate(zip(self.lines[:5], formatted_lines[:5]), 1):
                if orig != formatted:
                    print(f"  Line {i}: {repr(orig)} -> {repr(formatted)}")
        else:
            formatter.save_file(formatted_lines)
            print(f"✓ File formatted successfully")
        
        # Validate after formatting
        formatter.lines = formatted_lines
        is_valid = formatter.validate()
        if is_valid:
            print("✓ Formatted file passes validation")
        else:
            print("⚠ Formatted file still has issues")
            formatter.print_issues()


if __name__ == '__main__':
    main()
