#!/usr/bin/env python3
"""
Windows-specific build and upload script for PyPI distribution.
Uses the 'py' launcher instead of 'python'.
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path

def run_command(cmd, check=True):
    """Run a command and return the result."""
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if check and result.returncode != 0:
        print(f"Error: {result.stderr}")
        sys.exit(1)
    return result

def clean_build():
    """Clean previous build artifacts."""
    print("Cleaning previous build artifacts...")
    dirs_to_clean = ["build", "dist", "*.egg-info"]
    for dir_name in dirs_to_clean:
        if os.path.exists(dir_name):
            shutil.rmtree(dir_name)
            print(f"Removed {dir_name}")

def build_package():
    """Build the package."""
    print("Building package...")
    run_command("py -m build")

def check_package():
    """Check the package for common issues."""
    print("Checking package...")
    run_command("py -m twine check dist/*")

def upload_to_testpypi():
    """Upload to TestPyPI."""
    print("Uploading to TestPyPI...")
    run_command("py -m twine upload --repository testpypi dist/*")

def upload_to_pypi():
    """Upload to PyPI."""
    print("Uploading to PyPI...")
    run_command("py -m twine upload dist/*")

def main():
    """Main function."""
    if len(sys.argv) < 2:
        print("Usage: py build_windows.py [clean|build|check|test|upload]")
        print("  clean  - Clean build artifacts")
        print("  build  - Build the package")
        print("  check  - Check the package")
        print("  test   - Upload to TestPyPI")
        print("  upload - Upload to PyPI")
        sys.exit(1)

    command = sys.argv[1]

    if command == "clean":
        clean_build()
    elif command == "build":
        clean_build()
        build_package()
    elif command == "check":
        build_package()
        check_package()
    elif command == "test":
        clean_build()
        build_package()
        check_package()
        upload_to_testpypi()
    elif command == "upload":
        clean_build()
        build_package()
        check_package()
        upload_to_pypi()
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)

if __name__ == "__main__":
    main() 