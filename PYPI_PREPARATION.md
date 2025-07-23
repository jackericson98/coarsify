# PyPI Preparation Guide

This guide will help you prepare and upload the coarsify package to PyPI.

## Prerequisites

1. **PyPI Account**: Create an account on [PyPI](https://pypi.org/account/register/)
2. **TestPyPI Account**: Create an account on [TestPyPI](https://test.pypi.org/account/register/)
3. **API Tokens**: Generate API tokens for both PyPI and TestPyPI
4. **Development Tools**: Install required development packages

## Setup

### 1. Install Development Dependencies

```bash
pip install -r requirements-dev.txt
```

### 2. Configure PyPI Credentials

Create a `~/.pypirc` file with your credentials:

```ini
[distutils]
index-servers =
    pypi
    testpypi

[pypi]
repository = https://upload.pypi.org/legacy/
username = __token__
password = pypi-your-api-token-here

[testpypi]
repository = https://test.pypi.org/legacy/
username = __token__
password = pypi-your-test-api-token-here
```

### 3. Update Package Information

Before uploading, update the following files with your information:

- `setup.py`: Update author, email, and URLs
- `pyproject.toml`: Update author, email, and URLs
- `LICENSE`: Update copyright holder
- `README.md`: Update GitHub URLs and contact information

## Building and Testing

### 1. Clean Previous Builds

**Windows:**
```cmd
py build_windows.py clean
```

**Linux/Mac:**
```bash
python build_and_upload.py clean
```

### 2. Build Package

**Windows:**
```cmd
py build_windows.py build
```

**Linux/Mac:**
```bash
python build_and_upload.py build
```

### 3. Check Package

**Windows:**
```cmd
py build_windows.py check
```

**Linux/Mac:**
```bash
python build_and_upload.py check
```

This will verify that your package can be built and installed correctly.

## Testing on TestPyPI

### 1. Upload to TestPyPI

**Windows:**
```cmd
py build_windows.py test
```

**Linux/Mac:**
```bash
python build_and_upload.py test
```

### 2. Test Installation

```bash
pip install --index-url https://test.pypi.org/simple/ coarsify
```

### 3. Test Functionality

```bash
coarsify --help
coarsify-gui
```

## Uploading to PyPI

### 1. Final Check

**Windows:**
```cmd
py build_windows.py check
```

**Linux/Mac:**
```bash
python build_and_upload.py check
```

### 2. Upload to PyPI

**Windows:**
```cmd
py build_windows.py upload
```

**Linux/Mac:**
```bash
python build_and_upload.py upload
```

## Manual Commands

If you prefer to run commands manually:

**Windows:**
```cmd
# Clean
rmdir /s /q build dist *.egg-info

# Build
py -m build

# Check
py -m twine check dist/*

# Upload to TestPyPI
py -m twine upload --repository testpypi dist/*

# Upload to PyPI
py -m twine upload dist/*
```

**Linux/Mac:**
```bash
# Clean
rm -rf build/ dist/ *.egg-info/

# Build
python -m build

# Check
twine check dist/*

# Upload to TestPyPI
twine upload --repository testpypi dist/*

# Upload to PyPI
twine upload dist/*
```

## Version Management

The package uses `setuptools_scm` for automatic versioning based on git tags.

### Creating a New Release

1. **Update Version**: Create a new git tag
   ```bash
   git tag v1.0.1
   git push origin v1.0.1
   ```

2. **Build and Upload**: Use the build script
   ```bash
   python build_and_upload.py upload
   ```

## Troubleshooting

### Common Issues

1. **Authentication Errors**: Check your `~/.pypirc` file and API tokens
2. **Package Name Conflicts**: Ensure the package name is unique on PyPI
3. **Build Errors**: Check that all dependencies are properly specified
4. **Import Errors**: Verify that all modules are properly included in the package

### Getting Help

- [PyPI Help](https://pypi.org/help/)
- [Python Packaging User Guide](https://packaging.python.org/)
- [setuptools Documentation](https://setuptools.pypa.io/)

## Security Notes

- Never commit your `~/.pypirc` file to version control
- Use API tokens instead of passwords
- Regularly rotate your API tokens
- Test on TestPyPI before uploading to PyPI

## Next Steps

After successful upload:

1. **Verify Installation**: Test installation from PyPI
2. **Update Documentation**: Update README with PyPI installation instructions
3. **Create Release Notes**: Document changes in the new version
4. **Announce Release**: Share the release with your community 