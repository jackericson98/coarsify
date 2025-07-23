from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    with open("README.md", "r", encoding="utf-8") as fh:
        return fh.read()

# Read version from the version file
def get_version():
    try:
        version_file = os.path.join("coarsify", "src", "version.py")
        with open(version_file, "r") as f:
            exec(f.read())
        version = locals()["__version__"]
        if version == "unknown":
            return "1.0.2"
        return version
    except:
        return "1.0.2"

setup(
    name="coarsify",
    version=get_version(),
    author="John Ericson",
    author_email="jackericson98@gmail.com",
    description="A Python tool for coarse-graining molecular structures",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/jackericson98/coarsify",
    project_urls={
        "Bug Tracker": "https://github.com/jackericson98/coarsify/issues",
        "Documentation": "https://github.com/jackericson98/coarsify#readme",
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "coarsify"},
    packages=find_packages(where="coarsify"),
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.24.0",
        "pandas>=2.0.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=22.0",
            "flake8>=4.0",
            "mypy>=0.950",
        ],
        "docs": [
            "sphinx>=4.0",
            "sphinx-rtd-theme>=1.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "coarsify=coarsify.__main__:main",
            "coarsify-gui=coarsify.src.gui.gui:main",
        ],
    },
    include_package_data=True,
    package_data={
        "coarsify": [
            "data/*",
            "src/data/*",
        ],
    },
    keywords="molecular dynamics, coarse-graining, protein structure, chemistry, physics",
    zip_safe=False,
) 