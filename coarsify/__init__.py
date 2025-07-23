"""
Coarsify - A Python tool for coarse-graining molecular structures.

This package provides tools for converting detailed molecular structures
into simplified coarse-grained representations suitable for molecular
dynamics simulations and structural analysis.
"""

from .src.version import __version__
from .src.system.system import System as run


__all__ = [
    "run",
    "__version__",
]
