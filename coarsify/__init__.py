"""
Coarsify - A Python tool for coarse-graining molecular structures.

This package provides tools for converting detailed molecular structures
into simplified coarse-grained representations suitable for molecular
dynamics simulations and structural analysis.
"""

from .src.version import __version__

# Import main classes for easy access
try:
    from .src.system.system import System
    from .src.gui.GUI import settings_gui
except ImportError:
    # Handle case where dependencies aren't available
    pass

__all__ = [
    "System",
    "settings_gui",
    "__version__",
]