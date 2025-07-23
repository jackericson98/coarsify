#!/usr/bin/env python3
"""
Main entry point for the coarsify command-line interface.
"""

import os
import sys

# Add the parent directory to sys.path
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

# Add everything in /api/ to the module search path.
__path__ = [os.path.dirname(__file__), os.path.join(os.path.dirname(__file__), "api")]

def main():
    """Main function for the coarsify CLI."""
    try:
        from coarsify.src.system.system import System
        app = System()
        # Suppress the System object output by not returning it
        return 0
    except ImportError as e:
        print(f"Error importing coarsify modules: {e}")
        print("Please ensure coarsify is properly installed.")
        return 1
    except Exception as e:
        print(f"Error running coarsify: {e}")
        return 1

if __name__ == '__main__':
    main()
