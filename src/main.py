#!/usr/bin/env python3
"""
Main entry point for the Synthetic Image Generation application.
"""

import sys
import os
from pathlib import Path

# Add src directory to path
sys.path.insert(0, os.path.dirname(__file__))

from gui.main_window import main

if __name__ == "__main__":
    main()