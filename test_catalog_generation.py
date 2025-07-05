#!/usr/bin/env python3
"""
Test Catalog Generation

Simple test to verify that the catalog documentation generation works correctly.
"""

import sys
import os
sys.path.append('src')

from image_generator import ImageGenerator
from config_manager import ConfigManager
import numpy as np

def test_catalog_generation():
    """Test the catalog generation functionality."""
    print("=== TESTING CATALOG GENERATION ===\n")
    
    # Initialize components
    config_manager = ConfigManager()
    image_generator = ImageGenerator(config_manager)
    
    # Load configurations
    if not image_generator.load_configurations():
        print("Error: Could not load configurations")
        return False
    
    # Create mock source data that mimics what comes from GAIA queries
    mock_source_data = {
        'source_id': [1391924687894823296, 1296666402078585984, 1391687949297517824],
        'ra': [227.038063, 227.042, 227.046],
        'dec': [39.97024, 39.972, 39.974],
        'U': [11.345, 11.153, 11.930],
        'B': [11.282, 11.126, 11.475],
        'V': [10.649, 10.581, 10.686],
        'R': [10.332, 10.254, 10.256],
        'I': [9.929, 9.943, 9.878],
        'U_error': [0.025, 0.010, 0.008],
        'B_error': [0.007, 0.002, 0.002],
        'V_error': [0.009, 0.002, 0.001],
        'R_error': [0.010, 0.001, 0.001],
        'I_error': [0.003, 0.0003, 0.0003],
        'method': 'spectroscopic',
        'spectroscopic_bands': {
            0: ['U', 'B', 'V', 'R', 'I'],  # All bands spectroscopic for first source
            1: ['B', 'V', 'R', 'I'],       # U is polynomial for second source
            2: []                          # All bands polynomial for third source
        }
    }
    
    # Set mock source data
    image_generator.source_data = mock_source_data
    
    # Test catalog generation (no longer band-specific)
    test_fits_filename = "V*_TZ_Boo_20250105_120000_V.fits"
    
    print(f"Testing catalog generation for: {test_fits_filename}")
    print(f"Expected catalog name: V*_TZ_Boo_20250105_120000.catalog (unified)")
    
    try:
        catalog_path = image_generator.generate_catalog_documentation(test_fits_filename)
        
        if catalog_path:
            print(f"✓ Catalog generated successfully: {catalog_path}")
            
            # Verify the file exists and read a few lines
            if os.path.exists(catalog_path):
                with open(catalog_path, 'r') as f:
                    lines = f.readlines()
                
                print(f"✓ Catalog file created with {len(lines)} lines")
                print(f"✓ File size: {os.path.getsize(catalog_path)} bytes")
                
                # Show first few lines
                print("\nFirst 10 lines of catalog:")
                print("-" * 50)
                for i, line in enumerate(lines[:10]):
                    print(f"{i+1:2d}: {line.rstrip()}")
                
                # Show a data line
                print("\nSample data lines:")
                print("-" * 50)
                for i, line in enumerate(lines):
                    if not line.startswith('#') and not line.startswith('-') and line.strip():
                        print(f"Data: {line.rstrip()}")
                        if i > 15:  # Show first few data lines
                            break
                
                return True
            else:
                print(f"✗ Catalog file was not created: {catalog_path}")
                return False
        else:
            print("✗ Catalog generation returned None")
            return False
            
    except Exception as e:
        print(f"✗ Error during catalog generation: {e}")
        return False

def main():
    """Main test function."""
    success = test_catalog_generation()
    
    if success:
        print("\n🎉 CATALOG GENERATION TEST PASSED!")
        print("The catalog documentation system is working correctly.")
    else:
        print("\n❌ CATALOG GENERATION TEST FAILED!")
        print("There was an issue with the catalog documentation system.")
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())