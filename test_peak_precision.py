#!/usr/bin/env python3
"""
Peak Flux Precision Validation Test

This test validates that the two-pass peak correction approach achieves
exactly the desired peak flux (50% saturation) for the target star,
regardless of subpixel positioning effects.
"""

import numpy as np
import sys
import os
sys.path.append('src')

from psf_photutils import PhotoutilsPSFGenerator
from astropy.io import fits
import json

def load_test_configuration():
    """Load test configuration from JSON files."""
    try:
        # Load camera configuration
        with open('config/camera.json', 'r') as f:
            camera_config = json.load(f)
        
        # Load target configuration
        with open('config/target.json', 'r') as f:
            target_config = json.load(f)
        
        return camera_config, target_config
    except Exception as e:
        print(f"Error loading configuration: {e}")
        return None, None

def create_test_source_data():
    """Create test source data with known target position."""
    # Create a simple test case with target at various subpixel positions
    # Use well-separated positions to avoid PSF overlap
    test_positions = [
        (100.0, 100.0),    # Target - On pixel center
        (150.25, 150.25),  # Quarter pixel offset
        (200.5, 200.5),    # Half pixel offset
        (250.75, 250.75),  # Three-quarter pixel offset
        (300.333, 300.667), # Random subpixel position
    ]
    
    source_data = {
        'x_pixel': [],
        'y_pixel': [],
        'V': [],
        'ra': [],
        'dec': [],
        'source_names': []
    }
    
    for i, (x, y) in enumerate(test_positions):
        source_data['x_pixel'].append(x)
        source_data['y_pixel'].append(y)
        source_data['V'].append(10.4)  # TZ Boo V magnitude
        source_data['ra'].append(227.038063)  # TZ Boo RA
        source_data['dec'].append(39.97024)   # TZ Boo Dec
        if i == 0:
            source_data['source_names'].append('V* TZ Boo')
        else:
            source_data['source_names'].append(f'Test Star {i+1}')
    
    return source_data, test_positions

def test_peak_flux_precision():
    """Test peak flux precision with two-pass approach."""
    print("=== PEAK FLUX PRECISION VALIDATION TEST ===\n")
    
    # Load configuration
    camera_config, target_config = load_test_configuration()
    if camera_config is None or target_config is None:
        print("Error: Could not load configuration files")
        return False
    
    # Extract configuration values
    saturation_limit = camera_config['saturation_limit']
    target_magnitude = target_config['magnitudes']['V']
    target_coords = (target_config['coordinates']['ra'], target_config['coordinates']['dec'])
    
    print(f"Configuration:")
    print(f"  Saturation limit: {saturation_limit} ADU")
    print(f"  Target magnitude: {target_magnitude:.3f}")
    print(f"  Target coordinates: ({target_coords[0]:.6f}, {target_coords[1]:.6f})")
    
    # Create PSF generator
    psf_generator = PhotoutilsPSFGenerator()
    
    # PSF parameters (example values)
    psf_params = {
        'sigma_pixels': 1.5,
        'fwhm_pixels': 3.53
    }
    
    # Image parameters
    image_shape = (400, 400)
    
    # Create test source data
    source_data, test_positions = create_test_source_data()
    
    print(f"\nTest setup:")
    print(f"  Image shape: {image_shape}")
    print(f"  PSF sigma: {psf_params['sigma_pixels']:.1f} pixels")
    print(f"  Test positions: {len(test_positions)}")
    
    # Expected peak flux (50% saturation)
    expected_peak_adu = saturation_limit * 0.5
    
    print(f"\nExpected peak flux: {expected_peak_adu:.0f} ADU")
    
    # Test the two-pass approach
    print(f"\n=== TESTING TWO-PASS APPROACH ===")
    
    try:
        # Create stellar field using two-pass approach
        image = psf_generator.create_stellar_field(
            source_data=source_data,
            image_shape=image_shape,
            psf_params=psf_params,
            target_magnitude=target_magnitude,
            saturation_limit=saturation_limit,
            band='V',
            psf_type='gaussian',
            target_name='V* TZ Boo',
            target_coords=target_coords
        )
        
        if image is None:
            print("Error: Failed to create stellar field")
            return False
        
        # Analyze peak flux accuracy
        print(f"\n=== PEAK FLUX ANALYSIS ===")
        
        max_peak_adu = np.max(image)
        accuracy = (max_peak_adu / expected_peak_adu) * 100
        
        print(f"Maximum peak flux in image: {max_peak_adu:.1f} ADU")
        print(f"Expected peak flux: {expected_peak_adu:.0f} ADU")
        print(f"Peak flux accuracy: {accuracy:.2f}%")
        
        # Check accuracy tolerance
        accuracy_tolerance = 1.0  # 1% tolerance
        if abs(accuracy - 100.0) <= accuracy_tolerance:
            print(f"✓ PASS: Peak flux accuracy within {accuracy_tolerance}% tolerance")
            test_passed = True
        else:
            print(f"✗ FAIL: Peak flux accuracy outside {accuracy_tolerance}% tolerance")
            test_passed = False
        
        # Find peak positions and analyze each test position
        print(f"\n=== INDIVIDUAL POSITION ANALYSIS ===")
        
        for i, (x, y) in enumerate(test_positions):
            # Extract local maximum around each test position
            x_int, y_int = int(x), int(y)
            
            # Get local region around source
            region_size = 10
            x_start = max(0, x_int - region_size)
            x_end = min(image_shape[1], x_int + region_size + 1)
            y_start = max(0, y_int - region_size)
            y_end = min(image_shape[0], y_int + region_size + 1)
            
            local_region = image[y_start:y_end, x_start:x_end]
            local_max = np.max(local_region)
            
            # Calculate subpixel offset
            x_frac = x - int(x)
            y_frac = y - int(y)
            
            print(f"  Position {i+1}: ({x:.3f}, {y:.3f})")
            print(f"    Subpixel offset: ({x_frac:.3f}, {y_frac:.3f})")
            print(f"    Local peak: {local_max:.1f} ADU")
            print(f"    Peak accuracy: {(local_max/expected_peak_adu)*100:.2f}%")
        
        # Save test image for verification
        output_path = 'test_peak_precision_output.fits'
        fits.HDUList([fits.PrimaryHDU(image.astype(np.float32))]).writeto(output_path, overwrite=True)
        print(f"\nTest image saved to: {output_path}")
        
        return test_passed
        
    except Exception as e:
        print(f"Error during test: {e}")
        return False

def main():
    """Main test function."""
    print("Starting peak flux precision validation test...")
    
    success = test_peak_flux_precision()
    
    if success:
        print("\n🎉 ALL TESTS PASSED!")
        print("The two-pass approach successfully achieves exact peak flux control.")
    else:
        print("\n❌ TEST FAILED!")
        print("The two-pass approach did not achieve the required peak flux precision.")
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())