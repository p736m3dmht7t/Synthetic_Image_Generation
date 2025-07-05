#!/usr/bin/env python3
"""
Subpixel Position Precision Test

This test specifically validates that the correction factor properly accounts
for subpixel positioning effects by testing with challenging subpixel offsets.
"""

import numpy as np
import sys
import os
sys.path.append('src')

from psf_photutils import PhotoutilsPSFGenerator
from astropy.io import fits
import json

def test_subpixel_correction_factor():
    """Test correction factor calculation with various subpixel positions."""
    print("=== SUBPIXEL CORRECTION FACTOR TEST ===\n")
    
    # Load configuration
    with open('config/camera.json', 'r') as f:
        camera_config = json.load(f)
    with open('config/target.json', 'r') as f:
        target_config = json.load(f)
    
    saturation_limit = camera_config['saturation_limit']
    target_magnitude = target_config['magnitudes']['V']
    
    # PSF parameters
    psf_params = {
        'sigma_pixels': 1.5,
        'fwhm_pixels': 3.53
    }
    
    # Test various challenging subpixel positions (centered in a 50x50 image)
    test_positions = [
        (25.0, 25.0),       # Pixel center (baseline)
        (25.5, 25.5),       # Half pixel offset (worst case for Gaussian)
        (25.25, 25.75),     # Asymmetric offset
        (25.1, 25.9),       # Small offset
        (25.8, 25.2),       # Large offset
        (25.333, 25.667),   # Random positions
    ]
    
    psf_generator = PhotoutilsPSFGenerator()
    
    print(f"Configuration:")
    print(f"  Target magnitude: {target_magnitude:.3f}")
    print(f"  Saturation limit: {saturation_limit} ADU")
    print(f"  Expected peak: {saturation_limit * 0.5:.0f} ADU")
    print(f"  PSF sigma: {psf_params['sigma_pixels']:.1f} pixels")
    
    print(f"\nTesting {len(test_positions)} subpixel positions:")
    
    results = []
    
    for i, (x, y) in enumerate(test_positions):
        print(f"\n--- Position {i+1}: ({x:.3f}, {y:.3f}) ---")
        
        # Calculate correction factor for this position
        correction_factor = psf_generator.calculate_peak_correction_factor(
            x, y, target_magnitude, saturation_limit, psf_params, 'gaussian'
        )
        
        # Create test source data with single source at this position
        source_data = {
            'x_pixel': [x],
            'y_pixel': [y],
            'V': [target_magnitude],
            'ra': [227.038063],
            'dec': [39.97024],
            'source_names': ['V* TZ Boo']
        }
        
        # Create stellar field with two-pass approach
        image = psf_generator.create_stellar_field(
            source_data=source_data,
            image_shape=(50, 50),
            psf_params=psf_params,
            target_magnitude=target_magnitude,
            saturation_limit=saturation_limit,
            band='V',
            psf_type='gaussian',
            target_name='V* TZ Boo'
        )
        
        # Analyze actual peak achieved
        actual_peak = np.max(image)
        desired_peak = saturation_limit * 0.5
        accuracy = (actual_peak / desired_peak) * 100
        
        results.append({
            'position': (x, y),
            'subpixel_offset': (x - int(x), y - int(y)),
            'correction_factor': correction_factor,
            'actual_peak': actual_peak,
            'accuracy': accuracy
        })
        
        print(f"Result: {actual_peak:.1f} ADU (accuracy: {accuracy:.2f}%)")
    
    # Summary analysis
    print(f"\n=== SUMMARY ANALYSIS ===")
    print(f"{'Position':<12} {'Offset':<12} {'CorrFactor':<12} {'Peak':<8} {'Accuracy'}")
    print("-" * 60)
    
    all_passed = True
    tolerance = 2.0  # 2% tolerance for subpixel effects
    
    for result in results:
        x, y = result['position']
        x_off, y_off = result['subpixel_offset']
        cf = result['correction_factor']
        peak = result['actual_peak']
        acc = result['accuracy']
        
        status = "✓" if abs(acc - 100.0) <= tolerance else "✗"
        if abs(acc - 100.0) > tolerance:
            all_passed = False
        
        print(f"({x:4.1f},{y:4.1f})  ({x_off:4.2f},{y_off:4.2f})  {cf:8.4f}     {peak:6.0f}   {acc:6.2f}% {status}")
    
    # Statistical analysis
    correction_factors = [r['correction_factor'] for r in results]
    accuracies = [r['accuracy'] for r in results]
    
    print(f"\nStatistics:")
    print(f"  Correction factor range: {min(correction_factors):.4f} - {max(correction_factors):.4f}")
    print(f"  Accuracy range: {min(accuracies):.2f}% - {max(accuracies):.2f}%")
    print(f"  Mean accuracy: {np.mean(accuracies):.2f}%")
    print(f"  Std accuracy: {np.std(accuracies):.2f}%")
    
    if all_passed:
        print(f"\n🎉 ALL SUBPIXEL TESTS PASSED!")
        print(f"The correction factor successfully compensates for subpixel positioning effects.")
    else:
        print(f"\n❌ SOME TESTS FAILED!")
        print(f"The correction factor may need refinement for certain subpixel positions.")
    
    return all_passed

def main():
    """Main test function."""
    success = test_subpixel_correction_factor()
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())