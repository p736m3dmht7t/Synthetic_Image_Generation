#!/usr/bin/env python3
"""
FITS Header Extraction Utility

This script extracts FITS headers from user-provided example files and creates
a JSON template for use in synthetic image generation.
"""

import json
import sys
from pathlib import Path
from astropy.io import fits
import argparse


def extract_fits_header(fits_file_path, output_json_path=None):
    """
    Extract FITS header from a file and save as JSON template.
    
    Args:
        fits_file_path (str): Path to the FITS file
        output_json_path (str): Path for output JSON file (optional)
    
    Returns:
        dict: Extracted header as dictionary
    """
    try:
        # Open the FITS file
        with fits.open(fits_file_path) as hdul:
            # Get the primary header
            header = hdul[0].header
            
            # Convert header to dictionary, preserving key order
            header_dict = {}
            for key in header.keys():
                if key == 'COMMENT' or key == 'HISTORY':
                    # Handle multiple COMMENT/HISTORY entries
                    if key not in header_dict:
                        header_dict[key] = []
                    header_dict[key].extend(header[key])
                else:
                    value = header[key]
                    # Convert numpy types to Python types for JSON serialization
                    if hasattr(value, 'item'):
                        value = value.item()
                    header_dict[key] = value
            
            # If output path not specified, create default
            if output_json_path is None:
                output_json_path = Path("config/fits_header.json")
            
            # Ensure output directory exists
            output_path = Path(output_json_path)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Save to JSON file
            with open(output_path, 'w') as f:
                json.dump(header_dict, f, indent=2, default=str)
            
            print(f"FITS header extracted successfully!")
            print(f"Input file: {fits_file_path}")
            print(f"Output file: {output_path}")
            print(f"Extracted {len(header_dict)} header keywords")
            
            return header_dict
            
    except FileNotFoundError:
        print(f"Error: FITS file not found: {fits_file_path}")
        return None
    except Exception as e:
        print(f"Error extracting FITS header: {str(e)}")
        return None


def main():
    """Command line interface for FITS header extraction."""
    parser = argparse.ArgumentParser(
        description="Extract FITS header and save as JSON template"
    )
    parser.add_argument(
        "fits_file", 
        help="Path to the FITS file to extract header from"
    )
    parser.add_argument(
        "-o", "--output", 
        help="Output JSON file path (default: config/fits_header.json)"
    )
    
    args = parser.parse_args()
    
    # Extract the header
    result = extract_fits_header(args.fits_file, args.output)
    
    if result is None:
        sys.exit(1)


if __name__ == "__main__":
    main()