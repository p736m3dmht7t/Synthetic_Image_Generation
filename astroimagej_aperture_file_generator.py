#!/usr/bin/env python3
"""
AstroImageJ Aperture File Generator

Reads AstroImageJ .radec files and applies GAIA spectroscopic photometry
to generate accurate U, B, V, R, I magnitudes and uncertainties.

Usage:
    python astroimagej_aperture_file_generator.py input_file.radec
"""

import sys
import os
import re
import numpy as np
from pathlib import Path
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.gaia import Gaia
from gaiaxpy import generate, PhotometricSystem
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=UserWarning)


class AstroImageJProcessor:
    """Process AstroImageJ .radec files with GAIA photometry."""
    
    def __init__(self):
        self.photometric_system = PhotometricSystem.JKC_Std
        self.search_radius = 10.0  # arcseconds for GAIA source matching
        
    def parse_sexagesimal_coordinates(self, ra_str, dec_str):
        """
        Parse sexagesimal coordinates to decimal degrees.
        
        Parameters
        ----------
        ra_str : str
            RA in format "HH:MM:SS.SSS"
        dec_str : str
            Dec in format "+DD:MM:SS.SS" or "-DD:MM:SS.SS"
            
        Returns
        -------
        tuple
            (ra_deg, dec_deg) in decimal degrees
        """
        try:
            # Create SkyCoord object from sexagesimal strings
            coord = SkyCoord(ra=ra_str, dec=dec_str, unit=(u.hourangle, u.deg))
            return coord.ra.deg, coord.dec.deg
        except Exception as e:
            print(f"Error parsing coordinates {ra_str}, {dec_str}: {e}")
            return None, None
    
    def decimal_to_sexagesimal(self, ra_deg, dec_deg):
        """
        Convert decimal coordinates back to sexagesimal format.
        
        Parameters
        ----------
        ra_deg : float
            RA in decimal degrees
        dec_deg : float
            Dec in decimal degrees
            
        Returns
        -------
        tuple
            (ra_str, dec_str) in AstroImageJ format
        """
        coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
        ra_str = coord.ra.to_string(unit=u.hour, sep=':', precision=3)
        dec_str = coord.dec.to_string(unit=u.deg, sep=':', precision=2, alwayssign=True)
        return ra_str, dec_str
    
    def parse_radec_file(self, filepath):
        """
        Parse AstroImageJ .radec file.
        
        Parameters
        ----------
        filepath : str
            Path to input .radec file
            
        Returns
        -------
        list
            List of dictionaries with parsed data
        """
        sources = []
        
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
            
            for line_num, line in enumerate(lines, 1):
                line = line.strip()
                
                # Skip comments and empty lines
                if line.startswith('#') or not line:
                    continue
                
                # Parse comma-separated values
                parts = [p.strip() for p in line.split(',')]
                if len(parts) >= 5:
                    ra_str = parts[0]
                    dec_str = parts[1]
                    ref_star = int(parts[2]) if parts[2] != '' else None
                    centroid = int(parts[3]) if parts[3] != '' else None
                    magnitude = float(parts[4]) if parts[4] != '' else 99.999
                    
                    # Convert to decimal coordinates
                    ra_deg, dec_deg = self.parse_sexagesimal_coordinates(ra_str, dec_str)
                    
                    if ra_deg is not None and dec_deg is not None:
                        sources.append({
                            'line_num': line_num,
                            'ra_orig_str': ra_str,
                            'dec_orig_str': dec_str,
                            'ra_deg': ra_deg,
                            'dec_deg': dec_deg,
                            'ref_star': ref_star,
                            'centroid': centroid,
                            'orig_magnitude': magnitude
                        })
                    else:
                        print(f"Warning: Could not parse coordinates on line {line_num}: {line}")
                        
        except Exception as e:
            print(f"Error reading file {filepath}: {e}")
            return []
        
        return sources
    
    def find_gaia_source(self, ra_deg, dec_deg):
        """
        Find nearest GAIA source for given coordinates.
        
        Parameters
        ----------
        ra_deg : float
            RA in decimal degrees
        dec_deg : float
            Dec in decimal degrees
            
        Returns
        -------
        dict or None
            GAIA source info or None if not found
        """
        try:
            coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)
            
            # Query GAIA with search radius
            result = Gaia.query_object_async(
                coordinate=coord,
                width=self.search_radius * u.arcsec,
                height=self.search_radius * u.arcsec
            )
            
            if result and len(result) > 0:
                # Return the nearest source
                source = result[0]
                return {
                    'source_id': int(source['source_id']),
                    'ra': float(source['ra']),
                    'dec': float(source['dec']),
                    'g_mag': float(source['phot_g_mean_mag']) if source['phot_g_mean_mag'] else None
                }
            
        except Exception as e:
            print(f"Error querying GAIA for {ra_deg:.6f}, {dec_deg:.6f}: {e}")
        
        return None
    
    def get_gaia_photometry(self, source_ids):
        """
        Get GAIA spectroscopic photometry for source IDs.
        
        Parameters
        ----------
        source_ids : list
            List of GAIA source IDs
            
        Returns
        -------
        dict
            Dictionary mapping source_id to photometry data
        """
        if not source_ids:
            return {}
        
        try:
            # Generate synthetic photometry
            synthetic_phot = generate(source_ids, photometric_system=self.photometric_system)
            
            if synthetic_phot is None or len(synthetic_phot) == 0:
                return {}
            
            # Parse results
            photometry = {}
            for _, row in synthetic_phot.iterrows():
                source_id = int(row['source_id'])
                
                # Extract magnitudes and errors
                mags = {}
                errors = {}
                
                for band in ['U', 'B', 'V', 'R', 'I']:
                    mag_col = f'JkcStd_mag_{band}'
                    flux_col = f'JkcStd_flux_{band}'
                    flux_err_col = f'JkcStd_flux_error_{band}'
                    
                    # Get magnitude
                    if mag_col in row and not np.isnan(row[mag_col]):
                        mags[band] = float(row[mag_col])
                    else:
                        mags[band] = 99.999
                    
                    # Calculate magnitude error from flux and flux error
                    if (flux_col in row and flux_err_col in row and 
                        not np.isnan(row[flux_col]) and not np.isnan(row[flux_err_col]) and
                        row[flux_col] > 0 and row[flux_err_col] > 0):
                        snr = row[flux_col] / row[flux_err_col]
                        mag_error = 2.5 / np.log(10) / snr
                        errors[band] = float(mag_error)
                    else:
                        errors[band] = 0.999
                
                photometry[source_id] = {
                    'magnitudes': mags,
                    'errors': errors
                }
            
            return photometry
            
        except Exception as e:
            print(f"Error getting GAIA photometry: {e}")
            return {}
    
    def format_output_line(self, source, band, magnitude, error=None):
        """
        Format a single line for AstroImageJ output.
        
        Parameters
        ----------
        source : dict
            Source information
        band : str
            Photometric band (U, B, V, R, I)
        magnitude : float
            Magnitude value
        error : float, optional
            Magnitude error
            
        Returns
        -------
        str
            Formatted line for AstroImageJ
        """
        ra_str = source.get('ra_found_str', source['ra_orig_str'])
        dec_str = source.get('dec_found_str', source['dec_orig_str'])
        ref_star = source['ref_star'] if source['ref_star'] is not None else ''
        centroid = source['centroid'] if source['centroid'] is not None else ''
        
        # Format magnitude with appropriate precision
        if magnitude > 90:
            mag_str = "99.999"
        else:
            mag_str = f"{magnitude:.3f}"
        
        return f"{ra_str}, {dec_str}, {ref_star}, {centroid}, {mag_str}"
    
    def write_band_file(self, sources, band, output_dir, base_filename):
        """
        Write output file for a specific band.
        
        Parameters
        ----------
        sources : list
            List of source dictionaries
        band : str
            Photometric band (U, B, V, R, I)
        output_dir : str
            Output directory
        base_filename : str
            Base filename without extension
        """
        # Handle filename with existing band suffix
        # Replace existing _X suffix with _band, or add _band if none exists
        if base_filename.endswith('_V') or base_filename.endswith('_B') or base_filename.endswith('_R') or base_filename.endswith('_I') or base_filename.endswith('_U'):
            # Remove the last 2 characters (_X) and replace with new band
            output_filename = f"{base_filename[:-2]}_{band}.radec"
        else:
            # No existing band suffix, just add it
            output_filename = f"{base_filename}_{band}.radec"
        
        output_path = os.path.join(output_dir, output_filename)
        
        try:
            with open(output_path, 'w') as f:
                # Write header
                f.write("#RA in decimal or sexagesimal HOURS\n")
                f.write("#Dec in decimal or sexagesimal DEGREES\n")
                f.write("#Ref Star=0,1,missing (0=target star, 1=ref star, missing->first ap=target, others=ref)\n")
                f.write("#Centroid=0,1,missing (0=do not centroid, 1=centroid, missing=centroid)\n")
                f.write("#Apparent Magnitude or missing (value = apparent magnitude, or value > 99 or missing = no mag info)\n")
                f.write("#Add one comma separated line per aperture in the following format:\n")
                f.write("#RA, Dec, Ref Star, Centroid, Magnitude\n")
                
                # Write source data
                for source in sources:
                    if 'photometry' in source and band in source['photometry']['magnitudes']:
                        magnitude = source['photometry']['magnitudes'][band]
                    else:
                        magnitude = 99.999
                    
                    line = self.format_output_line(source, band, magnitude)
                    f.write(line + "\n")
            
            print(f"Created {band}-band file: {output_path}")
            
        except Exception as e:
            print(f"Error writing {band}-band file: {e}")
    
    def process_file(self, input_filepath, output_dir=None):
        """
        Process AstroImageJ .radec file and generate output files.
        
        Parameters
        ----------
        input_filepath : str
            Path to input .radec file
        output_dir : str, optional
            Output directory (default: ./output)
        """
        if output_dir is None:
            output_dir = "output"
        
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Parse input file
        print(f"Reading {input_filepath}...")
        sources = self.parse_radec_file(input_filepath)
        
        if not sources:
            print("No valid sources found in input file.")
            return
        
        print(f"Found {len(sources)} sources")
        
        # Find GAIA sources and collect source IDs
        print("Matching sources with GAIA catalog...")
        gaia_source_ids = []
        source_map = {}  # Maps GAIA source_id to source index
        
        for i, source in enumerate(sources):
            gaia_source = self.find_gaia_source(source['ra_deg'], source['dec_deg'])
            
            if gaia_source:
                source_id = gaia_source['source_id']
                gaia_source_ids.append(source_id)
                source_map[source_id] = i
                
                # Update source with GAIA coordinates
                sources[i]['gaia_source'] = gaia_source
                ra_found_str, dec_found_str = self.decimal_to_sexagesimal(
                    gaia_source['ra'], gaia_source['dec']
                )
                sources[i]['ra_found_str'] = ra_found_str
                sources[i]['dec_found_str'] = dec_found_str
                
                # Print coordinate comparison
                print(f"Coordinates {source['ra_orig_str']} {source['dec_orig_str']} -> "
                      f"{ra_found_str} {dec_found_str}")
            else:
                print(f"No GAIA match for {source['ra_orig_str']} {source['dec_orig_str']}")
        
        # Get photometry for all matched sources
        if gaia_source_ids:
            print(f"Getting photometry for {len(gaia_source_ids)} GAIA sources...")
            photometry_data = self.get_gaia_photometry(gaia_source_ids)
            
            # Add photometry to sources
            for source_id, phot in photometry_data.items():
                if source_id in source_map:
                    sources[source_map[source_id]]['photometry'] = phot
        
        # Generate output files for each band
        base_filename = Path(input_filepath).stem
        print(f"Generating output files...")
        
        for band in ['U', 'B', 'V', 'R', 'I']:
            self.write_band_file(sources, band, output_dir, base_filename)
        
        print("Processing complete!")


def main():
    """Main function."""
    if len(sys.argv) != 2:
        print("Usage: python astroimagej_aperture_file_generator.py input_file.radec")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    
    processor = AstroImageJProcessor()
    processor.process_file(input_file)


if __name__ == "__main__":
    main()