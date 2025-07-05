"""
Image Generation Engine

Main engine for generating synthetic multi-band astronomical images.
Coordinates all components to create realistic stellar field images.
"""

import numpy as np
import os
from pathlib import Path
import json

from catalog_query import CatalogQuery
from coordinates import CoordinateTransform
from optics import TelescopeOptics
try:
    from psf_photutils import PhotoutilsPSFGenerator
except ImportError:
    from .psf_photutils import PhotoutilsPSFGenerator


class ImageGenerator:
    """Main image generation engine for synthetic astronomical images."""
    
    def __init__(self, config_manager):
        """
        Initialize the image generator.
        
        Args:
            config_manager: ConfigManager instance
        """
        self.config_manager = config_manager
        self.catalog_query = CatalogQuery()
        self.coordinate_transform = CoordinateTransform()
        self.telescope_optics = TelescopeOptics()
        self.psf_generator = PhotoutilsPSFGenerator()
        
        # Generation parameters
        self.target_config = None
        self.telescope_config = None
        self.camera_config = None
        self.fits_header_config = None
        
        # Generated data
        self.source_data = None
        self.psf_params = None
        self.images = {}
        
    def load_configurations(self):
        """Load all configuration files."""
        try:
            self.target_config = self.config_manager.load_target_config()
            self.telescope_config = self.config_manager.load_telescope_config()
            self.camera_config = self.config_manager.load_camera_config()
            self.fits_header_config = self.config_manager.load_fits_header_config()
            
            print("Configurations loaded successfully")
            return True
            
        except Exception as e:
            print(f"Error loading configurations: {str(e)}")
            return False
    
    def validate_configurations(self):
        """Validate that all required configurations are present and valid."""
        try:
            issues = []
            
            # Check target configuration
            if not self.target_config or not self.target_config.get('target_name'):
                issues.append("Target name is required")
            
            target_ra = self.target_config.get('coordinates', {}).get('ra', 0)
            target_dec = self.target_config.get('coordinates', {}).get('dec', 0)
            if target_ra == 0 and target_dec == 0:
                issues.append("Target coordinates are required")
            
            target_v_mag = self.target_config.get('magnitudes', {}).get('V')
            if target_v_mag is None or target_v_mag == 0:
                issues.append("Target V magnitude is required")
            
            # Check telescope configuration
            if not self.telescope_config or not self.telescope_config.get('telescope_name'):
                issues.append("Telescope name is required")
            
            focal_length = self.telescope_config.get('focal_length_mm', 0)
            aperture = self.telescope_config.get('aperture_mm', 0)
            if focal_length == 0 or aperture == 0:
                issues.append("Telescope focal length and aperture are required")
            
            # Check camera configuration
            if not self.camera_config or not self.camera_config.get('camera_name'):
                issues.append("Camera name is required")
            
            pixel_size = self.camera_config.get('pixel_size_microns', 0)
            sensor_dims = self.camera_config.get('sensor_dimensions', {})
            width_pixels = sensor_dims.get('width_pixels', 0)
            height_pixels = sensor_dims.get('height_pixels', 0)
            
            if pixel_size == 0 or width_pixels == 0 or height_pixels == 0:
                issues.append("Camera pixel size and sensor dimensions are required")
            
            if issues:
                print("Configuration validation failed:")
                for issue in issues:
                    print(f"  - {issue}")
                return False
            
            print("Configuration validation passed")
            return True
            
        except Exception as e:
            print(f"Error validating configurations: {str(e)}")
            return False
    
    def setup_coordinate_system(self):
        """Set up the coordinate transformation system."""
        try:
            # Get target coordinates
            target_coords = self.target_config.get('coordinates', {})
            target_ra = target_coords.get('ra', 0)
            target_dec = target_coords.get('dec', 0)
            
            # Get camera parameters
            rotation_deg = self.camera_config.get('rotation_degrees', 0)
            pixel_size_microns = self.camera_config.get('pixel_size_microns', 0)
            
            # Get sensor dimensions
            sensor_dims = self.camera_config.get('sensor_dimensions', {})
            width_pixels = sensor_dims.get('width_pixels', 0)
            height_pixels = sensor_dims.get('height_pixels', 0)
            
            # Get telescope focal length
            focal_length_mm = self.telescope_config.get('focal_length_mm', 0)
            
            # Calculate plate scale
            plate_scale_arcsec_pixel = 206265 * pixel_size_microns / (focal_length_mm * 1000)
            
            # Set up coordinate transformation
            self.coordinate_transform.setup_transform(
                target_ra, target_dec, rotation_deg,
                plate_scale_arcsec_pixel, width_pixels, height_pixels
            )
            
            print("Coordinate system initialized")
            print(f"Target: RA={target_ra:.4f}°, Dec={target_dec:.4f}°")
            print(f"Plate scale: {plate_scale_arcsec_pixel:.2f} arcsec/pixel")
            
            return True
            
        except Exception as e:
            print(f"Error setting up coordinate system: {str(e)}")
            return False
    
    def calculate_psf_parameters(self):
        """Calculate PSF parameters for the current setup."""
        try:
            self.psf_params = self.telescope_optics.calculate_psf_parameters(
                self.telescope_config, self.camera_config
            )
            
            if self.psf_params is None:
                print("Error: Could not calculate PSF parameters")
                return False
            
            print("PSF parameters calculated")
            print(f"FWHM: {self.psf_params['fwhm_arcsec']:.2f} arcsec ({self.psf_params['fwhm_pixels']:.2f} pixels)")
            print(f"Sigma: {self.psf_params['sigma_pixels']:.2f} pixels")
            
            return True
            
        except Exception as e:
            print(f"Error calculating PSF parameters: {str(e)}")
            return False
    
    def query_catalog_sources(self):
        """Query catalog for sources in the field of view."""
        try:
            # Get target parameters
            target_coords = self.target_config.get('coordinates', {})
            target_ra = target_coords.get('ra', 0)
            target_dec = target_coords.get('dec', 0)
            target_v_mag = self.target_config.get('magnitudes', {}).get('V', 0)
            
            # Query catalog
            print("Querying GAIA catalog...")
            self.source_data = self.catalog_query.get_sources_in_field(
                target_ra, target_dec, target_v_mag,
                self.telescope_config, self.camera_config
            )
            
            if self.source_data is None:
                print("Error: Could not query catalog")
                return False
            
            print(f"Retrieved {len(self.source_data.get('ra', []))} catalog sources")
            
            # Filter sources to field of view and add pixel coordinates
            print("Converting coordinates to pixels...")
            filtered_sources = self.coordinate_transform.filter_sources_in_field(self.source_data)
            
            if filtered_sources is None:
                print("Error: Could not filter sources to field")
                return False
            
            self.source_data = filtered_sources
            print(f"Filtered to {len(self.source_data.get('ra', []))} sources in field")
            
            return True
            
        except Exception as e:
            print(f"Error querying catalog sources: {str(e)}")
            return False
    
    def generate_band_image(self, band):
        """
        Generate synthetic image for a specific photometric band.
        
        Args:
            band (str): Photometric band ('V', 'B', 'R', 'I')
        
        Returns:
            numpy.ndarray: Generated image
        """
        try:
            # Get image dimensions
            sensor_dims = self.camera_config.get('sensor_dimensions', {})
            width_pixels = sensor_dims.get('width_pixels', 0)
            height_pixels = sensor_dims.get('height_pixels', 0)
            image_shape = (height_pixels, width_pixels)
            
            # Get target magnitude for this band
            target_magnitude = self.target_config.get('magnitudes', {}).get(band)
            if target_magnitude is None or target_magnitude == 0:
                print(f"Skipping {band} band: No target magnitude specified")
                return None
            
            # Get saturation limit
            saturation_limit = self.camera_config.get('saturation_limit', 65535)
            
            print(f"Generating {band} band image...")
            print(f"Target magnitude: {target_magnitude:.2f}")
            print(f"Image size: {width_pixels} × {height_pixels} pixels")
            
            # Create stellar field
            image = self.psf_generator.create_stellar_field(
                self.source_data, image_shape, self.psf_params,
                target_magnitude, saturation_limit, band
            )
            
            if image is None:
                print(f"Error: Could not generate {band} band image")
                return None
            
            # Add noise (optional)
            noise_params = {'read_noise_sigma': 1.0}
            noisy_image = self.psf_generator.add_noise(image, 'both', noise_params)
            
            # Convert to integer type for FITS output
            final_image = np.clip(noisy_image, 0, saturation_limit).astype(np.uint16)
            
            print(f"{band} band image generated successfully")
            
            return final_image
            
        except Exception as e:
            print(f"Error generating {band} band image: {str(e)}")
            return None
    
    def generate_all_band_images(self):
        """Generate images for all photometric bands (V, B, R, I) that have magnitudes."""
        try:
            bands = ['V', 'B', 'R', 'I']
            self.images = {}
            skipped_bands = []
            
            for band in bands:
                print(f"\n--- Generating {band} Band ---")
                image = self.generate_band_image(band)
                
                if image is not None:
                    self.images[band] = image
                    print(f"{band} band: {image.shape}, dtype={image.dtype}")
                else:
                    skipped_bands.append(band)
                    print(f"Skipped {band} band (no magnitude available)")
            
            if len(self.images) == 0:
                print("Error: No band images could be generated (no magnitudes available)")
                return False
            
            print(f"\nSuccessfully generated {len(self.images)} band images")
            if skipped_bands:
                print(f"Skipped bands due to missing magnitudes: {', '.join(skipped_bands)}")
            return True
            
        except Exception as e:
            print(f"Error generating band images: {str(e)}")
            return False
    
    def get_generation_summary(self):
        """Get a summary of the image generation process."""
        try:
            if not self.source_data or not self.images:
                return None
            
            summary = {
                'target_name': self.target_config.get('target_name', 'Unknown'),
                'target_coordinates': {
                    'ra': self.target_config.get('coordinates', {}).get('ra', 0),
                    'dec': self.target_config.get('coordinates', {}).get('dec', 0)
                },
                'telescope': self.telescope_config.get('telescope_name', 'Unknown'),
                'camera': self.camera_config.get('camera_name', 'Unknown'),
                'num_sources': len(self.source_data.get('ra', [])),
                'bands_generated': list(self.images.keys()),
                'image_dimensions': {
                    'width': self.camera_config.get('sensor_dimensions', {}).get('width_pixels', 0),
                    'height': self.camera_config.get('sensor_dimensions', {}).get('height_pixels', 0)
                },
                'psf_fwhm_arcsec': self.psf_params.get('fwhm_arcsec', 0) if self.psf_params else 0,
                'field_of_view': self.psf_params.get('field_of_view', {}) if self.psf_params else {}
            }
            
            return summary
            
        except Exception as e:
            print(f"Error creating generation summary: {str(e)}")
            return None
    
    def generate_catalog_documentation(self, fits_filename):
        """
        Generate a human-readable catalog file documenting sources in the image.
        
        Args:
            fits_filename (str): Name of corresponding FITS file (used to derive base name)
            
        Returns:
            str: Path to generated catalog file
        """
        try:
            # Create catalog filename by removing band suffix and adding .catalog
            # Example: "V*_TZ_Boo_20250105_120000_V.fits" -> "V*_TZ_Boo_20250105_120000.catalog"
            from pathlib import Path
            output_dir = Path("output")
            output_dir.mkdir(exist_ok=True)
            
            # Remove the band suffix (last _X before .fits) to create a unified catalog name
            base_name = fits_filename.replace('.fits', '')
            if base_name.endswith('_V') or base_name.endswith('_B') or base_name.endswith('_R') or base_name.endswith('_I') or base_name.endswith('_U'):
                base_name = base_name[:-2]  # Remove last 2 characters (_X)
            
            catalog_filename = str(output_dir / f"{base_name}.catalog")
            
            # Get source data
            if not self.source_data:
                print(f"Warning: No source data available for catalog documentation")
                return None
            
            # Prepare data arrays
            source_ids = self.source_data.get('source_id', [])
            ras = self.source_data.get('ra', [])
            decs = self.source_data.get('dec', [])
            
            # Get method indicators for all bands
            spectroscopic_bands = self.source_data.get('spectroscopic_bands', {})
            
            if len(source_ids) == 0:
                print(f"Warning: No sources to document")
                return None
            
            print(f"Generating catalog documentation: {catalog_filename}")
            
            with open(catalog_filename, 'w') as f:
                # Write header
                f.write(f"# Source Catalog for {base_name}\n")
                f.write(f"# Bands: U, B, V, R, I (Johnson-Cousins)\n")
                f.write(f"# Target: {self.target_config.get('target_name', 'Unknown')}\n")
                f.write(f"# Generated: {self.get_timestamp()}\n")
                f.write(f"# Total Sources: {len(source_ids)}\n")
                f.write("#\n")
                f.write("# Magnitudes marked with '*' were derived using polynomial method\n")
                f.write("# Magnitudes without '*' were derived using GAIA spectroscopic method\n")
                f.write("#\n")
                f.write("# Column Headers:\n")
                f.write("#   Source_ID: GAIA DR3 source identifier\n")
                f.write("#   RA: Right Ascension (HH:MM:SS.SSS)\n")
                f.write("#   Dec: Declination (+DD:MM:SS.SS)\n")
                f.write("#   U,B,V,R,I: Johnson-Cousins magnitudes and errors\n")
                f.write("#\n")
                
                # Write column headers with proper spacing
                f.write("Source_ID                    RA              Dec               ")
                f.write("U                    B                    V                    ")
                f.write("R                    I\n")
                f.write("-" * 162 + "\n")
                
                # Process each source
                for i in range(len(source_ids)):
                    # Format source ID as string with "GAIA " prefix
                    source_id_str = f'"GAIA {source_ids[i]}"'
                    
                    # Convert coordinates to sexagesimal
                    ra_str, dec_str = self.format_coordinates(ras[i], decs[i])
                    
                    # Get all band magnitudes and errors for this source
                    all_bands = ['U', 'B', 'V', 'R', 'I']
                    mag_strings = []
                    
                    for band_name in all_bands:
                        band_mag_list = self.source_data.get(band_name, [])
                        band_err_list = self.source_data.get(f'{band_name}_error', [])
                        
                        if i < len(band_mag_list) and band_mag_list[i] is not None:
                            mag = band_mag_list[i]
                            err = band_err_list[i] if i < len(band_err_list) and band_err_list[i] is not None else 0.999
                            
                            # Check if this magnitude was derived spectroscopically
                            is_spectroscopic = (i in spectroscopic_bands and 
                                              band_name in spectroscopic_bands[i])
                            
                            # Add asterisk for polynomial method (always use 1 character for alignment)
                            asterisk = "*" if not is_spectroscopic else " "
                            
                            # Format magnitude with proper alignment
                            if np.isfinite(mag) and mag < 90:
                                mag_str = f"{mag:6.3f}{asterisk} +/- {err:5.3f}"
                            else:
                                mag_str = "     --  +/-    --"
                        else:
                            mag_str = "     --  +/-    --"
                        
                        mag_strings.append(mag_str)
                    
                    # Write formatted line
                    line = (f"{source_id_str:<25} {ra_str:<15} {dec_str:<17} "
                           f"{mag_strings[0]:<20} {mag_strings[1]:<20} {mag_strings[2]:<20} "
                           f"{mag_strings[3]:<20} {mag_strings[4]:<20}\n")
                    f.write(line)
                
                # Write footer
                f.write("-" * 162 + "\n")
                f.write(f"# End of catalog ({len(source_ids)} sources)\n")
            
            print(f"Catalog documentation saved: {catalog_filename}")
            return catalog_filename
            
        except Exception as e:
            print(f"Error generating catalog documentation: {e}")
            return None
    
    def format_coordinates(self, ra_deg, dec_deg):
        """
        Format coordinates to sexagesimal strings.
        
        Args:
            ra_deg (float): RA in decimal degrees
            dec_deg (float): Dec in decimal degrees
            
        Returns:
            tuple: (ra_str, dec_str) in HH:MM:SS.SSS, +DD:MM:SS.SS format
        """
        try:
            # Convert RA from degrees to hours
            ra_hours = ra_deg / 15.0
            
            # RA formatting
            ra_h = int(ra_hours)
            ra_m = int((ra_hours - ra_h) * 60)
            ra_s = ((ra_hours - ra_h) * 60 - ra_m) * 60
            ra_str = f"{ra_h:02d}:{ra_m:02d}:{ra_s:06.3f}"
            
            # Dec formatting
            dec_sign = "+" if dec_deg >= 0 else "-"
            dec_abs = abs(dec_deg)
            dec_d = int(dec_abs)
            dec_m = int((dec_abs - dec_d) * 60)
            dec_s = ((dec_abs - dec_d) * 60 - dec_m) * 60
            dec_str = f"{dec_sign}{dec_d:02d}:{dec_m:02d}:{dec_s:05.2f}"
            
            return ra_str, dec_str
            
        except Exception as e:
            print(f"Error formatting coordinates: {e}")
            return "00:00:00.000", "+00:00:00.00"
    
    def get_timestamp(self):
        """Get current timestamp string."""
        from datetime import datetime
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    def generate_images(self):
        """
        Main method to generate synthetic images.
        
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            print("=== Starting Image Generation ===")
            
            # Load and validate configurations
            if not self.load_configurations():
                return False
            
            if not self.validate_configurations():
                return False
            
            # Set up coordinate system
            if not self.setup_coordinate_system():
                return False
            
            # Calculate PSF parameters
            if not self.calculate_psf_parameters():
                return False
            
            # Query catalog sources
            if not self.query_catalog_sources():
                return False
            
            # Generate all band images
            if not self.generate_all_band_images():
                return False
            
            print("\n=== Image Generation Complete ===")
            
            # Print summary
            summary = self.get_generation_summary()
            if summary:
                print(f"Target: {summary['target_name']}")
                print(f"Sources: {summary['num_sources']}")
                print(f"Bands: {', '.join(summary['bands_generated'])}")
                print(f"Image size: {summary['image_dimensions']['width']} × {summary['image_dimensions']['height']}")
                print(f"PSF FWHM: {summary['psf_fwhm_arcsec']:.2f} arcsec")
            
            return True
            
        except Exception as e:
            print(f"Error in image generation: {str(e)}")
            return False