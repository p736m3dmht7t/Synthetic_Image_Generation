"""
FITS File Writer

Handles writing synthetic images to FITS files with proper headers
and metadata for astronomical use.
"""

import numpy as np
import json
from datetime import datetime
from pathlib import Path
from astropy.io import fits
from astropy.time import Time
import warnings

# Suppress FITS warnings
warnings.filterwarnings('ignore', category=fits.verify.VerifyWarning)


class FITSWriter:
    """Writes synthetic images to FITS files with proper headers."""
    
    def __init__(self, config_manager):
        """
        Initialize the FITS writer.
        
        Args:
            config_manager: ConfigManager instance
        """
        self.config_manager = config_manager
        self.base_header = None
        self.output_dir = Path("output")
        self.output_dir.mkdir(exist_ok=True)
    
    def load_header_template(self):
        """Load FITS header template from configuration."""
        try:
            self.base_header = self.config_manager.load_fits_header_config()
            
            if not self.base_header:
                print("Warning: No FITS header template found, using minimal header")
                self.base_header = self._create_minimal_header()
            
            print("FITS header template loaded")
            return True
            
        except Exception as e:
            print(f"Error loading FITS header template: {str(e)}")
            self.base_header = self._create_minimal_header()
            return False
    
    def _create_minimal_header(self):
        """Create minimal FITS header if no template available."""
        return {
            'SIMPLE': True,
            'BITPIX': 16,
            'NAXIS': 2,
            'EXTEND': True,
            'IMAGETYP': 'LIGHT',
            'EQUINOX': 2000.0,
            'CTYPE1': 'RA---TAN',
            'CTYPE2': 'DEC--TAN',
            'CUNIT1': 'deg',
            'CUNIT2': 'deg'
        }
    
    def update_header_for_target(self, header_dict, target_config, telescope_config, 
                                camera_config, band, psf_params=None):
        """
        Update header with target-specific information.
        
        Args:
            header_dict (dict): Base header dictionary
            target_config (dict): Target configuration
            telescope_config (dict): Telescope configuration
            camera_config (dict): Camera configuration
            band (str): Photometric band
            psf_params (dict): PSF parameters (optional)
        
        Returns:
            dict: Updated header dictionary
        """
        try:
            # Create a copy of the base header
            updated_header = header_dict.copy()
            
            # Update basic image parameters
            sensor_dims = camera_config.get('sensor_dimensions', {})
            updated_header['NAXIS1'] = sensor_dims.get('width_pixels', 0)
            updated_header['NAXIS2'] = sensor_dims.get('height_pixels', 0)
            
            # Target information
            target_coords = target_config.get('coordinates', {})
            updated_header['OBJECT'] = target_config.get('target_name', 'SYNTHETIC')
            updated_header['RA'] = target_coords.get('ra', 0.0)
            updated_header['DEC'] = target_coords.get('dec', 0.0)
            
            # Convert coordinates to string format if available
            ra_str = target_coords.get('ra_str', '')
            dec_str = target_coords.get('dec_str', '')
            if ra_str:
                updated_header['OBJCTRA'] = ra_str
            if dec_str:
                updated_header['OBJCTDEC'] = dec_str
            
            # Telescope information
            updated_header['TELESCOP'] = telescope_config.get('telescope_name', 'SYNTHETIC')
            updated_header['FOCALLEN'] = telescope_config.get('focal_length_mm', 0.0)
            updated_header['FOCRATIO'] = telescope_config.get('focal_ratio', 0.0)
            
            # Camera information
            updated_header['INSTRUME'] = camera_config.get('camera_name', 'SYNTHETIC')
            updated_header['XPIXSZ'] = camera_config.get('pixel_size_microns', 0.0)
            updated_header['YPIXSZ'] = camera_config.get('pixel_size_microns', 0.0)
            
            # Filter/band information
            updated_header['FILTER'] = band
            
            # Observation time (current time for synthetic data)
            now = Time.now()
            updated_header['DATE-OBS'] = now.iso
            updated_header['DATE-AVG'] = now.iso
            
            # Exposure time (use from template or default)
            if 'EXPTIME' not in updated_header:
                updated_header['EXPTIME'] = 100.0
            if 'EXPOSURE' not in updated_header:
                updated_header['EXPOSURE'] = updated_header['EXPTIME']
            
            # WCS information
            if psf_params and 'field_of_view' in psf_params:
                fov = psf_params['field_of_view']
                plate_scale = fov.get('plate_scale_arcsec_pixel', 0)
                
                # Reference pixel (center of image)
                updated_header['CRPIX1'] = sensor_dims.get('width_pixels', 0) / 2.0 + 0.5
                updated_header['CRPIX2'] = sensor_dims.get('height_pixels', 0) / 2.0 + 0.5
                
                # Reference coordinates (target position)
                updated_header['CRVAL1'] = target_coords.get('ra', 0.0)
                updated_header['CRVAL2'] = target_coords.get('dec', 0.0)
                
                # Pixel scale in degrees/pixel
                scale_deg_pixel = plate_scale / 3600.0
                
                # Simple CD matrix (no rotation for now)
                rotation_deg = camera_config.get('rotation_degrees', 0.0)
                cos_rot = np.cos(np.radians(rotation_deg))
                sin_rot = np.sin(np.radians(rotation_deg))
                
                updated_header['CD1_1'] = -scale_deg_pixel * cos_rot
                updated_header['CD1_2'] = scale_deg_pixel * sin_rot
                updated_header['CD2_1'] = scale_deg_pixel * sin_rot
                updated_header['CD2_2'] = scale_deg_pixel * cos_rot
                
                # Legacy keywords for compatibility
                updated_header['CDELT1'] = scale_deg_pixel
                updated_header['CDELT2'] = scale_deg_pixel
                updated_header['CROTA2'] = rotation_deg
            
            # Add synthetic data identifier
            updated_header['ORIGIN'] = 'Synthetic Image Generation'
            updated_header['SWCREATE'] = 'Python Synthetic Imaging'
            
            # Add generation timestamp
            updated_header['COMMENT'] = f'Synthetic {band} band image generated on {datetime.now().isoformat()}'
            
            # Target magnitude for this band
            target_mags = target_config.get('magnitudes', {})
            if band in target_mags:
                updated_header[f'TARG_{band}'] = target_mags[band]
            
            # PSF information
            if psf_params:
                updated_header['SEEING'] = psf_params.get('fwhm_arcsec', 0.0)
                updated_header['FWHM'] = psf_params.get('fwhm_pixels', 0.0)
            
            return updated_header
            
        except Exception as e:
            print(f"Error updating header for target: {str(e)}")
            return header_dict
    
    def create_fits_header(self, header_dict):
        """
        Convert header dictionary to FITS header object.
        
        Args:
            header_dict (dict): Header key-value pairs
        
        Returns:
            astropy.io.fits.Header: FITS header object
        """
        try:
            header = fits.Header()
            
            # Add standard keywords first
            standard_order = ['SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'EXTEND']
            
            for key in standard_order:
                if key in header_dict:
                    header[key] = header_dict[key]
            
            # Add remaining keywords
            for key, value in header_dict.items():
                if key not in standard_order:
                    if key == 'COMMENT':
                        # Handle multiple comments
                        if isinstance(value, list):
                            for comment in value:
                                header['COMMENT'] = str(comment)
                        else:
                            header['COMMENT'] = str(value)
                    elif key == 'HISTORY':
                        # Handle multiple history entries
                        if isinstance(value, list):
                            for history in value:
                                header['HISTORY'] = str(history)
                        else:
                            header['HISTORY'] = str(value)
                    else:
                        # Regular keyword
                        try:
                            header[key] = value
                        except Exception as e:
                            print(f"Warning: Could not add header keyword {key}: {e}")
            
            return header
            
        except Exception as e:
            print(f"Error creating FITS header: {str(e)}")
            return fits.Header()
    
    def write_fits_file(self, image_data, filename, header_dict):
        """
        Write image data to FITS file.
        
        Args:
            image_data (numpy.ndarray): Image data
            filename (str): Output filename
            header_dict (dict): Header information
        
        Returns:
            bool: True if successful
        """
        try:
            # Create FITS header
            fits_header = self.create_fits_header(header_dict)
            
            # Ensure image data is correct type
            if image_data.dtype != np.uint16:
                image_data = np.clip(image_data, 0, 65535).astype(np.uint16)
            
            # Create HDU
            hdu = fits.PrimaryHDU(data=image_data, header=fits_header)
            
            # Write file
            output_path = self.output_dir / filename
            hdu.writeto(output_path, overwrite=True)
            
            print(f"FITS file written: {output_path}")
            return True
            
        except Exception as e:
            print(f"Error writing FITS file {filename}: {str(e)}")
            return False
    
    def write_band_images(self, images, target_config, telescope_config, 
                         camera_config, psf_params=None):
        """
        Write all band images to FITS files.
        
        Args:
            images (dict): Dictionary of band -> image data
            target_config (dict): Target configuration
            telescope_config (dict): Telescope configuration
            camera_config (dict): Camera configuration
            psf_params (dict): PSF parameters (optional)
        
        Returns:
            dict: Dictionary of band -> output filename
        """
        try:
            # Load header template
            if not self.load_header_template():
                print("Warning: Using minimal header template")
            
            # Generate base filename
            target_name = target_config.get('target_name', 'SYNTHETIC').replace(' ', '_')
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            
            output_files = {}
            
            for band, image_data in images.items():
                # Create filename
                filename = f"{target_name}_{timestamp}_{band}.fits"
                
                # Update header for this band and target
                updated_header = self.update_header_for_target(
                    self.base_header, target_config, telescope_config,
                    camera_config, band, psf_params
                )
                
                # Write FITS file
                if self.write_fits_file(image_data, filename, updated_header):
                    output_files[band] = filename
                else:
                    print(f"Failed to write {band} band FITS file")
            
            if output_files:
                print(f"Successfully wrote {len(output_files)} FITS files")
                for band, filename in output_files.items():
                    print(f"  {band} band: {filename}")
            
            return output_files
            
        except Exception as e:
            print(f"Error writing band images: {str(e)}")
            return {}
    
    def get_output_directory(self):
        """Get the output directory path."""
        return self.output_dir
    
    def set_output_directory(self, directory):
        """
        Set the output directory.
        
        Args:
            directory (str or Path): Output directory path
        """
        self.output_dir = Path(directory)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Output directory set to: {self.output_dir}")
    
    def validate_fits_file(self, filename):
        """
        Validate a FITS file by attempting to read it.
        
        Args:
            filename (str): FITS filename
        
        Returns:
            dict: Validation results
        """
        try:
            filepath = self.output_dir / filename
            
            if not filepath.exists():
                return {'valid': False, 'error': 'File not found'}
            
            # Try to open and read the file
            with fits.open(filepath) as hdul:
                hdu = hdul[0]
                
                validation = {
                    'valid': True,
                    'filename': filename,
                    'dimensions': hdu.data.shape if hdu.data is not None else None,
                    'dtype': str(hdu.data.dtype) if hdu.data is not None else None,
                    'header_keywords': len(hdu.header),
                    'object': hdu.header.get('OBJECT', 'Unknown'),
                    'filter': hdu.header.get('FILTER', 'Unknown'),
                    'exptime': hdu.header.get('EXPTIME', 0),
                    'file_size_mb': filepath.stat().st_size / (1024 * 1024)
                }
                
                return validation
                
        except Exception as e:
            return {'valid': False, 'error': str(e)}
    
    def get_header_summary(self, header_dict):
        """
        Get a summary of header contents.
        
        Args:
            header_dict (dict): Header dictionary
        
        Returns:
            dict: Header summary
        """
        try:
            summary = {
                'total_keywords': len(header_dict),
                'object': header_dict.get('OBJECT', 'Unknown'),
                'ra': header_dict.get('RA', 0),
                'dec': header_dict.get('DEC', 0),
                'filter': header_dict.get('FILTER', 'Unknown'),
                'telescope': header_dict.get('TELESCOP', 'Unknown'),
                'instrument': header_dict.get('INSTRUME', 'Unknown'),
                'exposure_time': header_dict.get('EXPTIME', 0),
                'image_type': header_dict.get('IMAGETYP', 'Unknown')
            }
            
            return summary
            
        except Exception as e:
            print(f"Error creating header summary: {str(e)}")
            return {}