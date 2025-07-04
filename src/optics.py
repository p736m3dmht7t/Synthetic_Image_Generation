"""
Optics and Telescope Calculations

Handles telescope optics calculations including field of view, seeing,
and point spread function parameters.
"""

import numpy as np
import math


class TelescopeOptics:
    """Handles telescope optics and seeing calculations."""
    
    def __init__(self):
        """Initialize the telescope optics calculator."""
        pass
    
    def calculate_airy_disk_diameter(self, aperture_mm, wavelength_nm=550):
        """
        Calculate Airy disk diameter using the Airy disk formula.
        
        Args:
            aperture_mm (float): Telescope aperture in mm
            wavelength_nm (float): Wavelength in nanometers (default 550nm for V band)
        
        Returns:
            float: Airy disk diameter in arcseconds
        """
        try:
            # Convert wavelength to meters
            wavelength_m = wavelength_nm * 1e-9
            
            # Convert aperture to meters
            aperture_m = aperture_mm / 1000.0
            
            # Airy disk angular diameter (first dark ring)
            # θ = 2.44 * λ / D (in radians)
            airy_diameter_rad = 2.44 * wavelength_m / aperture_m
            
            # Convert to arcseconds
            airy_diameter_arcsec = airy_diameter_rad * 206265
            
            return airy_diameter_arcsec
            
        except Exception as e:
            print(f"Error calculating Airy disk diameter: {str(e)}")
            return None
    
    def calculate_dawes_limit(self, aperture_mm):
        """
        Calculate Dawes limit for double star resolution.
        
        Args:
            aperture_mm (float): Telescope aperture in mm
        
        Returns:
            float: Dawes limit in arcseconds
        """
        try:
            # Dawes limit formula: θ = 4.56 / D (where D is in mm)
            dawes_limit_arcsec = 4.56 / aperture_mm
            
            return dawes_limit_arcsec
            
        except Exception as e:
            print(f"Error calculating Dawes limit: {str(e)}")
            return None
    
    def calculate_theoretical_seeing(self, aperture_mm, wavelength_nm=550):
        """
        Calculate theoretical seeing limited by telescope optics.
        
        Args:
            aperture_mm (float): Telescope aperture in mm
            wavelength_nm (float): Wavelength in nanometers
        
        Returns:
            dict: Theoretical seeing parameters
        """
        try:
            # Calculate Airy disk diameter
            airy_diameter = self.calculate_airy_disk_diameter(aperture_mm, wavelength_nm)
            
            # Calculate Dawes limit
            dawes_limit = self.calculate_dawes_limit(aperture_mm)
            
            # Theoretical FWHM (approximately Airy disk diameter / 2)
            theoretical_fwhm = airy_diameter / 2.0
            
            seeing_data = {
                'airy_disk_diameter_arcsec': airy_diameter,
                'dawes_limit_arcsec': dawes_limit,
                'theoretical_fwhm_arcsec': theoretical_fwhm,
                'wavelength_nm': wavelength_nm
            }
            
            return seeing_data
            
        except Exception as e:
            print(f"Error calculating theoretical seeing: {str(e)}")
            return None
    
    def calculate_combined_seeing(self, aperture_mm, atmospheric_seeing_arcsec, wavelength_nm=550):
        """
        Calculate combined seeing from telescope optics and atmosphere.
        
        Args:
            aperture_mm (float): Telescope aperture in mm
            atmospheric_seeing_arcsec (float): Atmospheric seeing in arcseconds
            wavelength_nm (float): Wavelength in nanometers
        
        Returns:
            dict: Combined seeing parameters
        """
        try:
            # Get theoretical seeing
            theoretical = self.calculate_theoretical_seeing(aperture_mm, wavelength_nm)
            
            if theoretical is None:
                return None
            
            # Combined seeing using quadrature addition
            # FWHM_combined = sqrt(FWHM_telescope^2 + FWHM_atmosphere^2)
            combined_fwhm = math.sqrt(
                theoretical['theoretical_fwhm_arcsec']**2 + 
                atmospheric_seeing_arcsec**2
            )
            
            # For most cases, atmospheric seeing dominates
            limiting_factor = "atmospheric" if atmospheric_seeing_arcsec > theoretical['theoretical_fwhm_arcsec'] else "telescope"
            
            combined_data = {
                'telescope_fwhm_arcsec': theoretical['theoretical_fwhm_arcsec'],
                'atmospheric_fwhm_arcsec': atmospheric_seeing_arcsec,
                'combined_fwhm_arcsec': combined_fwhm,
                'limiting_factor': limiting_factor,
                'airy_disk_diameter_arcsec': theoretical['airy_disk_diameter_arcsec'],
                'dawes_limit_arcsec': theoretical['dawes_limit_arcsec']
            }
            
            return combined_data
            
        except Exception as e:
            print(f"Error calculating combined seeing: {str(e)}")
            return None
    
    def calculate_field_of_view(self, focal_length_mm, pixel_size_microns, width_pixels, height_pixels):
        """
        Calculate field of view for given telescope and camera parameters.
        
        Args:
            focal_length_mm (float): Telescope focal length in mm
            pixel_size_microns (float): Camera pixel size in microns
            width_pixels (int): Sensor width in pixels
            height_pixels (int): Sensor height in pixels
        
        Returns:
            dict: Field of view parameters
        """
        try:
            # Convert pixel size to mm
            pixel_size_mm = pixel_size_microns / 1000.0
            
            # Calculate sensor dimensions in mm
            sensor_width_mm = width_pixels * pixel_size_mm
            sensor_height_mm = height_pixels * pixel_size_mm
            
            # Calculate field of view in radians
            fov_width_rad = 2 * math.atan(sensor_width_mm / (2 * focal_length_mm))
            fov_height_rad = 2 * math.atan(sensor_height_mm / (2 * focal_length_mm))
            
            # Convert to degrees and arcminutes
            fov_width_deg = math.degrees(fov_width_rad)
            fov_height_deg = math.degrees(fov_height_rad)
            fov_width_arcmin = fov_width_deg * 60
            fov_height_arcmin = fov_height_deg * 60
            
            # Calculate diagonal FOV
            diagonal_arcmin = math.sqrt(fov_width_arcmin**2 + fov_height_arcmin**2)
            
            # Calculate plate scale (arcsec/pixel)
            plate_scale_arcsec_pixel = 206265 * pixel_size_microns / (focal_length_mm * 1000)
            
            # Calculate pixel scales
            pixel_scale_arcmin_pixel = plate_scale_arcsec_pixel / 60.0
            pixel_scale_deg_pixel = pixel_scale_arcmin_pixel / 60.0
            
            fov_data = {
                'width_deg': fov_width_deg,
                'height_deg': fov_height_deg,
                'width_arcmin': fov_width_arcmin,
                'height_arcmin': fov_height_arcmin,
                'diagonal_arcmin': diagonal_arcmin,
                'plate_scale_arcsec_pixel': plate_scale_arcsec_pixel,
                'pixel_scale_arcmin_pixel': pixel_scale_arcmin_pixel,
                'pixel_scale_deg_pixel': pixel_scale_deg_pixel
            }
            
            return fov_data
            
        except Exception as e:
            print(f"Error calculating field of view: {str(e)}")
            return None
    
    def calculate_psf_parameters(self, telescope_config, camera_config):
        """
        Calculate point spread function parameters for image generation.
        
        Args:
            telescope_config (dict): Telescope configuration
            camera_config (dict): Camera configuration
        
        Returns:
            dict: PSF parameters
        """
        try:
            # Extract telescope parameters
            aperture_mm = telescope_config.get('aperture_mm', 0)
            focal_length_mm = telescope_config.get('focal_length_mm', 0)
            atmospheric_seeing = telescope_config.get('seeing', {}).get('atmospheric_arcsec', 1.5)
            
            # Extract camera parameters
            pixel_size_microns = camera_config.get('pixel_size_microns', 0)
            
            # Calculate combined seeing
            combined_seeing = self.calculate_combined_seeing(aperture_mm, atmospheric_seeing)
            
            if combined_seeing is None:
                return None
            
            # Calculate field of view
            sensor_dims = camera_config.get('sensor_dimensions', {})
            fov = self.calculate_field_of_view(
                focal_length_mm,
                pixel_size_microns,
                sensor_dims.get('width_pixels', 0),
                sensor_dims.get('height_pixels', 0)
            )
            
            if fov is None:
                return None
            
            # Calculate PSF size in pixels
            fwhm_arcsec = combined_seeing['combined_fwhm_arcsec']
            fwhm_pixels = fwhm_arcsec / fov['plate_scale_arcsec_pixel']
            
            # Calculate Gaussian sigma (FWHM = 2.355 * sigma)
            sigma_pixels = fwhm_pixels / 2.355
            
            psf_data = {
                'fwhm_arcsec': fwhm_arcsec,
                'fwhm_pixels': fwhm_pixels,
                'sigma_pixels': sigma_pixels,
                'telescope_limited': combined_seeing['limiting_factor'] == 'telescope',
                'atmospheric_limited': combined_seeing['limiting_factor'] == 'atmospheric',
                'combined_seeing': combined_seeing,
                'field_of_view': fov
            }
            
            return psf_data
            
        except Exception as e:
            print(f"Error calculating PSF parameters: {str(e)}")
            return None
    
    def update_telescope_config_with_calculations(self, telescope_config):
        """
        Update telescope configuration with calculated seeing values.
        
        Args:
            telescope_config (dict): Telescope configuration
        
        Returns:
            dict: Updated telescope configuration
        """
        try:
            # Extract parameters
            aperture_mm = telescope_config.get('aperture_mm', 0)
            atmospheric_seeing = telescope_config.get('seeing', {}).get('atmospheric_arcsec', 1.5)
            
            # Calculate theoretical seeing
            theoretical = self.calculate_theoretical_seeing(aperture_mm)
            
            if theoretical is None:
                return telescope_config
            
            # Calculate combined seeing
            combined = self.calculate_combined_seeing(aperture_mm, atmospheric_seeing)
            
            if combined is None:
                return telescope_config
            
            # Update configuration
            updated_config = telescope_config.copy()
            
            # Update seeing section
            if 'seeing' not in updated_config:
                updated_config['seeing'] = {}
            
            updated_config['seeing']['theoretical_arcsec'] = theoretical['theoretical_fwhm_arcsec']
            updated_config['seeing']['combined_fwhm_arcsec'] = combined['combined_fwhm_arcsec']
            updated_config['seeing']['airy_disk_diameter_arcsec'] = theoretical['airy_disk_diameter_arcsec']
            updated_config['seeing']['dawes_limit_arcsec'] = theoretical['dawes_limit_arcsec']
            updated_config['seeing']['limiting_factor'] = combined['limiting_factor']
            
            return updated_config
            
        except Exception as e:
            print(f"Error updating telescope config: {str(e)}")
            return telescope_config