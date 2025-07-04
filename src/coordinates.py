"""
Coordinate Transformations

Handles coordinate transformations between RA/Dec and pixel coordinates,
including camera rotation and field of view calculations.
"""

import numpy as np
import math
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
import warnings

# Suppress warnings
warnings.filterwarnings('ignore')


class CoordinateTransform:
    """Handles coordinate transformations for synthetic image generation."""
    
    def __init__(self):
        """Initialize the coordinate transformation system."""
        self.wcs = None
        self.image_center_ra = None
        self.image_center_dec = None
        self.rotation_angle = 0.0
        self.plate_scale = None
        self.image_dimensions = None
    
    def setup_transform(self, center_ra_deg, center_dec_deg, rotation_deg, 
                       plate_scale_arcsec_pixel, width_pixels, height_pixels):
        """
        Set up coordinate transformation parameters.
        
        Args:
            center_ra_deg (float): Image center RA in degrees
            center_dec_deg (float): Image center Dec in degrees
            rotation_deg (float): Camera rotation in degrees (Y-axis relative to north)
            plate_scale_arcsec_pixel (float): Plate scale in arcsec/pixel
            width_pixels (int): Image width in pixels
            height_pixels (int): Image height in pixels
        """
        self.image_center_ra = center_ra_deg
        self.image_center_dec = center_dec_deg
        self.rotation_angle = rotation_deg
        self.plate_scale = plate_scale_arcsec_pixel
        self.image_dimensions = (width_pixels, height_pixels)
        
        # Create simple WCS for transformations
        self._create_wcs()
    
    def _create_wcs(self):
        """Create a WCS object for coordinate transformations."""
        try:
            self.wcs = WCS(naxis=2)
            
            # Set up basic WCS parameters
            width_pixels, height_pixels = self.image_dimensions
            
            # Reference pixel (center of image)
            self.wcs.wcs.crpix = [width_pixels/2 + 0.5, height_pixels/2 + 0.5]
            
            # Reference coordinates (center of field)
            self.wcs.wcs.crval = [self.image_center_ra, self.image_center_dec]
            
            # Coordinate types
            self.wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
            
            # Pixel scale in degrees/pixel
            scale_deg_pixel = self.plate_scale / 3600.0
            
            # Create CD matrix including rotation
            # Standard convention: rotation is measured from north towards east
            cos_rot = math.cos(math.radians(self.rotation_angle))
            sin_rot = math.sin(math.radians(self.rotation_angle))
            
            # CD matrix for tangent projection with rotation
            self.wcs.wcs.cd = np.array([
                [-scale_deg_pixel * cos_rot, scale_deg_pixel * sin_rot],
                [scale_deg_pixel * sin_rot, scale_deg_pixel * cos_rot]
            ])
            
        except Exception as e:
            print(f"Error creating WCS: {str(e)}")
            self.wcs = None
    
    def radec_to_pixel(self, ra_deg, dec_deg):
        """
        Convert RA/Dec coordinates to pixel coordinates.
        
        Args:
            ra_deg (float or array): Right ascension in degrees
            dec_deg (float or array): Declination in degrees
        
        Returns:
            tuple: (x_pixels, y_pixels) coordinates
        """
        if self.wcs is None:
            print("Error: WCS not initialized")
            return None, None
        
        try:
            # Handle both single values and arrays
            ra_array = np.atleast_1d(ra_deg)
            dec_array = np.atleast_1d(dec_deg)
            
            # Convert to pixel coordinates
            sky_coords = SkyCoord(ra_array, dec_array, unit='deg')
            x_pixels, y_pixels = self.wcs.world_to_pixel(sky_coords)
            
            # Return scalars if input was scalar
            if np.isscalar(ra_deg):
                return float(x_pixels), float(y_pixels)
            else:
                return x_pixels, y_pixels
                
        except Exception as e:
            print(f"Error converting RA/Dec to pixels: {str(e)}")
            return None, None
    
    def pixel_to_radec(self, x_pixels, y_pixels):
        """
        Convert pixel coordinates to RA/Dec coordinates.
        
        Args:
            x_pixels (float or array): X pixel coordinates
            y_pixels (float or array): Y pixel coordinates
        
        Returns:
            tuple: (ra_deg, dec_deg) coordinates
        """
        if self.wcs is None:
            print("Error: WCS not initialized")
            return None, None
        
        try:
            # Handle both single values and arrays
            x_array = np.atleast_1d(x_pixels)
            y_array = np.atleast_1d(y_pixels)
            
            # Convert to world coordinates
            sky_coords = self.wcs.pixel_to_world(x_array, y_array)
            
            ra_deg = sky_coords.ra.degree
            dec_deg = sky_coords.dec.degree
            
            # Return scalars if input was scalar
            if np.isscalar(x_pixels):
                return float(ra_deg), float(dec_deg)
            else:
                return ra_deg, dec_deg
                
        except Exception as e:
            print(f"Error converting pixels to RA/Dec: {str(e)}")
            return None, None
    
    def calculate_field_corners(self):
        """
        Calculate RA/Dec coordinates of the field corners.
        
        Returns:
            dict: Corner coordinates and field boundaries
        """
        if self.wcs is None or self.image_dimensions is None:
            print("Error: WCS or image dimensions not initialized")
            return None
        
        try:
            width_pixels, height_pixels = self.image_dimensions
            
            # Define corner pixel coordinates
            corners_pixels = [
                (0, 0),  # Bottom-left
                (width_pixels, 0),  # Bottom-right
                (width_pixels, height_pixels),  # Top-right
                (0, height_pixels)  # Top-left
            ]
            
            corners_radec = []
            for x_pix, y_pix in corners_pixels:
                ra, dec = self.pixel_to_radec(x_pix, y_pix)
                if ra is not None and dec is not None:
                    corners_radec.append((ra, dec))
            
            if len(corners_radec) != 4:
                print("Error: Could not calculate all corner coordinates")
                return None
            
            # Calculate field boundaries
            ra_values = [corner[0] for corner in corners_radec]
            dec_values = [corner[1] for corner in corners_radec]
            
            # Handle RA wraparound at 0/360 degrees
            ra_min = min(ra_values)
            ra_max = max(ra_values)
            
            # Check for RA wraparound
            if ra_max - ra_min > 180:
                # Likely wraparound case
                ra_adjusted = [(ra if ra > 180 else ra + 360) for ra in ra_values]
                ra_min = min(ra_adjusted) % 360
                ra_max = max(ra_adjusted) % 360
            
            field_data = {
                'corners_radec': corners_radec,
                'ra_min': ra_min,
                'ra_max': ra_max,
                'dec_min': min(dec_values),
                'dec_max': max(dec_values),
                'center_ra': self.image_center_ra,
                'center_dec': self.image_center_dec
            }
            
            return field_data
            
        except Exception as e:
            print(f"Error calculating field corners: {str(e)}")
            return None
    
    def is_source_in_field(self, ra_deg, dec_deg, margin_pixels=0):
        """
        Check if a source is within the field of view.
        
        Args:
            ra_deg (float or array): Source RA in degrees
            dec_deg (float or array): Source Dec in degrees
            margin_pixels (float): Margin in pixels around field edge
        
        Returns:
            bool or array: True if source is in field
        """
        if self.image_dimensions is None:
            return False
        
        try:
            # Convert to pixel coordinates
            x_pixels, y_pixels = self.radec_to_pixel(ra_deg, dec_deg)
            
            if x_pixels is None or y_pixels is None:
                return False
            
            # Check if within image bounds
            width_pixels, height_pixels = self.image_dimensions
            
            in_field = ((x_pixels >= -margin_pixels) & 
                       (x_pixels <= width_pixels + margin_pixels) &
                       (y_pixels >= -margin_pixels) & 
                       (y_pixels <= height_pixels + margin_pixels))
            
            return in_field
            
        except Exception as e:
            print(f"Error checking if source is in field: {str(e)}")
            return False
    
    def filter_sources_in_field(self, source_data, margin_pixels=0):
        """
        Filter catalog sources to only those within the field of view.
        
        Args:
            source_data (dict): Source data with 'ra' and 'dec' keys
            margin_pixels (float): Margin in pixels around field edge
        
        Returns:
            dict: Filtered source data with pixel coordinates added
        """
        try:
            ra_deg = source_data['ra']
            dec_deg = source_data['dec']
            
            # Convert to pixel coordinates
            x_pixels, y_pixels = self.radec_to_pixel(ra_deg, dec_deg)
            
            if x_pixels is None or y_pixels is None:
                print("Error: Could not convert coordinates to pixels")
                return None
            
            # Check which sources are in field
            in_field_mask = self.is_source_in_field(ra_deg, dec_deg, margin_pixels)
            
            # Filter all arrays in source_data
            filtered_data = {}
            for key, values in source_data.items():
                if isinstance(values, (list, np.ndarray)) and len(values) == len(ra_deg):
                    filtered_data[key] = np.array(values)[in_field_mask]
                else:
                    filtered_data[key] = values
            
            # Add pixel coordinates
            filtered_data['x_pixels'] = x_pixels[in_field_mask]
            filtered_data['y_pixels'] = y_pixels[in_field_mask]
            
            print(f"Filtered {len(ra_deg)} sources to {len(filtered_data['ra'])} sources in field")
            
            return filtered_data
            
        except Exception as e:
            print(f"Error filtering sources in field: {str(e)}")
            return None
    
    def calculate_angular_separation(self, ra1_deg, dec1_deg, ra2_deg, dec2_deg):
        """
        Calculate angular separation between two points on the sky.
        
        Args:
            ra1_deg, dec1_deg (float): First position in degrees
            ra2_deg, dec2_deg (float): Second position in degrees
        
        Returns:
            float: Angular separation in degrees
        """
        try:
            coord1 = SkyCoord(ra1_deg, dec1_deg, unit='deg')
            coord2 = SkyCoord(ra2_deg, dec2_deg, unit='deg')
            
            separation = coord1.separation(coord2)
            return separation.degree
            
        except Exception as e:
            print(f"Error calculating angular separation: {str(e)}")
            return None
    
    def get_transform_summary(self):
        """
        Get a summary of the current coordinate transformation setup.
        
        Returns:
            dict: Summary of transformation parameters
        """
        if self.image_dimensions is None:
            return None
        
        width_pixels, height_pixels = self.image_dimensions
        
        summary = {
            'center_ra_deg': self.image_center_ra,
            'center_dec_deg': self.image_center_dec,
            'rotation_deg': self.rotation_angle,
            'plate_scale_arcsec_pixel': self.plate_scale,
            'image_width_pixels': width_pixels,
            'image_height_pixels': height_pixels,
            'field_width_arcmin': (width_pixels * self.plate_scale) / 60.0,
            'field_height_arcmin': (height_pixels * self.plate_scale) / 60.0,
            'wcs_initialized': self.wcs is not None
        }
        
        return summary