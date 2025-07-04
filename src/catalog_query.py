"""
Catalog Query Module

Handles queries to the GAIA DR3 catalog and magnitude conversions
for synthetic image generation.
"""

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import warnings

# Suppress astroquery warnings
warnings.filterwarnings('ignore', category=UserWarning, module='astroquery')


class CatalogQuery:
    """Handles astronomical catalog queries and magnitude conversions."""
    
    def __init__(self):
        """Initialize the catalog query system."""
        self.gaia_table = None
        self.converted_magnitudes = None
    
    def query_gaia_region(self, ra_deg, dec_deg, radius_arcmin, limiting_magnitude=16.0):
        """
        Query GAIA DR3 catalog for sources in a circular region.
        
        Args:
            ra_deg (float): Right ascension in degrees
            dec_deg (float): Declination in degrees
            radius_arcmin (float): Search radius in arcminutes
            limiting_magnitude (float): Limiting magnitude for query
        
        Returns:
            astropy.table.Table: GAIA query results
        """
        try:
            # Create coordinate object
            coord = SkyCoord(ra=ra_deg*u.degree, dec=dec_deg*u.degree, frame='icrs')
            
            # Convert radius to degrees
            radius_deg = radius_arcmin / 60.0
            
            # Build GAIA query
            query = f"""
            SELECT 
                source_id,
                ra, dec,
                pmra, pmdec,
                parallax,
                phot_g_mean_mag,
                phot_bp_mean_mag,
                phot_rp_mean_mag,
                bp_rp,
                g_rp
            FROM gaiadr3.gaia_source
            WHERE 
                CONTAINS(
                    POINT('ICRS', ra, dec),
                    CIRCLE('ICRS', {ra_deg}, {dec_deg}, {radius_deg})
                ) = 1
                AND phot_g_mean_mag < {limiting_magnitude}
                AND phot_g_mean_mag IS NOT NULL
                AND phot_bp_mean_mag IS NOT NULL
                AND phot_rp_mean_mag IS NOT NULL
            ORDER BY phot_g_mean_mag ASC
            """
            
            print(f"Querying GAIA DR3 catalog...")
            print(f"Center: RA={ra_deg:.4f}°, Dec={dec_deg:.4f}°")
            print(f"Radius: {radius_arcmin:.2f} arcmin")
            print(f"Limiting magnitude: {limiting_magnitude}")
            
            # Execute query
            job = Gaia.launch_job(query)
            self.gaia_table = job.get_results()
            
            print(f"Found {len(self.gaia_table)} sources")
            
            return self.gaia_table
            
        except Exception as e:
            print(f"Error querying GAIA catalog: {str(e)}")
            return None
    
    def query_gaia_rectangular(self, ra_deg, dec_deg, width_arcmin, height_arcmin, limiting_magnitude=16.0):
        """
        Query GAIA DR3 catalog for sources in a rectangular region.
        
        Args:
            ra_deg (float): Center right ascension in degrees
            dec_deg (float): Center declination in degrees
            width_arcmin (float): Width of search region in arcminutes
            height_arcmin (float): Height of search region in arcminutes
            limiting_magnitude (float): Limiting magnitude for query
        
        Returns:
            astropy.table.Table: GAIA query results
        """
        try:
            # Convert to degrees
            width_deg = width_arcmin / 60.0
            height_deg = height_arcmin / 60.0
            
            # Calculate bounds
            ra_min = ra_deg - width_deg / 2.0
            ra_max = ra_deg + width_deg / 2.0
            dec_min = dec_deg - height_deg / 2.0
            dec_max = dec_deg + height_deg / 2.0
            
            # Handle RA wraparound at 0/360 degrees
            if ra_min < 0:
                ra_min += 360
            if ra_max > 360:
                ra_max -= 360
            
            # Build GAIA query
            if ra_min > ra_max:  # RA wraparound case
                ra_condition = f"(ra >= {ra_min} OR ra <= {ra_max})"
            else:
                ra_condition = f"(ra >= {ra_min} AND ra <= {ra_max})"
            
            query = f"""
            SELECT 
                source_id,
                ra, dec,
                pmra, pmdec,
                parallax,
                phot_g_mean_mag,
                phot_bp_mean_mag,
                phot_rp_mean_mag,
                bp_rp,
                g_rp
            FROM gaiadr3.gaia_source
            WHERE 
                {ra_condition}
                AND dec >= {dec_min}
                AND dec <= {dec_max}
                AND phot_g_mean_mag < {limiting_magnitude}
                AND phot_g_mean_mag IS NOT NULL
                AND phot_bp_mean_mag IS NOT NULL
                AND phot_rp_mean_mag IS NOT NULL
            ORDER BY phot_g_mean_mag ASC
            """
            
            print(f"Querying GAIA DR3 catalog (rectangular)...")
            print(f"Center: RA={ra_deg:.4f}°, Dec={dec_deg:.4f}°")
            print(f"Size: {width_arcmin:.2f}' × {height_arcmin:.2f}'")
            print(f"Limiting magnitude: {limiting_magnitude}")
            
            # Execute query
            job = Gaia.launch_job(query)
            self.gaia_table = job.get_results()
            
            print(f"Found {len(self.gaia_table)} sources")
            
            return self.gaia_table
            
        except Exception as e:
            print(f"Error querying GAIA catalog: {str(e)}")
            return None
    
    def convert_gaia_to_johnson_cousins(self, gaia_table=None):
        """
        Convert GAIA photometry to Johnson V, B and Cousins R, I magnitudes.
        
        Uses empirical transformations from literature.
        
        Args:
            gaia_table (astropy.table.Table): GAIA catalog results
        
        Returns:
            dict: Dictionary with V, B, R, I magnitude arrays
        """
        if gaia_table is None:
            gaia_table = self.gaia_table
        
        if gaia_table is None:
            print("Error: No GAIA catalog data available")
            return None
        
        try:
            # Extract GAIA magnitudes
            g_mag = np.array(gaia_table['phot_g_mean_mag'])
            bp_mag = np.array(gaia_table['phot_bp_mean_mag'])
            rp_mag = np.array(gaia_table['phot_rp_mean_mag'])
            
            # Calculate color indices
            bp_rp = bp_mag - rp_mag
            g_rp = g_mag - rp_mag
            
            # Transformation equations (simplified empirical relations)
            # These are approximate transformations - more sophisticated ones exist
            
            # Johnson V magnitude
            V = g_mag - 0.01760 - 0.006860 * bp_rp - 0.1732 * bp_rp**2
            
            # Johnson B magnitude  
            B = g_mag + 0.3130 + 0.2271 * bp_rp + 0.01397 * bp_rp**2
            
            # Cousins R magnitude
            R = g_mag - 0.4980 - 0.0916 * bp_rp - 0.0594 * bp_rp**2
            
            # Cousins I magnitude
            I = g_mag - 0.7597 - 0.3346 * bp_rp - 0.0424 * bp_rp**2
            
            # Handle invalid values
            valid_mask = np.isfinite(V) & np.isfinite(B) & np.isfinite(R) & np.isfinite(I)
            
            self.converted_magnitudes = {
                'V': V[valid_mask],
                'B': B[valid_mask],
                'R': R[valid_mask], 
                'I': I[valid_mask],
                'ra': np.array(gaia_table['ra'])[valid_mask],
                'dec': np.array(gaia_table['dec'])[valid_mask],
                'source_id': np.array(gaia_table['source_id'])[valid_mask],
                'valid_mask': valid_mask
            }
            
            print(f"Converted magnitudes for {len(V[valid_mask])} sources")
            print(f"V magnitude range: {np.min(V[valid_mask]):.2f} to {np.max(V[valid_mask]):.2f}")
            print(f"B magnitude range: {np.min(B[valid_mask]):.2f} to {np.max(B[valid_mask]):.2f}")
            print(f"R magnitude range: {np.min(R[valid_mask]):.2f} to {np.max(R[valid_mask]):.2f}")
            print(f"I magnitude range: {np.min(I[valid_mask]):.2f} to {np.max(I[valid_mask]):.2f}")
            
            return self.converted_magnitudes
            
        except Exception as e:
            print(f"Error converting magnitudes: {str(e)}")
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
            fov_width_rad = 2 * np.arctan(sensor_width_mm / (2 * focal_length_mm))
            fov_height_rad = 2 * np.arctan(sensor_height_mm / (2 * focal_length_mm))
            
            # Convert to degrees and arcminutes
            fov_width_deg = np.degrees(fov_width_rad)
            fov_height_deg = np.degrees(fov_height_rad)
            fov_width_arcmin = fov_width_deg * 60
            fov_height_arcmin = fov_height_deg * 60
            
            # Calculate diagonal FOV
            diagonal_arcmin = np.sqrt(fov_width_arcmin**2 + fov_height_arcmin**2)
            
            # Calculate plate scale (arcsec/pixel)
            plate_scale_arcsec_pixel = 206265 * pixel_size_microns / (focal_length_mm * 1000)
            
            fov_data = {
                'width_deg': fov_width_deg,
                'height_deg': fov_height_deg,
                'width_arcmin': fov_width_arcmin,
                'height_arcmin': fov_height_arcmin,
                'diagonal_arcmin': diagonal_arcmin,
                'plate_scale_arcsec_pixel': plate_scale_arcsec_pixel
            }
            
            print(f"Field of view: {fov_width_arcmin:.2f}' × {fov_height_arcmin:.2f}'")
            print(f"Diagonal FOV: {diagonal_arcmin:.2f}'")
            print(f"Plate scale: {plate_scale_arcsec_pixel:.2f} arcsec/pixel")
            
            return fov_data
            
        except Exception as e:
            print(f"Error calculating field of view: {str(e)}")
            return None
    
    def get_sources_in_field(self, target_ra_deg, target_dec_deg, target_magnitude, 
                            telescope_config, camera_config, magnitude_offset=5.0):
        """
        Get all sources in the field of view for a given target and setup.
        
        Args:
            target_ra_deg (float): Target RA in degrees
            target_dec_deg (float): Target Dec in degrees
            target_magnitude (float): Target V magnitude
            telescope_config (dict): Telescope configuration
            camera_config (dict): Camera configuration
            magnitude_offset (float): Limiting magnitude offset from target
        
        Returns:
            dict: Sources with converted magnitudes
        """
        try:
            # Calculate field of view
            fov = self.calculate_field_of_view(
                telescope_config['focal_length_mm'],
                camera_config['pixel_size_microns'],
                camera_config['sensor_dimensions']['width_pixels'],
                camera_config['sensor_dimensions']['height_pixels']
            )
            
            if fov is None:
                return None
            
            # Set limiting magnitude
            limiting_magnitude = target_magnitude + magnitude_offset
            
            # Query catalog
            catalog = self.query_gaia_rectangular(
                target_ra_deg, target_dec_deg,
                fov['width_arcmin'], fov['height_arcmin'],
                limiting_magnitude
            )
            
            if catalog is None:
                return None
            
            # Convert magnitudes
            converted_mags = self.convert_gaia_to_johnson_cousins(catalog)
            
            if converted_mags is None:
                return None
            
            # Add field of view information
            converted_mags['field_of_view'] = fov
            converted_mags['limiting_magnitude'] = limiting_magnitude
            converted_mags['target_ra'] = target_ra_deg
            converted_mags['target_dec'] = target_dec_deg
            
            return converted_mags
            
        except Exception as e:
            print(f"Error getting sources in field: {str(e)}")
            return None
    
    def lookup_target(self, target_name):
        """
        Lookup target information from Simbad database.
        
        Args:
            target_name (str): Target name to lookup
        
        Returns:
            dict: Target information including coordinates and magnitudes
        """
        try:
            # Configure Simbad to return the data we need
            Simbad.add_votable_fields('B', 'V', 'R', 'I')
            
            print(f"Looking up target: {target_name}")
            
            # Query Simbad
            result_table = Simbad.query_object(target_name)
            
            if result_table is None or len(result_table) == 0:
                print(f"No results found for target: {target_name}")
                return None
            
            # Extract the first result
            target_data = result_table[0]
            
            # Get coordinates (handle direct degree format)
            ra_deg = float(target_data['ra'])
            dec_deg = float(target_data['dec'])
            coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg)
            
            # Extract magnitudes (handle missing values)
            def get_magnitude(flux_field):
                try:
                    value = target_data[flux_field]
                    if value is None or np.ma.is_masked(value):
                        return None
                    return float(value)
                except:
                    return None
            
            # Format coordinates with proper precision
            ra_deg_formatted = float(f"{coord.ra.deg:.6g}")  # 6 significant digits
            dec_deg_formatted = float(f"{coord.dec.deg:+.6g}")  # 6 significant digits with sign
            
            # Format RA string with 2 decimal places for seconds
            ra_str_formatted = coord.ra.to_string(unit=u.hourangle, sep=':', precision=2)
            
            # Format Dec string with 1 decimal place for seconds and always show sign
            dec_str_formatted = coord.dec.to_string(unit=u.deg, sep=':', precision=1, alwayssign=True)
            
            # Format magnitudes to 3 decimal places
            def format_magnitude(mag):
                if mag is None:
                    return None
                return round(float(mag), 3)
            
            target_info = {
                'name': target_data['main_id'].decode('utf-8') if isinstance(target_data['main_id'], bytes) else str(target_data['main_id']),
                'ra_deg': ra_deg_formatted,
                'dec_deg': dec_deg_formatted,
                'ra_str': ra_str_formatted,
                'dec_str': dec_str_formatted,
                'magnitudes': {
                    'B': format_magnitude(get_magnitude('B')),
                    'V': format_magnitude(get_magnitude('V')),
                    'R': format_magnitude(get_magnitude('R')),
                    'I': format_magnitude(get_magnitude('I'))
                }
            }
            
            print(f"Target found: {target_info['name']}")
            print(f"Coordinates: RA={target_info['ra_deg']:.6f}°, Dec={target_info['dec_deg']:.6f}°")
            print(f"Magnitudes: B={target_info['magnitudes']['B']}, V={target_info['magnitudes']['V']}, R={target_info['magnitudes']['R']}, I={target_info['magnitudes']['I']}")
            
            return target_info
            
        except Exception as e:
            print(f"Error looking up target {target_name}: {str(e)}")
            return None