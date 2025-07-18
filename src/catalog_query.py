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
try:
    from .gaia_photometry import GaiaPhotometry
except ImportError:
    from gaia_photometry import GaiaPhotometry

# Suppress astroquery warnings
warnings.filterwarnings('ignore', category=UserWarning, module='astroquery')


class CatalogQuery:
    """Handles astronomical catalog queries and magnitude conversions."""
    
    def __init__(self):
        """Initialize the catalog query system."""
        self.gaia_table = None
        self.converted_magnitudes = None
        self.gaia_photometry = GaiaPhotometry()
        self.use_gaia_spectroscopy = True
    
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
    
    def convert_gaia_to_johnson_cousins_spectroscopic(self, gaia_table=None):
        """
        Convert GAIA photometry to Johnson-Kron-Cousins magnitudes using high-resolution spectra.
        
        Uses gaiaxpy to generate synthetic photometry from BP/RP spectra.
        Falls back to polynomial method if spectroscopic data is unavailable.
        
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
            # Extract source IDs for spectroscopic processing
            source_ids = gaia_table['source_id'].tolist()
            
            # Attempt GAIA spectroscopic conversion
            magnitudes_dict, errors_dict = self.gaia_photometry.process_source_list(source_ids)
            
            if magnitudes_dict:
                # Convert to arrays for consistency with original interface
                V_list, B_list, R_list, I_list = [], [], [], []
                ra_list, dec_list, source_id_list = [], [], []
                
                for i, source_id in enumerate(source_ids):
                    if source_id in magnitudes_dict:
                        mags = magnitudes_dict[source_id]
                        
                        # Only include sources with valid V, B, R, I magnitudes
                        if all(np.isfinite(mags[band]) for band in ['V', 'B', 'R', 'I']):
                            V_list.append(mags['V'])
                            B_list.append(mags['B'])
                            R_list.append(mags['R'])
                            I_list.append(mags['I'])
                            ra_list.append(gaia_table['ra'][i])
                            dec_list.append(gaia_table['dec'][i])
                            source_id_list.append(source_id)
                
                if V_list:
                    V = np.array(V_list)
                    B = np.array(B_list)
                    R = np.array(R_list)
                    I = np.array(I_list)
                    
                    # Track which magnitudes were derived spectroscopically
                    spectroscopic_bands = {i: ['U', 'B', 'V', 'R', 'I'] for i in range(len(V_list))}
                    
                    self.converted_magnitudes = {
                        'V': V,
                        'B': B,
                        'R': R,
                        'I': I,
                        'U': V,  # U magnitude (placeholder - would need specific conversion)
                        'ra': np.array(ra_list),
                        'dec': np.array(dec_list),
                        'source_id': np.array(source_id_list),
                        'method': 'spectroscopic',
                        'spectroscopic_bands': spectroscopic_bands,
                        'V_error': np.full(len(V), 0.01),  # Default error estimate
                        'B_error': np.full(len(B), 0.01),
                        'R_error': np.full(len(R), 0.01),
                        'I_error': np.full(len(I), 0.01),
                        'U_error': np.full(len(V), 0.02),  # Higher error for U
                    }
                    
                    print(f"Spectroscopic conversion successful for {len(V)} sources")
                    print(f"V magnitude range: {np.min(V):.2f} to {np.max(V):.2f}")
                    print(f"B magnitude range: {np.min(B):.2f} to {np.max(B):.2f}")
                    print(f"R magnitude range: {np.min(R):.2f} to {np.max(R):.2f}")
                    print(f"I magnitude range: {np.min(I):.2f} to {np.max(I):.2f}")
                    
                    return self.converted_magnitudes
            
            # Fallback to polynomial method if spectroscopic conversion fails
            print("Falling back to polynomial magnitude conversion")
            return self.convert_gaia_to_johnson_cousins_polynomial(gaia_table)
            
        except Exception as e:
            print(f"Error in spectroscopic conversion: {str(e)}")
            print("Falling back to polynomial magnitude conversion")
            return self.convert_gaia_to_johnson_cousins_polynomial(gaia_table)
    
    def convert_gaia_to_johnson_cousins_polynomial(self, gaia_table=None):
        """
        Convert GAIA photometry to Johnson V, B and Cousins R, I magnitudes using polynomial transforms.
        
        Uses empirical polynomial transformations (fallback method).
        
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
            
            # Transformation equations (polynomial approximations)
            # These are approximate transformations - spectroscopic method is more accurate
            
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
            
            # Add U magnitude approximation
            U = g_mag + 0.7840 + 0.5384 * bp_rp + 0.0755 * bp_rp**2
            
            # Track that all magnitudes used polynomial method (no spectroscopic data)
            polynomial_bands = {i: [] for i in range(len(V[valid_mask]))}  # Empty = all polynomial
            
            self.converted_magnitudes = {
                'V': V[valid_mask],
                'B': B[valid_mask],
                'R': R[valid_mask], 
                'I': I[valid_mask],
                'U': U[valid_mask],
                'ra': np.array(gaia_table['ra'])[valid_mask],
                'dec': np.array(gaia_table['dec'])[valid_mask],
                'source_id': np.array(gaia_table['source_id'])[valid_mask],
                'valid_mask': valid_mask,
                'method': 'polynomial',
                'spectroscopic_bands': polynomial_bands,
                'V_error': np.full(len(V[valid_mask]), 0.05),  # Higher error for polynomial
                'B_error': np.full(len(B[valid_mask]), 0.05),
                'R_error': np.full(len(R[valid_mask]), 0.05),
                'I_error': np.full(len(I[valid_mask]), 0.05),
                'U_error': np.full(len(U[valid_mask]), 0.10),  # Even higher error for U polynomial
            }
            
            print(f"Polynomial conversion for {len(V[valid_mask])} sources")
            print(f"V magnitude range: {np.min(V[valid_mask]):.2f} to {np.max(V[valid_mask]):.2f}")
            print(f"B magnitude range: {np.min(B[valid_mask]):.2f} to {np.max(B[valid_mask]):.2f}")
            print(f"R magnitude range: {np.min(R[valid_mask]):.2f} to {np.max(R[valid_mask]):.2f}")
            print(f"I magnitude range: {np.min(I[valid_mask]):.2f} to {np.max(I[valid_mask]):.2f}")
            
            return self.converted_magnitudes
            
        except Exception as e:
            print(f"Error converting magnitudes: {str(e)}")
            return None
    
    def convert_gaia_to_johnson_cousins(self, gaia_table=None):
        """
        Convert GAIA photometry to Johnson V, B and Cousins R, I magnitudes.
        
        Automatically selects the best available method:
        1. Spectroscopic conversion using gaiaxpy (preferred)
        2. Polynomial conversion (fallback)
        
        Args:
            gaia_table (astropy.table.Table): GAIA catalog results
        
        Returns:
            dict: Dictionary with V, B, R, I magnitude arrays
        """
        if self.use_gaia_spectroscopy:
            return self.convert_gaia_to_johnson_cousins_spectroscopic(gaia_table)
        else:
            return self.convert_gaia_to_johnson_cousins_polynomial(gaia_table)
    
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
            
            # Preserve full coordinate precision (no rounding)
            ra_deg_formatted = coord.ra.deg  # Full precision
            dec_deg_formatted = coord.dec.deg  # Full precision
            
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
    
    def lookup_target_by_coordinates(self, ra_deg, dec_deg, search_radius_arcsec=30.0):
        """
        Lookup target information from Simbad database by coordinates.
        
        Args:
            ra_deg (float): Right ascension in degrees
            dec_deg (float): Declination in degrees
            search_radius_arcsec (float): Search radius in arcseconds
        
        Returns:
            dict: Target information including coordinates and magnitudes, or None if not found
        """
        try:
            # Configure Simbad to return the data we need
            Simbad.add_votable_fields('B', 'V', 'R', 'I')
            
            # Create coordinate object
            coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
            
            print(f"Looking up target at coordinates: RA={ra_deg:.6f}°, Dec={dec_deg:.6f}°")
            print(f"Search radius: {search_radius_arcsec:.1f} arcsec")
            
            # Query Simbad around the coordinates
            result_table = Simbad.query_region(coord, radius=search_radius_arcsec*u.arcsec)
            
            if result_table is None or len(result_table) == 0:
                print(f"No objects found at coordinates RA={ra_deg:.6f}°, Dec={dec_deg:.6f}°")
                return None
            
            # Calculate distances for all results and find the closest
            min_distance = float('inf')
            closest_target = None
            
            for row in result_table:
                # Get coordinates from each result
                obj_ra_deg = float(row['ra'])
                obj_dec_deg = float(row['dec'])
                obj_coord = SkyCoord(ra=obj_ra_deg*u.deg, dec=obj_dec_deg*u.deg)
                
                # Calculate separation
                separation = coord.separation(obj_coord).arcsec
                
                if separation < min_distance:
                    min_distance = separation
                    closest_target = row
            
            # Use the closest target
            target_data = closest_target
            print(f"Selected closest object at {min_distance:.2f} arcsec separation")
            
            # Get coordinates from Simbad result (handle direct degree format)
            simbad_ra_deg = float(target_data['ra'])
            simbad_dec_deg = float(target_data['dec'])
            simbad_coord = SkyCoord(ra=simbad_ra_deg*u.deg, dec=simbad_dec_deg*u.deg)
            
            # Extract magnitudes (handle missing values)
            def get_magnitude(flux_field):
                try:
                    value = target_data[flux_field]
                    if value is None or np.ma.is_masked(value):
                        return None
                    return float(value)
                except:
                    return None
            
            # Format RA string with 2 decimal places for seconds
            ra_str_formatted = simbad_coord.ra.to_string(unit=u.hourangle, sep=':', precision=2)
            
            # Format Dec string with 1 decimal place for seconds and always show sign
            dec_str_formatted = simbad_coord.dec.to_string(unit=u.deg, sep=':', precision=1, alwayssign=True)
            
            # Format magnitudes to 3 decimal places
            def format_magnitude(mag):
                if mag is None:
                    return None
                return round(float(mag), 3)
            
            target_info = {
                'name': target_data['main_id'].decode('utf-8') if isinstance(target_data['main_id'], bytes) else str(target_data['main_id']),
                'ra_deg': simbad_ra_deg,
                'dec_deg': simbad_dec_deg,
                'ra_str': ra_str_formatted,
                'dec_str': dec_str_formatted,
                'magnitudes': {
                    'B': format_magnitude(get_magnitude('B')),
                    'V': format_magnitude(get_magnitude('V')),
                    'R': format_magnitude(get_magnitude('R')),
                    'I': format_magnitude(get_magnitude('I'))
                }
            }
            
            # Calculate separation from search coordinates
            search_coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg)
            separation = search_coord.separation(simbad_coord).arcsec
            
            print(f"Target found: {target_info['name']}")
            print(f"Simbad coordinates: RA={target_info['ra_deg']:.6f}°, Dec={target_info['dec_deg']:.6f}°")
            print(f"Separation from search: {separation:.2f} arcsec")
            print(f"Magnitudes: B={target_info['magnitudes']['B']}, V={target_info['magnitudes']['V']}, R={target_info['magnitudes']['R']}, I={target_info['magnitudes']['I']}")
            
            return target_info
            
        except Exception as e:
            print(f"Error looking up coordinates RA={ra_deg:.6f}°, Dec={dec_deg:.6f}°: {str(e)}")
            return None

    def find_gaia_source_for_target(self, ra_deg, dec_deg, search_radius_arcsec=10.0):
        """
        Find GAIA source_id and DR3 magnitudes for a target at given coordinates.
        
        Args:
            ra_deg (float): Right ascension in degrees
            dec_deg (float): Declination in degrees  
            search_radius_arcsec (float): Search radius in arcseconds (default 10.0)
            
        Returns:
            dict or None: GAIA source information with source_id and DR3 magnitudes
        """
        try:
            # Create coordinate object
            coord = SkyCoord(ra=ra_deg*u.degree, dec=dec_deg*u.degree, frame='icrs')
            
            # Convert radius to degrees
            radius_deg = search_radius_arcsec / 3600.0
            
            # Build GAIA query for nearby sources
            query = f"""
            SELECT TOP 5
                source_id,
                ra, dec,
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
                AND phot_g_mean_mag IS NOT NULL
                AND phot_bp_mean_mag IS NOT NULL
                AND phot_rp_mean_mag IS NOT NULL
            ORDER BY phot_g_mean_mag ASC
            """
            
            print(f"Searching for GAIA source near RA={ra_deg:.6f}°, Dec={dec_deg:.6f}°")
            print(f"Search radius: {search_radius_arcsec:.1f} arcsec")
            
            # Execute query
            job = Gaia.launch_job(query)
            results = job.get_results()
            
            if len(results) == 0:
                print("No GAIA sources found within search radius")
                return None
            
            # Find closest source
            best_source = None
            best_separation = float('inf')
            
            for source in results:
                source_coord = SkyCoord(ra=source['ra']*u.degree, dec=source['dec']*u.degree, frame='icrs')
                separation = coord.separation(source_coord).arcsecond
                
                if separation < best_separation:
                    best_separation = separation
                    best_source = source
            
            if best_source is None:
                print("No suitable GAIA source found")
                return None
                
            # Get GAIA DR3 magnitudes using spectroscopic method with polynomial fallback
            source_id = int(best_source['source_id'])
            magnitudes, errors = self.gaia_photometry.process_source_list([source_id])
            
            # Track which method was used for marking with asterisk
            used_polynomial = False
            
            if magnitudes is None or len(magnitudes) == 0:
                print(f"Spectroscopic photometry failed for source_id {source_id}, falling back to polynomial method")
                # Fallback to polynomial method
                mag_data, used_polynomial = self._get_polynomial_magnitudes_for_source(best_source)
                if mag_data is None:
                    print(f"Failed to get GAIA DR3 magnitudes for source_id {source_id}")
                    return None
            else:
                mag_data = magnitudes[source_id]
            
            gaia_info = {
                'source_id': source_id,
                'ra_deg': float(best_source['ra']),
                'dec_deg': float(best_source['dec']),
                'separation_arcsec': best_separation,
                'magnitudes': {
                    'B': round(mag_data['B'], 3) if mag_data['B'] is not None else None,
                    'V': round(mag_data['V'], 3) if mag_data['V'] is not None else None,
                    'R': round(mag_data['R'], 3) if mag_data['R'] is not None else None,
                    'I': round(mag_data['I'], 3) if mag_data['I'] is not None else None
                },
                'used_polynomial': used_polynomial
            }
            
            print(f"Found GAIA source_id: {source_id}")
            print(f"GAIA position: RA={gaia_info['ra_deg']:.6f}°, Dec={gaia_info['dec_deg']:.6f}°")
            print(f"Separation from target: {best_separation:.2f} arcsec")
            print(f"GAIA DR3 magnitudes: B={gaia_info['magnitudes']['B']}, V={gaia_info['magnitudes']['V']}, R={gaia_info['magnitudes']['R']}, I={gaia_info['magnitudes']['I']}")
            
            return gaia_info
            
        except Exception as e:
            print(f"Error finding GAIA source for target: {str(e)}")
            return None
    
    def _get_polynomial_magnitudes_for_source(self, source):
        """
        Calculate polynomial-based magnitudes for a single GAIA source.
        
        Args:
            source: GAIA source record with photometry
            
        Returns:
            tuple: (magnitude_dict, used_polynomial_flag) or (None, False) if failed
        """
        try:
            # Extract GAIA photometry
            g_mag = float(source['phot_g_mean_mag'])
            bp_mag = float(source['phot_bp_mean_mag']) if source['phot_bp_mean_mag'] is not None else None
            rp_mag = float(source['phot_rp_mean_mag']) if source['phot_rp_mean_mag'] is not None else None
            
            # Calculate BP-RP color
            if bp_mag is not None and rp_mag is not None:
                bp_rp = bp_mag - rp_mag
            else:
                print("Missing BP or RP photometry for polynomial conversion")
                return None, False
            
            # Check for valid color range (typical range is -0.5 to 4.0)
            if not (-1.0 <= bp_rp <= 5.0):
                print(f"BP-RP color {bp_rp:.3f} outside valid range for polynomial conversion")
                return None, False
            
            # Polynomial transformation equations (same as in convert_gaia_to_johnson_cousins_polynomial)
            V = g_mag - 0.01760 - 0.006860 * bp_rp - 0.1732 * bp_rp**2
            B = g_mag + 0.3130 + 0.2271 * bp_rp + 0.01397 * bp_rp**2
            R = g_mag - 0.4980 - 0.0916 * bp_rp - 0.0594 * bp_rp**2
            I = g_mag - 0.7597 - 0.1311 * bp_rp - 0.0753 * bp_rp**2
            
            print(f"Polynomial conversion: V={V:.3f}, B={B:.3f}, R={R:.3f}, I={I:.3f}")
            
            return {
                'B': B if np.isfinite(B) else None,
                'V': V if np.isfinite(V) else None,
                'R': R if np.isfinite(R) else None,
                'I': I if np.isfinite(I) else None
            }, True
            
        except Exception as e:
            print(f"Error in polynomial magnitude calculation: {str(e)}")
            return None, False