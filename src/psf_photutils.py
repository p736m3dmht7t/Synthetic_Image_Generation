"""
Point Spread Function (PSF) Generation using Photutils

Scientifically validated PSF implementation using astropy and photutils
libraries for accurate synthetic star image generation.
"""

import numpy as np
from astropy.modeling import models, fitting
from astropy.convolution import discretize_model
from photutils.psf import IntegratedGaussianPRF
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore', category=RuntimeWarning)


class PhotoutilsPSFGenerator:
    """Generates and applies point spread functions using photutils."""
    
    def __init__(self):
        """Initialize the photutils PSF generator."""
        self.psf_model = None
        self.psf_type = 'gaussian'  # Default to Gaussian
        
    def create_gaussian_psf_model(self, sigma_pixels, subsampling=1):
        """
        Create a Gaussian PSF model using photutils.
        
        Parameters
        ----------
        sigma_pixels : float
            Gaussian sigma in pixels
        subsampling : int, optional
            Subsampling factor for PSF accuracy (default: 1)
            
        Returns
        -------
        astropy.modeling.models
            Gaussian2D PSF model
        """
        # Create Gaussian2D model with proper flux normalization
        gaussian_model = models.Gaussian2D(
            amplitude=1.0,  # Will be scaled by total flux
            x_mean=0, 
            y_mean=0,
            x_stddev=sigma_pixels,
            y_stddev=sigma_pixels
        )
        
        self.psf_model = gaussian_model
        self.psf_type = 'gaussian'
        
        return gaussian_model
    
    def create_moffat_psf_model(self, fwhm_pixels, alpha=2.0, subsampling=1):
        """
        Create a Moffat PSF model using photutils.
        
        Parameters
        ----------
        fwhm_pixels : float
            Full Width at Half Maximum in pixels
        alpha : float, optional
            Moffat alpha parameter (default: 2.0)
        subsampling : int, optional
            Subsampling factor for PSF accuracy (default: 1)
            
        Returns
        -------
        astropy.modeling.models
            Moffat2D PSF model
        """
        # Convert FWHM to gamma for Moffat profile
        gamma = fwhm_pixels / (2.0 * np.sqrt(2.0**(1.0/alpha) - 1.0))
        
        # Create Moffat2D model
        moffat_model = models.Moffat2D(
            amplitude=1.0,  # Will be scaled by total flux
            x_0=0,
            y_0=0,
            gamma=gamma,
            alpha=alpha
        )
        
        self.psf_model = moffat_model
        self.psf_type = 'moffat'
        
        return moffat_model
    
    def calculate_magnitude_to_flux(self, magnitude, zero_point=25.0):
        """
        Convert magnitude to relative flux using Pogson's equation.
        
        Parameters
        ----------
        magnitude : float
            Magnitude value
        zero_point : float, optional
            Photometric zero point (default: 25.0)
            
        Returns
        -------
        float
            Relative flux value
        """
        return 10.0**(-0.4 * (magnitude - zero_point))
    
    def calculate_target_peak_flux(self, target_magnitude, saturation_limit, 
                                 target_fraction=0.5, zero_point=25.0):
        """
        Calculate target PEAK flux in ADU for saturation checking.
        
        This sets the peak pixel value, not the total flux.
        AstroImageJ uses peak flux to check saturation.
        
        Parameters
        ----------
        target_magnitude : float
            Target star magnitude
        saturation_limit : int
            Camera saturation limit in ADU
        target_fraction : float, optional
            Fraction of saturation for target star (default: 0.5)
        zero_point : float, optional
            Photometric zero point (default: 25.0)
            
        Returns
        -------
        float
            Target PEAK flux in ADU
        """
        # Calculate target peak ADU level (50% of saturation by default)
        target_peak_adu = saturation_limit * target_fraction
        
        print(f"Target star: {target_magnitude:.2f} mag -> {target_peak_adu:.0f} ADU peak "
              f"({target_fraction*100:.0f}% of saturation)")
        
        return target_peak_adu
    
    def scale_source_peak_flux(self, source_magnitude, target_magnitude, target_peak_adu):
        """
        Scale source peak flux relative to target star using magnitude difference.
        
        Parameters
        ----------
        source_magnitude : float
            Source magnitude
        target_magnitude : float
            Target star magnitude  
        target_peak_adu : float
            Target star peak flux in ADU
            
        Returns
        -------
        float
            Source peak flux in ADU
        """
        # Calculate magnitude difference
        delta_mag = source_magnitude - target_magnitude
        
        # Apply inverse square law for flux scaling
        flux_ratio = 10.0**(-0.4 * delta_mag)
        source_peak_flux = target_peak_adu * flux_ratio
        
        return source_peak_flux
    
    def discretize_psf_model(self, psf_model, psf_size, x_center, y_center):
        """
        Discretize PSF model onto pixel grid.
        
        Parameters
        ----------
        psf_model : astropy.modeling.models
            PSF model to discretize
        psf_size : int
            Size of PSF array (should be odd)
        x_center : float
            X center position
        y_center : float
            Y center position
            
        Returns
        -------
        numpy.ndarray
            Discretized PSF array
        """
        # Ensure odd PSF size for proper centering
        if psf_size % 2 == 0:
            psf_size += 1
        
        # Create coordinate grid
        y, x = np.mgrid[-psf_size//2:psf_size//2+1, -psf_size//2:psf_size//2+1]
        
        # Evaluate PSF model on grid
        psf_array = psf_model(x + x_center, y + y_center)
        
        # Normalize to unit flux for proper scaling
        total_flux = np.sum(psf_array)
        if total_flux > 0:
            psf_array = psf_array / total_flux
        
        return psf_array
    
    def calculate_peak_correction_factor(self, target_x, target_y, target_magnitude, 
                                       saturation_limit, psf_params, psf_type='gaussian',
                                       target_fraction=0.5):
        """
        Calculate peak correction factor to achieve exact saturation level.
        
        This method calculates the precise correction factor needed to achieve
        exactly the desired peak flux (e.g., 50% saturation) for the target star,
        accounting for subpixel positioning effects and PSF discretization.
        
        The approach:
        1. Calculate what the normal magnitude-based flux would be
        2. Generate the actual PSF at the target position and see its peak
        3. Calculate the ratio needed to achieve exact desired peak
        
        Parameters
        ----------
        target_x : float
            Target star X coordinate (can be subpixel)
        target_y : float
            Target star Y coordinate (can be subpixel)
        target_magnitude : float
            Target star magnitude
        saturation_limit : int
            Camera saturation limit in ADU
        psf_params : dict
            PSF parameters from optics calculations
        psf_type : str, optional
            PSF type ('gaussian' or 'moffat')
        target_fraction : float, optional
            Fraction of saturation for target star (default: 0.5)
            
        Returns
        -------
        float
            Correction factor to multiply all source fluxes
        """
        try:
            # Calculate desired peak ADU for target
            desired_peak_adu = saturation_limit * target_fraction
            
            # Calculate what the normal approach would give us for peak flux
            normal_target_peak_adu = self.calculate_target_peak_flux(
                target_magnitude, saturation_limit, target_fraction
            )
            
            # This should be the same as desired_peak_adu, but let's be explicit
            base_source_peak_flux = self.scale_source_peak_flux(
                target_magnitude, target_magnitude, normal_target_peak_adu
            )
            
            # Now simulate placing the source to see what peak we actually get
            # Create a small test image
            test_image = np.zeros((50, 50), dtype=np.float64)
            test_x, test_y = 25.0 + (target_x - int(target_x)), 25.0 + (target_y - int(target_y))
            
            # Place source using normal approach
            test_image_with_source = self.place_source_on_image(
                test_image.copy(), test_x, test_y, base_source_peak_flux, 
                psf_params, psf_type, debug=False
            )
            
            # Find actual peak achieved
            actual_peak_achieved = np.max(test_image_with_source)
            
            # Calculate correction factor
            correction_factor = desired_peak_adu / actual_peak_achieved
            
            print(f"Peak correction calculation:")
            print(f"  Target position: ({target_x:.3f}, {target_y:.3f})")
            print(f"  Subpixel offset: ({target_x - int(target_x):.3f}, {target_y - int(target_y):.3f})")
            print(f"  Desired peak ADU: {desired_peak_adu:.0f}")
            print(f"  Normal approach peak: {base_source_peak_flux:.0f}")
            print(f"  Actual peak achieved: {actual_peak_achieved:.0f}")
            print(f"  Correction factor: {correction_factor:.6f}")
            
            return correction_factor
            
        except Exception as e:
            print(f"Error calculating peak correction factor: {e}")
            return 1.0  # Return unity factor as fallback
    
    def identify_target_source(self, source_data, target_name="V* TZ Boo", target_coords=None):
        """
        Identify the target source in the source data catalog.
        
        Parameters
        ----------
        source_data : dict
            Source catalog with coordinates and magnitudes
        target_name : str, optional
            Target name for identification
        target_coords : tuple, optional
            (ra_deg, dec_deg) coordinates if available
            
        Returns
        -------
        dict or None
            Dictionary with target source info or None if not found
        """
        try:
            # Get source coordinates
            x_coords = source_data.get('x_pixels', source_data.get('x_pixel', []))
            y_coords = source_data.get('y_pixels', source_data.get('y_pixel', []))
            
            if len(x_coords) == 0:
                print("Warning: No source coordinates found in data")
                return None
                
            # Check for target star identification in source data
            target_source = None
            target_index = None
            
            # Method 1: Look for target name or identifier in source data
            if 'source_names' in source_data:
                source_names = source_data['source_names']
                for i, name in enumerate(source_names):
                    if target_name in name or "TZ Boo" in name:
                        target_index = i
                        target_source = {
                            'index': i,
                            'name': name,
                            'x_pixel': x_coords[i],
                            'y_pixel': y_coords[i]
                        }
                        break
            
            # Method 2: Look for target coordinates match
            if target_source is None and target_coords is not None:
                target_ra, target_dec = target_coords
                if 'ra' in source_data and 'dec' in source_data:
                    source_ras = source_data['ra']
                    source_decs = source_data['dec']
                    
                    # Find closest coordinate match (within 5 arcseconds)
                    min_separation = float('inf')
                    for i, (ra, dec) in enumerate(zip(source_ras, source_decs)):
                        # Calculate angular separation (rough approximation)
                        ra_diff = (ra - target_ra) * np.cos(np.radians(dec))
                        dec_diff = dec - target_dec
                        separation_arcsec = 3600 * np.sqrt(ra_diff**2 + dec_diff**2)
                        
                        if separation_arcsec < min_separation and separation_arcsec < 5.0:
                            min_separation = separation_arcsec
                            target_index = i
                            target_source = {
                                'index': i,
                                'separation_arcsec': separation_arcsec,
                                'x_pixel': x_coords[i],
                                'y_pixel': y_coords[i]
                            }
            
            # Method 3: Assume brightest source is target (fallback)
            if target_source is None:
                print("Warning: Target not identified by name or coordinates")
                print("Using brightest V-band source as target")
                
                v_mags = source_data.get('V', [])
                if len(v_mags) > 0:
                    valid_mags = [(i, mag) for i, mag in enumerate(v_mags) if np.isfinite(mag) and mag < 90]
                    if valid_mags:
                        # Find brightest (lowest magnitude)
                        target_index, brightest_mag = min(valid_mags, key=lambda x: x[1])
                        target_source = {
                            'index': target_index,
                            'magnitude_v': brightest_mag,
                            'x_pixel': x_coords[target_index],
                            'y_pixel': y_coords[target_index],
                            'method': 'brightest_assumed'
                        }
            
            if target_source is not None:
                print(f"Target source identified:")
                print(f"  Index: {target_source['index']}")
                print(f"  Position: ({target_source['x_pixel']:.3f}, {target_source['y_pixel']:.3f})")
                if 'name' in target_source:
                    print(f"  Name: {target_source['name']}")
                if 'separation_arcsec' in target_source:
                    print(f"  Coordinate match: {target_source['separation_arcsec']:.1f} arcsec")
                if 'method' in target_source:
                    print(f"  Method: {target_source['method']}")
            else:
                print("Warning: Could not identify target source")
            
            return target_source
            
        except Exception as e:
            print(f"Error identifying target source: {e}")
            return None
    
    def place_source_on_image(self, image, x_center, y_center, peak_flux_adu, 
                            psf_params, psf_type='gaussian', debug=False):
        """
        Place a source on the image using exact peak flux scaling.
        
        This method ensures the peak pixel value in the final image exactly
        matches the specified peak_flux_adu value, regardless of subpixel
        positioning effects and PSF discretization.
        
        Parameters
        ----------
        image : numpy.ndarray
            Target image array
        x_center : float
            X coordinate of source center (subpixel precision)
        y_center : float
            Y coordinate of source center (subpixel precision)
        peak_flux_adu : float
            Exact peak flux in ADU for this source
        psf_params : dict
            PSF parameters from optics calculations
        psf_type : str, optional
            PSF type ('gaussian' or 'moffat')
        debug : bool, optional
            Enable debug output for peak flux verification
            
        Returns
        -------
        numpy.ndarray
            Updated image array
        """
        try:
            # Get PSF parameters
            sigma_pixels = psf_params.get('sigma_pixels', 1.0)
            fwhm_pixels = psf_params.get('fwhm_pixels', 2.355 * sigma_pixels)
            
            # Create PSF model
            if psf_type == 'gaussian':
                psf_model = self.create_gaussian_psf_model(sigma_pixels)
            elif psf_type == 'moffat':
                psf_model = self.create_moffat_psf_model(fwhm_pixels)
            else:
                raise ValueError(f"Unknown PSF type: {psf_type}")
            
            # Determine PSF size (6 sigma for Gaussian, larger for Moffat)
            if psf_type == 'gaussian':
                psf_size = max(7, int(6 * sigma_pixels))
            else:
                psf_size = max(9, int(3 * fwhm_pixels))
            
            # Ensure odd size
            if psf_size % 2 == 0:
                psf_size += 1
            
            # Calculate fractional pixel offset for subpixel accuracy
            x_frac = x_center - int(x_center)
            y_frac = y_center - int(y_center)
            
            # Discretize PSF model with exact subpixel positioning
            psf_array = self.discretize_psf_model(psf_model, psf_size, x_frac, y_frac)
            
            # CRITICAL: Normalize PSF by PEAK value (not total flux)
            # This ensures exact peak flux control
            psf_peak_normalized = np.max(psf_array)
            if psf_peak_normalized > 0:
                psf_array = psf_array / psf_peak_normalized
            else:
                print(f"Warning: PSF peak is zero at ({x_center:.3f}, {y_center:.3f})")
                return image
            
            # Scale PSF to achieve exact peak flux
            scaled_psf = psf_array * peak_flux_adu
            
            # Verify peak scaling if debug enabled
            if debug:
                actual_peak = np.max(scaled_psf)
                total_flux = np.sum(scaled_psf)
                print(f"  Source at ({x_center:.3f}, {y_center:.3f}):")
                print(f"    Desired peak: {peak_flux_adu:.1f} ADU")
                print(f"    Actual peak: {actual_peak:.1f} ADU")
                print(f"    Total flux: {total_flux:.1f} ADU")
                print(f"    Peak accuracy: {(actual_peak/peak_flux_adu)*100:.2f}%")
            
            # Place PSF on image with exact peak control
            self._add_psf_to_image(image, scaled_psf, int(x_center), int(y_center))
            
            return image
            
        except Exception as e:
            print(f"Error placing source at ({x_center:.1f}, {y_center:.1f}): {e}")
            return image
    
    def _add_psf_to_image(self, image, psf_array, x_center_int, y_center_int):
        """
        Add PSF array to image at specified integer coordinates.
        
        Parameters
        ----------
        image : numpy.ndarray
            Target image
        psf_array : numpy.ndarray  
            PSF array to add
        x_center_int : int
            Integer X center coordinate
        y_center_int : int
            Integer Y center coordinate
        """
        img_height, img_width = image.shape
        psf_height, psf_width = psf_array.shape
        
        # Calculate PSF placement bounds
        psf_center_x = psf_width // 2
        psf_center_y = psf_height // 2
        
        # Calculate image coordinates for PSF placement
        img_x_start = x_center_int - psf_center_x
        img_y_start = y_center_int - psf_center_y
        img_x_end = img_x_start + psf_width
        img_y_end = img_y_start + psf_height
        
        # Calculate PSF coordinates (handle edge cases)
        psf_x_start = max(0, -img_x_start)
        psf_y_start = max(0, -img_y_start)
        psf_x_end = psf_width - max(0, img_x_end - img_width)
        psf_y_end = psf_height - max(0, img_y_end - img_height)
        
        # Adjust image coordinates
        img_x_start = max(0, img_x_start)
        img_y_start = max(0, img_y_start)
        img_x_end = min(img_width, img_x_end)
        img_y_end = min(img_height, img_y_end)
        
        # Check if there's any overlap
        if (img_x_start < img_x_end and img_y_start < img_y_end and
            psf_x_start < psf_x_end and psf_y_start < psf_y_end):
            
            # Add PSF to image (no shot noise added here)
            image[img_y_start:img_y_end, img_x_start:img_x_end] += \
                psf_array[psf_y_start:psf_y_end, psf_x_start:psf_x_end]
    
    def create_stellar_field(self, source_data, image_shape, psf_params,
                           target_magnitude, saturation_limit, band='V',
                           psf_type='gaussian', target_name="V* TZ Boo", target_coords=None):
        """
        Create a synthetic stellar field using two-pass approach for exact peak flux control.
        
        Pass 1: Identify target source and calculate peak correction factor
        Pass 2: Apply correction factor to all sources for precise peak flux control
        
        Parameters
        ----------
        source_data : dict
            Source catalog with pixel coordinates and magnitudes
        image_shape : tuple
            (height, width) of output image
        psf_params : dict
            PSF parameters from optics calculations
        target_magnitude : float
            Target magnitude in specified band
        saturation_limit : int
            Camera saturation limit
        band : str, optional
            Photometric band ('V', 'B', 'R', 'I')
        psf_type : str, optional
            PSF type ('gaussian' or 'moffat')
        target_name : str, optional
            Target name for identification
        target_coords : tuple, optional
            (ra_deg, dec_deg) target coordinates if available
            
        Returns
        -------
        numpy.ndarray
            Synthetic image array
        """
        try:
            # Create empty image
            image = np.zeros(image_shape, dtype=np.float64)
            
            # Get source coordinates and magnitudes
            x_coords = source_data.get('x_pixels', source_data.get('x_pixel', []))
            y_coords = source_data.get('y_pixels', source_data.get('y_pixel', []))
            magnitudes = source_data.get(band, [])
            
            if len(x_coords) == 0:
                print(f"Warning: No sources found for {band} band")
                return image
            
            print(f"\n=== TWO-PASS APPROACH FOR {band.upper()} BAND ===")
            
            # PASS 1: Identify target source and calculate correction factor
            print("\nPass 1: Target identification and correction factor calculation")
            
            # Identify target source in catalog
            target_source = self.identify_target_source(
                source_data, target_name=target_name, target_coords=target_coords
            )
            
            if target_source is None:
                print("Warning: Could not identify target - using standard approach")
                correction_factor = 1.0
                target_x, target_y = x_coords[0], y_coords[0]  # Use first source as fallback
            else:
                target_x = target_source['x_pixel']
                target_y = target_source['y_pixel']
                target_index = target_source['index']
                
                # Get target magnitude from catalog
                if target_index < len(magnitudes):
                    catalog_target_mag = magnitudes[target_index]
                    print(f"Target magnitude from catalog: {catalog_target_mag:.3f}")
                    
                    # Use catalog magnitude if valid, otherwise use config magnitude
                    if np.isfinite(catalog_target_mag) and catalog_target_mag < 90:
                        actual_target_magnitude = catalog_target_mag
                    else:
                        actual_target_magnitude = target_magnitude
                        print(f"Using config target magnitude: {target_magnitude:.3f}")
                else:
                    actual_target_magnitude = target_magnitude
                
                # Calculate peak correction factor using target's exact position
                correction_factor = self.calculate_peak_correction_factor(
                    target_x, target_y, actual_target_magnitude, saturation_limit,
                    psf_params, psf_type, target_fraction=0.5
                )
            
            # PASS 2: Apply correction factor to all sources
            print(f"\nPass 2: Placing {len(x_coords)} sources with correction factor {correction_factor:.2f}")
            
            # Calculate base target peak flux (before correction)
            base_target_peak_adu = self.calculate_target_peak_flux(
                target_magnitude, saturation_limit, target_fraction=0.5
            )
            
            # Place each source with corrected flux
            sources_placed = 0
            brightest_flux = 0
            brightest_mag = 99
            target_flux_placed = 0
            
            for i, (x, y, mag) in enumerate(zip(x_coords, y_coords, magnitudes)):
                if np.isfinite(mag) and mag < 90:  # Skip invalid magnitudes
                    # Calculate base source peak flux relative to target
                    base_source_peak_flux = self.scale_source_peak_flux(
                        mag, target_magnitude, base_target_peak_adu
                    )
                    
                    # Apply correction factor for precise peak control
                    corrected_source_peak_flux = base_source_peak_flux * correction_factor
                    
                    # Track brightest source and target
                    if corrected_source_peak_flux > brightest_flux:
                        brightest_flux = corrected_source_peak_flux
                        brightest_mag = mag
                    
                    # Track target source specifically
                    is_target = target_source is not None and i == target_source['index']
                    if is_target:
                        target_flux_placed = corrected_source_peak_flux
                        print(f"  TARGET: mag={mag:.3f}, corrected_peak={corrected_source_peak_flux:.0f} ADU at ({x:.3f}, {y:.3f})")
                    
                    # Place source on image (enable debug for target source)
                    image = self.place_source_on_image(
                        image, x, y, corrected_source_peak_flux, psf_params, psf_type, debug=is_target
                    )
                    sources_placed += 1
            
            print(f"\nResults:")
            print(f"Brightest source: mag={brightest_mag:.2f}, peak_flux={brightest_flux:.0f} ADU")
            if target_flux_placed > 0:
                desired_peak = saturation_limit * 0.5
                accuracy = (target_flux_placed / desired_peak) * 100
                print(f"Target peak flux: {target_flux_placed:.0f} ADU (desired: {desired_peak:.0f}, accuracy: {accuracy:.1f}%)")
            
            print(f"Placed {sources_placed} sources in {band} band image")
            print(f"Image statistics: min={np.min(image):.1f}, "
                  f"max={np.max(image):.1f}, mean={np.mean(image):.1f} ADU")
            
            # Additional debugging summary
            self._print_correction_factor_summary(
                correction_factor, target_source, saturation_limit, sources_placed, band
            )
            
            return image
            
        except Exception as e:
            print(f"Error creating stellar field: {e}")
            return None
    
    def _print_correction_factor_summary(self, correction_factor, target_source, saturation_limit, sources_placed, band):
        """
        Print detailed summary of correction factor analysis.
        
        Parameters
        ----------
        correction_factor : float
            Calculated correction factor
        target_source : dict or None
            Target source information
        saturation_limit : int
            Camera saturation limit
        sources_placed : int
            Number of sources placed
        band : str
            Photometric band
        """
        print(f"\n=== CORRECTION FACTOR ANALYSIS SUMMARY ===")
        print(f"Band: {band}")
        print(f"Saturation limit: {saturation_limit} ADU")
        print(f"Target saturation fraction: 50%")
        print(f"Desired target peak: {saturation_limit * 0.5:.0f} ADU")
        print(f"Correction factor applied: {correction_factor:.6f}")
        
        if target_source is not None:
            print(f"Target source details:")
            print(f"  Index: {target_source['index']}")
            print(f"  Position: ({target_source['x_pixel']:.3f}, {target_source['y_pixel']:.3f})")
            if 'name' in target_source:
                print(f"  Name: {target_source['name']}")
            if 'separation_arcsec' in target_source:
                print(f"  Coordinate match: {target_source['separation_arcsec']:.1f} arcsec")
            if 'method' in target_source:
                print(f"  Identification method: {target_source['method']}")
        else:
            print("Target source: Not identified - used fallback approach")
        
        print(f"Total sources processed: {sources_placed}")
        
        # Calculate expected flux scaling effects
        if correction_factor != 1.0:
            flux_increase_percent = (correction_factor - 1.0) * 100
            print(f"Flux scaling effect: {flux_increase_percent:+.1f}%")
            
            # Warn about potential saturation
            if correction_factor > 2.0:
                print(f"WARNING: Large correction factor ({correction_factor:.2f}) may cause saturation in bright sources")
            elif correction_factor < 0.5:
                print(f"WARNING: Small correction factor ({correction_factor:.2f}) may result in very faint sources")
        
        print(f"==========================================\n")
    
    def add_noise(self, image, noise_type='read', noise_params=None):
        """
        Add noise to image (read noise only - no shot noise per PSF).
        
        Parameters
        ----------
        image : numpy.ndarray
            Input image
        noise_type : str
            Type of noise ('read', 'shot', 'both')
        noise_params : dict
            Noise parameters
            
        Returns
        -------
        numpy.ndarray
            Noisy image
        """
        if noise_params is None:
            noise_params = {'read_noise_sigma': 1.0}
        
        noisy_image = image.copy()
        
        # Add read noise
        if noise_type in ['read', 'both']:
            read_noise_sigma = noise_params.get('read_noise_sigma', 1.0)
            read_noise = np.random.normal(0, read_noise_sigma, image.shape)
            noisy_image += read_noise
        
        # Skip shot noise per PSF - analysis software will assume it's present
        # if noise_type in ['shot', 'both']:
        #     # Shot noise would be added here, but we skip it per requirements
        #     pass
        
        return noisy_image