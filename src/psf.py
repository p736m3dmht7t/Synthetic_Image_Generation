"""
Point Spread Function (PSF) Generation

Handles the creation and application of point spread functions for
synthetic star image generation.
"""

import numpy as np
import math
from scipy import ndimage
from scipy.stats import multivariate_normal


class PSFGenerator:
    """Generates and applies point spread functions for stellar sources."""
    
    def __init__(self):
        """Initialize the PSF generator."""
        pass
    
    def generate_gaussian_psf(self, sigma_pixels, psf_size=None):
        """
        Generate a 2D Gaussian PSF.
        
        Args:
            sigma_pixels (float): Gaussian sigma in pixels
            psf_size (int): Size of PSF array (default: 6*sigma, minimum 7)
        
        Returns:
            numpy.ndarray: 2D PSF array, normalized to sum = 1
        """
        try:
            # Determine PSF size
            if psf_size is None:
                psf_size = max(7, int(6 * sigma_pixels))
                # Ensure odd size for proper centering
                if psf_size % 2 == 0:
                    psf_size += 1
            
            # Create coordinate arrays
            center = psf_size // 2
            x = np.arange(psf_size) - center
            y = np.arange(psf_size) - center
            X, Y = np.meshgrid(x, y)
            
            # Generate 2D Gaussian
            psf = np.exp(-(X**2 + Y**2) / (2 * sigma_pixels**2))
            
            # Normalize to sum = 1
            psf = psf / np.sum(psf)
            
            return psf
            
        except Exception as e:
            print(f"Error generating Gaussian PSF: {str(e)}")
            return None
    
    def generate_moffat_psf(self, alpha, beta, psf_size=None):
        """
        Generate a Moffat PSF profile (often better for atmospheric seeing).
        
        Args:
            alpha (float): Scale parameter in pixels
            beta (float): Power law index (typically 2-5 for atmospheric seeing)
            psf_size (int): Size of PSF array
        
        Returns:
            numpy.ndarray: 2D Moffat PSF array, normalized to sum = 1
        """
        try:
            # Determine PSF size
            if psf_size is None:
                psf_size = max(7, int(8 * alpha))
                if psf_size % 2 == 0:
                    psf_size += 1
            
            # Create coordinate arrays
            center = psf_size // 2
            x = np.arange(psf_size) - center
            y = np.arange(psf_size) - center
            X, Y = np.meshgrid(x, y)
            
            # Calculate radius
            r = np.sqrt(X**2 + Y**2)
            
            # Generate Moffat profile
            psf = (1 + (r / alpha)**2)**(-beta)
            
            # Normalize to sum = 1
            psf = psf / np.sum(psf)
            
            return psf
            
        except Exception as e:
            print(f"Error generating Moffat PSF: {str(e)}")
            return None
    
    def fwhm_to_sigma(self, fwhm_pixels):
        """
        Convert FWHM to Gaussian sigma.
        
        Args:
            fwhm_pixels (float): FWHM in pixels
        
        Returns:
            float: Gaussian sigma in pixels
        """
        return fwhm_pixels / (2 * math.sqrt(2 * math.log(2)))
    
    def sigma_to_fwhm(self, sigma_pixels):
        """
        Convert Gaussian sigma to FWHM.
        
        Args:
            sigma_pixels (float): Gaussian sigma in pixels
        
        Returns:
            float: FWHM in pixels
        """
        return sigma_pixels * (2 * math.sqrt(2 * math.log(2)))
    
    def magnitude_to_flux(self, magnitude, zero_point=25.0):
        """
        Convert magnitude to relative flux.
        
        Args:
            magnitude (float or array): Magnitude
            zero_point (float): Zero point magnitude
        
        Returns:
            float or array: Relative flux
        """
        return 10**((zero_point - magnitude) / 2.5)
    
    def calculate_target_adu(self, target_magnitude, saturation_limit, target_fraction=0.5):
        """
        Calculate target ADU value for the main target.
        Sets target star to specified fraction of saturation regardless of magnitude.
        
        Args:
            target_magnitude (float): Target magnitude (for reference, not used in calculation)
            saturation_limit (int): Camera saturation limit
            target_fraction (float): Fraction of saturation for target (default 0.5)
        
        Returns:
            float: Target ADU value
        """
        target_adu = saturation_limit * target_fraction
        print(f"Target magnitude {target_magnitude:.2f} -> {target_adu:.0f} ADU ({target_fraction*100:.0f}% of saturation)")
        return target_adu
    
    def place_psf_on_image(self, image, psf, x_center, y_center, peak_adu):
        """
        Place a PSF at specified coordinates on an image.
        
        Args:
            image (numpy.ndarray): Target image array
            psf (numpy.ndarray): PSF array (normalized to sum=1)
            x_center (float): X coordinate of PSF center
            y_center (float): Y coordinate of PSF center
            peak_adu (float): Desired peak pixel value in ADU
        
        Returns:
            numpy.ndarray: Updated image array
        """
        try:
            # Get image and PSF dimensions
            img_height, img_width = image.shape
            psf_height, psf_width = psf.shape
            
            # Calculate PSF placement bounds
            psf_center_x = psf_width // 2
            psf_center_y = psf_height // 2
            
            # Calculate image coordinates for PSF placement
            img_x_start = int(x_center - psf_center_x)
            img_y_start = int(y_center - psf_center_y)
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
                
                # Calculate scaling factor to achieve desired peak pixel value
                psf_peak = np.max(psf)
                if psf_peak > 0:
                    scale_factor = peak_adu / psf_peak
                else:
                    scale_factor = 0
                
                # Add scaled PSF to image
                scaled_psf = scale_factor * psf[psf_y_start:psf_y_end, psf_x_start:psf_x_end]
                image[img_y_start:img_y_end, img_x_start:img_x_end] += scaled_psf
            
            return image
            
        except Exception as e:
            print(f"Error placing PSF on image: {str(e)}")
            return image
    
    def add_noise(self, image, noise_type='poisson', noise_params=None):
        """
        Add noise to the image.
        
        Args:
            image (numpy.ndarray): Input image
            noise_type (str): Type of noise ('poisson', 'gaussian', 'both')
            noise_params (dict): Noise parameters
        
        Returns:
            numpy.ndarray: Noisy image
        """
        try:
            noisy_image = image.copy()
            
            if noise_params is None:
                noise_params = {}
            
            if noise_type in ['poisson', 'both']:
                # Poisson noise (shot noise)
                # Ensure non-negative values for Poisson
                positive_image = np.maximum(noisy_image, 0)
                poisson_noise = np.random.poisson(positive_image) - positive_image
                noisy_image += poisson_noise
            
            if noise_type in ['gaussian', 'both']:
                # Gaussian read noise
                read_noise_sigma = noise_params.get('read_noise_sigma', 1.0)
                gaussian_noise = np.random.normal(0, read_noise_sigma, image.shape)
                noisy_image += gaussian_noise
            
            return noisy_image
            
        except Exception as e:
            print(f"Error adding noise: {str(e)}")
            return image
    
    def apply_saturation(self, image, saturation_limit):
        """
        Apply saturation limit to image.
        
        Args:
            image (numpy.ndarray): Input image
            saturation_limit (float): Maximum ADU value
        
        Returns:
            numpy.ndarray: Saturated image
        """
        return np.clip(image, 0, saturation_limit)
    
    def create_stellar_field(self, source_data, image_shape, psf_params, 
                           target_magnitude, saturation_limit, band='V'):
        """
        Create a synthetic stellar field image.
        
        Args:
            source_data (dict): Source catalog with pixel coordinates and magnitudes
            image_shape (tuple): (height, width) of output image
            psf_params (dict): PSF parameters from optics calculations
            target_magnitude (float): Target magnitude in specified band
            saturation_limit (int): Camera saturation limit
            band (str): Photometric band ('V', 'B', 'R', 'I')
        
        Returns:
            numpy.ndarray: Synthetic image array
        """
        try:
            # Create empty image
            image = np.zeros(image_shape, dtype=np.float64)
            
            # Get PSF parameters
            sigma_pixels = psf_params.get('sigma_pixels', 1.0)
            
            # Generate PSF
            psf = self.generate_gaussian_psf(sigma_pixels)
            if psf is None:
                print("Error: Could not generate PSF")
                return None
            
            # Calculate target flux (50% of saturation)
            target_adu = self.calculate_target_adu(target_magnitude, saturation_limit)
            
            # First, place the target star at image center with correct brightness
            center_x = image_shape[1] // 2
            center_y = image_shape[0] // 2
            print(f"Placing target star (mag {target_magnitude:.2f}) at center ({center_x}, {center_y}) with ADU {target_adu:.0f}")
            image = self.place_psf_on_image(image, psf, center_x, center_y, target_adu)
            
            # Get source coordinates and magnitudes
            x_pixels = source_data.get('x_pixels', [])
            y_pixels = source_data.get('y_pixels', [])
            magnitudes = source_data.get(band, [])
            
            if len(x_pixels) == 0:
                print(f"Warning: No sources found for band {band}")
                return image
            
            print(f"Placing {len(x_pixels)} catalog sources in {band} band image")
            
            # Calculate peak ADU for catalog sources relative to target
            # Brighter sources (lower magnitude) get higher ADU values
            magnitudes = np.array(magnitudes)
            magnitude_diff = magnitudes - target_magnitude  # Positive for dimmer sources
            relative_brightness = 10**(-magnitude_diff / 2.5)  # Standard magnitude scale
            source_peak_adu = target_adu * relative_brightness
            
            # Place each source
            for i, (x, y, peak_adu) in enumerate(zip(x_pixels, y_pixels, source_peak_adu)):
                # Apply saturation limit to individual source
                peak_adu_limited = min(peak_adu, saturation_limit)
                
                # Place PSF
                image = self.place_psf_on_image(image, psf, x, y, peak_adu_limited)
                
                # Progress indicator for large catalogs
                if (i + 1) % 100 == 0:
                    print(f"Placed {i + 1}/{len(x_pixels)} sources")
            
            # Apply final saturation
            image = self.apply_saturation(image, saturation_limit)
            
            print(f"Created {band} band image with {len(x_pixels)} sources")
            print(f"Image statistics: min={np.min(image):.1f}, max={np.max(image):.1f}, mean={np.mean(image):.1f}")
            
            return image
            
        except Exception as e:
            print(f"Error creating stellar field: {str(e)}")
            return None
    
    def get_psf_info(self, psf_params):
        """
        Get PSF information summary.
        
        Args:
            psf_params (dict): PSF parameters
        
        Returns:
            dict: PSF information
        """
        try:
            fwhm_pixels = psf_params.get('fwhm_pixels', 0)
            fwhm_arcsec = psf_params.get('fwhm_arcsec', 0)
            sigma_pixels = psf_params.get('sigma_pixels', 0)
            
            info = {
                'fwhm_pixels': fwhm_pixels,
                'fwhm_arcsec': fwhm_arcsec,
                'sigma_pixels': sigma_pixels,
                'psf_size_pixels': max(7, int(6 * sigma_pixels)),
                'telescope_limited': psf_params.get('telescope_limited', False),
                'atmospheric_limited': psf_params.get('atmospheric_limited', True)
            }
            
            return info
            
        except Exception as e:
            print(f"Error getting PSF info: {str(e)}")
            return None