import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astropy import units as u
from gaiaxpy import generate, PhotometricSystem
import warnings
from typing import List, Optional, Dict, Tuple


class GaiaPhotometry:
    """
    GAIA-based photometry using high-resolution BP/RP spectra.
    
    This class provides more accurate Johnson-Kron-Cousins photometry
    compared to polynomial transformations by using GAIA's spectroscopic data.
    """
    
    def __init__(self):
        self.photometric_system = PhotometricSystem.JKC_Std
        
    def get_source_ids_from_coordinates(self, ra: float, dec: float, 
                                      width: float = 1.0, height: float = 1.0) -> List[int]:
        """
        Get GAIA source IDs from coordinates.
        
        Parameters
        ----------
        ra : float
            Right ascension in degrees
        dec : float
            Declination in degrees
        width : float, optional
            Search width in arcseconds (default: 1.0)
        height : float, optional
            Search height in arcseconds (default: 1.0)
            
        Returns
        -------
        List[int]
            List of GAIA source IDs
        """
        coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')
        width_qty = u.Quantity(width, u.arcsec)
        height_qty = u.Quantity(height, u.arcsec)
        
        try:
            result = Gaia.query_object_async(
                coordinate=coord, 
                width=width_qty, 
                height=height_qty
            )
            return result['source_id'].tolist()
        except Exception as e:
            warnings.warn(f"Failed to query GAIA coordinates: {e}")
            return []
    
    def generate_synthetic_photometry(self, source_ids: List[int]) -> Optional[pd.DataFrame]:
        """
        Generate synthetic photometry for multiple source IDs.
        
        Parameters
        ----------
        source_ids : List[int]
            List of GAIA source IDs
            
        Returns
        -------
        pd.DataFrame or None
            DataFrame with synthetic photometry results, or None if failed
        """
        if not source_ids:
            return None
            
        try:
            # Generate synthetic photometry using gaiaxpy
            synthetic_photometry = generate(
                source_ids, 
                photometric_system=self.photometric_system
            )
            
            # Check if any sources have missing BP/RP spectra
            if synthetic_photometry is not None and len(synthetic_photometry) > 0:
                missing_spectra = len(source_ids) - len(synthetic_photometry)
                if missing_spectra > 0:
                    warnings.warn(f"{missing_spectra} sources lack BP/RP spectra for synthetic photometry")
                    
            return synthetic_photometry
        except Exception as e:
            # Common error cases: sources without BP/RP spectra, network issues, etc.
            if "No valid sources" in str(e) or "No BP/RP spectra" in str(e):
                warnings.warn(f"No sources have BP/RP spectra available for synthetic photometry")
            else:
                warnings.warn(f"Failed to generate synthetic photometry: {e}")
            return None
    
    def extract_johnson_cousins_magnitudes(self, photometry_df: pd.DataFrame) -> Dict[int, Dict[str, float]]:
        """
        Extract Johnson-Kron-Cousins magnitudes from gaiaxpy results.
        
        Parameters
        ----------
        photometry_df : pd.DataFrame
            DataFrame from gaiaxpy generate() function
            
        Returns
        -------
        Dict[int, Dict[str, float]]
            Dictionary mapping source_id to magnitude dictionary
        """
        result = {}
        
        for _, row in photometry_df.iterrows():
            source_id = int(row['source_id'])
            magnitudes = {
                'U': row.get('JkcStd_mag_U', np.nan),
                'B': row.get('JkcStd_mag_B', np.nan),
                'V': row.get('JkcStd_mag_V', np.nan),
                'R': row.get('JkcStd_mag_R', np.nan),
                'I': row.get('JkcStd_mag_I', np.nan)
            }
            result[source_id] = magnitudes
            
        return result
    
    def get_magnitude_errors(self, photometry_df: pd.DataFrame) -> Dict[int, Dict[str, float]]:
        """
        Calculate magnitude errors from flux and flux_error columns.
        
        Parameters
        ----------
        photometry_df : pd.DataFrame
            DataFrame from gaiaxpy generate() function
            
        Returns
        -------
        Dict[int, Dict[str, float]]
            Dictionary mapping source_id to magnitude error dictionary
        """
        result = {}
        
        for _, row in photometry_df.iterrows():
            source_id = int(row['source_id'])
            mag_errors = {}
            
            for band in ['U', 'B', 'V', 'R', 'I']:
                flux_col = f'JkcStd_flux_{band}'
                flux_error_col = f'JkcStd_flux_error_{band}'
                
                if flux_col in row and flux_error_col in row:
                    flux = row[flux_col]
                    flux_error = row[flux_error_col]
                    
                    if flux > 0 and flux_error > 0:
                        snr = flux / flux_error
                        mag_error = 2.5 / np.log(10) / snr
                        mag_errors[band] = mag_error
                    else:
                        mag_errors[band] = np.nan
                else:
                    mag_errors[band] = np.nan
                    
            result[source_id] = mag_errors
            
        return result
    
    def process_source_list(self, source_ids: List[int]) -> Tuple[Dict[int, Dict[str, float]], 
                                                                Dict[int, Dict[str, float]]]:
        """
        Process a list of source IDs to get magnitudes and errors.
        
        Parameters
        ----------
        source_ids : List[int]
            List of GAIA source IDs
            
        Returns
        -------
        Tuple[Dict, Dict]
            Tuple of (magnitudes_dict, errors_dict)
        """
        photometry_df = self.generate_synthetic_photometry(source_ids)
        
        if photometry_df is None:
            return {}, {}
            
        magnitudes = self.extract_johnson_cousins_magnitudes(photometry_df)
        errors = self.get_magnitude_errors(photometry_df)
        
        return magnitudes, errors