"""
Configuration Manager

Handles loading, saving, and validation of JSON configuration files
for the synthetic image generation project.
"""

import json
import os
from pathlib import Path
from typing import Dict, Any, Optional


class ConfigManager:
    """Manages configuration files for the synthetic image generation project."""
    
    def __init__(self, config_dir: str = "config"):
        """
        Initialize the configuration manager.
        
        Args:
            config_dir (str): Directory containing configuration files
        """
        self.config_dir = Path(config_dir)
        self.config_dir.mkdir(parents=True, exist_ok=True)
        
        # Define configuration file paths
        self.target_file = self.config_dir / "target.json"
        self.telescope_file = self.config_dir / "telescope.json"
        self.camera_file = self.config_dir / "camera.json"
        self.fits_header_file = self.config_dir / "fits_header.json"
    
    def load_target_config(self) -> Dict[str, Any]:
        """Load target configuration from JSON file."""
        return self._load_json_file(self.target_file, self._get_default_target_config())
    
    def load_telescope_config(self) -> Dict[str, Any]:
        """Load telescope configuration from JSON file."""
        return self._load_json_file(self.telescope_file, self._get_default_telescope_config())
    
    def load_camera_config(self) -> Dict[str, Any]:
        """Load camera configuration from JSON file."""
        return self._load_json_file(self.camera_file, self._get_default_camera_config())
    
    def load_fits_header_config(self) -> Dict[str, Any]:
        """Load FITS header configuration from JSON file."""
        return self._load_json_file(self.fits_header_file, {})
    
    def save_target_config(self, config: Dict[str, Any]) -> bool:
        """Save target configuration to JSON file."""
        if self._validate_target_config(config):
            return self._save_json_file(self.target_file, config)
        return False
    
    def save_telescope_config(self, config: Dict[str, Any]) -> bool:
        """Save telescope configuration to JSON file."""
        if self._validate_telescope_config(config):
            return self._save_json_file(self.telescope_file, config)
        return False
    
    def save_camera_config(self, config: Dict[str, Any]) -> bool:
        """Save camera configuration to JSON file."""
        if self._validate_camera_config(config):
            return self._save_json_file(self.camera_file, config)
        return False
    
    def save_fits_header_config(self, config: Dict[str, Any]) -> bool:
        """Save FITS header configuration to JSON file."""
        return self._save_json_file(self.fits_header_file, config)
    
    def _load_json_file(self, file_path: Path, default_config: Dict[str, Any]) -> Dict[str, Any]:
        """Load JSON configuration file with fallback to default."""
        try:
            if file_path.exists():
                with open(file_path, 'r') as f:
                    return json.load(f)
            else:
                # Create default config file if it doesn't exist
                self._save_json_file(file_path, default_config)
                return default_config
        except (json.JSONDecodeError, IOError) as e:
            print(f"Error loading {file_path}: {e}")
            return default_config
    
    def _save_json_file(self, file_path: Path, config: Dict[str, Any]) -> bool:
        """Save configuration to JSON file."""
        try:
            with open(file_path, 'w') as f:
                json.dump(config, f, indent=2, default=str)
            return True
        except IOError as e:
            print(f"Error saving {file_path}: {e}")
            return False
    
    def _get_default_target_config(self) -> Dict[str, Any]:
        """Get default target configuration."""
        return {
            "target_name": "",
            "coordinates": {
                "ra": 0.0,
                "dec": 0.0,
                "ra_str": "",
                "dec_str": ""
            },
            "magnitudes": {
                "V": 0.0,
                "B": 0.0,
                "R": 0.0,
                "I": 0.0
            },
            "notes": ""
        }
    
    def _get_default_telescope_config(self) -> Dict[str, Any]:
        """Get default telescope configuration."""
        return {
            "telescope_name": "",
            "aperture_mm": 0.0,
            "focal_length_mm": 0.0,
            "focal_ratio": 0.0,
            "seeing": {
                "theoretical_arcsec": 0.0,
                "atmospheric_arcsec": 1.5,
                "combined_fwhm_arcsec": 0.0
            },
            "notes": ""
        }
    
    def _get_default_camera_config(self) -> Dict[str, Any]:
        """Get default camera configuration."""
        return {
            "camera_name": "",
            "pixel_size_microns": 0.0,
            "sensor_dimensions": {
                "width_pixels": 0,
                "height_pixels": 0,
                "width_mm": 0.0,
                "height_mm": 0.0
            },
            "rotation_degrees": 0.0,
            "saturation_limit": 65535,
            "notes": ""
        }
    
    def _validate_target_config(self, config: Dict[str, Any]) -> bool:
        """Validate target configuration structure."""
        try:
            required_keys = ["target_name", "coordinates", "magnitudes"]
            for key in required_keys:
                if key not in config:
                    print(f"Missing required key in target config: {key}")
                    return False
            
            # Validate coordinates
            coord_keys = ["ra", "dec", "ra_str", "dec_str"]
            for key in coord_keys:
                if key not in config["coordinates"]:
                    print(f"Missing coordinate key: {key}")
                    return False
            
            # Validate magnitudes
            mag_keys = ["V", "B", "R", "I"]
            for key in mag_keys:
                if key not in config["magnitudes"]:
                    print(f"Missing magnitude key: {key}")
                    return False
            
            return True
        except Exception as e:
            print(f"Error validating target config: {e}")
            return False
    
    def _validate_telescope_config(self, config: Dict[str, Any]) -> bool:
        """Validate telescope configuration structure."""
        try:
            required_keys = ["telescope_name", "aperture_mm", "focal_length_mm", "seeing"]
            for key in required_keys:
                if key not in config:
                    print(f"Missing required key in telescope config: {key}")
                    return False
            
            # Validate seeing
            seeing_keys = ["theoretical_arcsec", "atmospheric_arcsec", "combined_fwhm_arcsec"]
            for key in seeing_keys:
                if key not in config["seeing"]:
                    print(f"Missing seeing key: {key}")
                    return False
            
            return True
        except Exception as e:
            print(f"Error validating telescope config: {e}")
            return False
    
    def _validate_camera_config(self, config: Dict[str, Any]) -> bool:
        """Validate camera configuration structure."""
        try:
            required_keys = ["camera_name", "pixel_size_microns", "sensor_dimensions", "saturation_limit"]
            for key in required_keys:
                if key not in config:
                    print(f"Missing required key in camera config: {key}")
                    return False
            
            # Validate sensor dimensions
            sensor_keys = ["width_pixels", "height_pixels", "width_mm", "height_mm"]
            for key in sensor_keys:
                if key not in config["sensor_dimensions"]:
                    print(f"Missing sensor dimension key: {key}")
                    return False
            
            return True
        except Exception as e:
            print(f"Error validating camera config: {e}")
            return False