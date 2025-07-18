"""
Image Visualization System

Handles display and analysis of synthetic astronomical images
using matplotlib for interactive visualization.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from astropy.io import fits
from astropy.visualization import ZScaleInterval, ImageNormalize
from astropy.stats import sigma_clipped_stats
import tkinter as tk
from tkinter import ttk
import warnings

# Suppress matplotlib warnings
warnings.filterwarnings('ignore', category=UserWarning)


class ImageVisualizer:
    """Handles visualization of synthetic astronomical images."""
    
    def __init__(self):
        """Initialize the image visualizer."""
        self.images = {}
        self.current_band = None
        self.figure = None
        self.canvas = None
        self.source_data = None
        self.show_sources = False
        
        # Visualization parameters
        self.stretch = 'zscale'
        self.colormap = 'gray'
        self.contrast = 1.0
        self.brightness = 0.0
        self.invert_display = False  # Toggle for inverted display
    
    def load_images(self, images_dict):
        """
        Load images for visualization.
        
        Args:
            images_dict (dict): Dictionary of band -> image array
        """
        self.images = images_dict.copy()
        if self.images and self.current_band is None:
            self.current_band = list(self.images.keys())[0]
        
        print(f"Loaded {len(self.images)} images for visualization")
    
    def load_fits_images(self, fits_files):
        """
        Load images from FITS files.
        
        Args:
            fits_files (dict): Dictionary of band -> filename
        """
        try:
            self.images = {}
            
            for band, filename in fits_files.items():
                try:
                    with fits.open(filename) as hdul:
                        self.images[band] = hdul[0].data.astype(np.float64)
                    print(f"Loaded {band} band from {filename}")
                except Exception as e:
                    print(f"Error loading {filename}: {str(e)}")
            
            if self.images and self.current_band is None:
                self.current_band = list(self.images.keys())[0]
                
            return len(self.images) > 0
            
        except Exception as e:
            print(f"Error loading FITS images: {str(e)}")
            return False
    
    def set_source_data(self, source_data):
        """
        Set source data for overlay display.
        
        Args:
            source_data (dict): Source catalog with pixel coordinates
        """
        self.source_data = source_data
        print(f"Source data loaded: {len(source_data.get('x_pixels', []))} sources")
    
    def calculate_image_statistics(self, image):
        """
        Calculate image statistics for display scaling.
        
        Args:
            image (numpy.ndarray): Image array
        
        Returns:
            dict: Image statistics
        """
        try:
            # Use sigma clipping to get robust statistics
            mean, median, std = sigma_clipped_stats(image, sigma=3.0)
            
            stats = {
                'mean': mean,
                'median': median,
                'std': std,
                'min': np.min(image),
                'max': np.max(image),
                'shape': image.shape
            }
            
            return stats
            
        except Exception as e:
            print(f"Error calculating image statistics: {str(e)}")
            return {}
    
    def normalize_image(self, image, stretch='zscale'):
        """
        Normalize image for display.
        
        Args:
            image (numpy.ndarray): Image array
            stretch (str): Stretch algorithm ('zscale', 'linear', 'log', 'sqrt')
        
        Returns:
            numpy.ndarray: Normalized image
        """
        try:
            if stretch == 'zscale':
                # Use ZScale algorithm (common in astronomy)
                try:
                    interval = ZScaleInterval()
                    vmin, vmax = interval.get_limits(image)
                    
                    # Check if ZScale failed or returned unreasonable results
                    img_min, img_max = np.min(image), np.max(image)
                    if (vmin == vmax or not np.isfinite(vmin) or not np.isfinite(vmax) or 
                        vmax < img_max * 0.1):  # ZScale max is less than 10% of actual max
                        print(f"ZScale failed (returned [{vmin:.1f}, {vmax:.1f}] for image range [{img_min:.1f}, {img_max:.1f}]), using linear stretch")
                        vmin, vmax = np.percentile(image, [0.1, 99.9])
                        if vmin == vmax:
                            vmin, vmax = img_min, img_max
                    
                    normalized = np.clip((image - vmin) / (vmax - vmin), 0, 1)
                except Exception as e:
                    print(f"ZScale error ({e}), using linear stretch")
                    img_min, img_max = np.min(image), np.max(image)
                    vmin, vmax = np.percentile(image, [0.1, 99.9])
                    if vmin == vmax:
                        vmin, vmax = img_min, img_max
                    normalized = np.clip((image - vmin) / (vmax - vmin), 0, 1)
                
            elif stretch == 'linear':
                # Linear stretch - for astronomical images, use smart scaling
                img_min, img_max = np.min(image), np.max(image)
                
                # Check if this looks like an astronomical image (high dynamic range)
                dynamic_range = img_max / (img_min + 1e-10)  # Avoid division by zero
                
                if dynamic_range > 1000:  # High dynamic range suggests astronomical data
                    # Use a small percentile for vmin and full max for vmax
                    vmin = np.percentile(image, 0.1)  # Skip true zeros
                    vmax = img_max
                else:
                    # Normal percentile stretch for regular images
                    vmin, vmax = np.percentile(image, [1, 99])
                
                if vmin == vmax:  # Handle constant images
                    normalized = np.ones_like(image) * 0.5
                else:
                    normalized = np.clip((image - vmin) / (vmax - vmin), 0, 1)
                
            elif stretch == 'log':
                # Logarithmic stretch
                image_positive = image - np.min(image) + 1
                log_image = np.log10(image_positive)
                vmin, vmax = np.percentile(log_image, [1, 99])
                if vmin == vmax:
                    normalized = np.ones_like(image) * 0.5
                else:
                    normalized = np.clip((log_image - vmin) / (vmax - vmin), 0, 1)
                
            elif stretch == 'sqrt':
                # Square root stretch - for astronomical images, use smart scaling
                img_min, img_max = np.min(image), np.max(image)
                dynamic_range = img_max / (img_min + 1e-10)
                
                image_positive = image - img_min
                sqrt_image = np.sqrt(image_positive)
                
                if dynamic_range > 1000:  # Astronomical image
                    vmin = np.sqrt(np.percentile(image_positive, 0.1))
                    vmax = np.sqrt(img_max - img_min)
                else:
                    vmin, vmax = np.percentile(sqrt_image, [1, 99])
                
                if vmin == vmax:
                    normalized = np.ones_like(image) * 0.5
                else:
                    normalized = np.clip((sqrt_image - vmin) / (vmax - vmin), 0, 1)
                
            else:
                # Default to linear
                vmin, vmax = np.percentile(image, [1, 99])
                if vmin == vmax:
                    normalized = np.ones_like(image) * 0.5
                else:
                    normalized = np.clip((image - vmin) / (vmax - vmin), 0, 1)
            
            # Apply contrast and brightness adjustments
            normalized = np.clip(normalized * self.contrast + self.brightness, 0, 1)
            
            # Apply inversion if requested (for display only)
            if self.invert_display:
                normalized = 1.0 - normalized
            
            return normalized
            
        except Exception as e:
            print(f"Error normalizing image: {str(e)}")
            # Emergency fallback - simple linear normalization
            try:
                img_min, img_max = np.min(image), np.max(image)
                if img_min == img_max:
                    return np.ones_like(image) * 0.5
                else:
                    return np.clip((image - img_min) / (img_max - img_min), 0, 1)
            except:
                return image
    
    def create_matplotlib_figure(self, parent_frame, figsize=(8, 6)):
        """
        Create matplotlib figure embedded in tkinter frame.
        
        Args:
            parent_frame: Tkinter parent widget
            figsize (tuple): Figure size in inches
        
        Returns:
            tuple: (figure, canvas) objects
        """
        try:
            # Create figure
            self.figure = Figure(figsize=figsize, dpi=100)
            self.figure.patch.set_facecolor('white')
            
            # Create canvas
            self.canvas = FigureCanvasTkAgg(self.figure, parent_frame)
            self.canvas.get_tk_widget().pack(fill="both", expand=True)
            
            return self.figure, self.canvas
            
        except Exception as e:
            print(f"Error creating matplotlib figure: {str(e)}")
            return None, None
    
    def plot_image(self, band=None, show_sources=None):
        """
        Plot image with optional source overlay.
        
        Args:
            band (str): Photometric band to display
            show_sources (bool): Whether to show source overlay
        """
        try:
            
            if not self.images:
                print("No images loaded for visualization")
                return
            
            if band is None:
                band = self.current_band
            if band not in self.images:
                print(f"Band {band} not available")
                return
            
            if show_sources is not None:
                self.show_sources = show_sources
            
            # Just clear the figure - don't destroy/recreate canvas
            if self.figure:
                self.figure.clear()
            
            # Get image
            image = self.images[band]
            
            # Normalize image
            normalized_image = self.normalize_image(image, self.stretch)
            
            # Create subplot
            ax = self.figure.add_subplot(111)
            
            # Enhance display for sparse star fields
            display_image = normalized_image.copy()
            bright_pixels = normalized_image[normalized_image > 0.1]   # Values > 10% of max
            
            # For sparse astronomical images, enhance contrast dramatically
            if len(bright_pixels) < normalized_image.size * 0.001:  # Less than 0.1% bright pixels
                # More aggressive enhancement for tiny stars
                gamma = 0.2  # Aggressive gamma for visibility
                display_image = np.power(display_image, gamma)
                
                # Scale the non-zero values more aggressively
                nonzero_mask = display_image > 0
                if np.any(nonzero_mask):
                    # Get the 95th percentile of non-zero values
                    top_5_percent = np.percentile(display_image[nonzero_mask], 95)
                    if top_5_percent > 0:
                        display_image = np.clip(display_image / (top_5_percent * 0.5), 0, 1)
            
            im = ax.imshow(display_image, cmap=self.colormap, origin='lower', aspect='auto', interpolation='nearest')
            
            # Add source overlay
            if self.show_sources and self.source_data:
                self.add_source_overlay(ax, band)
            
            # Set title and labels
            stats = self.calculate_image_statistics(image)
            title = f"{band} Band - Shape: {stats.get('shape', 'Unknown')}"
            title += f" - Mean: {stats.get('mean', 0):.1f} ± {stats.get('std', 0):.1f}"
            ax.set_title(title, fontsize=10)
            ax.set_xlabel('X (pixels)')
            ax.set_ylabel('Y (pixels)')
            
            # Add colorbar
            try:
                cbar = self.figure.colorbar(im, ax=ax, shrink=0.8)
                cbar.set_label('ADU', rotation=270, labelpad=15)
            except:
                pass  # Skip colorbar if it fails
            
            # Tight layout
            self.figure.tight_layout()
            
            # Update canvas
            if self.canvas:
                self.canvas.draw()
            else:
                print("WARNING: No canvas to update!")
            
            self.current_band = band
            
        except Exception as e:
            print(f"Error plotting image: {str(e)}")
    
    def add_source_overlay(self, ax, band):
        """
        Add source overlay to the plot.
        
        Args:
            ax: Matplotlib axes object
            band (str): Current photometric band
        """
        try:
            if not self.source_data:
                return
            
            # Get source positions
            x_pixels = self.source_data.get('x_pixels', [])
            y_pixels = self.source_data.get('y_pixels', [])
            magnitudes = self.source_data.get(band, [])
            
            
            if len(x_pixels) == 0:
                return
            
            # Handle case where magnitudes may not match positions
            if len(magnitudes) > 0 and len(magnitudes) == len(x_pixels):
                # Scale marker size inversely with magnitude (brighter = larger)
                max_mag = np.max(magnitudes)
                min_mag = np.min(magnitudes)
                mag_range = max_mag - min_mag
                
                if mag_range > 0:
                    # Normalize magnitudes to marker sizes (10-100 pixels)
                    normalized_mags = (max_mag - magnitudes) / mag_range
                    marker_sizes = 10 + normalized_mags * 90
                else:
                    marker_sizes = [20] * len(magnitudes)
                
                # Plot sources as circles
                ax.scatter(x_pixels, y_pixels, s=marker_sizes, 
                          facecolors='none', edgecolors='red', 
                          alpha=0.7, linewidths=1)
            else:
                # Fallback: show all positions with uniform size
                marker_sizes = [20] * len(x_pixels)
                ax.scatter(x_pixels, y_pixels, s=marker_sizes, 
                          facecolors='none', edgecolors='red', 
                          alpha=0.7, linewidths=1)
            
            # Add text for brightest sources
            if len(magnitudes) > 0 and len(magnitudes) == len(x_pixels):
                # Show labels for brightest 10 sources
                bright_indices = np.argsort(magnitudes)[:min(10, len(magnitudes))]
                
                for idx in bright_indices:
                    ax.annotate(f'{magnitudes[idx]:.1f}', 
                              (x_pixels[idx], y_pixels[idx]),
                              xytext=(5, 5), textcoords='offset points',
                              fontsize=8, color='yellow',
                              bbox=dict(boxstyle='round,pad=0.2', 
                                       facecolor='black', alpha=0.7))
            
            print(f"Added overlay for {len(x_pixels)} sources")
            
        except Exception as e:
            print(f"Error adding source overlay: {str(e)}")
    
    def create_comparison_plot(self):
        """Create side-by-side comparison of all bands."""
        try:
            if not self.images:
                print("No images loaded for comparison")
                return
            
            bands = list(self.images.keys())
            n_bands = len(bands)
            
            if n_bands == 0:
                return
            
            # Clear previous plot
            if self.figure:
                self.figure.clear()
            
            # Create subplots
            if n_bands == 1:
                subplot_layout = (1, 1)
            elif n_bands == 2:
                subplot_layout = (1, 2)
            elif n_bands <= 4:
                subplot_layout = (2, 2)
            else:
                subplot_layout = (2, 3)  # Up to 6 bands
            
            for i, band in enumerate(bands[:6]):  # Limit to 6 bands max
                ax = self.figure.add_subplot(subplot_layout[0], subplot_layout[1], i+1)
                
                # Get and normalize image
                image = self.images[band]
                normalized_image = self.normalize_image(image, self.stretch)
                
                # Display image
                im = ax.imshow(normalized_image, cmap=self.colormap, origin='lower')
                
                # Add source overlay if enabled
                if self.show_sources and self.source_data:
                    self.add_source_overlay(ax, band)
                
                # Set title
                stats = self.calculate_image_statistics(image)
                ax.set_title(f"{band} Band (μ={stats.get('mean', 0):.0f})", fontsize=10)
                ax.set_xlabel('X (pixels)', fontsize=8)
                ax.set_ylabel('Y (pixels)', fontsize=8)
                
                # Smaller tick labels
                ax.tick_params(labelsize=8)
            
            # Tight layout
            self.figure.tight_layout()
            
            # Update canvas
            if self.canvas:
                self.canvas.draw()
            
        except Exception as e:
            print(f"Error creating comparison plot: {str(e)}")
    
    def save_plot(self, filename, dpi=150):
        """
        Save current plot to file.
        
        Args:
            filename (str): Output filename
            dpi (int): Resolution in DPI
        """
        try:
            if self.figure:
                self.figure.savefig(filename, dpi=dpi, bbox_inches='tight')
                print(f"Plot saved to {filename}")
                return True
            else:
                print("No figure to save")
                return False
                
        except Exception as e:
            print(f"Error saving plot: {str(e)}")
            return False
    
    def update_display_settings(self, stretch=None, colormap=None, 
                              contrast=None, brightness=None, invert_display=None):
        """
        Update display settings and refresh plot.
        
        Args:
            stretch (str): Stretch algorithm
            colormap (str): Colormap name
            contrast (float): Contrast adjustment
            brightness (float): Brightness adjustment
            invert_display (bool): Whether to invert display
        """
        if stretch is not None:
            self.stretch = stretch
        if colormap is not None:
            self.colormap = colormap
        if contrast is not None:
            self.contrast = contrast
        if brightness is not None:
            self.brightness = brightness
        if invert_display is not None:
            self.invert_display = invert_display
        
        # Refresh current plot
        if self.current_band:
            self.plot_image(self.current_band)
    
    def get_available_bands(self):
        """Get list of available bands."""
        return list(self.images.keys())
    
    def get_image_info(self, band=None):
        """
        Get information about current or specified image.
        
        Args:
            band (str): Band to get info for (default: current)
        
        Returns:
            dict: Image information
        """
        if band is None:
            band = self.current_band
        
        if band not in self.images:
            return {}
        
        image = self.images[band]
        stats = self.calculate_image_statistics(image)
        
        info = {
            'band': band,
            'shape': image.shape,
            'dtype': str(image.dtype),
            'statistics': stats,
            'total_pixels': image.size,
            'memory_mb': image.nbytes / (1024 * 1024)
        }
        
        return info