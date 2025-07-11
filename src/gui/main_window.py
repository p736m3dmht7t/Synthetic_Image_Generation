"""
Main GUI Window for Synthetic Image Generation

Provides tabbed interface for configuration management and image generation.
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import sys
import os
from pathlib import Path

# Add src directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from config_manager import ConfigManager
from fits_header_extractor import extract_fits_header


class MainWindow:
    """Main application window with tabbed interface."""
    
    def __init__(self):
        """Initialize the main window."""
        self.root = tk.Tk()
        self.root.title("Synthetic Image Generation for Astro Photometry")
        self.root.geometry("800x600")
        
        # Initialize configuration manager
        self.config_manager = ConfigManager()
        
        # Configure TTK styles for button feedback
        self.setup_button_styles()
        
        # Create main interface
        self.create_menu()
        self.create_main_interface()
        
        # Load initial configurations
        self.load_all_configs()
    
    def setup_button_styles(self):
        """Set up TTK styles for visual feedback."""
        style = ttk.Style()
        # Configure success button style (green)
        style.configure("Success.TButton", background="lightgreen")
    
    def create_menu(self):
        """Create the main menu bar."""
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)
        
        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Load Configurations", command=self.load_all_configs)
        file_menu.add_command(label="Save All Configurations", command=self.save_all_configs)
        file_menu.add_separator()
        file_menu.add_command(label="Import FITS Header", command=self.import_fits_header)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)
        
        # Generate menu
        generate_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Generate", menu=generate_menu)
        generate_menu.add_command(label="Generate Images", command=self.generate_images)
    
    def create_main_interface(self):
        """Create the main tabbed interface."""
        # Create notebook for tabs
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Create tabs
        self.create_target_tab()
        self.create_telescope_tab()
        self.create_camera_tab()
        self.create_fits_header_tab()
        self.create_generation_tab()
    
    def create_target_tab(self):
        """Create the target configuration tab."""
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="Target")
        
        # Target name
        ttk.Label(frame, text="Target Name:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.target_name = tk.StringVar()
        target_entry = ttk.Entry(frame, textvariable=self.target_name, width=30)
        target_entry.grid(row=0, column=1, padx=5, pady=5)
        target_entry.bind('<Return>', self.on_target_name_enter)
        
        # Search button for target name (using tk.Button for better color support)
        self.search_target_btn = tk.Button(frame, text="Search", command=self.lookup_target, 
                                         relief="raised", borderwidth=1)
        self.search_target_btn.grid(row=0, column=2, padx=5, pady=5)
        
        # Coordinates section
        ttk.Label(frame, text="Coordinates", font=("Arial", 10, "bold")).grid(row=1, column=0, columnspan=2, sticky="w", padx=5, pady=(15,5))
        
        ttk.Label(frame, text="RA (degrees):").grid(row=2, column=0, sticky="w", padx=5, pady=5)
        self.target_ra = tk.StringVar()
        ra_entry = ttk.Entry(frame, textvariable=self.target_ra, width=15)
        ra_entry.grid(row=2, column=1, sticky="w", padx=5, pady=5)
        ra_entry.bind('<Return>', self.on_coordinates_enter)
        
        # Search button for coordinates (using tk.Button for better color support)
        self.search_coord_btn = tk.Button(frame, text="Search", command=self.search_by_coordinates,
                                        relief="raised", borderwidth=1)
        self.search_coord_btn.grid(row=2, column=2, padx=5, pady=5)
        
        ttk.Label(frame, text="Dec (degrees):").grid(row=3, column=0, sticky="w", padx=5, pady=5)
        self.target_dec = tk.StringVar()
        dec_entry = ttk.Entry(frame, textvariable=self.target_dec, width=15)
        dec_entry.grid(row=3, column=1, sticky="w", padx=5, pady=5)
        dec_entry.bind('<Return>', self.on_coordinates_enter)
        
        ttk.Label(frame, text="RA (HH:MM:SS):").grid(row=4, column=0, sticky="w", padx=5, pady=5)
        self.target_ra_str = tk.StringVar()
        ra_str_entry = ttk.Entry(frame, textvariable=self.target_ra_str, width=15)
        ra_str_entry.grid(row=4, column=1, sticky="w", padx=5, pady=5)
        ra_str_entry.bind('<Return>', self.on_hexagesimal_enter)
        
        # Search button for hexagesimal coordinates (using tk.Button for better color support)
        self.search_hex_btn = tk.Button(frame, text="Search", command=self.search_by_hexagesimal,
                                       relief="raised", borderwidth=1)
        self.search_hex_btn.grid(row=4, column=2, padx=5, pady=5)
        
        ttk.Label(frame, text="Dec (+DD:MM:SS):").grid(row=5, column=0, sticky="w", padx=5, pady=5)
        self.target_dec_str = tk.StringVar()
        dec_str_entry = ttk.Entry(frame, textvariable=self.target_dec_str, width=15)
        dec_str_entry.grid(row=5, column=1, sticky="w", padx=5, pady=5)
        dec_str_entry.bind('<Return>', self.on_hexagesimal_enter)
        
        # Magnitudes section
        ttk.Label(frame, text="Magnitudes", font=("Arial", 10, "bold")).grid(row=6, column=0, columnspan=4, sticky="w", padx=5, pady=(15,5))
        
        # Column headers
        ttk.Label(frame, text="Band", font=("Arial", 9, "bold")).grid(row=7, column=0, sticky="w", padx=5, pady=2)
        ttk.Label(frame, text="Simbad", font=("Arial", 9, "bold")).grid(row=7, column=1, sticky="w", padx=5, pady=2)
        ttk.Label(frame, text="GAIA DR3", font=("Arial", 9, "bold")).grid(row=7, column=2, sticky="w", padx=5, pady=2)
        ttk.Label(frame, text="Source ID", font=("Arial", 9, "bold")).grid(row=7, column=3, sticky="w", padx=5, pady=2)
        
        # B magnitude row
        ttk.Label(frame, text="B magnitude:").grid(row=8, column=0, sticky="w", padx=5, pady=5)
        self.target_mag_b = tk.StringVar()
        ttk.Entry(frame, textvariable=self.target_mag_b, width=10).grid(row=8, column=1, sticky="w", padx=5, pady=5)
        self.target_gaia_mag_b = tk.StringVar()
        gaia_b_entry = ttk.Entry(frame, textvariable=self.target_gaia_mag_b, width=10, state="readonly")
        gaia_b_entry.grid(row=8, column=2, sticky="w", padx=5, pady=5)
        
        # V magnitude row
        ttk.Label(frame, text="V magnitude:").grid(row=9, column=0, sticky="w", padx=5, pady=5)
        self.target_mag_v = tk.StringVar()
        ttk.Entry(frame, textvariable=self.target_mag_v, width=10).grid(row=9, column=1, sticky="w", padx=5, pady=5)
        self.target_gaia_mag_v = tk.StringVar()
        gaia_v_entry = ttk.Entry(frame, textvariable=self.target_gaia_mag_v, width=10, state="readonly")
        gaia_v_entry.grid(row=9, column=2, sticky="w", padx=5, pady=5)
        
        # R magnitude row
        ttk.Label(frame, text="R magnitude:").grid(row=10, column=0, sticky="w", padx=5, pady=5)
        self.target_mag_r = tk.StringVar()
        ttk.Entry(frame, textvariable=self.target_mag_r, width=10).grid(row=10, column=1, sticky="w", padx=5, pady=5)
        self.target_gaia_mag_r = tk.StringVar()
        gaia_r_entry = ttk.Entry(frame, textvariable=self.target_gaia_mag_r, width=10, state="readonly")
        gaia_r_entry.grid(row=10, column=2, sticky="w", padx=5, pady=5)
        
        # I magnitude row
        ttk.Label(frame, text="I magnitude:").grid(row=11, column=0, sticky="w", padx=5, pady=5)
        self.target_mag_i = tk.StringVar()
        ttk.Entry(frame, textvariable=self.target_mag_i, width=10).grid(row=11, column=1, sticky="w", padx=5, pady=5)
        self.target_gaia_mag_i = tk.StringVar()
        gaia_i_entry = ttk.Entry(frame, textvariable=self.target_gaia_mag_i, width=10, state="readonly")
        gaia_i_entry.grid(row=11, column=2, sticky="w", padx=5, pady=5)
        
        # GAIA source ID (spans all magnitude rows)
        self.target_gaia_source_id = tk.StringVar()
        gaia_id_entry = ttk.Entry(frame, textvariable=self.target_gaia_source_id, width=15, state="readonly")
        gaia_id_entry.grid(row=8, column=3, rowspan=4, sticky="new", padx=5, pady=5)
        
        # Notes
        ttk.Label(frame, text="Notes:").grid(row=12, column=0, sticky="nw", padx=5, pady=5)
        self.target_notes = tk.Text(frame, height=3, width=40)
        self.target_notes.grid(row=12, column=1, columnspan=2, padx=5, pady=5)
        
        # Save button
        ttk.Button(frame, text="Save Target Config", command=self.save_target_config).grid(row=13, column=0, columnspan=4, pady=10)
    
    def lookup_target(self):
        """Lookup target information from Simbad and populate fields."""
        target_name = self.target_name.get().strip()
        
        if not target_name:
            messagebox.showwarning("Warning", "Please enter a target name first.")
            return
        
        success = False
        try:
            # Import and create catalog query instance
            from catalog_query import CatalogQuery
            catalog_query = CatalogQuery()
            
            # Show progress
            self.root.config(cursor="wait")
            self.set_button_progress(self.search_target_btn, "Simbad...")
            
            # Lookup target
            target_info = catalog_query.lookup_target(target_name)
            
            if target_info:
                # Populate coordinate fields with exactly 6 decimal places
                ra_str = f"{target_info['ra_deg']:.6f}"
                dec_str = f"{target_info['dec_deg']:+.6f}"  # Always show sign for Dec
                self.target_ra.set(ra_str)
                self.target_dec.set(dec_str)
                self.target_ra_str.set(target_info['ra_str'])
                self.target_dec_str.set(target_info['dec_str'])
                
                # Populate Simbad magnitude fields (only if values are available) with proper formatting
                if target_info['magnitudes']['B'] is not None:
                    self.target_mag_b.set(f"{target_info['magnitudes']['B']:.3f}")
                else:
                    self.target_mag_b.set("")
                if target_info['magnitudes']['V'] is not None:
                    self.target_mag_v.set(f"{target_info['magnitudes']['V']:.3f}")
                else:
                    self.target_mag_v.set("")
                if target_info['magnitudes']['R'] is not None:
                    self.target_mag_r.set(f"{target_info['magnitudes']['R']:.3f}")
                else:
                    self.target_mag_r.set("")  # Clear R magnitude if not available
                if target_info['magnitudes']['I'] is not None:
                    self.target_mag_i.set(f"{target_info['magnitudes']['I']:.3f}")
                else:
                    self.target_mag_i.set("")
                
                # Update target name to canonical form
                self.target_name.set(target_info['name'])
                
                # Show GAIA progress phase
                self.set_button_progress(self.search_target_btn, "GAIA DR3...")
                
                # Now lookup GAIA information for the target
                gaia_info = catalog_query.find_gaia_source_for_target(
                    target_info['ra_deg'], target_info['dec_deg']
                )
                
                if gaia_info:
                    # Check if polynomial method was used
                    asterisk = "*" if gaia_info.get('used_polynomial', False) else ""
                    
                    # Populate GAIA magnitude fields with asterisk marking if polynomial
                    if gaia_info['magnitudes']['B'] is not None:
                        self.target_gaia_mag_b.set(f"{gaia_info['magnitudes']['B']:.3f}{asterisk}")
                    else:
                        self.target_gaia_mag_b.set("")
                    if gaia_info['magnitudes']['V'] is not None:
                        self.target_gaia_mag_v.set(f"{gaia_info['magnitudes']['V']:.3f}{asterisk}")
                    else:
                        self.target_gaia_mag_v.set("")
                    if gaia_info['magnitudes']['R'] is not None:
                        self.target_gaia_mag_r.set(f"{gaia_info['magnitudes']['R']:.3f}{asterisk}")
                    else:
                        self.target_gaia_mag_r.set("")
                    if gaia_info['magnitudes']['I'] is not None:
                        self.target_gaia_mag_i.set(f"{gaia_info['magnitudes']['I']:.3f}{asterisk}")
                    else:
                        self.target_gaia_mag_i.set("")
                    
                    # Set GAIA source ID
                    self.target_gaia_source_id.set(str(gaia_info['source_id']))
                    
                    # Show success immediately
                    self.set_button_success(self.search_target_btn)
                    success = True
                else:
                    # Clear GAIA fields if no GAIA source found
                    self.target_gaia_mag_b.set("")
                    self.target_gaia_mag_v.set("")
                    self.target_gaia_mag_r.set("")
                    self.target_gaia_mag_i.set("")
                    self.target_gaia_source_id.set("")
                    print("Warning: No GAIA source found for target")
                    # Still show success since Simbad lookup succeeded
                    self.set_button_success(self.search_target_btn)
                    success = True
                
            else:
                self.set_button_failure(self.search_target_btn, "Failed: Simbad")
                messagebox.showerror("Error", f"Target '{target_name}' not found in Simbad database.")
                success = False
                
        except Exception as e:
            self.set_button_failure(self.search_target_btn, "Failed")
            messagebox.showerror("Error", f"Error looking up target: {str(e)}")
            success = False
        
        finally:
            # Reset cursor
            self.root.config(cursor="")
    
    def on_target_name_enter(self, event=None):
        """Handle Enter key press in target name field."""
        self.lookup_target()
    
    def on_coordinates_enter(self, event=None):
        """Handle Enter key press in coordinate fields."""
        self.search_by_coordinates()
    
    def on_hexagesimal_enter(self, event=None):
        """Handle Enter key press in hexagesimal coordinate fields."""
        self.search_by_hexagesimal()
    
    def search_by_coordinates(self):
        """Search for target by RA/Dec coordinates in Simbad database."""
        try:
            ra_str = self.target_ra.get().strip()
            dec_str = self.target_dec.get().strip()
            
            if not ra_str or not dec_str:
                messagebox.showwarning("Warning", "Please enter both RA and Dec coordinates.")
                return
            
            # Validate coordinates
            ra_val = float(ra_str)
            dec_val = float(dec_str)
            
            # Basic validation ranges
            if not (0 <= ra_val <= 360):
                messagebox.showerror("Error", "RA must be between 0 and 360 degrees.")
                return
            if not (-90 <= dec_val <= 90):
                messagebox.showerror("Error", "Dec must be between -90 and +90 degrees.")
                return
            
            # Import and create catalog query instance
            from catalog_query import CatalogQuery
            catalog_query = CatalogQuery()
            
            # Show progress
            self.root.config(cursor="wait")
            self.set_button_progress(self.search_coord_btn, "Simbad...")
            
            # Lookup target by coordinates
            target_info = catalog_query.lookup_target_by_coordinates(ra_val, dec_val)
            
            if target_info:
                # Update all fields with Simbad data
                self.target_name.set(target_info['name'])
                
                # Update coordinate fields with Simbad precision
                ra_str = f"{target_info['ra_deg']:.6f}"
                dec_str = f"{target_info['dec_deg']:+.6f}"  # Always show sign for Dec
                self.target_ra.set(ra_str)
                self.target_dec.set(dec_str)
                self.target_ra_str.set(target_info['ra_str'])
                self.target_dec_str.set(target_info['dec_str'])
                
                # Update Simbad magnitude fields with proper formatting
                if target_info['magnitudes']['B'] is not None:
                    self.target_mag_b.set(f"{target_info['magnitudes']['B']:.3f}")
                else:
                    self.target_mag_b.set("")
                if target_info['magnitudes']['V'] is not None:
                    self.target_mag_v.set(f"{target_info['magnitudes']['V']:.3f}")
                else:
                    self.target_mag_v.set("")
                if target_info['magnitudes']['R'] is not None:
                    self.target_mag_r.set(f"{target_info['magnitudes']['R']:.3f}")
                else:
                    self.target_mag_r.set("")
                if target_info['magnitudes']['I'] is not None:
                    self.target_mag_i.set(f"{target_info['magnitudes']['I']:.3f}")
                else:
                    self.target_mag_i.set("")
                
                # Show GAIA progress phase
                self.set_button_progress(self.search_coord_btn, "GAIA DR3...")
                
                # Now lookup GAIA information for the target
                gaia_info = catalog_query.find_gaia_source_for_target(
                    target_info['ra_deg'], target_info['dec_deg']
                )
                
                if gaia_info:
                    # Check if polynomial method was used
                    asterisk = "*" if gaia_info.get('used_polynomial', False) else ""
                    
                    # Populate GAIA magnitude fields with asterisk marking if polynomial
                    if gaia_info['magnitudes']['B'] is not None:
                        self.target_gaia_mag_b.set(f"{gaia_info['magnitudes']['B']:.3f}{asterisk}")
                    else:
                        self.target_gaia_mag_b.set("")
                    if gaia_info['magnitudes']['V'] is not None:
                        self.target_gaia_mag_v.set(f"{gaia_info['magnitudes']['V']:.3f}{asterisk}")
                    else:
                        self.target_gaia_mag_v.set("")
                    if gaia_info['magnitudes']['R'] is not None:
                        self.target_gaia_mag_r.set(f"{gaia_info['magnitudes']['R']:.3f}{asterisk}")
                    else:
                        self.target_gaia_mag_r.set("")
                    if gaia_info['magnitudes']['I'] is not None:
                        self.target_gaia_mag_i.set(f"{gaia_info['magnitudes']['I']:.3f}{asterisk}")
                    else:
                        self.target_gaia_mag_i.set("")
                    
                    # Set GAIA source ID
                    self.target_gaia_source_id.set(str(gaia_info['source_id']))
                    
                    # Show success immediately
                    self.set_button_success(self.search_coord_btn)
                else:
                    # Clear GAIA fields if no GAIA source found
                    self.target_gaia_mag_b.set("")
                    self.target_gaia_mag_v.set("")
                    self.target_gaia_mag_r.set("")
                    self.target_gaia_mag_i.set("")
                    self.target_gaia_source_id.set("")
                    print("Warning: No GAIA source found for target")
                    # Still show success since Simbad lookup succeeded
                    self.set_button_success(self.search_coord_btn)
            else:
                self.set_button_failure(self.search_coord_btn, "Failed: Simbad")
                messagebox.showinfo("No Results", f"No objects found at coordinates:\nRA: {ra_val:.6f}°\nDec: {dec_val:+.6f}°")
            
        except ValueError:
            self.set_button_failure(self.search_coord_btn, "Failed")
            messagebox.showerror("Error", "Invalid coordinate format. Please enter numeric values.")
        except Exception as e:
            self.set_button_failure(self.search_coord_btn, "Failed")
            messagebox.showerror("Error", f"Error searching coordinates: {str(e)}")
        
        finally:
            # Reset cursor
            self.root.config(cursor="")
    
    def search_by_hexagesimal(self):
        """Search for target by RA/Dec hexagesimal coordinates in Simbad database."""
        try:
            ra_str = self.target_ra_str.get().strip()
            dec_str = self.target_dec_str.get().strip()
            
            if not ra_str or not dec_str:
                messagebox.showwarning("Warning", "Please enter both RA and Dec in HH:MM:SS format.")
                return
            
            # Parse hexagesimal coordinates to degrees
            try:
                from astropy.coordinates import SkyCoord
                from astropy import units as u
                
                # Parse RA and Dec strings - SkyCoord can handle various formats
                coord = SkyCoord(ra=ra_str, dec=dec_str, unit=(u.hourangle, u.deg))
                ra_deg = coord.ra.deg
                dec_deg = coord.dec.deg
                
            except Exception as parse_error:
                messagebox.showerror("Error", f"Invalid coordinate format. Please use HH:MM:SS for RA and +DD:MM:SS for Dec.\n\nExample: RA: 22:29:10.26, Dec: +58:24:54.7\n\nError: {str(parse_error)}")
                return
            
            # Import and create catalog query instance
            from catalog_query import CatalogQuery
            catalog_query = CatalogQuery()
            
            # Show progress
            self.root.config(cursor="wait")
            self.set_button_progress(self.search_hex_btn, "Simbad...")
            
            # Lookup target by converted coordinates
            target_info = catalog_query.lookup_target_by_coordinates(ra_deg, dec_deg)
            
            if target_info:
                # Update all fields with Simbad data
                self.target_name.set(target_info['name'])
                
                # Update coordinate fields with Simbad precision
                ra_deg_str = f"{target_info['ra_deg']:.6f}"
                dec_deg_str = f"{target_info['dec_deg']:+.6f}"  # Always show sign for Dec
                self.target_ra.set(ra_deg_str)
                self.target_dec.set(dec_deg_str)
                self.target_ra_str.set(target_info['ra_str'])
                self.target_dec_str.set(target_info['dec_str'])
                
                # Update magnitude fields with proper formatting
                if target_info['magnitudes']['B'] is not None:
                    self.target_mag_b.set(f"{target_info['magnitudes']['B']:.3f}")
                else:
                    self.target_mag_b.set("")
                if target_info['magnitudes']['V'] is not None:
                    self.target_mag_v.set(f"{target_info['magnitudes']['V']:.3f}")
                else:
                    self.target_mag_v.set("")
                if target_info['magnitudes']['R'] is not None:
                    self.target_mag_r.set(f"{target_info['magnitudes']['R']:.3f}")
                else:
                    self.target_mag_r.set("")
                if target_info['magnitudes']['I'] is not None:
                    self.target_mag_i.set(f"{target_info['magnitudes']['I']:.3f}")
                else:
                    self.target_mag_i.set("")
                
                # Show GAIA progress phase
                self.set_button_progress(self.search_hex_btn, "GAIA DR3...")
                
                # Now lookup GAIA information for the target
                gaia_info = catalog_query.find_gaia_source_for_target(
                    target_info['ra_deg'], target_info['dec_deg']
                )
                
                if gaia_info:
                    # Check if polynomial method was used
                    asterisk = "*" if gaia_info.get('used_polynomial', False) else ""
                    
                    # Populate GAIA magnitude fields with asterisk marking if polynomial
                    if gaia_info['magnitudes']['B'] is not None:
                        self.target_gaia_mag_b.set(f"{gaia_info['magnitudes']['B']:.3f}{asterisk}")
                    else:
                        self.target_gaia_mag_b.set("")
                    if gaia_info['magnitudes']['V'] is not None:
                        self.target_gaia_mag_v.set(f"{gaia_info['magnitudes']['V']:.3f}{asterisk}")
                    else:
                        self.target_gaia_mag_v.set("")
                    if gaia_info['magnitudes']['R'] is not None:
                        self.target_gaia_mag_r.set(f"{gaia_info['magnitudes']['R']:.3f}{asterisk}")
                    else:
                        self.target_gaia_mag_r.set("")
                    if gaia_info['magnitudes']['I'] is not None:
                        self.target_gaia_mag_i.set(f"{gaia_info['magnitudes']['I']:.3f}{asterisk}")
                    else:
                        self.target_gaia_mag_i.set("")
                    
                    # Set GAIA source ID
                    self.target_gaia_source_id.set(str(gaia_info['source_id']))
                    
                    # Show success immediately
                    self.set_button_success(self.search_hex_btn)
                else:
                    # Clear GAIA fields if no GAIA source found
                    self.target_gaia_mag_b.set("")
                    self.target_gaia_mag_v.set("")
                    self.target_gaia_mag_r.set("")
                    self.target_gaia_mag_i.set("")
                    self.target_gaia_source_id.set("")
                    print("Warning: No GAIA source found for target")
                    # Still show success since Simbad lookup succeeded
                    self.set_button_success(self.search_hex_btn)
            else:
                self.set_button_failure(self.search_hex_btn, "Failed: Simbad")
                messagebox.showinfo("No Results", f"No objects found at coordinates:\nRA: {ra_str}\nDec: {dec_str}")
            
        except Exception as e:
            self.set_button_failure(self.search_hex_btn, "Failed")
            messagebox.showerror("Error", f"Error searching hexagesimal coordinates: {str(e)}")
        
        finally:
            # Reset cursor
            self.root.config(cursor="")
    
    def set_button_progress(self, button, phase_text):
        """Set button to show progress phase."""
        button.configure(text=phase_text, fg="blue", font=("Arial", 9, "bold"),
                        relief="solid", borderwidth=2)
        self.root.update()
    
    def set_button_success(self, button):
        """Set button to green for success feedback."""
        # Change to success appearance with checkmark
        button.configure(text="✓ Success", fg="green", font=("Arial", 9, "bold"),
                        relief="solid", borderwidth=2)
        # Ensure cursor is normal and force update immediately
        self.root.config(cursor="")
        self.root.update()
        # Schedule reset after 3 seconds - always reset to "Search"
        self.root.after(3000, lambda: self.reset_button_color(button, "Search"))
    
    def set_button_failure(self, button, error_text):
        """Set button to show temporary failure state."""
        # Store original text
        original_text = button.cget('text')
        # Change to failure appearance
        button.configure(text=f"✗ {error_text}", fg="red", font=("Arial", 9, "bold"),
                        relief="solid", borderwidth=2)
        # Force update immediately
        self.root.update()
        # Schedule reset after 3 seconds
        self.root.after(3000, lambda: self.reset_button_color(button, original_text))
    
    def reset_button_color(self, button, original_text="Search"):
        """Reset specific button to default color."""
        try:
            # Reset to default button appearance
            button.configure(text=original_text, fg="black", font=("Arial", 9, "normal"),
                           relief="raised", borderwidth=1)
        except Exception as e:
            pass
    
    def reset_button_colors(self):
        """Reset all search buttons to default color."""
        try:
            self.reset_button_color(self.search_target_btn, "Search")
            self.reset_button_color(self.search_coord_btn, "Search")
            self.reset_button_color(self.search_hex_btn, "Search")
        except:
            pass  # Ignore if buttons don't exist yet
    
    def create_telescope_tab(self):
        """Create the telescope configuration tab."""
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="Telescope")
        
        # Telescope name
        ttk.Label(frame, text="Telescope Name:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.telescope_name = tk.StringVar()
        ttk.Entry(frame, textvariable=self.telescope_name, width=30).grid(row=0, column=1, padx=5, pady=5)
        
        # Telescope specs
        ttk.Label(frame, text="Aperture (mm):").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.telescope_aperture = tk.DoubleVar()
        aperture_entry = ttk.Entry(frame, textvariable=self.telescope_aperture, width=15)
        aperture_entry.grid(row=1, column=1, sticky="w", padx=5, pady=5)
        aperture_entry.bind('<KeyRelease>', self.on_telescope_param_change)
        aperture_entry.bind('<FocusOut>', self.on_telescope_param_change)
        
        ttk.Label(frame, text="Focal Length (mm):").grid(row=2, column=0, sticky="w", padx=5, pady=5)
        self.telescope_focal_length = tk.DoubleVar()
        focal_length_entry = ttk.Entry(frame, textvariable=self.telescope_focal_length, width=15)
        focal_length_entry.grid(row=2, column=1, sticky="w", padx=5, pady=5)
        focal_length_entry.bind('<KeyRelease>', self.on_telescope_param_change)
        focal_length_entry.bind('<FocusOut>', self.on_telescope_param_change)
        
        ttk.Label(frame, text="Focal Ratio:").grid(row=3, column=0, sticky="w", padx=5, pady=5)
        self.telescope_focal_ratio = tk.DoubleVar()
        focal_ratio_entry = ttk.Entry(frame, textvariable=self.telescope_focal_ratio, width=15)
        focal_ratio_entry.grid(row=3, column=1, sticky="w", padx=5, pady=5)
        focal_ratio_entry.bind('<KeyRelease>', self.on_telescope_param_change)
        focal_ratio_entry.bind('<FocusOut>', self.on_telescope_param_change)
        
        # Seeing section
        ttk.Label(frame, text="Seeing", font=("Arial", 10, "bold")).grid(row=4, column=0, columnspan=2, sticky="w", padx=5, pady=(15,5))
        
        ttk.Label(frame, text="Theoretical (arcsec):").grid(row=5, column=0, sticky="w", padx=5, pady=5)
        self.telescope_theoretical_seeing = tk.DoubleVar()
        theoretical_entry = ttk.Entry(frame, textvariable=self.telescope_theoretical_seeing, width=15)
        theoretical_entry.grid(row=5, column=1, sticky="w", padx=5, pady=5)
        theoretical_entry.bind('<KeyRelease>', self.on_seeing_change)
        theoretical_entry.bind('<FocusOut>', self.on_seeing_change)
        
        ttk.Label(frame, text="Atmospheric (arcsec):").grid(row=6, column=0, sticky="w", padx=5, pady=5)
        self.telescope_atmospheric_seeing = tk.DoubleVar()
        atmospheric_entry = ttk.Entry(frame, textvariable=self.telescope_atmospheric_seeing, width=15)
        atmospheric_entry.grid(row=6, column=1, sticky="w", padx=5, pady=5)
        atmospheric_entry.bind('<KeyRelease>', self.on_seeing_change)
        atmospheric_entry.bind('<FocusOut>', self.on_seeing_change)
        
        ttk.Label(frame, text="Combined FWHM (arcsec):").grid(row=7, column=0, sticky="w", padx=5, pady=5)
        self.telescope_combined_fwhm = tk.DoubleVar()
        combined_entry = ttk.Entry(frame, textvariable=self.telescope_combined_fwhm, width=15)
        combined_entry.grid(row=7, column=1, sticky="w", padx=5, pady=5)
        
        # Notes
        ttk.Label(frame, text="Notes:").grid(row=8, column=0, sticky="nw", padx=5, pady=5)
        self.telescope_notes = tk.Text(frame, height=3, width=40)
        self.telescope_notes.grid(row=8, column=1, padx=5, pady=5)
        
        # Save button
        ttk.Button(frame, text="Save Telescope Config", command=self.save_telescope_config).grid(row=9, column=0, columnspan=2, pady=10)
    
    def create_camera_tab(self):
        """Create the camera configuration tab."""
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="Camera")
        
        # Camera name
        ttk.Label(frame, text="Camera Name:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.camera_name = tk.StringVar()
        ttk.Entry(frame, textvariable=self.camera_name, width=30).grid(row=0, column=1, padx=5, pady=5)
        
        # Camera specs
        ttk.Label(frame, text="Pixel Size (microns):").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.camera_pixel_size = tk.DoubleVar()
        pixel_size_entry = ttk.Entry(frame, textvariable=self.camera_pixel_size, width=15)
        pixel_size_entry.grid(row=1, column=1, sticky="w", padx=5, pady=5)
        pixel_size_entry.bind('<KeyRelease>', self.on_camera_param_change)
        pixel_size_entry.bind('<FocusOut>', self.on_camera_param_change)
        
        # Sensor dimensions
        ttk.Label(frame, text="Sensor Dimensions", font=("Arial", 10, "bold")).grid(row=2, column=0, columnspan=2, sticky="w", padx=5, pady=(15,5))
        
        ttk.Label(frame, text="Width (pixels):").grid(row=3, column=0, sticky="w", padx=5, pady=5)
        self.camera_width_pixels = tk.IntVar()
        width_pixels_entry = ttk.Entry(frame, textvariable=self.camera_width_pixels, width=15)
        width_pixels_entry.grid(row=3, column=1, sticky="w", padx=5, pady=5)
        width_pixels_entry.bind('<KeyRelease>', self.on_camera_param_change)
        width_pixels_entry.bind('<FocusOut>', self.on_camera_param_change)
        
        ttk.Label(frame, text="Height (pixels):").grid(row=4, column=0, sticky="w", padx=5, pady=5)
        self.camera_height_pixels = tk.IntVar()
        height_pixels_entry = ttk.Entry(frame, textvariable=self.camera_height_pixels, width=15)
        height_pixels_entry.grid(row=4, column=1, sticky="w", padx=5, pady=5)
        height_pixels_entry.bind('<KeyRelease>', self.on_camera_param_change)
        height_pixels_entry.bind('<FocusOut>', self.on_camera_param_change)
        
        ttk.Label(frame, text="Width (mm):").grid(row=5, column=0, sticky="w", padx=5, pady=5)
        self.camera_width_mm = tk.DoubleVar()
        width_mm_entry = ttk.Entry(frame, textvariable=self.camera_width_mm, width=15)
        width_mm_entry.grid(row=5, column=1, sticky="w", padx=5, pady=5)
        width_mm_entry.bind('<KeyRelease>', self.on_camera_param_change)
        width_mm_entry.bind('<FocusOut>', self.on_camera_param_change)
        
        ttk.Label(frame, text="Height (mm):").grid(row=6, column=0, sticky="w", padx=5, pady=5)
        self.camera_height_mm = tk.DoubleVar()
        height_mm_entry = ttk.Entry(frame, textvariable=self.camera_height_mm, width=15)
        height_mm_entry.grid(row=6, column=1, sticky="w", padx=5, pady=5)
        height_mm_entry.bind('<KeyRelease>', self.on_camera_param_change)
        height_mm_entry.bind('<FocusOut>', self.on_camera_param_change)
        
        # Other settings
        ttk.Label(frame, text="Rotation (degrees):").grid(row=7, column=0, sticky="w", padx=5, pady=5)
        self.camera_rotation = tk.DoubleVar()
        ttk.Entry(frame, textvariable=self.camera_rotation, width=15).grid(row=7, column=1, sticky="w", padx=5, pady=5)
        
        ttk.Label(frame, text="Saturation Limit:").grid(row=8, column=0, sticky="w", padx=5, pady=5)
        self.camera_saturation = tk.IntVar()
        ttk.Entry(frame, textvariable=self.camera_saturation, width=15).grid(row=8, column=1, sticky="w", padx=5, pady=5)
        
        # Notes
        ttk.Label(frame, text="Notes:").grid(row=9, column=0, sticky="nw", padx=5, pady=5)
        self.camera_notes = tk.Text(frame, height=3, width=40)
        self.camera_notes.grid(row=9, column=1, padx=5, pady=5)
        
        # Save button
        ttk.Button(frame, text="Save Camera Config", command=self.save_camera_config).grid(row=10, column=0, columnspan=2, pady=10)
    
    def create_fits_header_tab(self):
        """Create the FITS header configuration tab."""
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="FITS Header")
        
        # Instructions
        ttk.Label(frame, text="FITS Header Template", font=("Arial", 12, "bold")).grid(row=0, column=0, columnspan=2, sticky="w", padx=5, pady=5)
        ttk.Label(frame, text="Import a FITS file to extract header template, then customize as needed.").grid(row=1, column=0, columnspan=2, sticky="w", padx=5, pady=5)
        
        # Import button
        ttk.Button(frame, text="Import FITS Header", command=self.import_fits_header).grid(row=2, column=0, padx=5, pady=10)
        
        # Header display/edit area
        ttk.Label(frame, text="Header Keywords:").grid(row=3, column=0, sticky="nw", padx=5, pady=5)
        
        # Create scrollable text area for header
        header_frame = ttk.Frame(frame)
        header_frame.grid(row=4, column=0, columnspan=2, sticky="nsew", padx=5, pady=5)
        
        self.fits_header_text = tk.Text(header_frame, height=15, width=80)
        header_scrollbar = ttk.Scrollbar(header_frame, orient="vertical", command=self.fits_header_text.yview)
        self.fits_header_text.configure(yscrollcommand=header_scrollbar.set)
        
        self.fits_header_text.grid(row=0, column=0, sticky="nsew")
        header_scrollbar.grid(row=0, column=1, sticky="ns")
        
        header_frame.grid_rowconfigure(0, weight=1)
        header_frame.grid_columnconfigure(0, weight=1)
        
        # Save button
        ttk.Button(frame, text="Save FITS Header Config", command=self.save_fits_header_config).grid(row=5, column=0, columnspan=2, pady=10)
        
        # Configure grid weights
        frame.grid_rowconfigure(4, weight=1)
        frame.grid_columnconfigure(0, weight=1)
    
    def create_generation_tab(self):
        """Create the image generation tab."""
        frame = ttk.Frame(self.notebook)
        self.notebook.add(frame, text="Generate")
        
        # Status
        ttk.Label(frame, text="Image Generation", font=("Arial", 12, "bold")).grid(row=0, column=0, columnspan=2, sticky="w", padx=5, pady=5)
        
        # Configuration status
        self.config_status = tk.StringVar()
        self.config_status.set("Configuration status: Not checked")
        ttk.Label(frame, textvariable=self.config_status).grid(row=1, column=0, columnspan=2, sticky="w", padx=5, pady=5)
        
        # Check config button
        ttk.Button(frame, text="Check Configuration", command=self.check_configuration).grid(row=2, column=0, padx=5, pady=10)
        
        # Generate button
        self.generate_button = ttk.Button(frame, text="Generate Images", command=self.generate_images)
        self.generate_button.grid(row=2, column=1, padx=5, pady=10)
        
        # Progress bar
        self.progress = ttk.Progressbar(frame, mode='indeterminate')
        self.progress.grid(row=3, column=0, columnspan=2, sticky="ew", padx=5, pady=5)
        
        # Status text
        self.status_text = tk.Text(frame, height=10, width=60)
        self.status_text.grid(row=4, column=0, columnspan=2, sticky="nsew", padx=5, pady=5)
        
        # Configure grid weights
        frame.grid_rowconfigure(4, weight=1)
        frame.grid_columnconfigure(1, weight=1)
    
    def create_image_display(self, images, source_data=None):
        """Create image display window."""
        try:
            from visualizer import ImageVisualizer
            
            # Create new window
            display_window = tk.Toplevel(self.root)
            display_window.title("Synthetic Images")
            display_window.geometry("1000x800")
            
            # Create notebook for different views
            notebook = ttk.Notebook(display_window)
            notebook.pack(fill="both", expand=True, padx=10, pady=10)
            
            # Individual band tabs - each gets its own visualizer
            band_visualizers = {}
            for band in images.keys():
                frame = ttk.Frame(notebook)
                notebook.add(frame, text=f"{band} Band")
                
                # Create controls frame
                controls_frame = ttk.Frame(frame)
                controls_frame.pack(fill="x", padx=5, pady=5)
                
                # Display controls
                ttk.Label(controls_frame, text="Stretch:").grid(row=0, column=0, padx=5)
                stretch_var = tk.StringVar(value="linear")  # Start with linear
                stretch_combo = ttk.Combobox(controls_frame, textvariable=stretch_var, 
                                           values=["linear", "zscale", "log", "sqrt"], width=10)
                stretch_combo.grid(row=0, column=1, padx=5)
                
                ttk.Label(controls_frame, text="Colormap:").grid(row=0, column=2, padx=5)
                cmap_var = tk.StringVar(value="gray")
                cmap_combo = ttk.Combobox(controls_frame, textvariable=cmap_var,
                                        values=["gray", "viridis", "plasma", "hot", "cool"], width=10)
                cmap_combo.grid(row=0, column=3, padx=5)
                
                # Source overlay checkbox
                show_sources_var = tk.BooleanVar()
                sources_check = ttk.Checkbutton(controls_frame, text="Show Sources", 
                                              variable=show_sources_var)
                sources_check.grid(row=0, column=4, padx=10)
                
                # Invert display checkbox
                invert_var = tk.BooleanVar()
                invert_check = ttk.Checkbutton(controls_frame, text="Invert", 
                                             variable=invert_var)
                invert_check.grid(row=0, column=5, padx=10)
                
                # Create individual visualizer for this band
                band_visualizer = ImageVisualizer()
                band_visualizer.load_images({band: images[band]})  # Only load this band's image
                if source_data:
                    band_visualizer.set_source_data(source_data)
                band_visualizers[band] = band_visualizer
                
                # Update function for this specific band
                def update_display(band=band, sv=stretch_var, cv=cmap_var, ssv=show_sources_var, iv=invert_var, viz=band_visualizer):
                    viz.update_display_settings(stretch=sv.get(), colormap=cv.get(), invert_display=iv.get())
                    viz.plot_image(band=band, show_sources=ssv.get())
                
                # Manual update button
                update_btn = ttk.Button(controls_frame, text="Update", 
                                      command=update_display)
                update_btn.grid(row=0, column=6, padx=10)
                
                # Create matplotlib figure for this band only
                figure, canvas = band_visualizer.create_matplotlib_figure(frame, figsize=(10, 8))
                
                # Initial plot with better settings
                band_visualizer.stretch = 'linear'  # Start with linear stretch
                band_visualizer.plot_image(band=band)
            
            # Comparison tab
            comp_frame = ttk.Frame(notebook)
            notebook.add(comp_frame, text="Comparison")
            
            # Comparison controls
            comp_controls_frame = ttk.Frame(comp_frame)
            comp_controls_frame.pack(fill="x", padx=5, pady=5)
            
            comp_show_sources_var = tk.BooleanVar()
            comp_sources_check = ttk.Checkbutton(comp_controls_frame, text="Show Sources", 
                                                variable=comp_show_sources_var)
            comp_sources_check.pack(side="left", padx=10)
            
            comp_invert_var = tk.BooleanVar()
            comp_invert_check = ttk.Checkbutton(comp_controls_frame, text="Invert", 
                                               variable=comp_invert_var)
            comp_invert_check.pack(side="left", padx=10)
            
            # Create comparison figure
            comp_visualizer = ImageVisualizer()
            comp_visualizer.load_images(images)
            if source_data:
                comp_visualizer.set_source_data(source_data)
            
            def update_comparison(viz=comp_visualizer, ssv=comp_show_sources_var, iv=comp_invert_var):
                viz.show_sources = ssv.get()
                viz.invert_display = iv.get()
                viz.create_comparison_plot()
            
            comp_update_btn = ttk.Button(comp_controls_frame, text="Update Comparison", 
                                       command=update_comparison)
            comp_update_btn.pack(side="left", padx=10)
            
            comp_figure, comp_canvas = comp_visualizer.create_matplotlib_figure(comp_frame, figsize=(12, 8))
            comp_visualizer.create_comparison_plot()
            
            # Store references to prevent garbage collection
            display_window.band_visualizers = band_visualizers
            display_window.comp_visualizer = comp_visualizer
            
        except Exception as e:
            print(f"Error creating image display: {str(e)}")
            messagebox.showerror("Error", f"Failed to create image display: {str(e)}")
    
    def load_all_configs(self):
        """Load all configuration files."""
        try:
            # Load target config
            target_config = self.config_manager.load_target_config()
            self.target_name.set(target_config.get("target_name", ""))
            coords = target_config.get("coordinates", {})
            # Format RA and Dec with exactly 6 decimal places
            ra_val = coords.get("ra", 0.0)
            dec_val = coords.get("dec", 0.0)
            self.target_ra.set(f"{ra_val:.6f}")
            self.target_dec.set(f"{dec_val:+.6f}")  # Always show sign for Dec
            self.target_ra_str.set(coords.get("ra_str", ""))
            self.target_dec_str.set(coords.get("dec_str", ""))
            mags = target_config.get("magnitudes", {})
            # Format magnitudes to 3 decimal places as strings, handle None values
            b_val = mags.get("B", None)
            v_val = mags.get("V", None)
            r_val = mags.get("R", None)
            i_val = mags.get("I", None)
            self.target_mag_b.set(f"{b_val:.3f}" if b_val is not None and b_val != 0.0 else "")
            self.target_mag_v.set(f"{v_val:.3f}" if v_val is not None and v_val != 0.0 else "")
            self.target_mag_r.set(f"{r_val:.3f}" if r_val is not None and r_val != 0.0 else "")
            self.target_mag_i.set(f"{i_val:.3f}" if i_val is not None and i_val != 0.0 else "")
            
            # Load GAIA data if available
            gaia_data = target_config.get("gaia_data", {})
            gaia_mags = gaia_data.get("magnitudes", {})
            gaia_b_val = gaia_mags.get("B", None)
            gaia_v_val = gaia_mags.get("V", None)
            gaia_r_val = gaia_mags.get("R", None)
            gaia_i_val = gaia_mags.get("I", None)
            self.target_gaia_mag_b.set(f"{gaia_b_val:.3f}" if gaia_b_val is not None and gaia_b_val != 0.0 else "")
            self.target_gaia_mag_v.set(f"{gaia_v_val:.3f}" if gaia_v_val is not None and gaia_v_val != 0.0 else "")
            self.target_gaia_mag_r.set(f"{gaia_r_val:.3f}" if gaia_r_val is not None and gaia_r_val != 0.0 else "")
            self.target_gaia_mag_i.set(f"{gaia_i_val:.3f}" if gaia_i_val is not None and gaia_i_val != 0.0 else "")
            self.target_gaia_source_id.set(gaia_data.get("source_id", ""))
            
            self.target_notes.delete(1.0, tk.END)
            self.target_notes.insert(1.0, target_config.get("notes", ""))
            
            # Load telescope config
            telescope_config = self.config_manager.load_telescope_config()
            self.telescope_name.set(telescope_config.get("telescope_name", ""))
            self.telescope_aperture.set(telescope_config.get("aperture_mm", 0.0))
            self.telescope_focal_length.set(telescope_config.get("focal_length_mm", 0.0))
            self.telescope_focal_ratio.set(telescope_config.get("focal_ratio", 0.0))
            seeing = telescope_config.get("seeing", {})
            self.telescope_theoretical_seeing.set(seeing.get("theoretical_arcsec", 0.0))
            self.telescope_atmospheric_seeing.set(seeing.get("atmospheric_arcsec", 1.5))
            self.telescope_combined_fwhm.set(seeing.get("combined_fwhm_arcsec", 0.0))
            self.telescope_notes.delete(1.0, tk.END)
            self.telescope_notes.insert(1.0, telescope_config.get("notes", ""))
            
            # Load camera config
            camera_config = self.config_manager.load_camera_config()
            self.camera_name.set(camera_config.get("camera_name", ""))
            self.camera_pixel_size.set(camera_config.get("pixel_size_microns", 0.0))
            sensor = camera_config.get("sensor_dimensions", {})
            self.camera_width_pixels.set(sensor.get("width_pixels", 0))
            self.camera_height_pixels.set(sensor.get("height_pixels", 0))
            self.camera_width_mm.set(sensor.get("width_mm", 0.0))
            self.camera_height_mm.set(sensor.get("height_mm", 0.0))
            self.camera_rotation.set(camera_config.get("rotation_degrees", 0.0))
            self.camera_saturation.set(camera_config.get("saturation_limit", 65535))
            self.camera_notes.delete(1.0, tk.END)
            self.camera_notes.insert(1.0, camera_config.get("notes", ""))
            
            # Load FITS header config
            fits_header_config = self.config_manager.load_fits_header_config()
            self.fits_header_text.delete(1.0, tk.END)
            if fits_header_config:
                import json
                header_json = json.dumps(fits_header_config, indent=2)
                self.fits_header_text.insert(1.0, header_json)
            
            # Configuration loaded successfully (removed excessive info popup)
            
        except Exception as e:
            messagebox.showerror("Error", f"Error loading configurations: {str(e)}")
    
    def save_all_configs(self):
        """Save all configuration files."""
        try:
            self.save_target_config()
            self.save_telescope_config()
            self.save_camera_config()
            self.save_fits_header_config()
            # All configurations saved successfully (removed excessive info popup)
        except Exception as e:
            messagebox.showerror("Error", f"Error saving configurations: {str(e)}")
    
    def save_target_config(self):
        """Save target configuration."""
        try:
            # Convert string values back to numbers for saving
            ra_str = self.target_ra.get().strip()
            dec_str = self.target_dec.get().strip()
            ra_val = float(ra_str) if ra_str else 0.0
            dec_val = float(dec_str) if dec_str else 0.0
            
            # Handle magnitude strings - convert to float if not empty, otherwise None
            b_str = self.target_mag_b.get().strip()
            v_str = self.target_mag_v.get().strip()
            r_str = self.target_mag_r.get().strip()
            i_str = self.target_mag_i.get().strip()
            
            b_val = float(b_str) if b_str else None
            v_val = float(v_str) if v_str else None
            r_val = float(r_str) if r_str else None
            i_val = float(i_str) if i_str else None
        except ValueError:
            messagebox.showerror("Error", "Invalid numeric values in coordinate or magnitude fields!")
            return
        
        # Handle GAIA magnitude strings - convert to float if not empty, otherwise None
        try:
            gaia_b_str = self.target_gaia_mag_b.get().strip()
            gaia_v_str = self.target_gaia_mag_v.get().strip()
            gaia_r_str = self.target_gaia_mag_r.get().strip()
            gaia_i_str = self.target_gaia_mag_i.get().strip()
            
            gaia_b_val = float(gaia_b_str) if gaia_b_str else None
            gaia_v_val = float(gaia_v_str) if gaia_v_str else None
            gaia_r_val = float(gaia_r_str) if gaia_r_str else None
            gaia_i_val = float(gaia_i_str) if gaia_i_str else None
            
            gaia_source_id = self.target_gaia_source_id.get().strip()
            gaia_source_id_val = gaia_source_id if gaia_source_id else None
        except ValueError:
            messagebox.showerror("Error", "Invalid numeric values in GAIA magnitude fields!")
            return
            
        config = {
            "target_name": self.target_name.get(),
            "coordinates": {
                "ra": ra_val,
                "dec": dec_val,
                "ra_str": self.target_ra_str.get(),
                "dec_str": self.target_dec_str.get()
            },
            "magnitudes": {
                "B": b_val,
                "V": v_val,
                "R": r_val,
                "I": i_val
            },
            "gaia_data": {
                "source_id": gaia_source_id_val,
                "magnitudes": {
                    "B": gaia_b_val,
                    "V": gaia_v_val,
                    "R": gaia_r_val,
                    "I": gaia_i_val
                }
            },
            "notes": self.target_notes.get(1.0, tk.END).strip()
        }
        
        if not self.config_manager.save_target_config(config):
            messagebox.showerror("Error", "Failed to save target configuration!")
    
    def save_telescope_config(self):
        """Save telescope configuration."""
        config = {
            "telescope_name": self.telescope_name.get(),
            "aperture_mm": self.telescope_aperture.get(),
            "focal_length_mm": self.telescope_focal_length.get(),
            "focal_ratio": self.telescope_focal_ratio.get(),
            "seeing": {
                "theoretical_arcsec": self.telescope_theoretical_seeing.get(),
                "atmospheric_arcsec": self.telescope_atmospheric_seeing.get(),
                "combined_fwhm_arcsec": self.telescope_combined_fwhm.get()
            },
            "notes": self.telescope_notes.get(1.0, tk.END).strip()
        }
        
        if not self.config_manager.save_telescope_config(config):
            messagebox.showerror("Error", "Failed to save telescope configuration!")
    
    def save_camera_config(self):
        """Save camera configuration."""
        config = {
            "camera_name": self.camera_name.get(),
            "pixel_size_microns": self.camera_pixel_size.get(),
            "sensor_dimensions": {
                "width_pixels": self.camera_width_pixels.get(),
                "height_pixels": self.camera_height_pixels.get(),
                "width_mm": self.camera_width_mm.get(),
                "height_mm": self.camera_height_mm.get()
            },
            "rotation_degrees": self.camera_rotation.get(),
            "saturation_limit": self.camera_saturation.get(),
            "notes": self.camera_notes.get(1.0, tk.END).strip()
        }
        
        if not self.config_manager.save_camera_config(config):
            messagebox.showerror("Error", "Failed to save camera configuration!")
    
    def save_fits_header_config(self):
        """Save FITS header configuration."""
        try:
            import json
            header_text = self.fits_header_text.get(1.0, tk.END).strip()
            if header_text:
                config = json.loads(header_text)
                if not self.config_manager.save_fits_header_config(config):
                    messagebox.showerror("Error", "Failed to save FITS header configuration!")
            else:
                messagebox.showwarning("Warning", "No FITS header data to save!")
        except json.JSONDecodeError as e:
            messagebox.showerror("Error", f"Invalid JSON format: {str(e)}")
        except Exception as e:
            messagebox.showerror("Error", f"Error saving FITS header: {str(e)}")
    
    def import_fits_header(self):
        """Import FITS header from a file."""
        file_path = filedialog.askopenfilename(
            title="Select FITS File",
            filetypes=[("FITS files", "*.fits"), ("All files", "*.*")]
        )
        
        if file_path:
            try:
                # Extract header using the utility
                header_dict = extract_fits_header(file_path)
                if header_dict:
                    # Display in text area
                    import json
                    header_json = json.dumps(header_dict, indent=2)
                    self.fits_header_text.delete(1.0, tk.END)
                    self.fits_header_text.insert(1.0, header_json)
                    messagebox.showinfo("Success", f"FITS header imported successfully from {file_path}")
                else:
                    messagebox.showerror("Error", "Failed to extract FITS header!")
            except Exception as e:
                messagebox.showerror("Error", f"Error importing FITS header: {str(e)}")
    
    def check_configuration(self):
        """Check if all configurations are valid."""
        # Simple validation - check if essential fields are filled
        issues = []
        
        if not self.target_name.get():
            issues.append("Target name is required")
        
        if not self.telescope_name.get():
            issues.append("Telescope name is required")
        
        if not self.camera_name.get():
            issues.append("Camera name is required")
        
        if issues:
            self.config_status.set(f"Configuration issues: {', '.join(issues)}")
            return False
        else:
            self.config_status.set("Configuration is valid!")
            return True
    
    def generate_images(self):
        """Generate synthetic images."""
        if not self.check_configuration():
            messagebox.showerror("Error", "Please fix configuration issues before generating images.")
            return
        
        try:
            # Import required modules
            from image_generator import ImageGenerator
            from fits_writer import FITSWriter
            from visualizer import ImageVisualizer
            
            # Clear status text
            self.status_text.delete(1.0, tk.END)
            self.status_text.insert(tk.END, "Starting image generation...\n")
            self.status_text.update()
            
            # Start progress indicator
            self.progress.config(mode='indeterminate')
            self.progress.start()
            
            # Save current configurations
            self.save_all_configs()
            
            # Create image generator
            generator = ImageGenerator(self.config_manager)
            
            # Generate images
            self.status_text.insert(tk.END, "Generating synthetic images...\n")
            self.status_text.update()
            
            success = generator.generate_images()
            
            if success:
                # Write FITS files
                self.status_text.insert(tk.END, "Writing FITS files...\n")
                self.status_text.update()
                
                fits_writer = FITSWriter(self.config_manager)
                output_files = fits_writer.write_band_images(
                    generator.images,
                    generator.target_config,
                    generator.telescope_config,
                    generator.camera_config,
                    generator.psf_params,
                    generator.source_data,
                    generator
                )
                
                if output_files:
                    self.status_text.insert(tk.END, f"Successfully generated {len(output_files)} FITS files:\n")
                    for band, filename in output_files.items():
                        self.status_text.insert(tk.END, f"  {band} band: {filename}\n")
                    
                    # Create visualization
                    self.status_text.insert(tk.END, "Creating visualization...\n")
                    self.status_text.update()
                    
                    self.create_image_display(generator.images, generator.source_data)
                    
                    self.status_text.insert(tk.END, "Image generation complete!\n")
                    messagebox.showinfo("Success", f"Generated {len(output_files)} synthetic images successfully!")
                else:
                    self.status_text.insert(tk.END, "Error: Failed to write FITS files\n")
                    messagebox.showerror("Error", "Failed to write FITS files")
            else:
                self.status_text.insert(tk.END, "Error: Image generation failed\n")
                messagebox.showerror("Error", "Image generation failed")
                
        except Exception as e:
            self.status_text.insert(tk.END, f"Error: {str(e)}\n")
            messagebox.showerror("Error", f"Image generation failed: {str(e)}")
        
        finally:
            # Stop progress indicator
            self.progress.stop()
            self.progress.config(mode='determinate')
    
    def on_telescope_param_change(self, event=None):
        """Handle changes to telescope parameters and perform automatic calculations."""
        try:
            # Get values, handling empty fields
            def get_field_value(var):
                try:
                    value = var.get()
                    return value if value > 0 else 0
                except (tk.TclError, ValueError):
                    return 0
            
            aperture = get_field_value(self.telescope_aperture)
            focal_length = get_field_value(self.telescope_focal_length)
            focal_ratio = get_field_value(self.telescope_focal_ratio)
            
            # Count how many fields have valid values
            valid_fields = sum([aperture > 0, focal_length > 0, focal_ratio > 0])
            
            # Only calculate if exactly 2 fields have values and 1 is empty/zero
            if valid_fields == 2:
                if aperture > 0 and focal_length > 0 and focal_ratio == 0:
                    # Calculate focal ratio: f/# = focal_length / aperture
                    calculated_ratio = focal_length / aperture
                    # Format to 1 decimal place
                    formatted_ratio = round(calculated_ratio, 1)
                    self.telescope_focal_ratio.set(formatted_ratio)
                elif aperture > 0 and focal_ratio > 0 and focal_length == 0:
                    # Calculate focal length: focal_length = aperture * f/#
                    calculated_focal_length = aperture * focal_ratio
                    # Format to 1 decimal place
                    formatted_fl = round(calculated_focal_length, 1)
                    self.telescope_focal_length.set(formatted_fl)
                elif focal_length > 0 and focal_ratio > 0 and aperture == 0:
                    # Calculate aperture: aperture = focal_length / f/#
                    calculated_aperture = focal_length / focal_ratio
                    # Format to 1 decimal place
                    formatted_aperture = round(calculated_aperture, 1)
                    self.telescope_aperture.set(formatted_aperture)
            
            # Calculate theoretical seeing if aperture is available
            if aperture > 0:
                # Dawes limit: resolution = 116 / aperture_mm arcsec
                theoretical_seeing = 116.0 / aperture
                self.telescope_theoretical_seeing.set(round(theoretical_seeing, 2))
                
                # Trigger combined FWHM calculation
                self.calculate_combined_fwhm()
                
        except (tk.TclError, ValueError, ZeroDivisionError):
            # Ignore invalid values during typing
            pass
    
    def on_seeing_change(self, event=None):
        """Handle changes to seeing parameters and calculate combined FWHM."""
        self.calculate_combined_fwhm()
    
    def calculate_combined_fwhm(self):
        """Calculate combined FWHM from theoretical and atmospheric seeing."""
        try:
            theoretical = self.telescope_theoretical_seeing.get()
            atmospheric = self.telescope_atmospheric_seeing.get()
            
            if theoretical > 0 and atmospheric > 0:
                # Combined FWHM: sqrt(theoretical^2 + atmospheric^2)
                combined = (theoretical**2 + atmospheric**2)**0.5
                self.telescope_combined_fwhm.set(round(combined, 2))
                
        except (tk.TclError, ValueError):
            # Ignore invalid values during typing
            pass
    
    def on_camera_param_change(self, event=None):
        """Handle changes to camera parameters and perform automatic calculations."""
        try:
            pixel_size_microns = self.camera_pixel_size.get()
            width_pixels = self.camera_width_pixels.get()
            height_pixels = self.camera_height_pixels.get()
            width_mm = self.camera_width_mm.get()
            height_mm = self.camera_height_mm.get()
            
            # Convert pixel size to mm for calculations
            pixel_size_mm = pixel_size_microns / 1000.0 if pixel_size_microns > 0 else 0
            
            # Calculate mm dimensions if pixel size and pixel dimensions are known
            if pixel_size_microns > 0 and width_pixels > 0 and width_mm == 0:
                calculated_width_mm = width_pixels * pixel_size_mm
                formatted_width_mm = round(calculated_width_mm, 3)
                self.camera_width_mm.set(formatted_width_mm)
                
            if pixel_size_microns > 0 and height_pixels > 0 and height_mm == 0:
                calculated_height_mm = height_pixels * pixel_size_mm
                formatted_height_mm = round(calculated_height_mm, 3)
                self.camera_height_mm.set(formatted_height_mm)
            
            # Calculate pixel dimensions if pixel size and mm dimensions are known
            if pixel_size_microns > 0 and width_mm > 0 and width_pixels == 0:
                calculated_width_pixels = int(round(width_mm / pixel_size_mm))
                self.camera_width_pixels.set(calculated_width_pixels)
                
            if pixel_size_microns > 0 and height_mm > 0 and height_pixels == 0:
                calculated_height_pixels = int(round(height_mm / pixel_size_mm))
                self.camera_height_pixels.set(calculated_height_pixels)
            
            # Calculate pixel size if pixel and mm dimensions are known
            if width_pixels > 0 and width_mm > 0 and pixel_size_microns == 0:
                calculated_pixel_size_mm = width_mm / width_pixels
                calculated_pixel_size_microns = calculated_pixel_size_mm * 1000.0
                formatted_pixel_size = round(calculated_pixel_size_microns, 2)
                self.camera_pixel_size.set(formatted_pixel_size)
            elif height_pixels > 0 and height_mm > 0 and pixel_size_microns == 0:
                calculated_pixel_size_mm = height_mm / height_pixels
                calculated_pixel_size_microns = calculated_pixel_size_mm * 1000.0
                formatted_pixel_size = round(calculated_pixel_size_microns, 2)
                self.camera_pixel_size.set(formatted_pixel_size)
                
        except (tk.TclError, ValueError, ZeroDivisionError):
            # Ignore invalid values during typing
            pass
    
    def run(self):
        """Run the main application loop."""
        self.root.mainloop()


def main():
    """Main entry point for the GUI application."""
    app = MainWindow()
    app.run()


if __name__ == "__main__":
    main()