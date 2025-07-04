# Synthetic Image Generation for Astro Photometry

A Python application for generating synthetic astronomical images to support observation campaign planning and photometry research.

## Overview

This project creates realistic synthetic images of the sky for astro photometry applications. The system generates multi-band FITS images (V, B, R, I) using real catalog data from GAIA DR3, accounting for telescope optics, atmospheric seeing, and camera characteristics.

## Features

- **Multi-band Image Generation**: Creates synthetic images in Johnson V, B and Cousins R, I photometric bands
- **Real Catalog Integration**: Uses GAIA DR3 catalog for accurate star positions and magnitudes
- **Telescope Simulation**: Accounts for telescope optics, atmospheric seeing, and point spread function
- **GUI Configuration**: User-friendly interface for managing telescope, camera, and target configurations
- **FITS Header Management**: Extracts and customizes FITS headers from user-provided templates
- **Observation Planning**: Generates images at 50% full well capacity for optimal exposure planning

## Installation

### Prerequisites

- Python 3.13.5 or higher
- Virtual environment (recommended)

### Setup

1. Clone or download the project
2. Create and activate virtual environment:
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Running the Application

```bash
python src/main.py
```

This launches the GUI application for configuring and generating synthetic images.

### Configuration

The system uses JSON configuration files stored in the `config/` directory:

- **`target.json`** - Target specification (name, coordinates, magnitudes)
- **`telescope.json`** - Telescope parameters (aperture, focal length, seeing)
- **`camera.json`** - Camera specifications (pixel size, dimensions, rotation)
- **`fits_header.json`** - FITS header template extracted from example files

### Image Generation Process

1. **Configure Target**: Set target object name or RA/Dec coordinates
2. **Set Telescope Parameters**: Define aperture, focal length, and seeing conditions
3. **Configure Camera**: Specify pixel size, sensor dimensions, and rotation
4. **Import FITS Headers**: Extract headers from existing FITS files as templates
5. **Generate Images**: Create synthetic V, B, R, I band images

## Technical Details

### Magnitude System
- Uses Johnson V, B and Cousins R, I photometric bands
- Converts GAIA photometry to standard bands
- Default limiting magnitude: Target magnitude + 5

### Point Spread Function
- Combines telescope optics (Airy disk) and atmospheric seeing
- Uses Dawes limit and user-specified atmospheric conditions
- Gaussian PSF model with calculated FWHM

### Output Format
- 16-bit unsigned integer FITS files
- Separate file for each photometric band
- Customizable FITS headers based on user templates
- Sources brighter than saturation limit are clipped

## Project Structure

```
Synthetic_Image_Generation/
├── src/                    # Source code
│   ├── main.py            # Application entry point
│   ├── gui/               # GUI components
│   ├── catalog_query.py   # GAIA catalog interface
│   ├── config_manager.py  # Configuration management
│   ├── coordinates.py     # Coordinate transformations
│   ├── fits_header_extractor.py  # FITS header extraction
│   ├── fits_writer.py     # FITS file output
│   ├── image_generator.py # Core image generation
│   ├── optics.py          # Telescope optics calculations
│   ├── psf.py             # Point spread function
│   └── visualizer.py      # Image visualization
├── config/                # Configuration files
├── output/                # Generated FITS files
├── tasks/                 # Development planning
└── requirements.txt       # Python dependencies
```

## Dependencies

- **numpy** - Numerical computations
- **scipy** - Scientific computing and signal processing
- **astropy** - Astronomy utilities and FITS handling
- **astroquery** - GAIA catalog queries
- **matplotlib** - Visualization and GUI

## Development Status

This project is actively under development. See `tasks/todo.md` for the current development plan and progress tracking.

## Use Cases

- **Observation Planning**: Generate realistic images for exposure time estimation
- **Photometry Research**: Create controlled datasets for testing photometric algorithms
- **Education**: Demonstrate astronomical imaging principles and techniques
- **Instrument Testing**: Simulate observations for camera and telescope validation

## Contributing

This project follows a structured development approach with detailed task tracking. Please refer to the development documentation in `tasks/todo.md` for current priorities and implementation guidelines.

## License

[Add license information as appropriate]