# Synthetic Image Generation Project - Development Plan

## Project Overview
Creating a synthetic image generation system for astro photometry with GUI configuration management.

## Step-by-Step Development Plan

### Phase 1: Project Foundation
- [X] **Task 1.1**: Create organized directory structure
  - `src/` - Main source code
  - `config/` - Configuration files (JSON templates)
  - `output/` - Generated FITS files
  - `tests/` - Unit tests (future)

- [X] **Task 1.2**: Create requirements.txt with essential packages
  - `numpy` - Numerical computations
  - `astropy` - Astronomy utilities and FITS handling
  - `astroquery` - GAIA catalog queries
  - `matplotlib` - Visualization

- [X] **Task 1.3**: Install dependencies in virtual environment
  - Run `pip install -r requirements.txt`

### Phase 2: Configuration System
- [X] **Task 2.1**: Create FITS header extraction utility
  - `src/fits_header_extractor.py` - Python script using astropy
  - Reads any FITS file and extracts relevant keywords
  - Saves to fits_header.json template

- [X] **Task 2.2**: Create JSON configuration templates
  - `config/target.json` - Target specification (name, RA/Dec, magnitudes)
  - `config/telescope.json` - Telescope details (aperture, focal length, seeing)
  - `config/camera.json` - Camera specs (pixel size, dimensions, rotation, saturation)

- [X] **Task 2.3**: Create configuration manager module
  - `src/config_manager.py` - Load/save/validate JSON configurations

### Phase 3: GUI Application
- [X] **Task 3.1**: Create main GUI window structure
  - `src/gui/main_window.py` - Main application window
  - Tab-based interface for each configuration type

- [X] **Task 3.2**: Create configuration edit panels
  - `src/gui/target_panel.py` - Target configuration editor
  - `src/gui/telescope_panel.py` - Telescope configuration editor
  - `src/gui/camera_panel.py` - Camera configuration editor
  - `src/gui/fits_header_panel.py` - FITS header editor with import functionality

- [X] **Task 3.3**: Add FITS header import functionality
  - File browser to select example FITS file
  - Run extraction script and populate fits_header.json
  - Allow user customization of headers per target

### Phase 4: Core Astronomy Functionality
- [X] **Task 4.1**: Implement catalog query system
  - `src/catalog_query.py` - GAIA DR3 queries by RA/Dec field
  - Magnitude conversion (GAIA → V, B, R, I bands)
  - Limiting magnitude filtering

- [X] **Task 4.2**: Implement telescope optics calculations
  - `src/optics.py` - Field of view calculations
  - PSF/seeing calculations (Airy disk, Dawes limit)
  - Combined telescope + atmospheric seeing

- [X] **Task 4.3**: Implement coordinate transformations
  - `src/coordinates.py` - RA/Dec to pixel coordinates
  - Camera rotation handling
  - Field of view boundary calculations

### Phase 5: Image Generation
- [X] **Task 5.1**: Create PSF generation system
  - `src/psf.py` - Point spread function modeling
  - Gaussian PSF with calculated FWHM
  - Source brightness scaling

- [X] **Task 5.2**: Implement multi-band image synthesis
  - `src/image_generator.py` - Main image generation engine
  - Generate V, B, R, I band images
  - 50% target brightness scaling
  - Saturation clipping at maximum ADU

- [X] **Task 5.3**: Create source placement system
  - Place catalog sources in image coordinates
  - Apply appropriate PSF to each source
  - Handle magnitude-dependent brightness

### Phase 6: Output and Visualization
- [X] **Task 6.1**: Implement FITS file output
  - `src/fits_writer.py` - 16-bit unsigned integer FITS files
  - Apply customizable headers from user's template
  - Separate file for each band

- [X] **Task 6.2**: Create visualization system
  - `src/visualizer.py` - Interactive matplotlib display
  - Multi-band image comparison
  - Source identification overlay

- [X] **Task 6.3**: Integrate GUI with generation pipeline
  - Add "Generate Images" button to main GUI
  - Progress indicator for long operations
  - Display generated images in GUI

### Phase 7: Integration and Testing
- [X] **Task 7.1**: Create main application entry point
  - `src/main.py` - Launch GUI application

- [ ] **Task 7.2**: Test full workflow
  - Create sample configurations
  - Import FITS headers from example file
  - Generate test images
  - Verify FITS output quality

- [ ] **Task 7.3**: Add error handling and validation
  - Configuration validation
  - User input validation
  - Graceful error messages

## Current Status
- Virtual environment: ✅ Active (Python 3.13.5)
- Project directory: ✅ Created
- FITS example file: ✅ Available for header extraction

## Key Workflow Notes
- FITS header extraction is part of user workflow, not design-time
- Each target can have customized FITS headers
- GUI provides import functionality for FITS header templates
- User can modify headers as needed for each observation campaign

## Next Steps
Ready to begin Phase 1: Project Foundation

## Review Section
*To be filled in as tasks are completed*

### Completed Tasks Summary
*Will be updated as work progresses*

### Key Decisions Made
*Will document important architectural decisions*

### Known Issues
*Will track any issues discovered during development*