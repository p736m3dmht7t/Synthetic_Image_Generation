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

## Flux Calculation Coefficient Review

### Issue: Missing Literature References for GAIA-to-Johnson-Cousins Transformations

**Location**: `src/catalog_query.py:195-208`

**Problem**: The empirical transformation coefficients used to convert GAIA photometry (G, BP, RP) to Johnson-Cousins magnitudes (V, B, R, I) lack proper literature citations. The code contains only generic comments about "empirical transformations from literature" without specific references.

**Current Coefficients**:
- **V band**: V = G - 0.01760 - 0.006860×(BP-RP) - 0.1732×(BP-RP)²
- **B band**: B = G + 0.3130 + 0.2271×(BP-RP) + 0.01397×(BP-RP)²
- **R band**: R = G - 0.4980 - 0.0916×(BP-RP) - 0.0594×(BP-RP)²
- **I band**: I = G - 0.7597 - 0.3346×(BP-RP) - 0.0424×(BP-RP)²

**Required Actions**:
1. **Literature Search**: Identify the original paper(s) that derived these specific coefficients
2. **Coefficient Validation**: Verify coefficients against published transformations
3. **Documentation Update**: Add proper citations to code comments
4. **Alternative Transformations**: Research if more accurate/recent transformations exist
5. **Uncertainty Analysis**: Determine transformation accuracy and color range validity

**Potential Sources to Check**:
- GAIA Data Release documentation
- Photometric transformation papers in AJ, A&A, MNRAS
- Recent studies comparing GAIA to Johnson-Cousins systems
- ESA GAIA mission publications

**Impact**: Without proper validation, synthetic image magnitudes may have systematic errors affecting flux calculations throughout the pipeline.

### Issue: Custom PSF Implementation Instead of Scientific Libraries

**Location**: `src/psf.py:21-96`

**Problem**: The project implements PSF generation from scratch using basic NumPy operations rather than utilizing scientifically accepted and validated functions from astropy or photutils libraries.

**Current Implementation**:
- Manual Gaussian PSF: `psf = np.exp(-(X**2 + Y**2) / (2 * sigma_pixels**2))`
- Manual Moffat PSF: `psf = (1 + (r / alpha)**2)**(-beta)`
- Custom normalization and discretization algorithms

**Required Actions**:
1. **Library Integration**: Replace custom PSF functions with astropy/photutils equivalents
2. **Scientific Validation**: Use peer-reviewed PSF models from established astronomy libraries
3. **Enhanced Features**: Leverage advanced PSF capabilities (elliptical profiles, subpixel sampling)
4. **Performance Optimization**: Benefit from optimized library implementations
5. **Documentation**: Use well-documented, scientifically accepted PSF models

**Potential Libraries to Use**:
- `astropy.modeling.functional_models.Gaussian2D`
- `astropy.modeling.functional_models.Moffat2D`
- `photutils.psf` module for advanced PSF handling
- `astropy.convolution` for PSF convolution operations

**Impact**: Custom PSF implementations may lack scientific rigor, accuracy, and advanced features available in established astronomy libraries. Using validated scientific libraries ensures reproducible and accurate synthetic image generation.

### Issue: Missing Realistic Sky Background Implementation

**Location**: `src/image_generator.py`, `src/psf.py`, configuration files

**Problem**: The synthetic image generation system produces unrealistic images with pure black backgrounds (ADU = 0), completely lacking sky background modeling that characterizes real astronomical observations. This fundamental omission severely limits the system's utility for realistic observation planning and photometry simulation.

**Current State**:
- Images created with `np.zeros()` - pure black backgrounds
- No sky brightness implementation (mag/arcsec²)
- No Bortle scale or light pollution modeling
- No sky-related Poisson noise
- No background flux calibration
- Noise model only includes stellar photons + read noise

**Required Implementation**:

**1. Sky Configuration System**:
- Create `config/sky.json` for observation site characterization
- Include Bortle scale ratings (1-9)
- Band-specific sky brightness values (V, B, R, I mag/arcsec²)
- Typical values: Bortle 1 (22.0 mag/arcsec²), Bortle 4 (21.5), Bortle 7 (19.0)
- Seasonal/lunar phase variations
- Atmospheric extinction coefficients

**2. Sky Brightness Calculation Module**:
- `src/sky_background.py` - Dedicated sky modeling module
- Convert Bortle scale to quantitative sky brightness
- Band-dependent sky brightness lookup tables
- Altitude/airmass corrections
- Moon phase and position effects
- Zodiacal light modeling

**3. Background Flux Calibration**:
- Use same zero point (25.0 mag) and target star reference (50% saturation)
- Convert sky brightness (mag/arcsec²) to flux using Pogson's equation
- Scale to ADU per pixel using plate scale and pixel size
- Account for telescope aperture and exposure time
- Ensure consistent flux calibration with stellar sources

**4. Enhanced Noise Model**:
- Add sky background Poisson noise: `sqrt(sky_adu_per_pixel)`
- Total noise budget: `sqrt(star_adu + sky_adu + read_noise²)`
- Background-limited vs. read-limited noise regimes
- Proper signal-to-noise ratio calculations

**5. Image Generation Integration**:
- Modify `create_stellar_field()` to add uniform sky background
- Apply sky ADU to entire image before placing sources
- Include sky noise in total noise calculation
- Provide realistic detection limits and photometry accuracy

**6. Configuration Parameters**:
```json
{
  "sky_conditions": {
    "bortle_scale": 4,
    "sky_brightness_mag_per_arcsec2": {
      "V": 21.5,
      "B": 22.8, 
      "R": 20.8,
      "I": 19.9
    },
    "atmospheric_extinction": {
      "V": 0.12,
      "B": 0.20,
      "R": 0.08,
      "I": 0.05
    },
    "moon_phase": 0.0,
    "airmass": 1.2,
    "zodiacal_light": true
  }
}
```

**7. GUI Integration**:
- Add sky conditions panel to main GUI
- Bortle scale selection widget
- Sky brightness value display/editing
- Preview of sky background level
- Integration with existing configuration tabs

**Scientific Considerations**:
- Use established sky brightness measurements from literature
- Implement proper photometric calibration chain
- Account for telescope throughput and camera quantum efficiency
- Include atmospheric extinction effects
- Model realistic noise characteristics for different sky conditions

**Validation Requirements**:
- Compare synthetic images with real astronomical observations
- Verify limiting magnitude calculations
- Test photometry accuracy under various sky conditions
- Validate signal-to-noise ratios for faint sources

**Impact**: Without realistic sky background, synthetic images cannot accurately represent actual observing conditions, limiting their utility for observation planning, instrument testing, and photometric analysis. This implementation is critical for producing scientifically valid synthetic astronomical images.