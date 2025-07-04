# Synthetic Image Generation for Astro Photometry

## Project Overview

The objective of this project is to generate synthetic images of the sky for astro photometry. The images will be used to plan for an upcoming observation campaign.

## Dependencies

This project uses common astronomy and scientific Python packages:

- `numpy`
- `astropy`
- `astroquery`
- `matplotlib`

> **Note:** Other packages may be added as the project is developed. Claude should always ask before adding a new package and carefully consider why they are being added and what alternatives exist that may be more commonly used in the astronomy and variable star observing community.

## Configuration Files

The system uses JSON configuration files to store reusable settings:

### Target Configuration (`target.json`)
- Target specification by name or RA/Dec coordinates
- Target lookup and storage of significant details
- **Magnitudes:** Johnson V and B bands, Cousins R and I bands

### Telescope Configuration (`telescope.json`)
- Observing telescope details (aperture and focal length)
- **Seeing calculation:** Combined telescope optics and atmospheric seeing
  - Theoretical seeing estimated using Airy disk and Dawes limit formulas
  - User-specified atmospheric seeing conditions
  - Representative FWHM for point spread function

### Camera Configuration (`camera.json`)
- Camera specifications:
  - Pixel size (microns)
  - Sensor dimensions (x and y width)
  - Camera rotation (Y-axis relative to celestial north in degrees)
- **Saturation limit:** User-configurable (typically ~65,535 for 16-bit cameras)

### FITS Header Template (`fits_header.json`)
- Extracted from user-provided example FITS file using `src/fits_header_extractor.py`
- Ensures consistency and allows reuse across sessions
- **Note:** Target-specific keywords (RA, DEC, OBJECT, etc.) will need to be updated when user selects a new target

## Image Generation Process

### 1. Exposure and Limiting Magnitude
- Target brightness set to **50% of full well capacity**
- Default limiting magnitude: **Target magnitude + 5**
- **Saturation handling:** Sources brighter than saturation limit are clipped at maximum ADU

### 2. Multi-band Imaging
Synthetic images generated for each photometric band:
- **V** (Johnson V)
- **B** (Johnson B) 
- **R** (Cousins R)
- **I** (Cousins I)

### 3. Field of View Calculation
Calculated using standard astrophotography formulas:
- Telescope focal length
- Camera pixel size
- Sensor dimensions

### 4. FITS File Creation
- **Format:** 16-bit unsigned integer
- **Headers:** Applied from fits_header.json template
- **Output:** Separate FITS file for each band

## Data Processing Pipeline

### 1. Catalog Query
- **Source:** GAIA DR3 catalog
- **Region:** RA/Dec field of view
- **Limit:** Sources brighter than limiting magnitude
- **Magnitude conversion:** GAIA photometry → V, B, R, I bands

### 2. Point Spread Function (PSF)
- Accounts for both telescope optics and atmospheric seeing
- Uses combined FWHM from telescope configuration

### 3. Source Placement
- Synthetic representations of each light source
- Appropriate PSF applied to each source
- Magnitude-dependent brightness scaling

## Output

- **Storage:** Individual FITS files for each band
- **Display:** Interactive matplotlib visualization
- **Planning:** Ready for observation campaign use