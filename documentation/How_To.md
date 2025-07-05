*How To Generate Synthetic Photometry In Johnson Kron Cousins Standard Colors from GAIA High Resolution Spectra

This method requires several packages
python```
# Step 0: The required packages
from astroquery.gaia import Gaia
import astropy.units as u
from astropy.coordinates import SkyCoord
from gaiaxpy import generate, PhotometricSystem
```

Starting with a known RA and Dec, follow a two step process

First, find the GAIA Source ID.  Notice that the source_id must be a list:
python```
# Star 114 548 RA 22 41 36.833 Dec +00 59 05.80
RA = 15*(22 + 41/60 + 36.833/3600) * u.deg
Dec = (59/60 + 5.8/3600) * u.deg
coord = SkyCoord(ra=RA, dec=Dec, frame='icrs')
width = u.Quantity(1, u.arcsec)
height = u.Quantity(1, u.arcsec)
r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
```

Next, generate the synthetic photometry:
python```
# Step 2: Generate the photometry from the BP and RP spectra
synthetic_photometry = generate(source_id, photometric_system=PhotometricSystem.JKC_Std)
```

The columns returned are:
Index(['source_id', 'JkcStd_mag_U', 'JkcStd_mag_B', 'JkcStd_mag_V',
       'JkcStd_mag_R', 'JkcStd_mag_I', 'JkcStd_flux_U', 'JkcStd_flux_B',
       'JkcStd_flux_V', 'JkcStd_flux_R', 'JkcStd_flux_I',
       'JkcStd_flux_error_U', 'JkcStd_flux_error_B', 'JkcStd_flux_error_V',
       'JkcStd_flux_error_R', 'JkcStd_flux_error_I'],
      dtype='object')

From this object, we are interested in:
'JkcStd_mag_U'
'JkcStd_mag_B'
'JkcStd_mag_V'
'JkcStd_mag_R'
'JkcStd_mag_I'

The magnitude errors can also be calculated from the flux and flux error for each color.