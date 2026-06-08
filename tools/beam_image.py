import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
from matplotlib.patches import Circle

fits_file = 'mypipelinerun_ABELL2631_4-MFS-image.pb.fits'  # Replace with your file path
with fits.open(fits_file) as hdul:
    data = hdul[0].data.squeeze()

# Parameters
target_value = 0.3
tolerance = 0.01  # Accept values in [0.29, 0.31]

# Find pixels near target value
mask = np.abs(data - target_value) < tolerance

# Label connected regions
labeled = label(mask)
regions = regionprops(labeled)

# Plot image with circles
fig, ax = plt.subplots()
im=ax.imshow(data, origin='lower', cmap='gray')
contours = ax.contour(data, levels=[0.3], colors='red', linewidths=2)

plt.colorbar(im, ax=ax)
plt.show()
