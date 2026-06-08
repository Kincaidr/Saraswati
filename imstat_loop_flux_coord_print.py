import glob
from casatools import image
from casatasks import imstat
import numpy as np
import re
import csv

# Helper functions to convert HMS/DMS strings to degrees and back

def hms_to_deg(hms_str):
    """Convert RA string 'hh:mm:ss.sss' or 'hh.mm.ss.sss' to degrees."""
    hms_str = hms_str.replace('.', ':') if '.' in hms_str and ':' not in hms_str else hms_str
    parts = re.split('[:.]', hms_str)
    if len(parts) < 3:
        return None
    h, m, s = float(parts[0]), float(parts[1]), float(parts[2])
    return 15 * (h + m / 60 + s / 3600)

def dms_to_deg(dms_str):
    """Convert Dec string '+dd:mm:ss.sss' or '+dd.mm.ss.sss' to degrees."""
    dms_str = dms_str.strip()
    dms_str = dms_str.replace('.', ':') if '.' in dms_str and ':' not in dms_str else dms_str
    sign = 1
    if dms_str.startswith('-'):
        sign = -1
        dms_str = dms_str[1:]
    elif dms_str.startswith('+'):
        dms_str = dms_str[1:]
    parts = dms_str.split(':')
    if len(parts) < 3:
        return None
    d, m, s = float(parts[0]), float(parts[1]), float(parts[2])
    return sign * (d + m / 60 + s / 3600)

def deg_to_hms(ra_deg):
    """Convert RA degrees to hh:mm:ss.sss string"""
    total_seconds = ra_deg * 240
    h = int(total_seconds // 3600)
    m = int((total_seconds % 3600) // 60)
    s = total_seconds % 60
    return f"{h:02d}:{m:02d}:{s:06.3f}"

def deg_to_dms(dec_deg):
    """Convert Dec degrees to ±dd:mm:ss.sss string"""
    sign = '+' if dec_deg >= 0 else '-'
    dec_deg = abs(dec_deg)
    d = int(dec_deg)
    m = int((dec_deg - d) * 60)
    s = (dec_deg - d - m / 60) * 3600
    return f"{sign}{d:02d}:{m:02d}:{s:06.3f}"

# File and region file names
File = ['/home/kincaid/Desktop/Saraswati_codes/A2631/flux_scale/image_DI_Clustered.DeeperDeconv.AP4.int.restored.fits.img.rgd.fits']
#File= ['mypipelinerun_vlass_casa1.alpha','mypipelinerun_vlass_casa1_sc_simple.alpha','mypipelinerun_vlass_casa1_field_comb_peel.alpha']

region_file = 'bright_srcs2.crtf'  # replace with your actual region file name

def parse_box_region_center(region_line):
    pattern = r'box\s*\[\[\s*([^,]+),\s*([^\]]+)\],\s*\[\s*([^,]+),\s*([^\]]+)\]\]'
    m = re.search(pattern, region_line)
    if not m:
        return None, None
    ra1_str, dec1_str, ra2_str, dec2_str = m.group(1).strip(), m.group(2).strip(), m.group(3).strip(), m.group(4).strip()
    ra1_deg = hms_to_deg(ra1_str)
    ra2_deg = hms_to_deg(ra2_str)
    dec1_deg = dms_to_deg(dec1_str)
    dec2_deg = dms_to_deg(dec2_str)
    if None in (ra1_deg, ra2_deg, dec1_deg, dec2_deg):
        return None, None
    ra_center_deg = (ra1_deg + ra2_deg) / 2
    dec_center_deg = (dec1_deg + dec2_deg) / 2
    ra_center_str = deg_to_hms(ra_center_deg)
    dec_center_str = deg_to_dms(dec_center_deg)
    return ra_center_str, dec_center_str

with open(region_file) as f:
    region_lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]

n_regions = len(region_lines)
n_images = len(File)

flux_table = np.zeros((n_regions, n_images))
peak_table = np.zeros((n_regions, n_images))
ra_list = []
dec_list = []

for region_line in region_lines:
    ra, dec = parse_box_region_center(region_line)
    ra_list.append(ra)
    dec_list.append(dec)

for idx, (ra, dec) in enumerate(zip(ra_list, dec_list)):
    if ra is None or dec is None:
        print(f"Warning: Could not parse RA/Dec from region line {idx+1}: {region_lines[idx]}")

for i, region_line in enumerate(region_lines):
    print(f"Region {i+1}: {region_line}")
    for j, image_file in enumerate(File):
        try:
            stats = imstat(imagename=image_file, region=region_line)
            flux = stats['flux'][0] * 1000  # mJy
            peak = stats['max'][0] * 1000   # mJy/beam
            flux_table[i, j] = flux
            peak_table[i, j] = peak
        except Exception as e:
            print(f"Error with {image_file}, region {i+1}: {e}")
            #flux_table[i, j] = np.nan
            peak_table[i, j] = np.nan

print("\n### Flux Density Table (mJy):")
header = f"{'Region':<6} {'RA':<15} {'Dec':<15}" + "".join(f"{img:^25}" for img in File)
print(header)
for i in range(n_regions):
    ra_val = ra_list[i] if ra_list[i] is not None else 'N/A'
    dec_val = dec_list[i] if dec_list[i] is not None else 'N/A'
    row = f"{i+1:<6} {ra_val:<15} {dec_val:<15}" + \
          "".join(f"{flux_table[i, j]:>12.3f} mJy" for j in range(n_images))
    print(row)

print("\n### Peak Intensity Table (mJy/beam):")
print(header)
for i in range(n_regions):
    ra_val = ra_list[i] if ra_list[i] is not None else 'N/A'
    dec_val = dec_list[i] if dec_list[i] is not None else 'N/A'
    row = f"{i+1:<6} {ra_val:<15} {dec_val:<15}" + \
          "".join(f"{peak_table[i, j]:>12.3f} mJy" for j in range(n_images))
    print(row)

# --- Save flux_table to CSV ---
flux_csv_file = 'new_flux_values.csv'
with open(flux_csv_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    header_row = ['Region'] + [f"Flux_{img}" for img in File]
    writer.writerow(header_row)
    for i in range(n_regions):
        row = [i+1] + list(np.round(flux_table[i], 3))
        writer.writerow(row)

# --- Save peak_table to CSV ---
peak_csv_file = 'new_peak_values.csv'
with open(peak_csv_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    header_row = ['Region'] + [f"Peak_{img}" for img in File]
    writer.writerow(header_row)
    for i in range(n_regions):
        row = [i+1] + list(np.round(peak_table[i], 3))
        writer.writerow(row)

print(f"\nFlux values saved to {flux_csv_file}")
print(f"Peak values saved to {peak_csv_file}")

# --- Generate DS9 region file based on flux differences ---
ds9_region_file = 'flux_difference_circles.reg'

# Indexes of the two images to compare
imageA_idx = File.index('/home/kincaid/Desktop/Saraswati_codes/A2631/flux_scale/image_DI_Clustered.DeeperDeconv.AP4.int.restored.fits.img.rgd.fits')
imageB_idx = File.index('/home/kincaid/Desktop/Saraswati_codes/A2631/flux_scale/A2631_FIRST_large.fits.img.conv.fits')

flux_A = flux_table[:, imageA_idx]
flux_B = flux_table[:, imageB_idx]
flux_diff = np.abs(flux_A - flux_B)

# --- Map difference to radius (50" to 10") using percentile-based scaling ---
low, high = np.nanpercentile(flux_diff, [5, 95])
diff_range = high - low if high != low else 1.0
scaled = np.clip((flux_diff - low) / diff_range, 0, 1)
radii_arcsec = 50 - 40 * scaled  # inverse mapping
radii_arcsec = np.clip(radii_arcsec, 10, 50)  # clamp
radii_arcsec_int = np.round(radii_arcsec).astype(int)

# --- Write DS9 region file ---
ds9_region_file = 'flux_diff_circles.reg'
with open(ds9_region_file, 'w') as f:
    f.write('# Region file format: DS9 version 4.1\n')
    f.write('global color=red font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0\n')
    f.write('fk5\n')
    for ra, dec, radius in zip(ra_list, dec_list, radii_arcsec_int):
        if ra is None or dec is None:
            continue
        f.write(f'circle({ra},{dec},{radius}\")\n')

print(f"DS9 region file written to: {ds9_region_file}")
