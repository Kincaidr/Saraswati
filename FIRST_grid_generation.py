from astropy.coordinates import SkyCoord
import astropy.units as u

# Parameters
start_coord = SkyCoord('23:41:41.3736 +1:13:29.197', unit=(u.hourangle, u.deg))
n_ra_steps = 5  # rows down (you can change this)
n_dec_steps = 6
  # columns across (you can change this)
step_size_deg = 0.5 # 30 arcmin = 0.5 deg

coords = []

for ra_step in range(n_ra_steps):
    for dec_step in range(n_dec_steps):
        # For each position, RA decreases as you go "down" (East)
        # and Dec increases as you slide "across" (North-South)
        new_coord = SkyCoord(
            ra=start_coord.ra - ra_step * step_size_deg * u.deg,
            dec=start_coord.dec - dec_step * step_size_deg * u.deg
        )
        ra_str = new_coord.ra.to_string(unit=u.hour, sep=':', pad=True, precision=4)
        dec_str = new_coord.dec.to_string(sep=':', alwayssign=True, pad=True, precision=4)
        coords.append(f"{ra_str} {dec_str}")

# Write to file
with open('grid_coordinates.txt', 'w') as f:
    for c in coords:
        f.write(c + '\n')

print("Generated coordinates saved to 'grid_coordinates.txt'")
