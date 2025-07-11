def convovle_cutout(cutout_file):
    from astropy.convolution import convolve, Gaussian2DKernel
    from astropy.io import fits
    with fits.open(cutout_file) as hdul:
        data = hdul[0].data

    kernel_fwhm = np.sqrt(target_beam**2 - original_beam**2)
    kernel_stddev_pix = (kernel_fwhm / pix_size) / 2.3548  # FWHM to stddev
    kernel = Gaussian2DKernel(x_stddev=kernel_stddev_pix)
    convolved_data = convolve(data, kernel)
    convolved_filename = cutout_file.replace('.fits', '_convolved.fits')
    fits.writeto(convolved_filename, convolved_data, overwrite=True)
    print(f"Convolved cutout saved to {convolved_filename}")


if __name__ == "__main__":
    import numpy as np
    import sys

    if len(sys.argv) != 4:
        print("Usage: python convovle_cutout.py <cutout_file> <original_beam> <target_beam>")
        sys.exit(1)

    cutout_file = sys.argv[1]
    original_beam = float(sys.argv[2])
    target_beam = float(sys.argv[3])
    pix_size = 1.5

    convovle_cutout(cutout_file)