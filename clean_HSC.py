from astropy.table import Table

def clean_table(table):
    cat=Table.read(table)
    redshift=cat['zph']
    clean=cat['clean']
    breakpoint()
    mask1=clean==1
    mask2 = (redshift > 0) & (redshift < 4)
    cat=cat[mask1 & mask2]
    cat.write(path+'HSC_MeerKAT_combined_clean.fits', overwrite=True)

if __name__ == "__main__":
    path='/home/kincaid/Desktop/Saraswati_codes/catalogs/'
    table = path+'HSC_MeerKAT_combined.fits'  # Replace with your actual table file
    clean_table(table)
    print("Table cleaned and saved'")