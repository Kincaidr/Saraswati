import bdsf

def catalog_generation(fits_image, res_image=False):
    outfile=outname+'_srl.fits'
    img = bdsf.process_image(fits_image, rms_box=(40,40),rms_box_bright=(20,20),adaptive_thresh=150,thresh_isl=4.0,thresh_pix=5.0,
                 detection_image=fits_image,interactive=False,clobber=True,spectralindex_do = False,atrous_do = False, shapelet_do=False)
    img.write_catalog(outfile=outfile,format='fits', catalog_type='srl',clobber=True)
    print("Real catalog written")

if __name__ == "__main__":
    images=['image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored.0.fits','image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored.1.fits','image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored.2.fits']
    out_names=['image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored_998','image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored_1283','image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored_1569']
    for image,outname in zip(images,out_names):
        catalog_generation(image, outname)