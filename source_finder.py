import bdsf

def catalog_generation(fits_image, output_path, output_name, res_image=True):
    outfile=output_path+output_name+'_srl.fits'
    img = bdsf.process_image(fits_image, rms_box=(40,40),rms_box_bright=(20,20),adaptive_thresh=150,thresh_isl=3.0,thresh_pix=5.0,
                 detection_image=fits_image,interactive=False,clobber=True,spectralindex_do = False,atrous_do = False, shapelet_do=False)
    img.write_catalog(outfile=outfile,format='fits', catalog_type='srl',clobber=True)
    print("Real catalog written")

    if res_image:
        img.export_image(outfile=path+output_name+"_full_res_map.fits",clobber=True,img_type='gaus_resid')
        img.export_image(outfile=path+output_name+"_full_rms_map.fits",clobber=True,img_type='rms')
        print("residual image written")
    else:
        print("residual image not written")
    return(outfile)

if __name__ == "__main__":
    name='A2631'
    path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/'
    out_path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'
    #images=['image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored.0.fits','image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored.1.fits','image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored.2.fits']
    #images=['image_DI_Clustered.DeeperDeconv.AP.int.restored.fits']
    #images=['image_DI_Clustered.DeeperDeconv.AP4.int.restored.fits']
    images=['A2631_cutout_image.fits']
    #out_names=['image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored_998','image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored_1283','image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored_1569']
    out_names=name
    for i, image in enumerate(images):
        catalog_generation(path+image, out_path, out_names)
