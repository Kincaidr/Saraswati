import bdsf

def catalog_generation(fits_image, output_path, output_name, res_image=True):
    img = bdsf.process_image(fits_image, rms_box=(40,40),rms_box_bright=(20,20),adaptive_thresh=150,thresh_isl=3.0,thresh_pix=5.0,
                 detection_image=fits_image,interactive=False,clobber=True,spectralindex_do = False,atrous_do = False, shapelet_do=False)
    outfile=output_path+output_name+'_srl.fits'
    img.write_catalog(outfile=outfile,format='fits', catalog_type='srl',clobber=True)
    outfile=output_path+output_name+'_gaus_srl.fits'
    img.write_catalog(outfile=outfile,format='fits', catalog_type='gaul',clobber=True)
    print("Real catalog written")

    if res_image:
        img.export_image(outfile=path+output_name+"_res_map.fits",clobber=True,img_type='gaus_resid')
        img.export_image(outfile=path+output_name+"_rms_map.fits",clobber=True,img_type='rms')
        img.export_image(outfile=path+output_name+"_gaus_model.fits",clobber=True,img_type='gaus_model')
        img.export_image(outfile=path+output_name+"_island_mask.fits",clobber=True,img_type='island_mask')
        print("residual image written")
    else:
        print("residual image not written")
    return(outfile)

if __name__ == "__main__":
    name='Zwcl2341'  # Change this to the desired cluster name
    path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/images/'
    out_path='/home/kincaid/Desktop/Saraswati_codes/'+name+'/catalogs/'
    #images=['image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored.0.fits','image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored.1.fits','image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored.2.fits']
    #images=['Zwcl2341_FIRST.fits.img.conv.fits']
    images=['image_DI_Clustered.DeeperDeconv.AP.int.restored.pbcor.fits']# image_DI_Clustered.DeeperDeconv.AP4.int.restored.pbcor.fits
    #out_names=['image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored_998','image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored_1283','image_DI_Clustered.DeeperDeconv.AP4.cube.int.restored_1569']
    out_names=name
    for i, image in enumerate(images):
        catalog_generation(path+image, out_path, out_names)
