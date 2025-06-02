from PIL import Image
import matplotlib.pyplot as plt



def plot():
        path1='/home/kincaid/Desktop/Saraswati_codes/Zwcl2341/images/cutouts/'
        path2='/home/kincaid/Desktop/Saraswati_codes/A2631/images/cutouts/'
        optical_image1 = path2+'rgb_cutout_2056_2057.jpg'       
        optical_image2 = path1+'rgb_cutout_394_395.jpg'
        optical_image3 = path1+'rgb_cutout_189_190.jpg'
        optical_image4 = path1+'rgb_cutout_triple_source_1161_1162_1163.jpg'
        
        cutouts = [(optical_image1,5,2,1),
                (optical_image2,5,2,2),
                (optical_image3,5,2,3),
                (optical_image4,5,2,4),]
        
        
        fig = plt.figure(figsize=(8, 16))
        plt.subplots_adjust(wspace=0.05, hspace=0.05, right=0.85,bottom=0.05) 

        for  image,r,c,p in cutouts:
            print(r,c,p)
            ax = plt.subplot(r,c,p)
            img= Image.open(image)
            ax.imshow(img)
            ax.axis('off')
        plt.tight_layout()
        plt.subplots_adjust(wspace=0.05, hspace=0.05, right=0.85, bottom=0)
        plt.savefig('plots/optical_radio_cutouts.png', bbox_inches='tight', pad_inches=0,dpi=300)
        plt.show()


if __name__ == '__main__':
    output_path = '/home/kincaid/Desktop/Saraswati_codes/A2631/plots/'  
    image_path =  '/home/kincaid/Desktop/Saraswati_codes/'
    plot( )     