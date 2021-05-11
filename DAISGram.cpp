#include <iostream>
#include <string>

#include "dais_exc.h"
#include "tensor.h"
#include "libbmp.h"
#include "DAISGram.h"

using namespace std;

/**
 * Load a bitmap from file
 *
 * @param filename String containing the path of the file
 */

void DAISGram::load_image(string filename){
    BmpImg img = BmpImg();
    img.read(filename.c_str());
    const int h = img.get_height();
    const int w = img.get_width();

    data = Tensor(h, w, 3, 0.0);

    for(int i=0;i<img.get_height();i++){
        for(int j=0;j<img.get_width();j++){
            data(i,j,0) =(float) img.red_at(j,i);
            data(i,j,1) = (float) img.green_at(j,i);
            data(i,j,2) = (float) img.blue_at(j,i);
        }
    }
}


/**
 * Save a DAISGram object to a bitmap file.
 *
 * Data is clamped to 0,255 before saving it.
 *
 * @param filename String containing the path where to store the image.
 */

void DAISGram::save_image(string filename){

    data.clamp(0,255);

    BmpImg img = BmpImg(getCols(), getRows());

    img.init(getCols(), getRows());

    for(int i=0;i<getRows();i++){
        for(int j=0;j<getCols();j++){
            img.set_pixel(j,i,(unsigned char) data(i,j,0),(unsigned char) data(i,j,1),(unsigned char) data(i,j,2));
        }
    }

    img.write(filename);

}

        int DAISGram::getRows(){
            return data.rows();
        }

        int DAISGram::getCols(){
            return data.cols();
        }

        int DAISGram::getDepth(){
            return data.depth();
        }

/**
         * Create a grayscale version of the object
         *
         * A grayscale image is produced by substituting each pixel with its average on all the channel
         *
         * @return returns a new DAISGram containing the modified object
         */
        DAISGram DAISGram::grayscale(){

            int i_max = data.rows() * data.cols() * data.depth();
            float sum = 0;
            DAISGram result;

            for(int i=0; i < data.rows(); i++){
                for (int j = 0; j < data.cols(); j++){
                    for (int k = 0; k < data.depth(); k++){
                        sum += data(i, j, k);
                    }
                }
            }

            float average = sum / i_max;

            result.data.init(data.rows(), data.cols(), data.depth(), average);

            return result;

        }

        DAISGram DAISGram::brighten(float bright) {
            DAISGram result;
            //result.data= (this->data)+bright;
            //result.data(this->data.rows(),this->data.cols(),this->data.depth());
            //cout<<result.data;
            result.data.clamp(0,255);
            return result;
        }


    /**
     * Create a Warhol effect on the image
     * 
     * This function returns a composition of 4 different images in which the:
     * - top left is the original image
     * - top right is the original image in which the Red and Green channel are swapped
     * - bottom left is the original image in which the Blue and Green channel are swapped
     * - bottom right is the original image in which the Red and Blue channel are swapped
     *  
     * The output image is twice the dimensions of the original one.
     * 
     * @return returns a new DAISGram containing the modified object
     */
    DAISGram DAISGram::warhol(){
        Tensor verde, rosso, blu, originale;
        originale = Tensor(data);
        verde = Tensor(data);
        for(int i=0; i<data.rows(); i++){
            for(int j=0; j<data.cols(); j++){
                for(int k=0; k<data.depth(); k++){
                    if(data.rows() == 0){
                        verde(i, j, k) = originale(i+1, j, k);
                        blu(i, j, k) = originale(i+2, j, k);
                    }else if(data.rows() == 1){
                        rosso(i, j, k) = originale(i+1, j, k);
                        verde(i, j, k) = originale(i-1, j, k);
                    }else{
                        rosso(i, j, k) = originale(i-1, j, k);
                        blu(i, j, k) = originale(i-2, j, k);
                    }
                }
            }
        }
        DAISGram risultato;
        originale.concat(verde, 0);     //concatenazione delle due immagini una vicino all'altra
        rosso.concat(blu, 0);           //concatenazione delle due immagini una vicino all'altra
        originale.concat(rosso, 1);     //concatenazione delle due immagini una sotto l'altra

        risultato.data = originale;
        return risultato;
    };


/**
 * Generate Random Image
 *
 * Generate a random image from nois
 *
 * @param h height of the image
 * @param w width of the image
 * @param d number of channels
 * @return returns a new DAISGram containing the generated image.
 */
void DAISGram::generate_random(int h, int w, int d){
    data = Tensor(h,w,d,0.0);
    data.init_random(128,50);
    data.rescale(255);
}
