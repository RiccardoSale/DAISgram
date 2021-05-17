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

void DAISGram::load_image(string filename) {
    BmpImg img = BmpImg();

    img.read(filename.c_str());
    const int h = img.get_height();
    const int w = img.get_width();

    data = Tensor(h, w, 3, 0.0);

    for (int i = 0; i < img.get_height(); i++) {
        for (int j = 0; j < img.get_width(); j++) {
            data(i, j, 0) = (float) img.red_at(j, i);
            data(i, j, 1) = (float) img.green_at(j, i);
            data(i, j, 2) = (float) img.blue_at(j, i);
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

void DAISGram::save_image(string filename) {

    data.clamp(0, 255);

    BmpImg img = BmpImg(getCols(), getRows());

    img.init(getCols(), getRows());

    for (int i = 0; i < getRows(); i++) {
        for (int j = 0; j < getCols(); j++) {
            img.set_pixel(j, i, (unsigned char) data(i, j, 0), (unsigned char) data(i, j, 1),
                          (unsigned char) data(i, j, 2));
        }
    }

    img.write(filename);

}

int DAISGram::getRows() {
    return data.rows();
}

int DAISGram::getCols() {
    return data.cols();
}

int DAISGram::getDepth() {
    return data.depth();
}

/**
         * Create a grayscale version of the object
         *
         * A grayscale image is produced by substituting each pixel with its average on all the channel
         *
         * @return returns a new DAISGram containing the modified object
         */
DAISGram DAISGram::grayscale() {
    float sum{};
    DAISGram result;
    result.data.init(data.rows(), data.cols(), data.depth());
    for (int i = 0; i < data.rows(); i++) {
        for (int j = 0; j < data.cols(); j++) {
            for (int k = 0; k < 1; k++) {
                sum = (data(i, j, 0) + data(i, j, 1) + data(i, j, 2)) / 3;
                result.data(i, j, 0) = sum;
                result.data(i, j, 1) = sum;
                result.data(i, j, 2) = sum;
            }
        }
    }

    return result;
}

DAISGram DAISGram::brighten(float bright) {
    DAISGram result;
    result.data = (this->data) + bright;
    result.data.clamp(0, 255);
    return result;
}

DAISGram DAISGram::warhol() {
    Tensor rg, bg, rb,originale;
    DAISGram risultato;
    originale=data;
    rg=data;
    bg=data;
    rb=data;

    for (int z = 0; z < data.depth(); z++) {
        for (int x = 0; x < data.rows(); x++) {
            for (int y = 0; y < data.cols(); y++) {
                if (z == 0) {
                    rg(x, y, z) = originale(x, y, z+1);
                    rb(x, y, z) = originale(x, y, z+2);
                } else if (z  == 1) {
                    rg(x, y, z) = originale(x, y, z-1);
                    bg(x, y, z) = originale(x, y, z+1);
                } else {
                    bg(x, y, z) = originale(x, y, z-1);
                    rb(x, y, z) = originale(x, y, z-2);
                }
            }
        }
    }

    originale=originale.concat(rg, 1);
    bg=bg.concat(rb, 1);
    risultato.data = originale.concat(bg,0);
    return risultato;
}

//risultato ->alpha*data + rhs.data*(1-alpha)
DAISGram DAISGram::blend(const DAISGram & rhs, float alpha){
    Tensor supp;
    DAISGram risultato;
    risultato.data = data*(alpha) + (rhs.data*(1-alpha));
    return risultato;
};

DAISGram DAISGram::sharpen() {
    DAISGram res;
    Tensor f(3, 3, 3);
    int max = f.depth() * f.cols() * f.rows();
    int cf = 0;
    for (int i = 0; i < max; i++) {
        if (cf == 9)
            cf = 0;
        if (cf % 2 == 0) {
            if (cf == 4)
                f.at(i) = 5;
            else
                f.at(i) = 0;
        } else
            f.at(i) = -1;
        cf++;
    }
    cout<<f;
    res.data = this->data.convolve(f);
    res.data.clamp(0,255);
    return res;
}

/**
         * Emboss the image
         * 
         * This function makes the image embossed (a light 3D effect) by convolving it with an
         * embossing filter
         * 
         * filter[3][3]
         *    -2 -1  0
         *    -1  1  1
         *     0  1  2
         * 
         * Before returning the image, the corresponding tensor should be clamped in [0,255]
         *  
         * @return returns a new DAISGram containing the modified object
         */
DAISGram DAISGram::emboss(){
    DAISGram res;
    Tensor f(3, 3, 3);
    int max = f.depth() * f.cols() * f.rows();
    int cf = 0;
    for (int i = 0; i < max; i++) {
        if (cf == 9){
            cf = 0;
        }
        if(cf%2 == 0){
            if(cf == 0){
                f.at(i) = -2;
            }

            if(cf == 2 || cf == 6){
                f.at(i) = 0;
            }

            if(cf == 4){
                f.at(i) = 1;
            }

            if(cf == 8){
                f.at(i) = 2;
            }
        }else{
            if(cf == 1 || cf == 3){
                f.at(i) = -1;
            }else{// cf == 5 || cf == 7
                f.at(i) = 1;
            }
        }
    
        cf++;
    }
    cout<<f;
    res.data = this->data.convolve(f);
    res.data.clamp(0,255);
    return res;

}

/**
         * Edges of an image
         * 
         * This function extract the edges of an image by using the convolution 
         * operator and the following filter
         * 
         * 
         * filter[3][3]
         * -1  -1  -1
         * -1   8  -1
         * -1  -1  -1
         * 
         * Remeber to convert the image to grayscale before running the convolution.
         * 
         * Before returning the image, the corresponding tensor should be clamped in [0,255]
         *  
         * @return returns a new DAISGram containing the modified object
         */  
        DAISGram DAISGram::edge(){
            DAISGram res;
            Tensor f(3, 3, 3);
            int max = f.depth() * f.cols() * f.rows();
            int cf = 0;
            for(int i=0; i<max; i++){
                if(cf != 4){
                    f.at(i) = -1;
                }else{
                    f.at(i) = 8;
                }
                cf++;
            }
            cout<<f;
            this->grayscale();
            res.data = this->data.convolve(f);
            res.data.clamp(0,255);
            return res;
        };

/**
         * Smooth the image
         * 
         * This function remove the noise in an image using convolution and an average filter
         * of size h*h:
         * 
         * c = 1/(h*h)
         * 
         * filter[3][3]
         *    c c c
         *    c c c
         *    c c c
         *  
         * @param h the size of the filter
         * @return returns a new DAISGram containing the modified object
         */
DAISGram DAISGram::smooth(int h){
    DAISGram res;
    Tensor f(h, h, h);
    int max = f.depth() * f.cols() * f.rows();
    
    float c = (float) 1 / (h*h) ;

    for (int i = 0; i < max; i++) {
        f.at(i) = c;
    }

    cout<<f;
    res.data = this->data.convolve(f);
    
    return res;
}

/**
         * Green Screen
         * 
         * This function substitutes a pixel with the corresponding one in a background image 
         * if its colors are in the surrounding (+- threshold) of a given color (rgb).
         * 
         * (rgb - threshold) <= pixel <= (rgb + threshold)
         * 
         * 
         * @param bkg The second image used as background
         * @param rgb[] The color to substitute (rgb[0] = RED, rgb[1]=GREEN, rgb[2]=BLUE) 
         * @param threshold[] The threshold to add/remove for each color (threshold[0] = RED, threshold[1]=GREEN, threshold[2]=BLUE) 
         * @return returns a new DAISGram containing the result.
         */  
        DAISGram DAISGram::greenscreen(DAISGram & bkg, int rgb[], float threshold[]){};
        /*    DAISGram res;
            float soglia_min[3];
            float soglia_max[3];
            for(int i=0; i<2; i++){
                soglia_min[i] = rgb[i] - threshold[i];
                soglia_max[i] = rgb[i] + threshold[i];
            }
                
            for (int z = 0; z < data.depth(); z++) {
                for (int x = 0; x < data.rows(); x++) {
                    for (int y = 0; y < data.cols(); y++) {
                        if((data(x, y, z) >= soglia_min[z]) && (data(x, y, z) <= (soglia_max[z]))){
                            res.data(x, y, z) = bkg.data(x, y, z);
                        }
                    }
                }
            }
        };*/

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
void DAISGram::generate_random(int h, int w, int d) {
    data = Tensor(h, w, d, 0.0);
    data.init_random(128, 50);
    data.rescale(255);
}
