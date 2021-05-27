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
    float sum;
    DAISGram res;
    res.data.init(data.rows(), data.cols(), data.depth());
    for (int i = 0; i < data.rows(); i++) {
        for (int j = 0; j < data.cols(); j++) {
                sum = (data(i, j, 0) + data(i, j, 1) + data(i, j, 2)) / 3;
                res.data(i, j, 0) = res.data(i, j, 1) = res.data(i, j, 2) = sum;
        }
    }

    res.save_image("prova4.bmp");
    DAISGram prova;
    prova.load_image("prova4.bmp");

    DAISGram nuovo;
    nuovo.load_image("results/dais_gray.bmp");
    cout<<"gray"<<(nuovo.data==prova.data)<<"\n";

    return res;
}

DAISGram DAISGram::brighten(float bright) {
    DAISGram res;
    res.data = (this->data) + bright;
    res.data.clamp(0, 255);

    DAISGram nuovo;
    nuovo.load_image("results/dais_brighten_20.bmp");
    cout<<"bright"<<(nuovo.data==res.data)<<"\n";

    return res;
}

DAISGram DAISGram::warhol() {
    Tensor rg, bg, rb, original;
    DAISGram res;
    original = rg = bg = rb = data;
    for (int z = 0; z < data.depth(); z++) {
        for (int x = 0; x < data.rows(); x++) {
            for (int y = 0; y < data.cols(); y++) {
                if (z == 0) {
                    rg(x, y, z) = original(x, y, z + 1);
                    rb(x, y, z) = original(x, y, z + 2);
                } else if (z == 1) {
                    rg(x, y, z) = original(x, y, z - 1);
                    bg(x, y, z) = original(x, y, z + 1);
                } else {
                    bg(x, y, z) = original(x, y, z - 1);
                    rb(x, y, z) = original(x, y, z - 2);
                }
            }
        }
    }
    original = original.concat(rg, 1);
    bg = bg.concat(rb, 1);
    res.data = original.concat(bg, 0);

    DAISGram nuovo;
    nuovo.load_image("results/dais_warhol.bmp");
    cout<<"warhol"<<(nuovo.data==res.data)<<"\n";

    return res;
}


DAISGram DAISGram::blend(const DAISGram &rhs, float alpha) {
    Tensor supp;
    DAISGram res;
    res.data = (data * alpha) + (rhs.data * (1 - alpha));

    res.save_image("prova1.bmp");
    DAISGram prova;
    prova.load_image("prova1.bmp");

    DAISGram nuovo;
    nuovo.load_image("results/blend/blend_0.75.bmp");
    cout<<"blend"<<(prova.data==nuovo.data)<<"\n";

    return res;
}

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

    res.data = this->data.convolve(f);
    res.data.clamp(0, 255);


    DAISGram nuovo;
    nuovo.load_image("results/dais_sharp.bmp");
    cout<<"sharp"<<(nuovo.data==res.data)<<"\n";

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
DAISGram DAISGram::emboss() {
    DAISGram res;
    Tensor f(3, 3, 3);
    int max = f.depth() * f.cols() * f.rows();
    int cf = 0;
    for (int i = 0; i < max; i++) {
        if (cf == 9) {
            cf = 0;
        }
        if (cf % 2 == 0) {
            if (cf == 0) {
                f.at(i) = -2;
            }
            if (cf == 2 || cf == 6) {
                f.at(i) = 0;
            }
            if (cf == 4) {
                f.at(i) = 1;
            }
            if (cf == 8) {
                f.at(i) = 2;
            }
        }else {
            if (cf == 1 || cf == 3) {
                f.at(i) = -1;
            } else {
                f.at(i) = 1;
            }
        }
        cf++;
    }
    res.data = this->data.convolve(f);
    res.data.clamp(0, 255);
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
DAISGram DAISGram::edge() {
    DAISGram res;
    Tensor f(3, 3, 3);
    int max = f.depth() * f.cols() * f.rows();
    int cf = 0;
    for (int i = 0; i < max; i++) {
        if(cf==9)
            cf=0;
        if (cf != 4) {
            f.at(i) = -1;
        } else {
            f.at(i) = 8;
        }
        cf++;
    }

    *this=this->grayscale();
    res.data=this->data.convolve(f);
    res.data.clamp(0,255);


    res.save_image("prova2.bmp");
    DAISGram prova;
    prova.load_image("prova2.bmp");


    DAISGram nuovo;
    nuovo.load_image("results/dais_edge.bmp");
    cout<<"edge"<<(nuovo.data==prova.data)<<"\n";

    return res;
}

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
DAISGram DAISGram::smooth(int h) {
    DAISGram res;
    float c = (1.0 / (h * h));
    Tensor f(h,h,h,c);
    res.data = data.convolve(f);


    res.save_image("prova3.bmp");
    DAISGram prova;
    prova.load_image("prova3.bmp");



    DAISGram nuovo;
    nuovo.load_image("results/dais_smooth_3.bmp");
    cout<<"smooth"<<(nuovo.data==prova.data)<<"\n";

    return res;
}


DAISGram DAISGram::greenscreen(DAISGram & bkg, int rgb[], float threshold[]){
    DAISGram res;
    res.data=this->data;
    int val=bkg.getCols()*bkg.getRows();
    int s_min[3];
    int s_max[3];
    for(int i=0; i<3; i++){
        s_min[i] = rgb[i] - threshold[i];
        s_max[i] = rgb[i] + threshold[i];
    }
    int val2=val+val;
    for(int i=0;i<val;i++) {
        if(res.data.at(i)>=s_min[0] && res.data.at(i)<=s_max[0])
            if(res.data.at(i+val)>=s_min[1] && res.data.at(i+val)<=s_max[1])
                if(res.data.at(i+val2)>=s_min[2] && res.data.at(i+val2)<=s_max[2]) {
                    res.data.at(i)=bkg.data.at(i);
                    res.data.at(i+val)=bkg.data.at(i+val);
                    res.data.at(i+val2)=bkg.data.at(i+val2);
                }
    }
    DAISGram nuovo;
    nuovo.load_image("results/greenscreen/seba_flower.bmp");
    cout<<"greenscale"<<(nuovo.data==res.data)<<"\n";

    return res;
}

DAISGram DAISGram::equalize() { // DA TERMINARE !!
    DAISGram res;
    res.data.init(data.rows(),data.cols(),data.depth());
    int len = 256;
    int c;
    int add;
    for (int depth{}; depth < data.depth(); depth++) {
        float cdf[256]={};
        float sc = {};
        float save;
        add = data.rows()*data.cols() * depth;
        c=0;
        while (c < data.rows()*data.cols()) { //occorrenze fatte //ciclo su matrice
            cdf[(int)data.at(c + add)]++;
            c++;
        }
        c = 0;
        while (c < len) { //creo cdf pesato
            if(cdf[c]!=0) {
                save = cdf[c];
                cdf[c] = save + sc;
                sc += save;
            }
            c++;
        }
        c = 0;
        while (cdf[c] == 0) {
            c++;
        }
        int cdf_min = cdf[c];
        int den = (res.getRows() * res.getCols()) -(cdf_min);
        for (int i = 0; i < len; i++) { //normalizzato cdf
            cdf[i] = (int)((((cdf[i] - cdf_min) / (den) ) * (len - 1)));
        }
        for (int i = 0; i < res.data.rows()*res.data.cols(); i++) {
            res.data.at(i+add) = cdf[(int)data.at(i+add)];
        }
    }

    DAISGram nuovo;
    nuovo.load_image("results/dais_equalize.bmp");
    cout<<(nuovo.data==res.data)<<"equalize"<<"\n";

    return res;
}

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
