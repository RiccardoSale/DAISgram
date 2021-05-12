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


DAISGram blend(const DAISGram & rhs, float alpha=0.5){
    Tensor supp;
    DAISGram risultato;
    risultato.data = this.operators*(alpha);
    supp = rhs.operators*(1-alpha);
    risultato.data.operator+(supp);
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
void DAISGram::generate_random(int h, int w, int d) {
    data = Tensor(h, w, d, 0.0);
    data.init_random(128, 50);
    data.rescale(255);
}
