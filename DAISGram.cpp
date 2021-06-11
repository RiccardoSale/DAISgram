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
    return res;
}

DAISGram DAISGram::brighten(float bright) {
    DAISGram res;
    res.data = (this->data) + bright;
    res.data.clamp(0, 255);
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
    return res;
}


DAISGram DAISGram::blend(const DAISGram &rhs, float alpha) {
    Tensor supp;
    DAISGram res;
    res.data = (data * alpha) + (rhs.data * (1 - alpha));
    return res;
}

DAISGram DAISGram::sharpen() {
    DAISGram res;
    int arr[9]={0,-1,0,-1,5,-1,0,-1,0};
    Tensor f(3,3,3,arr);
    res.data = this->data.convolve(f);
    res.data.clamp(0, 255);
    return res;
}

DAISGram DAISGram::emboss() {
    DAISGram res;
    int arr[9]={-2,-1,0,-1,1,1,0,1,2};
    Tensor f(3,3,3,arr);
    res.data = this->data.convolve(f);
    res.data.clamp(0, 255);
    return res;
}

DAISGram DAISGram::edge() {
    DAISGram res;
    int arr[9]={-1,-1,-1,-1,8,-1,-1,-1,-1};
    Tensor f(3,3,3,arr);
    *this=this->grayscale();
    res.data=this->data.convolve(f);
    res.data.clamp(0,255);
    return res;
}

DAISGram DAISGram::smooth(int h) {
    if(h%2!=0) {
        DAISGram res;
        float c = (1.0 / (h * h));
        Tensor f(h, h, h, c);
        res.data = data.convolve(f);
        return res;
    }else throw (filter_odd_dimensions());
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
    return res;
}

DAISGram DAISGram::equalize() {
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
        while (cdf[c] == 0) c++;
        int cdf_min = cdf[c];
        int den = (res.getRows() * res.getCols()) -(cdf_min);
        for (int i = 0; i < len; i++) { //normalizzato cdf
            cdf[i] = (int)((((cdf[i] - cdf_min) / (den) ) * (len - 1)));
        }
        for (int i = 0; i < res.data.rows()*res.data.cols(); i++) {
            res.data.at(i+add) = cdf[(int)data.at(i+add)];
        }
    }

    return res;
}

void DAISGram::generate_random(int h, int w, int d) {
    data = Tensor(h, w, d, 0.0);
    data.init_random(128, 50);
    data.rescale(255);
}
