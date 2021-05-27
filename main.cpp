#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tensor.h"
#include "libbmp.h"
#include "DAISGram.h"

void show_help() {
    printf("*** DAISGram ***\n");
    printf("\targ 1: input file name (img1) \n");
    printf("\targ 2: input file name (img2) \n");
    printf("\targ 3: operazione da effettuare (gray, brighten, blend, sharp, edge, emboss, smooth, warhol, equalize, chromakey) \n");
    printf("\targ 4: output file name\n");
    printf("\targ 5: Diversi significati in funzione dell'operazione (default 3):\n"
           "\t\t- [smooth]: kernel size \n"
           "\t\t- [brighten]: valore bright per aumentare la luminosità \n"
           "\t\t\n");
    printf("\targ 6: Diversi significati in funzione dell'operazione (default 1.0):\n"
           "\t\t- [blend] parametro alpha per il blending di due immagini");
    printf("\n");
}

int main(int argc, char *argv[]) {

    DAISGram bright;
    bright.load_image("dais.bmp");
    bright=bright.brighten(20);
    bright.save_image("bright.bmp");

    DAISGram gray;
    gray.load_image("dais.bmp");
    gray=gray.grayscale();
    gray.save_image("gray.bmp");

    DAISGram war;
    war.load_image("dais.bmp");
    war = war.warhol();
    war.save_image("warhol.bmp");

    DAISGram blend;
    blend.load_image("blend_a.bmp");
    DAISGram blend2;
    blend2.load_image("blend_b.bmp");
    blend=blend.blend(blend2,0.75);
    blend.save_image("blend.bmp");


    DAISGram sharp;
    sharp.load_image("dais.bmp");
    sharp = sharp.sharpen();
    sharp.save_image("sharp.bmp");

    DAISGram edge;
    edge.load_image("dais.bmp");
    edge = edge.edge();
    edge.save_image("edge.bmp");

    DAISGram emboss;
    emboss.load_image("blend_a.bmp");
    emboss = emboss.emboss();
    emboss.save_image("fiore_emboss.bmp");

    DAISGram smooth;
    smooth.load_image("dais.bmp");
    smooth = smooth.smooth(3);
    smooth.save_image("smooth.bmp");

    DAISGram green;
    green.load_image("gs_4.bmp");
    DAISGram green2;
    green2.load_image("gs_4_bkg.bmp");
    int rgb[3]={226,225,220};
    float tresh[3]={50,50,50};
    green=green.greenscreen(green2,rgb,tresh);
    green.save_image("greenscale.bmp");

    DAISGram equalize;
    equalize.load_image("dais.bmp");
    equalize=equalize.equalize();
    equalize.save_image("equalizzata.bmp");
}

