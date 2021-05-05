#include <iostream>
#include "tensor.h"
int main() {
    Tensor prova(3,3,3);
    prova.init_random(40,20);
    //Tensor due=prova-prova;
    //due.stampa();
    prova.stampa();
    prova.clamp(10,60);
    prova.stampa();
    Tensor pad=prova.padding(3,2);
    pad.stampa();
}
