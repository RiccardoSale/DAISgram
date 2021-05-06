#include <iostream>
#include "tensor.h"
int main() {
    Tensor prova(2,3,3);
    prova.init_random(40,20); //vedere perche ultimo elem sempre = 0 problema nelle () ???
    prova.stampa();
    cout<<prova(1,2,2);// voglio riga "1" colonna "3" dim "3" //sembra andare
    cout<<"\n";
    /*
    Tensor subset=prova.subset(0,3,0,2,0,3);
    subset.stampa();
    //Tensor due=prova+prova;
    //due.stampa();
    */
    cout<<"\n";
    cout<<"PROVA CLAMP"<<"\n";
    prova.clamp(10,60);
    prova.stampa();
    cout<<"\n";
    cout<<"PROVA PADD"<<"\n";
    Tensor pad=prova.padding(3,4);
    pad.stampa();

}
