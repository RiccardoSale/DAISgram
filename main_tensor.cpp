#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tensor.h"

void show_help(){
    printf("*** Tensors Operations ***\n");
    printf("\targ 1: input file name (tensor1) \n");
    printf("\targ 2: input file name (tensor2) \n");
    printf("\targ 3: operazione da effettuare (+,-,d,x,convolve, concat)\n");
    printf("\targ 4: output file name\n");
    printf("\targ 5: parametro axis della concat\n");
    printf("\n");
}

int main() {
    Tensor rescale(3,3,3);
    rescale.init_random(60,5);
    cout<<rescale<<"\n"<<"\n";
    rescale.rescale();
    cout<<rescale<<"\n"<<"\n";
    cout<<"inizio"<<"\n";
    Tensor s1(1,1,1,10);
    Tensor t;
    //t=t+5;
    cout<<t;
    cout<<"fine"<<"\n";
    Tensor p(3, 3, 3);
    p.init_random(60,5);
    cout<<p(1,2,1)<<"ELEM";
    cout<<p.getMax(1)<<"MAX"<<"\n";
    cout << p;
    cout<<"\n"<<"somma"<<"\n";
    Tensor somma=p+p;
    cout<< somma;
    cout<<"\n"<<"clamp"<<"\n";
    somma.clamp(10,130);
    cout<<somma;
    cout<<"\n"<<"padding"<<"\n";
    Tensor pad =somma.padding(3,2);
    cout<<pad;
    cout<<"\n"<<"SUBSET"<<"\n";
    Tensor sub=p.subset(0,1,0,3,0,3);
    cout<<sub;
    sub.write_file("prova.txt");
    Tensor r;
    r.read_file("prova.txt");
    cout<<"\n"<<"\n"<<"STAMPA TENSORE FILE"<<"\n";
    cout<<r;
    cout<<"\n"<<"CONCAT"<<"\n";
    Tensor concat(90,90,3,10);
    Tensor c2(90,90,3,5);
    Tensor res=c2.concat(concat,0);
    res.showSize();

    /*PROVE PER TESTARE LE ECCEZIONI DEI TENSORI A NULL*/
    cout<<"\n"<<"ECCEZIONI DEI TENSORI A NULL"<<"\n";
    Tensor t1;
    Tensor t2(3,3,3, 9);
    //Tensor t3 (t1); --> copy constructor --> funziona

    //cout  << (t1 == t2) ; --> operatore == --> funziona
    //cout  << (t2 == t1) ; --> funziona

    //Tensor t3 = t1 - t2; --> operatore - tra tensori --> funziona
    //Tensor t3 = t2 - t1; --> funziona

    //Tensor t3 = t1 + t2; --> operatore + tra tensori --> funziona
    //Tensor t3 = t2 + t1; --> funziona

    //Tensor t3 = t1 * t2; --> operatore * tra tensori --> funziona
    //Tensor t3 = t2 * t1; --> funziona

    //Tensor t3 = t1 / t2; --> operatore / tra tensori --> funziona
    //Tensor t3 = t2 / t1; --> funziona

    float number = 6;

    //Tensor t3 = t1 - number; --> operatore - --> funziona

    //Tensor t3 = t1 + number; --> operatore + --> funziona

    //Tensor t3 = t1 * number; --> operatore * --> funziona

    //Tensor t3 = t1 / number; --> operatore / --> funziona

    //Tensor t3 = t1; --> operatore = --> funziona


    //Tensor t3 = t1.concat(t2, 0); --> metodo concat --> funziona
    //Tensor t3 = t2.concat(t1, 0); --> funziona
    
    
    //cout << t3;
}

//int main (int argc, char * argv[]) {

    //char * fn_in_1;  /* file 1 */
    //char * fn_in_2;  /* file 2 */
    //char * operation; /* operazione da eseguire */
    //char * fn_out; /* output file */

    //int axis = 0; /* axis for concat */

    /* variabili di appoggio per le computazioni */
    //Tensor a,b,out;

    //if(argc<4){
    //    show_help();
    //    return 0;
    //}

    //fn_in_1 = argv[1];  /* file 1 */
    //fn_in_2 = argv[2];  /* file 2 */
    //operation = argv[3]; /* operazione da eseguire */
    //fn_out = argv[4]; /* output file */

    //if(argc>5) {
    //    axis = atoi(argv[5]);
    //}

    //a.read_file(fn_in_1);
    //b.read_file(fn_in_2);
    /*
    if (strcmp(operation, "+") == 0) {
        out=a+b; 
    }else if(strcmp(operation, "-") == 0) {
        out=a-b; 
    }else if(strcmp(operation, "x") == 0) {
        out=a*b; 
    }else if(strcmp(operation, "convolve") == 0) {
        //out=a.convolve(b);
    }else if(strcmp(operation, "concat") == 0) {
        //out=a.concat(b,axis);
    }else if(strcmp(operation, "d") == 0) {
        out=a/b; 
    }else {
        throw(unknown_operation());
    }

    //out.write_file(fn_out);

    return 0; /* ciao a tutti!*/
//}
