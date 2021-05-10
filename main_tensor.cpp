#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tensor.h"

void show_help(){
    printf("*** Tensors Operations ***\n");
    printf("\targ 1: input file name (tensor1) \n");
    printf("\targ 2: input file name (tensor2) \n");
    printf("\targ 3: operazione da effettuare (+,-,/,*,convolve, concat)\n");
    printf("\targ 4: output file name\n");
    printf("\targ 5: parametro axis della concat\n");
    printf("\n");
};

int main() {
    cout<<"inizio"<<"\n";
    Tensor p(3, 3, 3);
    p.init_random(60,5);

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

    /*Proviamo la concat*/
    cout<<"\n"<<"CONCAT"<<"\n";
    Tensor p1(3, 4, 3);
    p1.init_random();
    cout<< p1 <<endl << endl;

    Tensor p2(3, 3, 3);
    p2.init(3, 3, 3, 50);
    cout<< p2 <<endl << endl;

    Tensor result = p1.concat(p2, 1);
    cout<<result <<endl;

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
    }else if(strcmp(operation, "*") == 0) {
        out=a*b; 
    }else if(strcmp(operation, "convolve") == 0) {
        //out=a.convolve(b);
    }else if(strcmp(operation, "concat") == 0) {
        //out=a.concat(b,axis);
    }else if(strcmp(operation, "/") == 0) {
        out=a/b; 
    }else {
        throw(unknown_operation());
    }

    //out.write_file(fn_out);

    return 0; /* ciao a tutti!*/
//}
