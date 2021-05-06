#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.h"

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;

/*
 int at(int i,int j,int k){
     return k*r*c+(i*c+j);
     //j*r salto tot righe
 };
 */



/**
* Class constructor
*
* Parameter-less class constructor
*/
Tensor::Tensor(){
    data = nullptr;
    this->r = 0;
    this->c = 0;
    this->d = 0;
}
/**
* Class constructor
*
* Creates a new tensor of size r*c*d initialized at value v
*
* @param r
* @param c
* @param d
* @param v
* @return new Tensor
*/
Tensor::Tensor(int r, int c, int d, float v = 0.0){
    this->r=r;
    this->c=c;
    this->d=d;
    int i_max{r*c*d};
    data=new float[i_max];
    for(int i=0;i<i_max;i++){
        data[i]=v;
    }
};
/**
 * Class distructor
 *
 * Cleanup the data when deallocated
 */
Tensor::~Tensor(){
    delete data;
};

void init_progressive();

/**
 * Operator oveloding ()
 *
 * if indexes are out of bound throw index_out_of_bound() exception
 *
 * @return the value at location [i][j][k]
 */
float Tensor::operator()(int i, int j, int k) const{
    return data[k*r*c+(i*c+j)];
    //TODO EXCEPTION
};

/**
 * Operator overloading ()
 *
 * Return the pointer to the location [i][j][k] such that the operator (i,j,k) can be used to
 * modify tensor data.
 *
 * If indexes are out of bound throw index_out_of_bound() exception
 *
 * @return the pointer to the location [i][j][k]
 */
float& Tensor::operator()(int i, int j, int k){
    float& res=data[k*r*c+(i*c+j)];
    return res;
    //TODO EXCEPTION
};

/**
 * Copy constructor
 *
 * This constructor copy the data from another Tensor
 *
 * @return the new Tensor
 */
Tensor::Tensor(const Tensor& that){
    r=that.r;
    c=that.c;
    d=that.d;
    int i_max{r*c*d};
    data=new float[i_max];
    for(int i=0;i<i_max;i++) {
        data[i] = that.data[i];
    }
};
/**
 * Operator overloading -
 *
 * It performs the point-wise difference between two Tensors.
 *
 * result(i,j,k)=this(i,j,k)-rhs(i,j,k)
 *
 * The two tensors must have the same size otherwise throw a dimension_mismatch()
 *
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator-(const Tensor &rhs){
    if(r == rhs.r && c == rhs.c && d == rhs.d){
        Tensor result(r, c, d);
        int i_max{r*c*d};
        for(int i=0;i<i_max;i++) {
            result.data[i] = data[i] - rhs.data[i];
        }

        return result;
    }else{
        throw(dimension_mismatch());
    }
};

/**
 * Operator overloading +
 *
 * It performs the point-wise sum between two Tensors.
 *
 * result(i,j,k)=this(i,j,k)+rhs(i,j,k)
 *
 * The two tensors must have the same size otherwise throw a dimension_mismatch()
 *
 * @return returns a new Tensor containing the result of the operation
*/
Tensor Tensor::operator +(const Tensor &rhs){
    if(r == rhs.r && c == rhs.c && d == rhs.d){
        Tensor result(r, c, d);
        int i_max{r*c*d};
        for(int i=0;i<i_max;i++) {
            result.data[i] = data[i] + rhs.data[i];
        }

        return result;
    }else{
        throw(dimension_mismatch());
    }
};

/**
 * Operator overloading *
 *
 * It performs the point-wise product between two Tensors.
 *
 * result(i,j,k)=this(i,j,k)*rhs(i,j,k)
 *
 * The two tensors must have the same size otherwise throw a dimension_mismatch()
 *
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator*(const Tensor &rhs){
    if(r == rhs.r && c == rhs.c && d == rhs.d){
        Tensor result(r, c, d);
        int i_max{r*c*d};
        for(int i=0;i<i_max;i++) {
            result.data[i] = data[i] * rhs.data[i];
        }

        return result;
    }else{
        throw(dimension_mismatch());
    }
};

/**
 * Operator overloading /
 *
 * It performs the point-wise division between two Tensors.
 *
 * result(i,j,k)=this(i,j,k)/rhs(i,j,k)
 *
 * The two tensors must have the same size otherwise throw a dimension_mismatch()
 *
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator /(const Tensor &rhs){
    if(r == rhs.r && c == rhs.c && d == rhs.d){
        Tensor result(r, c, d);
        int i_max{r*c*d};
        for(int i=0;i<i_max;i++) {
            result.data[i] = data[i] / rhs.data[i];
        }

        return result;
    }else{
        throw(dimension_mismatch());
    }
};

/**
 * Operator overloading -
 *
 * It performs the point-wise difference between a Tensor and a constant
 *
 * result(i,j,k)=this(i,j,k)-rhs
 *
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator-(const float &rhs){
    Tensor result(r, c, d);
    int i_max{r*c*d};
    for(int i=0;i<i_max;i++) {
        result.data[i] = data[i]-rhs;
    }

    return result;
};

/**
 * Operator overloading +
 *
 * It performs the point-wise sum between a Tensor and a constant
 *
 * result(i,j,k)=this(i,j,k)+rhs
 *
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator+(const float &rhs){
    Tensor result(r, c, d);
    int i_max{r*c*d};
    for(int i=0;i<i_max;i++) {
        result.data[i] = data[i] + rhs;
    }

    return result;
};

/**
 * Operator overloading *
 *
 * It performs the point-wise product between a Tensor and a constant
 *
 * result(i,j,k)=this(i,j,k)*rhs
 *
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator*(const float &rhs){
    Tensor result(r, c, d);
    int i_max{r*c*d};
    for(int i=0;i<i_max;i++) {
        result.data[i] = data[i] * rhs;
    }

    return result;
};

/**
 * Operator overloading / between a Tensor and a constant
 *
 * It performs the point-wise division between a Tensor and a constant
 *
 * result(i,j,k)=this(i,j,k)/rhs
 *
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator/(const float &rhs){
    Tensor result(r, c, d);
    int i_max{r*c*d};
    for(int i=0;i<i_max;i++) {
        result.data[i] = data[i] / rhs;
    }

    return result;
};

/**
 * Operator overloading = (assignment)
 *
 * Perform the assignment between this object and another
 *
 * @return a reference to the receiver object
 */

Tensor & Tensor::operator=(const Tensor &other){
    r=other.r;
    c=other.c;
    d=other.d;  //aggiorno
    return *this; //NON SICURO 100% di questa riga
};


/**
 * Random Initialization
 * 
 * Perform a random initialization of the tensor
 * 
 * @param mean The mean
 * @param std  Standard deviation
 */
void Tensor::init_random(float mean, float std){
    if(data){

        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean,std);

        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                for(int k=0;k<d;k++){
                    this->operator()(i,j,k)= distribution(generator);
                }
            }
        }    

    }else{
        throw(tensor_not_initialized());
    }
}

/**
 * Constant Initialization
 *
 * Perform the initialization of the tensor to a value v
 *
 * @param r The number of rows
 * @param c The number of columns
 * @param d The depth
 * @param v The initialization value
 */
void Tensor::init(int r, int c, int d, float v=0.0){};//aspettare documentazione maggiore

/**
 * Tensor Clamp
 *
 * Clamp the tensor such that the lower value becomes low and the higher one become high.
 *
 * @param low Lower value
 * @param high Higher value
 */
void Tensor::clamp(float low, float high){//POTREBBE ESSERE CHE DOBBIAMO IMPOSTARE SOLO IL MAX E IL MIN
    int i_max=r*c*d;
    for(int i=0;i<i_max;i++){
        float elem = data[i];
        if(elem < low)
            data[i] = low;
        else if(elem > high)
            data[i] = high;
    }
};

/**
 * Tensor Rescaling
 *
 * Rescale the value of the tensor following this rule:
 *
 * newvalue(i,j,k) = ((data(i,j,k)-min(k))/(max(k)-min(k)))*new_max
 *
 * where max(k) and min(k) are the maximum and minimum value in the k-th channel.
 *
 * new_max is the new value for the maximum
 *
 * @param new_max New maximum vale
 */
void Tensor:: rescale(float new_max=1.0){
    //trovo max value dim 1
    float max{data[0]};
    float min{data[0]};
    for(int z=0;z<d;z++) {
        for (int y = 0; y < c; y++) { // Da testare probabilmente scambiare c e r
            int i = this->at(0, y, z); 
            for (int x = 0; x < r; x++) {
                if (data[i + r] < min)
                    min = data[i + r];
                else if (data[i + r] > max)
                    max = data[i + r];
            }
        }//a questo punto dovrei aver trovato maggiore e minore per quella dimensione
        float diff = max - min;
        for (int y = 0; y < c; y++) {
            int i = this->at(0, y, z);
            for (int x = 0; x < r; x++) {
                data[i + r] = ((data[i + r] - min) / (diff)) * new_max;
            }
        }
    }
};

/**
 * Tensor padding
 *
 * Zero pad a tensor in height and width, the new tensor will have the following dimensions:
 *
 * (rows+2*pad_h) x (cols+2*pad_w) x (depth)
 *
 * @param pad_h the height padding
 * @param pad_w the width padding
 * @return the padded tensor
 */
Tensor Tensor::padding(int pad_h, int pad_w){
    int new_x = c+2*pad_w; //nuovo numero colonne
    int new_y = r+2*pad_h; //nuovo numero righe
    Tensor res(new_y, new_x, d);
    int i=0;
    int ii=0;
    for(int z=0; z<d; z++){ //scorro dimensioni
        for(int y=0; y<new_y; y++){
            for(int x=0; x<new_x; x++){
                if(x < pad_w || x >= c+pad_w || y < pad_h || y >= r+pad_h){
                    res.data[ii++]=0.0;
                }else{
                    res.data[ii++]=data[i++];
                }
            }
        }
    }
    return res;
};

/**
 * Subset a tensor
 *
 * retuns a part of the tensor having the following indices:
 * row_start <= i < row_end
 * col_start <= j < col_end
 * depth_start <= k < depth_end
 *
 * The right extrema is NOT included
 *
 * @param row_start
 * @param row_end
 * @param col_start
 * @param col_end
 * @param depth_start
 * @param depth_end
 * @return the subset of the original tensor
 */
Tensor Tensor::subset(unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end, unsigned int depth_start, unsigned int depth_end){
    Tensor res(row_end-row_start,col_end-col_start,depth_end-depth_start);
    cout<<res.r<<"dsda";
    cout<<res.c<<"dsd";
    cout<<res.d<<"dsds"<<"\n";
    int i{0};
    Tensor att= (*this);
    for(int z=res.d;z<depth_end;z++){
        for(int y=res.c;y<col_end;y++){
            for(int x=res.r;x<row_end;x++){
               cout<< att.data[(z*r*c)+(y*r+(x))]<<"||"<<"AAAAAAAAA";
               res.data[i++]=att.data[(z*r*c)+(y*r+(x))];
            }
        }
    }
    return res;
};

/** 
     * Concatenate 
     * 
     * The function concatenates two tensors along a give axis
     * 
     * Example: this is of size 10x5x6 and rhs is of 25x5x6
     * 
     * if concat on axis 0 (row) the result will be a new Tensor of size 35x5x6
     * 
     * if concat on axis 1 (columns) the operation will fail because the number 
     * of rows are different (10 and 25).
     * 
     * In order to perform the concatenation is mandatory that all the dimensions 
     * different from the axis should be equal, other wise throw concat_wrong_dimension(). 
     *  
     * @param rhs The tensor to concatenate with
     * @param axis The axis along which perform the concatenation 
     * @return a new Tensor containing the result of the concatenation
     */
Tensor concat(const Tensor &rhs, int axis=0);


    /**
     * Convolution
     *
     * This function performs the convolution of the Tensor with a filter.
     *
     * The filter f must have odd dimensions and same depth.
     *
     * Remeber to apply the padding before running the convolution
     *
     * @param f The filter
     * @return a new Tensor containing the result of the convolution
     */
    Tensor convolve(const Tensor &f);

    /* UTILITY */

    /**
     * Rows
     *
     * @return the number of rows in the tensor
     */
    int rows();

    /**
     * Cols
     *
     * @return the number of columns in the tensor
     */
    int cols();

    /**
     * Depth
     *
     * @return the depth of the tensor
     */
    int depth();
    
    /**
     * Get minimum
     *
     * Compute the minimum value considering a particular index in the third dimension
     *
     * @return the minimum of data( , , k)
     */
    float getMin(int k);

    /**
     * Get maximum
     *
     * Compute the maximum value considering a particular index in the third dimension
     *
     * @return the maximum of data( , , k)
     */
    float getMax(int k);

    /**
     * showSize
     *
     * shows the dimensions of the tensor on the standard output.
     *
     * The format is the following:
     * rows" x "colums" x "depth
     *
     */
    void showSize();
    
    /* IOSTREAM */

    /**
     * Operator overloading <<
     *
     * Use the overaloading of << to show the content of the tensor.
     *
     * You are free to chose the output format, btw we suggest you to show the tensor by layer.
     *
     * [..., ..., 0]
     * [..., ..., 1]
     * ...
     * [..., ..., k]
     */
    friend ostream& operator<< (ostream& stream, const Tensor & obj);

    /**
     * Reading from file
     *
     * Load the content of a tensor from a textual file.
     *
     * The file should have this structure: the first three lines provide the dimensions while
     * the following lines contains the actual data by channel.
     *
     * For example, a tensor of size 4x3x2 will have the following structure:
     * 4
     * 3
     * 2
     * data(0,0,0)
     * data(0,1,0)
     * data(0,2,0)
     * data(1,0,0)
     * data(1,1,0)
     * .
     * .
     * .
     * data(3,1,1)
     * data(3,2,1)
     *
     * if the file is not reachable throw unable_to_read_file()
     *
     * @param filename the filename where the tensor is stored
     */
    void read_file(string filename);

    /**
     * Write the tensor to a file
     *
     * Write the content of a tensor to a textual file.
     *
     * The file should have this structure: the first three lines provide the dimensions while
     * the following lines contains the actual data by channel.
     *
     * For example, a tensor of size 4x3x2 will have the following structure:
     * 4
     * 3
     * 2
     * data(0,0,0)
     * data(0,1,0)
     * data(0,2,0)
     * data(1,0,0)
     * data(1,1,0)
     * .
     * .
     * .
     * data(3,1,1)
     * data(3,2,1)
     *
     * if the file is not reachable throw unable_to_read_file()
     *
     */
    void write_file(string filename);

};
