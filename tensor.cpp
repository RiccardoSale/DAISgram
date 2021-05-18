#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>
#include "dais_exc.h"
#include "tensor.h"
#include <algorithm>

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;

/**
* Class constructor
*
* Parameter-less class constructor
*/
Tensor::Tensor() {
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
Tensor::Tensor(int r, int c, int d, float v) {
    this->r = r;
    this->c = c;
    this->d = d;
    int i_max{r * c * d};
    data = new float[i_max];
    for (int i = 0; i < i_max; i++) {
        data[i] = v;
    }
};

/**
 * Class distructor
 *
 * Cleanup the data when deallocated
 */
Tensor::~Tensor() {
    delete[] data;
};

/**
 * Operator oveloding ()
 *
 * if indexes are out of bound throw index_out_of_bound() exception
 *
 * @return the value at location [i][j][k]
 */
float Tensor::operator()(int i, int j, int k) const {
    if (i > r || j > c || k > d) {
        throw (index_out_of_bound());
    } else {
        return data[k * r * c + (i * c + j)];
    }
}

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
float &Tensor::operator()(int i, int j, int k) {
    if (i > r || j > c || k > d) {
        throw (index_out_of_bound());
    } else {
        float &res = data[k * r * c + (i * c + j)];//dim *riga_base *col_base + righe *col_base +col
        return res;
    }
}

float &Tensor::at(int i) const {
    float &res = data[i];
    return res;
}

/**
 * Copy constructor
 *
 * This constructor copy the data from another Tensor
 *
 * @return the new Tensor
 */
Tensor::Tensor(const Tensor &that) {
    if (that.data){
        r = that.r;
        c = that.c;
        d = that.d;
        int i_max{r * c * d};
        data = new float[i_max];
        for (int i = 0; i < i_max; i++) {
            data[i] = that.data[i];
        }
    }else{
        throw (tensor_not_initialized());
    }
}

bool Tensor::operator==(const Tensor& rhs) const {
    if(rhs.data && data){
        if(r == rhs.r && c == rhs.c && d == rhs.d){
            int max = this->rows() * this->cols() * this->depth();
            for (int i = 0; i < max; i++)
                if (abs(this->data[i] - rhs.data[i]) >= EPSILON)
                    return false;
            return true;
        }else{
            throw (dimension_mismatch());
        }
    }else{
        throw (tensor_not_initialized());
    }
}
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
Tensor Tensor::operator-(const Tensor &rhs) const {
    if(rhs.data && data){
        if (r == rhs.r && c == rhs.c && d == rhs.d) {
            Tensor result(r, c, d);
            int i_max{r * c * d};
            for (int i = 0; i < i_max; i++) {
                result.data[i] = data[i] - rhs.data[i];
            }
            return result;
        } else {
            throw (dimension_mismatch());
        }
    }else{
        throw (tensor_not_initialized());
    }
}

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
Tensor Tensor::operator+(const Tensor &rhs) const {
    if(rhs.data && data){
        if (r == rhs.r && c == rhs.c && d == rhs.d) {
            Tensor result(r, c, d);
            int i_max{r * c * d};
            for (int i = 0; i < i_max; i++) {
                result.data[i] = data[i] + rhs.data[i];
            }
            return result;
        }else{
            throw (dimension_mismatch());
        }
    }else{
        throw (tensor_not_initialized());
    }
}

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
Tensor Tensor::operator*(const Tensor &rhs) const {
    if(rhs.data && data){
        if (r == rhs.r && c == rhs.c && d == rhs.d) {
            Tensor result(r, c, d);
            int i_max{r * c * d};
            for (int i = 0; i < i_max; i++) {
                result.data[i] = data[i] * rhs.data[i];
            }
            return result;
        } else {
            throw (dimension_mismatch());
        }
    }else{
        throw (tensor_not_initialized());   
    }
}

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
Tensor Tensor::operator/(const Tensor &rhs) const {
    if(rhs.data && data){
        if (r == rhs.r && c == rhs.c && d == rhs.d) {
            Tensor result(r, c, d);
            int i_max{r * c * d};
            for (int i = 0; i < i_max; i++) {
                result.data[i] = data[i] / rhs.data[i];
            }
            return result;
        } else {
            throw (dimension_mismatch());
        }
    }else{
        throw (tensor_not_initialized());
    }
}

/**
 * Operator overloading -
 *
 * It performs the point-wise difference between a Tensor and a constant
 *
 * result(i,j,k)=this(i,j,k)-rhs
 *
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator-(const float &rhs) const {
    Tensor result(r, c, d);
    int i_max{r * c * d};
    for (int i = 0; i < i_max; i++) {
        result.data[i] = data[i] - rhs;
    }
    return result;
}

/**
 * Operator overloading +
 *
 * It performs the point-wise sum between a Tensor and a constant
 *
 * result(i,j,k)=this(i,j,k)+rhs
 *
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator+(const float &rhs) const {
    Tensor result(r, c, d);
    int i_max{r * c * d};
    for (int i = 0; i < i_max; i++) {
        result.data[i] = data[i] + rhs;
    }
    return result;
}

/**
 * Operator overloading *
 *
 * It performs the point-wise product between a Tensor and a constant
 *
 * result(i,j,k)=this(i,j,k)*rhs
 *
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator*(const float &rhs) const {
    Tensor result(r, c, d);
    int i_max{r * c * d};
    for (int i = 0; i < i_max; i++) {
        result.data[i] = data[i] * rhs;
    }
    return result;
}

/**
 * Operator overloading / between a Tensor and a constant
 *
 * It performs the point-wise division between a Tensor and a constant
 *
 * result(i,j,k)=this(i,j,k)/rhs
 *
 * @return returns a new Tensor containing the result of the operation
 */
Tensor Tensor::operator/(const float &rhs) const {
    Tensor result(r, c, d);
    int i_max{r * c * d};
    for (int i = 0; i < i_max; i++) {
        result.data[i] = data[i] / rhs;
    }
    return result;
}

/**
 * Operator overloading = (assignment)
 *
 * Perform the assignment between this object and another
 *
 * @return a reference to the receiver object
 */

Tensor &Tensor::operator=(const Tensor &other) {    //controllare cambiare
    if(other.data){
        d = other.d;
        r = other.r;
        c = other.c;
        int i_max = r * c * d;
        //controllare con valgrind se serve delete
        if (data)
            delete[] data;
        data = new float[i_max];
        for (int i = 0; i < i_max; i++) {
            data[i] = other.data[i];
        }
        return *this;
    }else{
        throw (tensor_not_initialized());
    }
}


/**
 * Random Initialization
 * 
 * Perform a random initialization of the tensor
 * 
 * @param mean The mean
 * @param std  Standard deviation
 */
void Tensor::init_random(float mean, float std) {
    if (data) {

        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean, std);

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                for (int k = 0; k < d; k++) {
                    this->operator()(i, j, k) = distribution(generator);
                }
            }
        }

    } else {
        throw (tensor_not_initialized());
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
void Tensor::init(int r, int c, int d, float v) {
    this->r = r;
    this->c = c;
    this->d = d;
    int i_max{r * c * d};
    //controllare con valgrind
    if (data)
        delete[] data;
    data = new float[i_max];
    for (int i = 0; i < i_max; i++) {
        data[i] = v;
    }
}//aspettare documentazione maggiore

/**
 * Tensor Clamp
 *
 * Clamp the tensor such that the lower value becomes low and the higher one become high.
 *
 * @param low Lower value
 * @param high Higher value
 */
void Tensor::clamp(float low, float high) {
    int i_max = r * c * d;
    for (int i = 0; i < i_max; i++) {
        float elem = data[i];
        if (elem < low)
            data[i] = low;
        else if (elem > high)
            data[i] = high;
    }
}

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
void Tensor::rescale(float new_max) {
    //trovo max value dim 1
    float max{data[0]};
    float min{data[0]};
    for (int z = 0; z < d; z++) {
        for (int y = 0; y < c; y++) { // Da testare probabilmente scambiare c e r
            int i = z * r * c + (y * c);
            for (int x = 0; x < r; x++) {
                if (data[i + r] < min)
                    min = data[i + r];
                else if (data[i + r] > max)
                    max = data[i + r];
            }
        }//a questo punto dovrei aver trovato maggiore e minore per quella dimensione
        float diff = max - min;
        for (int y = 0; y < c; y++) {
            int i = z * r * c + (y * c);
            for (int x = 0; x < r; x++) {
                data[i + r] = ((data[i + r] - min) / (diff)) * new_max; //vedere caso 0/0
            }
        }
    }
}

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
Tensor Tensor::padding(int pad_h, int pad_w) const {
    int new_x = c + 2 * pad_w; //nuovo numero colonne
    int new_y = r + 2 * pad_h; //nuovo numero righe
    Tensor res(new_y, new_x, d);
    int i = 0;
    int ii = 0;
    for (int z = 0; z < d; z++) { //scorro dimensioni
        for (int y = 0; y < new_y; y++) {
            for (int x = 0; x < new_x; x++) {
                if (x < pad_w || x >= c + pad_w || y < pad_h || y >= r + pad_h) {
                    res.data[ii++] = 0.0;
                } else {
                    res.data[ii++] = data[i++];
                }
            }
        }
    }
    return res;
}

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
Tensor Tensor::subset(unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end,
                      unsigned int depth_start, unsigned int depth_end) const {
    Tensor res(row_end - row_start, col_end - col_start, depth_end - depth_start);
    int i = 0;
    for (int z = depth_start; z < depth_end; z++) {
        for (int y = col_start; y < col_end; y++) {
            for (int x = row_start; x < row_end; x++) {
                res.data[i++] = data[(z * r * c) + y * c + x];
            }
        }
    }
    return res;
}

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
Tensor Tensor::concat(const Tensor &rhs, int axis) const {
    if(rhs.data && data){
        if (axis == 0) { //l'asse è sulle righe
            if (c == rhs.c && d == rhs.d) {
                Tensor result(r + rhs.r, c, d);
                //result.r = r + rhs.r;

                int my_pos = 0;
                int rhs_pos = 0;
                int res_pos = 0;

                int my_dimension = r * c;
                int rhs_dimension = rhs.r * rhs.c;

                for (int i = 0; i < d; i++) {
                    for (int j = 0; j < my_dimension; j++) {
                        result.data[res_pos++] = data[my_pos++];
                    }
                    for (int k = 0; k < rhs_dimension; k++) {
                        result.data[res_pos++] = rhs.data[rhs_pos++];
                    }
                }
                return result;
            } else {
                throw (concat_wrong_dimension());
            }
        } else if (axis == 1) {
            if (r == rhs.r && d == rhs.d) {
                Tensor result(r, c + rhs.c, d);

                int my_pos = 0;
                int rhs_pos = 0;
                int res_pos = 0;

                //j e k arrivano fino a c perchè ogni c colonne sono una riga
                //i arriva fino a r * d perchè lo devo fare per tutte r le righe e per d dimensioni

                for (int i = 0; i < d; i++) {

                    for (int s = 0; s < r; s++) {
                        for (int j = 0; j < c; j++) {
                            result.data[res_pos++] = data[my_pos++];
                        }

                        for (int k = 0; k < rhs.c; k++) {
                            result.data[res_pos++] = rhs.data[rhs_pos++];
                        }
                    }

                }
                return result;
            } else {
                throw (concat_wrong_dimension());
            }
        }
    }else{
        throw (tensor_not_initialized());
    }
}

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

Tensor Tensor::convolve(const Tensor &f) const {
    Tensor res;
    res.init(rows(), cols(), depth()); //risultato uguale di dim alla "matrice" originale
    int p = (f.rows() - 1) / 2; //formuletta del padding
    Tensor pad = padding(p, p); //tensore che ho con il padding fatto

    int counter;
    int cf = 0;
    int cd = 0;
    float sum = 0;
    int quadrati_riga = pad.cols() - f.cols();
    int quadrati_altezza = pad.rows() - f.rows(); //quanti quadrati devo fare in altezza

    for (int depth = 0; depth < res.depth(); depth++) {
        for (int qa = 0; qa <= quadrati_altezza; qa++) {
            for (int a = 0; a <= quadrati_riga; a++) { //quanti quadrati faccio per riga
                //questi due for fanno un quadrato  //ogni volta che faccio un quadrato sulla riga devo aumentare counter di partenza di 1
                for (int i = 0; i < f.rows(); i++) { //ciclo per fare "3 righe" //doppio for per fare dimensione filtro
                    counter = (pad.cols() * (i + qa)) + a + (depth * pad.rows() *pad.cols()); // vado all inzio della prossima riga e aggiungo a per spostarmi di uno a destra in base al quadrato che sto facendo
                    for (int j = 0; j < f.cols(); j++) {//ciclo per fare le "3 colonne" ogni riga  //
                        sum += pad.data[counter++] * f.data[cf++];
                    }
                    //qa =profondita del cubo quindi in quel caso devo saltare più righe e dopodichè aggiungo a per spostartmi a destra
                }
                //dopo i due for ho un risultato da caricare
                res.data[cd++] = sum;
                sum = 0;
                cf = 0;//azzero contatore filtro
            }
        }
    }
    return res;
}

/* UTILITY */

/**
 * Rows
 *
 * @return the number of rows in the tensor
 */
int Tensor::rows() const {
    return r;
}

/**
 * Cols
 *
 * @return the number of columns in the tensor
 */
int Tensor::cols() const {
    return c;
}

/**
 * Depth
 *
 * @return the depth of the tensor
 */
int Tensor::depth() const {
    return d;
}

/**
 * Get minimum
 *
 * Compute the minimum value considering a particular index in the third dimension
 *
 * @return the minimum of data( , , k)
 */
float Tensor::getMin(int k) const {
    float min = data[0];
    int stop = k * r * c + r * c;
    for (int i = k * r * c; i < stop; i++) {
        if (data[i] < min)
            min = data[i];
    }
    return min;
}

/**
 * Get maximum
 *
 * Compute the maximum value considering a particular index in the third dimension
 *
 * @return the maximum of data( , , k)
 */
float Tensor::getMax(int k) const {
    float max = data[0];
    int stop = k * r * c + r * c;
    for (int i = k * r * c; i < stop; i++) {
        if (data[i] > max)
            max = data[i];
    }
    return max;
}

/**
 * showSize
 *
 * shows the dimensions of the tensor on the standard output.
 *
 * The format is the following:
 * rows" x "colums" x "depth
 *
 */
void Tensor::showSize() const {
    cout << r << " x " << c << " x " << d;
}

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
//friend ostream& operator<< (ostream& stream, const Tensor & obj);

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
void Tensor::read_file(string filename) {
    ifstream f(filename, ifstream::in);
    if (f.is_open()) {
        f >> r;
        f >> c;
        f >> d;
        int max = r * c * d;
        data = new float[max];
        for (int i = 0; i < max; i++)
            f >> data[i];
    } else
        throw (unable_to_read_file());

}

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
void Tensor::write_file(string filename) {
    ofstream f(filename, ofstream::out);
    if (f.is_open()) {
        f << r << endl;
        f << c << endl;
        f << d << endl;
        for (int i = 0; i < r * c * d; i++)
            f << data[i] << endl;
    } else
        throw (unable_to_read_file());
}

