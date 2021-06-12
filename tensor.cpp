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

Tensor::Tensor() {
    data = nullptr;
    this->r = this->c = this->d = 0;
}

Tensor::Tensor(int r, int c, int d, float v) {
    this->r = r;
    this->c = c;
    this->d = d;
    int i_max{r * c * d};
    data = new float[i_max];
    for (int i = 0; i < i_max; i++)
        data[i] = v;
}

Tensor::Tensor(int r,int c,int d,int a[]) {
    this->r = r;
    this->c = c;
    this->d = d;
    int i_max{r * c * d};
    data = new float[i_max];
    int ii{};
    for (int i = 0; i < i_max; i++) {
        if(c==9) c=0;
        data[i] = a[ii++];
        c++;
    }
}

Tensor::~Tensor() {
    delete[] data;
}

float Tensor::operator()(int i, int j, int k) const {
    if(data) {
        if (i >= r || j >= c || k >= d || i<0 || j<0 || k<0) throw (index_out_of_bound());
        return data[k * r * c + (i * c + j)]; //k*r*c = salto tot dimensioni / i*c =salto tot "righe" +j= aggiungo colonne
    }else throw (tensor_not_initialized());
}

float &Tensor::operator()(int i, int j, int k) {
    if(data) {
        if (i >= r || j >= c || k >= d || i<0 || j<0 || k<0) throw (index_out_of_bound());
        float &res = data[k * r * c + (i * c + j)];
        return res;
    }else throw (tensor_not_initialized());
}

float &Tensor::at(int i) {
    if(data) {
        float &res = data[i];
        return res;
    }else throw (tensor_not_initialized());
}

float Tensor::at(int i)const {
    if(data) {
        return data[i];
    }else throw (tensor_not_initialized());
}

Tensor::Tensor(const Tensor &that) {
    if (that.data){
        r = that.r;
        c = that.c;
        d = that.d;
        int i_max{r * c * d};
        data = new float[i_max];
        for (int i = 0; i < i_max; i++)
            data[i] = that.data[i];
    }else throw (tensor_not_initialized());
}

bool Tensor::operator==(const Tensor& rhs) const {
    if(rhs.data && data){
        if(r == rhs.r && c == rhs.c && d == rhs.d){
            int max = r * c * d;
            for (int i = 0; i < max; i++)
                if (abs(data[i] - rhs.data[i]) >= EPSILON)
                    return false;
            return true;
        }else{
            throw (dimension_mismatch());
        }
    }else throw (tensor_not_initialized());
}

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
    }else throw (tensor_not_initialized());
}

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
    }else throw (tensor_not_initialized());
}

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
    }else throw (tensor_not_initialized());
}

Tensor Tensor::operator/(const Tensor &rhs) const {
    if(rhs.data && data){
        if (r == rhs.r && c == rhs.c && d == rhs.d) {
            Tensor result(r, c, d);
            int i_max{r * c * d};
            for (int i = 0; i < i_max; i++) {
                if(rhs.data[i]==0) throw (unknown_exception());
                result.data[i] = data[i] / rhs.data[i];
            }
            return result;
        } else {
            throw (dimension_mismatch());
        }
    }else throw (tensor_not_initialized());
}

Tensor Tensor::operator-(const float &rhs) const {
    if(data){
        Tensor result(r, c, d);
        int i_max{r * c * d};
        for (int i = 0; i < i_max; i++) {
            result.data[i] = data[i] - rhs;
        }
        return result;
    }else throw (tensor_not_initialized());
}

Tensor Tensor::operator+(const float &rhs) const {
    if(data){
        Tensor result(r, c, d);
        int i_max{r * c * d};
        for (int i = 0; i < i_max; i++) {
            result.data[i] = data[i] + rhs;
        }
        return result;
    }else throw (tensor_not_initialized());
}

Tensor Tensor::operator*(const float &rhs) const {
    if(data){
        Tensor result(r, c, d);
        int i_max{r * c * d};
        for (int i = 0; i < i_max; i++) {
            result.data[i] = data[i] * rhs;
        }
        return result;
    }else throw (tensor_not_initialized());
}

Tensor Tensor::operator/(const float &rhs) const {
    if(data){
        if(rhs==0) throw (unknown_exception());
        Tensor result(r, c, d);
        int i_max{r * c * d};
        for (int i = 0; i < i_max; i++) {
            result.data[i] = data[i] / rhs;
        }
        return result;
    }else throw (tensor_not_initialized());
}

Tensor &Tensor::operator=(const Tensor &other) {
        int i_max = other.r * other.c * other.d;
        if(other.d!=d || other.r!=r || other.c!=c) {
            d = other.d;
            r = other.r;
            c = other.c;
            if (data)
                delete[] data;
            data = new float[i_max];
        }
        for (int i = 0; i < i_max; i++)
            data[i] = other.data[i];
        return *this;
}

void Tensor::init_random(float mean, float std) {
    if(data) {
        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean, std);

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                for (int k = 0; k < d; k++) {
                    this->operator()(i, j, k) = distribution(generator);
                }
            }
        }
    }else throw (tensor_not_initialized());
}

void Tensor::init(int r, int c, int d, float v) {
    if(r<0 || c<0 || d<0) throw (index_out_of_bound());
    this->r = r;
    this->c = c;
    this->d = d;
    int i_max{r * c * d};
    if (data)
        delete[] data;
    data = new float[i_max];
    for (int i = 0; i < i_max; i++)
        data[i] = v;
}

void Tensor::clamp(float low, float high) {
    if(data) {
        int i_max = r * c * d;
        for (int i = 0; i < i_max; i++) {
            float elem = data[i];
            if (elem < low)
                data[i] = low;
            else if (elem > high)
                data[i] = high;
        }
    }else throw (tensor_not_initialized());
}

void Tensor::rescale(float new_max) {
    //trovo max value dim 1
    if(data) {
        if(new_max>255 || new_max<0)throw unknown_exception();
        float max{data[0]};
        float min{data[0]};
        float val;
        for (int z = 0; z < d; z++) { //ciclo dimensioni
            for (int y = 0; y < r; y++) {
                int i = z * r * c + (y * c);
                for (int x = 0; x < r; x++) {
                    val = data[i + x];
                    if (val < min)
                        min = val;
                    else if (val > max)
                        max = val;
                }
            }//a questo punto dovrei aver trovato maggiore e minore per quella dimensione
            float diff = max - min;
            if(diff==0) throw (unknown_exception());
            for (int y = 0; y < r; y++) {
                int i = z * r * c + (y * c);
                for (int x = 0; x < c; x++) {
                    if ((data[i + x] - min) == 0)
                        data[i + x] = 0;
                    else {
                        data[i + x] = ((data[i + x] - min) / (diff)) * new_max;
                    }
                }
            }
        }
    }else throw (tensor_not_initialized());
}

Tensor Tensor::padding(int pad_h, int pad_w) const {
    if(data) {
        if(pad_h <0 || pad_w<0) throw (unknown_exception());
        int new_x = c + 2 * pad_w; //nuovo numero colonne
        int new_y = r + 2 * pad_h; //nuovo numero righe
        Tensor res(new_y, new_x, d);
        int i = 0;
        int ii = 0;
        for (int z = 0; z < d; z++) { //scorro dimensioni
            for (int y = 0; y < new_y; y++) {
                for (int x = 0; x < new_x; x++) {
                    if (x < pad_w || x >= c + pad_w || y < pad_h ||
                        y >= r + pad_h) { //aggiungiamo 0 o valori precedentemente esistenti in base alla posizione
                        res.data[ii++] = 0.0;
                    } else {
                        res.data[ii++] = data[i++];
                    }
                }
            }
        }
        return res;
    }else throw (tensor_not_initialized());
}

Tensor Tensor::subset(unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end,unsigned int depth_start, unsigned int depth_end) const {
    if(data) {
        if(row_start <0 || row_end <0 || col_start<0 || col_end<0 || depth_start <0 || depth_end<0 || row_start>row_end || col_start>col_end || depth_start>depth_end) throw index_out_of_bound();
        Tensor res(row_end - row_start, col_end - col_start, depth_end - depth_start);
        int i = 0;
        for (int z = depth_start; z < depth_end; z++) {
            for (int y = col_start; y < col_end; y++) {
                for (int x = row_start; x < row_end; x++) { //partendo dagli indici corretti ci ricaviamo la subset
                    res.data[i++] = data[(z * r * c) + y * c + x]; //salvo nel nuovo tensore il valore trovato
                }
            }
        }
        return res;
    }else throw (tensor_not_initialized());
}

Tensor Tensor::concat(const Tensor &rhs, int axis) const {
    if(rhs.data && data){
        Tensor res;
        if(axis <0 || axis >1) throw unknown_exception();
        int my_pos = 0;
        int rhs_pos = 0;
        int res_pos = 0;
        if (axis == 0) { //l'asse è sulle righe
            if (c == rhs.c && d == rhs.d) {
                res.init(r+rhs.r,c,d);

                int my_dimension = r * c;
                int rhs_dimension = rhs.r * rhs.c;

                for (int i = 0; i < d; i++) {
                    for (int m = 0; m < my_dimension; m++) { //copio prima tutto da uno e poi tutto dall'altro dato che ho "aumentato" num righe
                        res.data[res_pos++] = data[my_pos++];
                    }
                    for (int rh = 0; rh < rhs_dimension; rh++) {
                        res.data[res_pos++] = rhs.data[rhs_pos++];
                    }
                }
            } else {
                throw (concat_wrong_dimension());
            }
        } else if (axis == 1) { //asse colonne
            if (r == rhs.r && d == rhs.d) {
                res.init(r,c+rhs.c,d);

                for (int i = 0; i < d; i++) { //lo faccio per tot dim
                    for (int nr = 0; nr < r; nr++) {
                        for (int my = 0; my < c; my++) { //inserisco la riga di this
                            res.data[res_pos++] = data[my_pos++];
                        }
                        for (int rh = 0; rh < rhs.c; rh++) { //inserisco la riga di rhs andando a formare cosi la nuova riga di res e ripeto tot volte
                            res.data[res_pos++] = rhs.data[rhs_pos++];
                        }
                    }
                }
            } else {
                throw (concat_wrong_dimension());
            }
        }
        return res;
    }else{
        throw (tensor_not_initialized());
    }
}

Tensor Tensor::convolve(const Tensor &f) const {
    if(data && f.data) {
        Tensor res;
        res.init(rows(), cols(), depth()); //risultato uguale di dim alla "matrice" originale
        int p = (f.rows() - 1) / 2; //formuletta del padding
        Tensor pad = padding(p, p); //tensore che ho con il padding fatto
        int counter;
        int cf = 0;
        int cd = 0;
        float sum = 0;
        int ncube_rows = pad.cols() - f.cols();//quanti quadrati per riga
        int ncube_high = pad.rows() - f.rows(); //quanti quadrati devo fare in altezza

        for (int depth = 0; depth < res.depth(); depth++) {
            for (int qa = 0; qa <= ncube_high; qa++) {
                for (int a = 0; a <= ncube_rows; a++) { //quanti quadrati faccio per riga
                    //questi due for fanno un quadrato  //ogni volta che faccio un quadrato sulla riga devo aumentare counter di partenza di 1
                    for (int i = 0;
                         i < f.rows(); i++) { //ciclo per fare "3 righe" //doppio for per fare dimensione filtro
                        counter = (pad.cols() * (i + qa)) + a + (depth * pad.rows() *
                                                                 pad.cols()); // vado all inzio della prossima riga e aggiungo a per spostarmi di uno a destra in base al quadrato che sto facendo
                        for (int j = 0; j < f.cols(); j++) {//ciclo per fare le "3 colonne" ogni riga  //
                            sum += pad.data[counter++] * f.data[cf++];
                        }//qa =profondita del cubo quindi in quel caso devo saltare più righe e dopodichè aggiungo a per spostartmi a destra
                    }
                    res.data[cd++] = sum;//dopo i due for ho un risultato da caricare
                    sum = 0;
                    cf = 0;//azzero contatore filtro
                }
            }
        }
        return res;
    } else throw tensor_not_initialized();
}

/* UTILITY */

int Tensor::rows() const {
    return r;
}

int Tensor::cols() const {
    return c;
}

int Tensor::depth() const {
    return d;
}

float Tensor::getMin(int k) const {
    if(data) {
        if(k<0 || k>=d) throw index_out_of_bound();
        float min = data[0];
        int stop = k * r * c + r * c;
        for (int i = k * r * c; i < stop; i++)
            if (data[i] < min) min = data[i];
        return min;
    }else throw tensor_not_initialized();
}

float Tensor::getMax(int k) const {
    if(data) {
        if(k<0 || k>=d) throw index_out_of_bound();
        float max = data[0];
        int stop = k * r * c + r * c;
        for (int i = k * r * c; i < stop; i++)
            if (data[i] > max) max = data[i];
        return max;
    }else throw tensor_not_initialized();
}

void Tensor::showSize() const {
    cout << r << " x " << c << " x " << d;
}

ostream& operator<<(ostream& stream, const Tensor & obj){
    if(obj.data) {
        int i_max = obj.r * obj.c * obj.d;
        for (int i = 0; i < i_max; i++) {
            if (i % (obj.c * obj.r) == 0 && i != 0)
                stream << "\n" << "dim" << "\n";
            else if (i % obj.c == 0 && i != 0)
                stream << "\n";
            stream << obj.data[i] << "||";
        }
        return stream;
    }else throw tensor_not_initialized();
}

void Tensor::read_file(string filename) {
    ifstream f(filename);
    if (f.is_open()) {
        f >> r;
        f >> c;
        f >> d;
        int max = r * c * d;
        if(data)
            delete []data ;
        data = new float[max];
        for (int i = 0; i < max; i++)
            f >> data[i];
        f.close();
    } else
        throw (unable_to_read_file());
}

void Tensor::write_file(string filename) {
    ofstream f(filename);
    if (f.is_open()) {
        f << r << endl;
        f << c << endl;
        f << d << endl;
        for (int i = 0; i < r * c * d; i++)
            f << data[i] << endl;
        f.close();
    } else
        throw (unable_to_read_file());
}