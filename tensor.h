#ifndef TENSOR_H
#define TENSOR_H

#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;

//IMPLEMENTARE INSERIMENTO ALL INTERO DELL ARRAY


class Tensor{
private:
    float * data;
    int x,y,z;
public:
    //stampa provvisoria
    void stampa(){
        for(int zz=0;zz<z;zz++){
            for(int yy=0;yy<y;yy++){
                for(int xx=0;xx<x;xx++){
                    std::cout<<data[zz*x*y+(yy*x+xx)]<<"||";
                }
                std::cout<<"\n";
            }
            std::cout<<"\n"<<"nuova dim"<<"\n";
        }
    }


    /**
     * Class constructor
     * 
     * Parameter-less class constructor 
     */
    Tensor(){
        data = nullptr;
        this->x = 0;
        this->y = 0;
        this->z = 0;
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
    Tensor(int r, int c, int d, float v = 0.0){
        x=r;
        y=c;
        z=d;
        data=new float[x*y*z];
        for(int zz=0;zz<z;zz++){
            for(int yy=0;yy<y;yy++){
                for(int xx=0;xx<x;xx++){
                    data[zz*x*y+(yy*x+xx)]=0.0;
                }
            }
        }
    };
    /**
     * Class distructor
     * 
     * Cleanup the data when deallocated
     */
    ~Tensor(){
        delete data;
    };

    /**
     * Operator oveloding ()
     * 
     * if indexes are out of bound throw index_out_of_bound() exception
     * 
     * @return the value at location [i][j][k]
     */
    float operator()(int i, int j, int k) const{
        return data[k*x*y+(j*x+i)];
        //TODO EXCEPTION
    };

    /**
     * Operator oveloding ()
     * 
     * Return the pointer to the location [i][j][k] such that the operator (i,j,k) can be used to 
     * modify tensor data.
     * 
     * If indexes are out of bound throw index_out_of_bound() exception
     * 
     * @return the pointer to the location [i][j][k]
     */
    float &operator()(int i, int j, int k){
        float& res=data[k*x*y+(j*x+i)];
        return res;
        //TODO EXCEPTION
    };
    int at(int i,int j,int k){
        return k*x*y+(j*x+i);
    };
    /**
     * Copy constructor
     * 
     * This constructor copy the data from another Tensor
     *      
     * @return the new Tensor
     */
    Tensor(const Tensor& that){
    	x=that.x;
        y=that.y;
        z=that.z;
        data=new float[x*y*z];
        
        for(int zz=0; zz<z; zz++){    
            for(int yy=0; yy<y; yy++){
                for(int i=0; i<x; i++){
                    data(i, yy, zz) = that(i, yy, zz);
                }
            }
        }
    };

    /**
     * Operator oveloding -
     * 
     * Perform the point-wise difference between two Tensors.
     * 
     * lhs(i,j,k)=lhs(i,j,k)-rhs(i,j,k)
     * 
     * 
     * @return lhs with the result of the operation (lhs is passed by copy, so is a new lhs ;) )
     */
    friend Tensor operator-(Tensor lhs, const Tensor &rhs){
        Tensor res(lhs.x,lhs.y,lhs.z,0.0);
        for(int zz=0;zz<lhs.z;zz++){
            for(int yy=0;yy<lhs.y;yy++){
                for(int xx=0;xx<lhs.x;xx++){
                    res(xx,yy,zz)=lhs(xx,yy,zz)-rhs(xx,yy,zz);
                }
            }
        }
        return res;
    };

    /**
     * Operator oveloding +
     * 
     * Perform the point-wise difference between two Tensors.
     * 
     * lhs(i,j,k)=lhs(i,j,k)+rhs(i,j,k)
     * 
     * 
     * @return lhs with the result of the operation (lhs is passed by copy, so is a new lhs ;) )
     */
    friend Tensor operator+(Tensor lhs, const Tensor &rhs){
        Tensor res(lhs.x,lhs.y,lhs.z,0.0);
        for(int zz=0;zz<lhs.z;zz++){
            for(int yy=0;yy<lhs.y;yy++){
                for(int xx=0;xx<lhs.x;xx++){
                    res(xx,yy,zz)=lhs(xx,yy,zz)+rhs(xx,yy,zz);
                }
            }
        }
        return res;
    };

    /**
     * Operator oveloding *
     * 
     * Perform the point-wise difference between two Tensors.
     * 
     * lhs(i,j,k)=lhs(i,j,k)*rhs(i,j,k)
     * 
     * 
     * @return lhs with the result of the operation (lhs is passed by copy, so is a new lhs ;) )
     */
    friend Tensor operator*(Tensor lhs, const Tensor &rhs){
        Tensor res(lhs.x,lhs.y,lhs.z,0.0);
        for(int zz=0;zz<lhs.z;zz++){
            for(int yy=0;yy<lhs.y;yy++){
                for(int xx=0;xx<lhs.x;xx++){
                    res(xx,yy,zz)=lhs(xx,yy,zz)*rhs(xx,yy,zz);
                }
            }
        }
        return res;
    };
    
    /**
     * Operator oveloding /
     * 
     * Perform the point-wise difference between two Tensors.
     * 
     * lhs(i,j,k)=lhs(i,j,k)/rhs(i,j,k)
     * 
     * 
     * @return lhs with the result of the operation (lhs is passed by copy, so is a new lhs ;) )
     */
    friend Tensor operator/(Tensor lhs, const Tensor &rhs){
        Tensor res(lhs.x,lhs.y,lhs.z,0.0);
        for(int zz=0;zz<lhs.z;zz++){
            for(int yy=0;yy<lhs.y;yy++){
                for(int xx=0;xx<lhs.x;xx++){
                    res(xx,yy,zz)=lhs(xx,yy,zz)/rhs(xx,yy,zz);
                }
            }
        }
        return res;
    };

    /**
     * Operator oveloding - between a Tensor and a constant
     * 
     * Perform the point-wise difference between two Tensors.
     * 
     * lhs(i,j,k)=lhs(i,j,k)-rhs
     * 
     * 
     * @return lhs with the result of the operation (lhs is passed by copy, so is a new lhs ;) )
     */
    friend Tensor operator-(Tensor lhs, const float &rhs){
        Tensor res(lhs.x,lhs.y,lhs.z,0.0);
        for(int zz=0;zz<lhs.z;zz++){
            for(int yy=0;yy<lhs.y;yy++){
                for(int xx=0;xx<lhs.x;xx++){
                    res(xx,yy,zz)=lhs(xx,yy,zz)/rhs;
                }
            }
        }
        return res;
    };

    /**
     * Operator oveloding + between a Tensor and a constant
     * 
     * Perform the point-wise difference between two Tensors.
     * 
     * lhs(i,j,k)=lhs(i,j,k)+rhs
     * 
     * 
     * @return lhs with the result of the operation (lhs is passed by copy, so is a new lhs ;) )
     */
    friend Tensor operator+(Tensor lhs, const float &rhs){
        Tensor res(lhs.x,lhs.y,lhs.z,0.0);
        for(int zz=0;zz<lhs.z;zz++){
            for(int yy=0;yy<lhs.y;yy++){
                for(int xx=0;xx<lhs.x;xx++){
                    res(xx,yy,zz)=lhs(xx,yy,zz)+rhs;
                }
            }
        }
        return res;
    };

    /**
     * Operator oveloding * between a Tensor and a constant
     * 
     * Perform the point-wise difference between two Tensors.
     * 
     * lhs(i,j,k)=lhs(i,j,k)*rhs
     * 
     * 
     * @return lhs with the result of the operation (lhs is passed by copy, so is a new lhs ;) )
     */
    friend Tensor operator*(Tensor lhs, const float &rhs){
        Tensor res(lhs.x,lhs.y,lhs.z,0.0);
        for(int zz=0;zz<lhs.z;zz++){
            for(int yy=0;yy<lhs.y;yy++){
                for(int xx=0;xx<lhs.x;xx++){
                    res(xx,yy,zz)=lhs(xx,yy,zz)*rhs;
                }
            }
        }
        return res;
    };
    
    /**
     * Operator oveloding / between a Tensor and a constant
     * 
     * Perform the point-wise difference between two Tensors.
     * 
     * lhs(i,j,k)=lhs(i,j,k)/rhs
     * 
     * 
     * @return lhs with the result of the operation (lhs is passed by copy, so is a new lhs ;) )
     */
    friend Tensor operator/(Tensor lhs, const float &rhs){
        Tensor res(lhs.x,lhs.y,lhs.z,0.0);
        for(int zz=0;zz<lhs.z;zz++){
            for(int yy=0;yy<lhs.y;yy++){
                for(int xx=0;xx<lhs.x;xx++){
                    res(xx,yy,zz)=lhs(xx,yy,zz)/rhs;
                }
            }
        }
        return res;
    };

    /**
     * Operator oveloding = (assignment) 
     * 
     * Perform the assignment between this object and another
     * 
     * @return a pointer to the receiver object
     */
    Tensor & operator=(const Tensor &other);

    /**
     * Random Initialization
     * 
     * Perform a random initialization of the tensor
     * 
     * @param mean The mean
     * @param std  Standard deviation
     */
    void init_random(float mean=1.0, float std=0.0){
        /*
        if(data){
            float y1;
            float y2;
            float num;
            for(int i=0;i<r;i++){
                for(int j=0;j<c;j++){
                    for(int k=0;k<d;k++){
                        y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
                        y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
                        num = cos(2*PI*y2)*sqrt(-2.*log(y1));
                        this->operator()(i,j,k)= mean + num*std;
                    }
                }
            }    

        }else{
            throw(tensor_not_initialized());
        }*/
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
    void init(int r, int c, int d, float v=0.0);

    /**
     * Tensor Clamp
     * 
     * Clamp the tensor such that the lower value becomes low and the higher one become high.
     * 
     * @param low Lower value
     * @param high Higher value 
     */
    void clamp(float low, float high){
    	for(int zz=0; zz<z; zz++)
            for(int yy=0; yy<y; yy++)
                for(int xx=0; xx<x; xx++){
                    float elem = data(xx, yy, zz);
                    if(elem < low)
                        data(xx, yy, zz) = low;
                    else if(elem > high)
                        data(xx, yy, zz) = high;
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
     * MAX(K) e MIN(K) corrispondono al massimo valore e al minimo valore presente in quella specifica dimensione
     * 
     * new_max is the new value for the maximum
     * 
     * @param new_max New maximum vale
     */
    void rescale(float new_max=1.0){
        //trovo max value dim 1
        float max{};
        float min{};
        for(int zz=0;zz<z;zz++){
            for(int yy=0;yy<y;yy++){
                int indice=this->at(0,yy,zz);
                for(int xx=0;xx<x;xx++){
                    data[indice+x]
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
    Tensor padding(int pad_h, int pad_w){
        Tensor res(x+2*pad_h, y+2*pad_w, z);
        int i=0;
        float new_x = x+2*pad_h;
        float new_y = x+2*pad_w;
        for(int zz=0; zz<z; zz++){
            for(int yy=0; yy<new_y; yy++){
                for(int xx=0; xx<new_x; xx++){
                    if(xx < pad_h || xx > x+pad_h || yy < pad_w || yy > y+pad_w){   //da controllare !!
                        res.data(xx, yy, zz) = 0.0;
                    }else{
                        res.data(xx, yy, zz) = data[i++];
                    }
                }
            }
        }
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
    Tensor subset(unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end, unsigned int depth_start, unsigned int depth_end);

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
     * shows the dimensions of the tensor
     */
    void showSize(){
        cout<<this->rows()<<" "<<this->cols()<<" "<<this->depth()<<endl;
    }
    
    /* IOSTREAM */

    /**
     * Operator overloading <<
     * 
     * Use the overaloading of << to show the content of the tensor.
     * 
     * The content is shown considering each data( , , k).
     * 
     */
    friend ostream& operator<< (ostream& stream, const Tensor & obj);

    /**
     * Reading from file
     * 
     * Load the content of a tensor from file.
     * 
     * The file should have this structure: the first three lines provide the dimensions while 
     * the following lines contains the actual data by channel.
     * 
     * Example:
     * number of rows
     * number of columns
     * depth
     * value(0,0,0)
     * value(0,0,1)
     * value(0,0,2)
     * .
     * .
     * .
     * value(0,0,k)
     * value(0,1,0)
     * 
     * if the file is not reachable throw unable_to_read_file()
     * 
     * @param filename the filename where the tensor is stored
     */
    void read_file(string filename);

    /**
     * Write the tensor to a file
     * 
     * Write the content of a tensor to a file.
     * 
     * The file should have this structure: the first three lines provide the dimensions while 
     * the following lines contains the actual data by channel.
     * 
     * Example:
     * number of rows
     * number of columns
     * depth
     * value(0,0,0)
     * value(0,0,1)
     * value(0,0,2)
     * .
     * .
     * .
     * value(0,0,k)
     * value(0,1,0)
     * 
     * if the file is not reachable throw unable_to_read_file()
     * 
     */
    void write_file(string filename);
};

#endif
