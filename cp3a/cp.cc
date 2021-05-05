#include "math.h"
#include <vector>
#include <stdlib.h>
#include <malloc.h>
#include <iostream>
/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
typedef double double4_t __attribute__ ((vector_size (4 * sizeof(double))));

static double4_t* double4_t_alloc(std::size_t n) {
    void* tmp = 0;
    if (posix_memalign(&tmp, sizeof(double4_t), sizeof(double4_t) * n)) {
        throw std::bad_alloc();
    }
    return (double4_t*)tmp;
}



double sum(double4_t a)
{
    double sum = 0;
    int len = sizeof(a)/sizeof(double4_t);
    for(int i = 0; i < len; i++)
    {
        sum += a[i];
    }
    return sum;
}

double4_t corr(double4_t a)
{
    double4_t res;
    double sum = 0;

    int len = sizeof(a)/sizeof(double4_t);
    for(int i = 0; i < len ; i++)
    {
        sum += a[i];
    }
    double mean = sum/len;
    std::cout << mean << std::endl;


    double stde = 0;
    for(int x = 0; x < len; x++)
    {
        stde+=pow(a[x]-mean,2);
    }

    for(int x = 0; x < len; x++)
    {
        res[x] = (a[x]-mean)/stde;
    }

    return res;
}


void correlate(int ny, int nx, const float *data, float *result) 
{
    //elements per vector, using doubles here which take 64 bits each (4x in total 256bit == vector registry size)
    constexpr int nb = 4;
    // vectors per input row
    int nvrow = nx/nb;

    double4_t* vd = double4_t_alloc(ny*nvrow);
    
    for(int y = 0 ; y < ny ; y++)
    {
        for(int k = 0 ; k < nvrow ; k++)
        {
            for(int x = 0 ; x < nb ; x++)
            {
                vd[nvrow*y+k][x] = static_cast<double>(data[y*nx + x]);
            }
        }    
    }

    std::cout << "Taalllaaaa!!!" << std::endl;


    for(int y = 0 ; y < ny ; y++)
    {
        for(int k = y; k < nvrow ; k++)
        {
            double4_t row = vd[y*nvrow];
            double4_t row2 = vd[k*nvrow];
            row = corr(row);
            row2 = corr(row2);
            std::cout << "Taalllaaaa!!!" << std::endl;
            double4_t pairwisemultip = row*row2;
            result[k+y*ny] = sum(pairwisemultip)/nx;
            std::cout << result[k+y*ny] << std::endl;
        }
    }
    std::free(vd);
}




/*
        double sum = 0;
        // Calculate the mean of the row
        double mean = sum/nx;
        
        // Normalize the row by subtracting mean from each of the elements
        // Normalize so that the sum of the squares is 1
        double stde = 0;
        for(x = 0; x < nx; x++)
        {
            stde+=pow(data[y*nx+x]-mean,2);
        }
        stde = sqrt(stde/nx);
        #pragma omp parallel for
        for(x = 0; x < nx ; x++)
        {
            mat.push_back((static_cast<double>(data[y*nx + x])-mean)/stde);
        }
    }


    #pragma omp parallel for  
    for(int c = 0; c < ny; c++)
    {
        for(int i = c; i < ny; i++)
        {
            double s = 0;
            double ss = 0;
            for(int j = 0; j < nx; j++)
            {
                double4_t row = vd[]
                ss += mat[j+c*nx]*mat[j+i*nx];
            }
            result[i+c*ny] = ss/nx;
        }   
    }
}

for(int i = 0; i < nx ; i++)
{
    double4_t row = vd[i];
    double4_t row2 = vt[i];
    row*row2
}

*/