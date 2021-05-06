#include "math.h"
#include <vector>
#include <stdlib.h>
#include <malloc.h>
#include <iostream>
#include <limits>
/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
typedef double double4_t __attribute__ ((vector_size (4 * sizeof(double))));
constexpr double infty = std::numeric_limits<double>::infinity();

constexpr double4_t d4infty {
    infty, infty, infty, infty,
};

static double4_t* double4_t_alloc(std::size_t n) {
    void* tmp = 0;
    if (posix_memalign(&tmp, sizeof(double4_t), sizeof(double4_t) * n)) {
        throw std::bad_alloc();
    }
    return (double4_t*)tmp;
}



double sum(double4_t a, int nx)
{
    double sum = 0;
    for(int i = 0; i < std::min(4,nx); i++)
    {
        sum += a[i];
        std::cout << a[i] << std::endl;
    }
    return sum;
}

void corr(double4_t* a, int nvrow, int nx)
{
    double summa = 0;
    int len = nvrow;

    // Calculate the sum of the vectors of a row == sum of the row
    for(int i = 0; i < nvrow ; i++)
    {
        summa += sum(a[i], nx);
    }


    double mean = summa/len;
    double stde = 0;
    double rowstde = 0;

    // Calculate stde of the row
    for(int x = 0; x < len; x++)
    {
        stde = sum((a[x]-mean)*(a[x]-mean), nx);
        rowstde += stde;
    }


    // Finally, modify the array
    for(int x = 0; x < len; x++)
    {
        a[x] = (a[x]-mean)/rowstde;
    }
}


void correlate(int ny, int nx, const float *data, float *result) 
{
    //elements per vector, using doubles here which take 64 bits each (4x in total 256bit == vector registry size)
    constexpr int nb = 4;
    // vectors per input row
    //float nrow = ceil(nx/nb);
    int nvrow = (ny+nb-1)/nb;
    std::cout << "Number of vectors per row: "<< nvrow << std::endl;

    double4_t* vd = double4_t_alloc(ny*nvrow);
    double4_t* vt = double4_t_alloc(ny*nvrow);

    for(int y = 0 ; y < ny ; ++y)
    {
        for(int k = 0 ; k < nvrow ; ++k)
        {
            for(int x = 0 ; x < std::min(nb,nx) ; ++x)
            {
                vd[nvrow*y+k][x] = static_cast<double>(data[y*nx + x]);
                vd[nvrow*y+k][x] = static_cast<double>(data[y*nx + x]);
            }
        }    
    }

    for(int i = 0 ; i < ny ; i++)
    {
        for(int j = i; j < ny ; j++)
        {
            double4_t* row = double4_t_alloc(nvrow);
            double4_t* row2 = double4_t_alloc(nvrow);
            for(int c = 0; c < nvrow; c++)
            {
                // Choose two rows to calculate corr for
                row[c] = vd[c+i*nvrow];
                row2[c] = vt[c+j*nvrow];
            }
            for(int rivi = 0; rivi < nvrow; rivi++)
            {
                corr(row, nvrow, nx);
                corr(row2, nvrow, nx);
                result[j+i*ny] += sum(row[rivi]*row2[rivi], nx);
            }
            result[j+i*ny] /= nx;
            std::free(row);
            std::free(row2);
        }
    }
    std::free(vt);
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