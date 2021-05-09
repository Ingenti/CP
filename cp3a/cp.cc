#include "math.h"
#include <vector>
#include <stdlib.h>
#include <malloc.h>
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
constexpr double4_t dnan {
    infty, infty, infty, infty
};

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
    for(int i = 0; i < 4; i++)
    {
        if((a[i])==infty){continue;}
        sum += a[i];
    }
    return sum;
}

void corr(double4_t* a, int nvrow, int nx)
{
    double summa = 0;
    // Calculate the sum of the vectors of a row == sum of the row
    for(int i = 0; i < nvrow ; i++)
    {
        summa += sum(a[i]);
    }


    double mean = summa/(nx);
    double rowstde = 0;

    // Calculate stde of the row
    for(int x = 0; x < nvrow; x++)
    {
        rowstde += sum((a[x]-mean)*(a[x]-mean));
    }
    rowstde = sqrt(rowstde/nx);
    // Finally, modify the array
    for(int x = 0; x < nvrow; x++)
    {
        a[x] = (a[x]-mean)/rowstde;
    }
}


void correlate(int ny, int nx, const float *data, float *result) 
{
    //elements per vector, using doubles here which take 64 bits each (4x in total 256bit == vector registry size)
    constexpr int nb = 4;
    // vectors per input row
    int nvrow = static_cast<int>((ceil(static_cast<double>(nx)/static_cast<double>(nb))));

    double4_t* vd = double4_t_alloc(ny*nvrow);
    double4_t* vt = double4_t_alloc(ny*nvrow);
    int jj = nx%nb;
    int mimi = std::min(nb,nx);

    #pragma omp parallel for
    for(int y = 0 ; y < ny ; ++y)
    {
        for(int k = 0 ; k < nvrow ; ++k)
        {
            vd[nvrow*y+k] = dnan;
            vt[nvrow*y+k] = dnan;
            for(int x = 0 ; x < mimi ; ++x)
            {
                vd[nvrow*y+k][x] = static_cast<double>(data[y*nx + x + (k*4)]);
                vt[nvrow*y+k][x] = static_cast<double>(data[k*nx + x + (y*4)]);

            }
            
            if(k == (nvrow-2) && jj!=0)
            {
                vd[nvrow*y+k+1] = dnan;
                for(int t = 0; t < jj; t++)
                {
                    vt[nvrow*y+k+1][t] = static_cast<double>(data[y*nx + t + (nvrow-1)*4]);
                    vd[nvrow*y+k+1][t] = static_cast<double>(data[y*nx + t + (nvrow-1)*4]);
                }
                k++;
            }
        }    
    }

    #pragma omp parallel for
    for(int i = 0 ; i < ny ; ++i)
    {
        for(int j = i; j < ny ; ++j)
        {
            double4_t* row = double4_t_alloc(nvrow);
            double4_t* row2 = double4_t_alloc(nvrow);

            for(int c = 0; c < nvrow; ++c)
            {
                // Choose two rows to calculate corr for
                row[c] = vd[c+i*nvrow];
                row2[c] = vt[c+j*nvrow];
            }
            corr(row, nvrow, nx);
            corr(row2, nvrow, nx);
            double s = 0;
            for(int dvec = 0; dvec < nvrow; ++dvec)
            {
                
                s += sum(row[dvec]*row2[dvec]);
            }
            result[j+i*ny] = s/nx;
            std::free(row);
            std::free(row2);
        }
    }
    std::free(vt);
    std::free(vd);
}