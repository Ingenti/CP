#include "math.h"
#include <vector>

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/



void correlate(int ny, int nx, const float *data, float *result) 
{

    std::vector<double> mat = {};
    int x = 0;
    for(int y = 0 ; y < ny ; y++)
    {
        double sum = 0;
        for(x = 0 ; x < nx ; x++)
        {
            sum += static_cast<double>(data[y*nx + x]);
        }

        // Calculate the mean of the row
        double mean = (sum/nx);
        
        // Normalize the row by subtracting mean from each of the elements
        // Normalize so that the sum of the squares is 1
        double stde = 0;
        for(x = 0; x < nx; x++)
        {
            stde+=pow(data[y*nx+x]-mean,2);
        }
        stde = sqrt(stde/nx);

        for(x = 0; x < nx ; x++)
        {
            mat.push_back((static_cast<double>(data[y*nx + x])-mean)/stde);
        }
    }

    for(int c = 0; c < ny; c+=2)
    {
        for(int i = c; i < ny; i+=2)
        {
            double s = 0;
            double ss = 0;
            double sss = 0;
            for(int j = 0; j < nx; j++)
            {
                ss += mat[j+c*nx]*mat[j+i*nx];
                s += mat[j+(c+1)*nx]*mat[j+(i+1)*nx];
                sss += mat[j+c*nx]*mat[j+(i+1)*nx];
            }
            result[i+c*ny] = ss/nx;
            result[(i+1)+(c+1)*ny] = s/nx;
            result[(i+1)+c*ny] = sss/nx;

        }   
    }
}
