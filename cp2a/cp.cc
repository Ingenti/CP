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

double calculate_mean(double *temp)
{
    double sum = 0;
    int counter = 0;
    for(counter;counter < sizeof(temp)/sizeof(temp[0]);counter++)
    {
        sum += temp[counter];
    }

    return (sum/counter);
}


void correlate(int ny, int nx, const float *data, float *result) 
{

    std::vector<double> mat;

    int y = 0;
    int x = 0;
    for(y = 0 ; y < ny ; y++)
    {
        // Save a row in temp
        double temp[sizeof(data)/(sizeof(data[0])*nx)] = {};
        
        // Cast to double precision and save row as temp
        for (x = 0 ; x < nx ; x++)
        {
            temp[x] = static_cast<double>(data[y*nx + x]);
        }

        // Calculate the mean of the row
        double mean = calculate_mean(temp);
        
        // Normalize the row by subtracting mean from each of the elements
        // Normalize so that the sum of the squares is 1
        // Sum of the elements == mean*nx
        double squaresum = pow(mean*nx,2);

        for(x = 0; x < nx ; x++)
        {
            // Save values to mat
            mat.push_back((temp[x] - mean)/squaresum);
        }
    }

    double div = sizeof(result)/sizeof(result[0]);
    #pragma omp parallel for
    for(int i = 0; i<ny-1; i++)
    {
        for(int j = 0; j < nx; j++)
        {
            result[i+j*ny] = mat[j+i*nx]*mat[j+(i+1)*nx]/(div);
        }
    }
}
    

int main()

{
    
    return 0;

}