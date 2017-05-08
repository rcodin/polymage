#include <iostream>
#include <fstream>

#include "laplacian.h"

#define NREPS 1
#define NRUNS 1

int main ()
{
    const int R = 2560, C = 1536;
    const float alpha = 0.142857, beta = 1.000000;
    
    std::ifstream img_file ("img.bin", std::ios::in | std::ios::binary);
    unsigned short int * img_colour = new unsigned short int[3*(R+4)*(C+4)];
    img_file.read ((char *)img_colour, 3*(R+4)*(C+4)*sizeof (unsigned short int));
    img_file.close();
    
    short unsigned int * laplacian = new short unsigned int[3*R*C];
    
    for (int i = 0; i < NREPS; i++)
    {
        double sum_time = 0;
        
        for (int j = 0; j < NRUNS; j++)
        {
            clock_t t1 = clock();
         
            pipeline_laplacian (C, R, alpha, beta, (void *)img_colour, (void *)laplacian);
            
            clock_t t2 = clock();
            
            //std::cout << "Time taken " << (((double)(t2-t1)))/(CLOCKS_PER_SEC/1000.0f) << " ms" << std::endl;
            sum_time += (t2-t1);
        }
        
        std::cout << "Average Time Taken "<< (((double)(sum_time))/(CLOCKS_PER_SEC/1000.0f))/NRUNS << " ms"<< std::endl;
    }
    
    delete img_colour;
    delete laplacian;
}
