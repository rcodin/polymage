#include <iostream>
#include <fstream>

#include "laplacian.h"
#include <opencv2/opencv.hpp>

using namespace cv;

int main ()
{
    const int R = 2560, C = 1536;
    const float alpha = 0.142857, beta = 1.000000;
    Mat mat = imread ("img.png"); // BGR
    std::cout << "Channels " << mat.channels () << " Cols " << mat.cols
        << " Rows " << mat.rows << std::endl;
    
    assert (mat.channels () == 3 && mat.cols == C && mat.rows == R);
    
    unsigned short int * img_colour = new unsigned short int[3*(R+4)*(C+4)];
    
    //Copy image to array and bring channel dimension outside 
    for (int i = 0; i < R; i++)
    {
        for (int j = 0; j < C; j++)
        {
            for (int c = 0; c < 3; c++)
                img_colour [c*(R+4)*(C+4)+(i+2)*(C+4) + (j+2)] = mat.at<Vec3b> (i, j)[c];
        }
    }
    
    short unsigned int * laplacian = new short unsigned int[3*R*C];
    
    for (int i = 0; i < 5; i++)
    {
        double sum_time = 0;
        
        for (int j = 0; j < 10; j++)
        {
            clock_t t1 = clock ();
         
            pipeline_laplacian (C, R, alpha, beta, (void *)img_colour, (void *)laplacian);
            
            clock_t t2 = clock ();
            
            //std::cout << "Time taken " << (((double)(t2-t1)))/(CLOCKS_PER_SEC/1000.0f) << " ms" << std::endl;
            sum_time += (t2-t1);
        }
        
        std::cout << "Average Time Taken "<< (((double)(sum_time))/(CLOCKS_PER_SEC/1000.0f))/500 << " ms"<< std::endl;
    }
    
    delete img_colour;
    delete laplacian;
}
