/*
 * GaussianFilter.cpp
 *
 *  Created on: Jun 14, 2017
 *      Author: nvduc
 */
#include <iostream>
#include <cmath>
#include <iomanip>
#include "GaussianFilter.h"

double** filter(int size, double sigma){
    double r, s = 2.0 * sigma * sigma;  // Assigning standard deviation to 1.0
    double sum = 0.0;   // Initialization of sun for normalization
    double**  kernel = new double*[size];
    for(int i=0; i < size; i++)
    	kernel[i] = new double[size];
    for (int x = -size/2; x <= size/2; x++) // Loop to generate 5x5 kernel
    {
        for(int y = -size/2; y <= size/2; y++)
        {
        		r = sqrt(x*x + y*y);
        		kernel[x+ size/2][y + size/2] = (exp(-(r*r)/s))/(M_PI * s);
        		sum += kernel[x + size/2][y + size/2];
        }
    }
    for(int i = 0; i < size; ++i) // Loop to normalize the kernel
        for(int j = 0; j < size; ++j)
            kernel[i][j] /= sum;
    return kernel;
}
double** applyFilter(double** img, double** kernel, int img_size, int ksize){
	// find center position of kernel (half of kernel size)
	int kCols = ksize;
	int kRows = ksize;
	int kCenterX = kCols / 2;
	int kCenterY = kRows / 2;
	int i,j;
	int rows = img_size + 2;
	int cols = img_size + 2;
	int m,n;
	int mm, nn;
	int ii,jj;
	//
	double** out = new double*[img_size];
	for(i=0; i < img_size; i++)
		out[i] = new double[img_size];
	// padding
	double** img_pad = new double*[rows];
	for(i=0; i < rows; i++){
		img_pad[i] = new double[cols];
		for(j=0; j < cols; j++){
			if(i == 0 || j == 0 || i==rows-1 || j == cols -1){
				img_pad[i][j] = 0;
			}else{
				img_pad[i][j] = img[i-1][j-1];
			}
		}
	}
	double out_tmp[rows][cols];
	// covolution
	for(i=0; i < rows; ++i)              // rows
	{
	    for(j=0; j < cols; ++j)          // columns
	    {
	    	out_tmp[i][j] = 0;
	    	// calcuate out[i][j]
	        for(m=0; m < kRows; ++m)     // kernel rows
	        {
	            mm = kRows - 1 - m;      // row index of flipped kernel

	            for(n=0; n < kCols; ++n) // kernel columns
	            {
	                nn = kCols - 1 - n;  // column index of flipped kernel

	                // index of input signal, used for checking boundary
	                ii = i + (m - kCenterY);//
	                jj = j + (n - kCenterX);
	                // ignore input samples which are out of bound
	                if( ii >= 0 && ii < rows && jj >= 0 && jj < cols )
	                    out_tmp[i][j] += img_pad[ii][jj] * kernel[mm][nn];
	            }
	        }
	    }
	}
// extract invisible probability
	for(i=0; i < img_size; i++)
		for(j=0; j < img_size; j++){
			out[i][j] = out_tmp[i+1][j+1];
			if(j==0)
				out[i][j] = out_tmp[i+1][j+1] + out_tmp[i+1][0];
			if(j == img_size-1)
				out[i][j] = out_tmp[i+1][j+1] + out_tmp[i+1][rows-1];

		}
	return out;
}

