/*
 * GaussianFilter.h
 *
 *  Created on: Jun 14, 2017
 *      Author: nvduc
 */

#ifndef CLIENT_GAUSSIANFILTER_H_
#define CLIENT_GAUSSIANFILTER_H_

double** filter(int size, double sigma);
double** applyFilter(double** img, double** gk, int img_size, int ksize);



#endif /* CLIENT_GAUSSIANFILTER_H_ */
