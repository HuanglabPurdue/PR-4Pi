/*!
 * \file kernel.h
 * \author <your name>
 * \date April 20, 2012
 * \brief  Put prototypes of all Cuda Kernels here.
 */
#include "Model_Localization.hpp"
//#ifndef max
//#define max(a,b) (((a) > (b)) ? (a) : (b))
//#endif
//#ifndef min
//#define min(a,b) (((a) < (b)) ? (a) : (b))
//#endif
#ifndef KERNEL_H
#define KERNEL_H

__global__ void kernel_calcPSFPixel(float *SampledPSF, int SizeX, int SizeY, float *PSFs, float *dPSFx, float *dPSFy, float *dPSFz,
	float * X, float *Y, float *Z, 
	float SampleSpacingXY, float SampleSpacingZ, float StartX, float StartY, float StartZ, int N_int, int PSFSizeOut, int NPSFs);

__global__ void kernel_calcPSFI(float *psf, float *psfI, float *I, float *bg, int Nfit, int PSFsize);

__global__ void kernel_Localization(float *ParamIn, float *ParamNext, float *Convergence, float *FirstDev, float *SecondDev,
	int Nfit, int N_int, int FitBoxCenter, float lambda, float SampleSpacingXY);

__global__ void kernel_getdev(float *data, float *gainR, float *PSF, float *dPSFx, float *dPSFy, float *dPSFz, float *I, float *bg, int Nfit, int PSFsize,
	float *FirstDev, float *SecondDev);

__global__ void kernel_calCRLB(float *ParamF, float *ParamVar, int Nfit);
	
__global__ void kernel_calFisherM(float *PSF, float *dPSFx, float *dPSFy, float *dPSFz, float *I, float *bg, float *gainR, float *ParamF,
	int Q, int Nfit, int PSFSize);

__global__ void kernel_Err(float *PSF, float *Data, float *gainR, float *Error, int Nfit, int PSFSize);

__device__ void gencoeff(float f1, float f2, float dx, float x1, float x2, float *a, float *b);

__device__ float evallinear(float a, float b, float x);

__device__ void fundev(float *data, float *gainR, float *psf, float *dpsf, float I,float Id, float bg, float *dL, float *dL2, int PSFsize);
#endif