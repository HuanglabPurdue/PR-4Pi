
#include "cuda_runtime.h"
#include "definitions.h"
#include "kernel.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


__global__ void kernel_calcPSFPixel(float *SampledPSF, int SizeX, int SizeY, float *PSFs, float *dPSFx, float *dPSFy, float *dPSFz, float * X, float *Y, float *Z,
	float SampleSpacingXY, float SampleSpacingZ, float StartX, float StartY, float StartZ, int N_int, int PSFSizeOut, int NPSFs)

{
	__shared__ float PSF_Row_Samples[64], dPSFx_Row_Samples[64], dPSFy_Row_Samples[64], dPSFz_Row_Samples[64];	// store the interpolated pixels in here
	__shared__ float PSF_Row_Sum[16], dPSFx_Row_Sum[16], dPSFy_Row_Sum[16], dPSFz_Row_Sum[16];		// This is the sum across rows and columns of interpolated pixels
	 
	__shared__ float theta[5];				// X, Y, Z, I, Bg

	int idX = threadIdx.x;
	int PixelNumber = threadIdx.y;
	float F2d[2], X1z[2], dF2d_y[2], dF2d_x[2];
	float a, b;	// coefficient for linear interpolation
	//Get XYZ
	if ((threadIdx.x == 0) && (threadIdx.y == 0)){
		theta[0] = X[blockIdx.y] * SampleSpacingXY * N_int;
		theta[1] = Y[blockIdx.y] * SampleSpacingXY * N_int;
		theta[2] = Z[blockIdx.y];
	}

	__syncthreads();

	float Y_thread = (N_int*PixelNumber + idX + 0.5)*SampleSpacingXY - theta[1];
	
	int YBaseIndex = floor((Y_thread - StartY) / SampleSpacingXY);

	float SampleSpacingXYInv = 1 / SampleSpacingXY;
	float SampleSpacingZInv = 1 / SampleSpacingZ;
	int idZ = round((theta[2] - StartZ) * SampleSpacingZInv);


	//intialize PSF_Row_Sum counter
	if (threadIdx.x == 0)
	{
		PSF_Row_Sum[threadIdx.y] = 0;
		dPSFx_Row_Sum[threadIdx.y] = 0;
		dPSFy_Row_Sum[threadIdx.y] = 0;
		dPSFz_Row_Sum[threadIdx.y] = 0;
	}
	for (int ii = 0; ii < N_int; ii++) //go right in row
	{
		float X_thread = (blockIdx.x*N_int + ii + 0.5)*SampleSpacingXY - theta[0];


		//for interpolation we need the four surrounding points
		
		int XBaseIndex = floor((X_thread - StartX) * SampleSpacingXYInv);

		//using the follwing notation:
		//X1    X2
		//   o
		//X3	X4

		//These are values of the sampled points.
		for (int nn = 0; nn < 2; nn++)
		{

			int tmp2 = SizeY*SizeX*(idZ+nn) + SizeY*XBaseIndex + YBaseIndex;
			float F1 = SampledPSF[tmp2];
			float F2 = SampledPSF[tmp2 + SizeX];
			float F3 = SampledPSF[tmp2 + 1];
			float F4 = SampledPSF[tmp2 + SizeX + 1];

			//These are locations of the sampled points
			float X1x = XBaseIndex*SampleSpacingXY + StartX;
			float X1y = YBaseIndex*SampleSpacingXY + StartY;
			X1z[nn] = (idZ + nn)*SampleSpacingZ + StartZ;
			//Bilinear interpolation
			gencoeff(F1, F2, SampleSpacingXY, X1x, X1x + SampleSpacingXY,&a,&b);
			float X1X2 = evallinear(a,b,X_thread);
			float dX1X2_x = a;
			gencoeff(F3, F4, SampleSpacingXY, X1x, X1x + SampleSpacingXY, &a, &b);
			float X3X4 = evallinear(a, b, X_thread);
			float dX3X4_x = a;
			gencoeff(X1X2, X3X4, SampleSpacingXY, X1y, X1y + SampleSpacingXY, &a, &b);
			F2d[nn] = evallinear(a, b, Y_thread);
			dF2d_y[nn] = a;
			gencoeff(dX1X2_x, dX3X4_x, SampleSpacingXY, X1y, X1y + SampleSpacingXY, &a, &b);
			dF2d_x[nn] = evallinear(a, b, Y_thread);
		}
		gencoeff(F2d[0], F2d[1], SampleSpacingZ, X1z[0], X1z[1], &a, &b);
		PSF_Row_Samples[N_int*PixelNumber + idX] = evallinear(a, b, theta[2]);
		dPSFz_Row_Samples[N_int*PixelNumber + idX] = a;

		gencoeff(dF2d_x[0], dF2d_x[1], SampleSpacingZ, X1z[0], X1z[1], &a, &b);
		dPSFx_Row_Samples[N_int*PixelNumber + idX] = evallinear(a, b, theta[2]);

		gencoeff(dF2d_y[0], dF2d_y[1], SampleSpacingZ, X1z[0], X1z[1], &a, &b);
		dPSFy_Row_Samples[N_int*PixelNumber + idX] = evallinear(a, b, theta[2]);

		
		
		__syncthreads();
		//now sum over the row
		if (threadIdx.x == 0) 
		for (int jj = 0; jj < N_int; jj++)
		{
			PSF_Row_Sum[threadIdx.y] += PSF_Row_Samples[N_int*PixelNumber + jj];
			dPSFx_Row_Sum[threadIdx.y] += dPSFx_Row_Samples[N_int*PixelNumber + jj];
			dPSFy_Row_Sum[threadIdx.y] += dPSFy_Row_Samples[N_int*PixelNumber + jj];
			dPSFz_Row_Sum[threadIdx.y] += dPSFz_Row_Samples[N_int*PixelNumber + jj];
		}
	}

	//now return value for each pixel

	__syncthreads();
	if (threadIdx.x == 0)
	{
		PSFs[PSFSizeOut*PSFSizeOut*blockIdx.y + PSFSizeOut*blockIdx.x + threadIdx.y] = PSF_Row_Sum[threadIdx.y];
		dPSFx[PSFSizeOut*PSFSizeOut*blockIdx.y + PSFSizeOut*blockIdx.x + threadIdx.y] = dPSFx_Row_Sum[threadIdx.y];
		dPSFy[PSFSizeOut*PSFSizeOut*blockIdx.y + PSFSizeOut*blockIdx.x + threadIdx.y] = dPSFy_Row_Sum[threadIdx.y];
		dPSFz[PSFSizeOut*PSFSizeOut*blockIdx.y + PSFSizeOut*blockIdx.x + threadIdx.y] = dPSFz_Row_Sum[threadIdx.y];
	}
	

}

__global__ void kernel_calcPSFI(float *psf, float *psfI, float *I, float *bg, int Nfit, int PSFsize)
{
	const int tx = threadIdx.x;
	const int bx = blockIdx.x;
	const int BlockSize = blockDim.x;

	//Prevent read/write past end of array
	int j = BlockSize*bx + tx;
	if ((bx*BlockSize + tx) >= Nfit) return;

	for (int i = 0; i < PSFsize; i++)
	{
		psfI[j*PSFsize + i] = psf[j*PSFsize + i] * I[j] + bg[j];

	}
}

__device__ void gencoeff(float f1, float f2, float dx, float x1, float x2,float *a, float *b)
{
	a[0] = (f2 - f1) / dx;
	b[0] = (x2*f1 - x1*f2) / dx;
}

__device__ float evallinear(float a, float b, float x)
{
	float f;
	f = x*a + b;
	return f;
}