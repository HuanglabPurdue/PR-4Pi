#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda_runtime.h>
#include <math.h>
#include "definitions.h"
#include <mex.h>

extern "C" void kernel_Localization_wrapper(dim3, dim3, float *ParamIn, float *ParamNext, float *Convergence, float *FirstDev, float *SecondDev,
	int Nfit, int N_int, int FitBoxCenter, float lambda, float SampleSpacingXY);

extern "C" void kernel_getdev_wrapper(dim3, dim3, float *data, float *gainR, float *PSF, float *dPSFx, float *dPSFy, float *dPSFz, float *I, float *bg, int Nfit, int PSFsize,
	float *FirstDev, float *SecondDev);

extern "C" void kernel_calCRLB_wrapper(dim3, dim3, float *ParamF, float *ParamVar, int Nfit);

extern "C" void kernel_calFisherM_wrapper(dim3, dim3, float *PSF, float *dPSFx, float *dPSFy, float *dPSFz, float *I, float *bg, float *gainR, float *ParamF,
	int Q, int Nfit, int PSFSize);

extern "C" void kernel_calcPSFPixel_wrapper(dim3, dim3, float *SampledPSF, int SizeX, int SizeY, float *PSFs, float *dPSFx, float *dPSFy, float *dPSFz,
	float * X, float *Y, float *Z,
	float SampleSpacingXY, float SampleSpacingZ, float StartX, float StartY, float StartZ, int N_int, int PSFSizeOut, int NPSFs);

extern "C" void kernel_calcPSFI_wrapper(dim3, dim3, float *psf, float *psfI, float *I, float *bg, int Nfit, int PSFsize);

extern "C" void kernel_Err_wrapper(dim3, dim3, float *PSF, float *Data, float *gainR, float *Error, int Nfit, int PSFSize);

class Model_Localization{
public:
	float * Data;
	float * ParamIn;
	float * ParamNext;
	float * ParamVar;
	float * Convergence;
	float * Error;
	float * SampledPSF;
	float * X;
	float * Y;
	float * Z;
	float * I;
	float * Bg;
	float * PSFs;
	float * GainR;

	float * d_Data;
	float * d_ParamIn;
	float * d_ParamNext;
	float * d_ParamVar;
	float * d_Convergence;
	float * d_Error;
	float * d_FirstDev;
	float * d_SecondDev;
	float * d_ParamF;
	float * d_SampledPSF;
	float * d_PSFs;
	float * d_PSFI;
	float * d_dPSFx;
	float * d_dPSFy;
	float * d_dPSFz;
	float * d_GainR;

	float * d_X;
	float * d_Y;
	float * d_Z;
	float * d_I;
	float * d_Bg;

	int Nparam;
	int Nfit;
	int PSFSizeOut = 16;
	int PSFSize = PSFSizeOut*PSFSizeOut;

	int SampledPSF_XYPixels;  // Size of sampled PSF
	int SampledPSF_ZPixels;

	int FitBoxCenter;
	float Lambda = 1.0; // damping factor Levenberg-Marquardt algorithm
	float Deltax;
	float Deltaz;
	int NPSFs;


	float SampleSpacingXY;
	float SampleSpacingZ;
	float StartX;
	float StartY;
	float StartZ;



	void setup();
	void sendParams(float *, float *, float *, float *, float *, float *, float *);
	void calPSF();
	void mLocalization();
	void calCRLB();
	void calErr();
	void cleanUp();
	void CudaError(const char *);


};