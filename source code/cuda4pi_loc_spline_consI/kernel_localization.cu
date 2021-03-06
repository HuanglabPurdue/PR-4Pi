#include "cuda_runtime.h"
#include "definitions.h"
#include "kernel.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


__global__ void kernel_Localization(float *ParamIn, float *ParamNext, float *Convergence, float *FirstDev, float *SecondDev,
	int Nfit, int N_int, int FitBoxsize, float lambda, float SampleSpacingXY)
{
	const int tx = threadIdx.x;
	const int bx = blockIdx.x;
	const int BlockSize = blockDim.x;

	//Prevent read/write past end of array
	int j = BlockSize*bx + tx;
	if ((bx*BlockSize + tx) >= Nfit) return;

	
	float stepLimit[NP] = {0.03f, 0.03f, 0.04f, 400, 2}; // x,y,z step limits are in micron
	float x0_next[NPL];
	float dL_pos = 0, dL2_pos = 0;
	float dL_I, dL2_I; // photon and background
	float step[NPL];
	float rate = 1/(1 + lambda);
	float tmp;
	int s, p, k;
	// x,y,z
	for (p = 0; p < NP; p++)
	{
		for (s = 0; s < 4; s++)
		{
			dL_pos += FirstDev[s*NP*Nfit + j*NP + p];
			dL2_pos += SecondDev[s*NP*Nfit + j*NP + p];
		}
		tmp = -1 * dL_pos / dL2_pos * rate;
		step[p] = fminf(fmaxf(tmp, -stepLimit[p]), stepLimit[p]);
	}

	

	x0_next[0] = ParamIn[NPL*j + 0] + step[0] * (-1 / SampleSpacingXY / N_int);
	x0_next[1] = ParamIn[NPL*j + 1] + step[1] * (-1 / SampleSpacingXY / N_int);
	for (k = 2; k < NPL; k++)
	{		
		x0_next[k] = ParamIn[NPL*j + k] + step[k];
	}

	
	x0_next[3] = (x0_next[3] <= 100 ? 100 : x0_next[3]); // intensity is not less than 100	
	x0_next[4] = (x0_next[4] <= 0 ? 0.01f : x0_next[4]);// bg is not less than 0

	x0_next[0] = fminf(fmaxf(x0_next[0], 4), FitBoxsize - 4);// xy shift is within fitting box
	x0_next[1] = fminf(fmaxf(x0_next[1], 4), FitBoxsize - 4);
	x0_next[2] = fminf(fmaxf(x0_next[2], -1.4), 1.4);//z position is within -1.4 to 1.4 um
	
	for (k = 0; k < NPL; k++) {
		ParamNext[NPL*j + k] = x0_next[k];
		Convergence[NPL*j + k] = x0_next[k] - ParamIn[NPL*j + k];
	}

}

__global__ void kernel_getdev(float *data, float *gainR, float *PSF, float *dPSFx, float *dPSFy, float *dPSFz, float *I, float *bg, int Nfit, int PSFsize,
	float *FirstDev, float *SecondDev)
{
	const int tx = threadIdx.x;
	const int bx = blockIdx.x;
	const int BlockSize = blockDim.x;
	
	//Prevent read/write past end of array
	int j = BlockSize*bx + tx;
	if ((bx*BlockSize + tx) >= Nfit) return;
	
	float dL[NP], dL2[NP];
	float psfI;
	int k, i;
	for (k = 0; k < NP; k++)
	{
		dL[k] = 0;
		dL2[k] = 0;
	}
	fundev(&data[j*PSFsize], &gainR[j*PSFsize], &PSF[j*PSFsize], &dPSFx[j*PSFsize], I[j], I[j], bg[j], &dL[0], &dL2[0], PSFsize);
	fundev(&data[j*PSFsize], &gainR[j*PSFsize], &PSF[j*PSFsize], &dPSFy[j*PSFsize], I[j], I[j], bg[j], &dL[1], &dL2[1], PSFsize);
	fundev(&data[j*PSFsize], &gainR[j*PSFsize], &PSF[j*PSFsize], &dPSFz[j*PSFsize], I[j], I[j], bg[j], &dL[2], &dL2[2], PSFsize);
	fundev(&data[j*PSFsize], &gainR[j*PSFsize], &PSF[j*PSFsize], &PSF[j*PSFsize], I[j], 1.0, bg[j], &dL[3], &dL2[3], PSFsize);
	for (int i = 0; i < PSFsize; i++)
	{
		psfI = PSF[j*PSFsize + i] * I[j] + bg[j] + gainR[j*PSFsize + i];
		dL[4] += (data[j*PSFsize + i] + gainR[j*PSFsize + i]) / psfI - 1;
		dL2[4] += -1 * (data[j*PSFsize + i] + gainR[j*PSFsize + i]) / psfI / psfI;

	}

	for (int k = 0; k < NP; k++)
	{
		FirstDev[NP * j + k] = dL[k];
		SecondDev[NP * j + k] = dL2[k];
	}
}

__device__ void fundev(float *data, float *gainR, float *psf, float *dpsf, float I, float Id, float bg, float *dL, float *dL2, int PSFsize)
{
	float psfI;
	for (int i = 0; i < PSFsize; i++)
	{
		psfI = psf[i] * I + bg + gainR[i];
		dL[0] += ((data[i] + gainR[i]) / psfI - 1) * dpsf[i] * Id;
		dL2[0] += -1 * Id * Id * dpsf[i] * dpsf[i] * (data[i] + gainR[i]) / psfI / psfI;
	}
}