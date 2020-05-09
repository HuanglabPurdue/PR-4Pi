#include "Model_Localization.hpp"

void Model_Localization::setup(){

	//allocate memory on device
	cudaMalloc(&d_PSFs, Nfit*PSFSize*sizeof(float));
	cudaMalloc(&d_PSFI, Nfit*PSFSize*sizeof(float));
	cudaMalloc(&d_dPSFx, Nfit*PSFSize*sizeof(float));
	cudaMalloc(&d_dPSFy, Nfit*PSFSize*sizeof(float));
	cudaMalloc(&d_dPSFz, Nfit*PSFSize*sizeof(float));
	cudaMalloc(&d_Data, 4*Nfit*PSFSize*sizeof(float));
	cudaMalloc(&d_FirstDev, 4 * NP*Nfit*sizeof(float));
	cudaMalloc(&d_SecondDev, 4 * NP*Nfit*sizeof(float));
	cudaMalloc(&d_ParamIn, Nfit*NPL*sizeof(float));
	cudaMalloc(&d_ParamVar, Nfit*NPL*sizeof(float));
	cudaMalloc(&d_ParamF, 4*Nfit*NPL*NPL*sizeof(float));
	cudaMalloc(&d_ParamNext, Nfit*NPL*sizeof(float));
	cudaMalloc(&d_Convergence, Nfit*NPL*sizeof(float));
	cudaMalloc(&d_Error, Nfit * 2 * sizeof(float));
	cudaMalloc(&d_GainR, 4 * Nfit*PSFSize*sizeof(float));

	cudaMalloc(&d_SamplePSF, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float));
	cudaMalloc(&d_SamplePSFx, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float));
	cudaMalloc(&d_SamplePSFy, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float));
	cudaMalloc(&d_SamplePSFz, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float));
	cudaMalloc(&d_SamplePSFxy, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float));
	cudaMalloc(&d_SamplePSFxz, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float));
	cudaMalloc(&d_SamplePSFyz, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float));
	cudaMalloc(&d_SamplePSFxyz, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float));
	cudaMalloc(&d_X, Nfit*sizeof(float));
	cudaMalloc(&d_Y, Nfit*sizeof(float));
	cudaMalloc(&d_Z, Nfit*sizeof(float));
	cudaMalloc(&d_I, 4*Nfit*sizeof(float));
	cudaMalloc(&d_Bg, 4*Nfit*sizeof(float));
	cudaMalloc(&d_transx, 4 * Nfit*sizeof(float));
	cudaMalloc(&d_transy, 4 * Nfit*sizeof(float));
	CudaError("Malloc");

	//copy from host to device
	cudaMemcpy(d_transx, transx, 4 * Nfit*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_transy, transy, 4 * Nfit*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Data, Data, 4*Nfit*PSFSize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_GainR, GainR, 4 * Nfit*PSFSize*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_SamplePSF, SamplePSF, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_SamplePSFx, SamplePSFx, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_SamplePSFy, SamplePSFy, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_SamplePSFz, SamplePSFz, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_SamplePSFxy, SamplePSFxy, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_SamplePSFxz, SamplePSFxz, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_SamplePSFyz, SamplePSFyz, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_SamplePSFxyz, SamplePSFxyz, 4 * SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels*sizeof(float), cudaMemcpyHostToDevice);
	CudaError("Memcpy");

}

void Model_Localization::sendParams(float *X, float *Y, float *Z, float *I, float *Bg, float *ParamIn, float *Error){
	//copy from host to device
	cudaMemcpy(d_X, X, Nfit*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Y, Y, Nfit*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Z, Z, Nfit*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_I, I, 4*Nfit*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Bg, Bg, 4*Nfit*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_ParamIn, ParamIn, Nfit*NPL*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Error, Error, 2 * Nfit*sizeof(float), cudaMemcpyHostToDevice);
	CudaError("Memcpy");
}

void Model_Localization::cleanUp(){
	//free all memory we allocated. 
	cudaFree(d_PSFs);
	cudaFree(d_PSFI);
	cudaFree(d_dPSFx);
	cudaFree(d_dPSFy);
	cudaFree(d_dPSFz);
	cudaFree(d_Data);
	cudaFree(d_ParamIn);
	cudaFree(d_ParamNext);
	cudaFree(d_ParamVar);
	cudaFree(d_ParamF);
	cudaFree(d_FirstDev);
	cudaFree(d_SecondDev);
	cudaFree(d_Convergence);
	cudaFree(d_Error);
	cudaFree(d_X);
	cudaFree(d_Y);
	cudaFree(d_Z);
	cudaFree(d_I);
	cudaFree(d_Bg);
	cudaFree(d_transx);
	cudaFree(d_transy);
	cudaFree(d_SamplePSF);
	cudaFree(d_SamplePSFx);
	cudaFree(d_SamplePSFy);
	cudaFree(d_SamplePSFz);
	cudaFree(d_SamplePSFxy);
	cudaFree(d_SamplePSFxz);
	cudaFree(d_SamplePSFyz);
	cudaFree(d_SamplePSFxyz);
	cudaFree(d_GainR);
	CudaError("Cuda Free");

}

void Model_Localization::mLocalization(){
	//PSF model generation
	int N_int = 4;		//There are N_int*N_int number of interpoloted points. 
	int BlockSize = BSZ;
	int NBlock = Nfit / BlockSize + 1;
	dim3 dimBlock(N_int, PSFSizeOut);
	dim3 dimGrid(PSFSizeOut, Nfit);
	dim3 dimBlock2(BlockSize);
	dim3 dimGrid2(NBlock);

	CudaError("mLocalization::setup Dims");
	//mexPrintf("%d %d %d \n", N_int, PSFSizeOut, NPSFs);

	//kernel calls
	int Samplepsf_Pixels = SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels;
	for (int s = 0; s < 4; s++)
	{
		kernel_calcPSFPixel_wrapper(dimGrid, dimBlock, &d_SamplePSF[s*Samplepsf_Pixels], &d_SamplePSFx[s*Samplepsf_Pixels], &d_SamplePSFy[s*Samplepsf_Pixels], &d_SamplePSFz[s*Samplepsf_Pixels],
			&d_SamplePSFxy[s*Samplepsf_Pixels], &d_SamplePSFxz[s*Samplepsf_Pixels], &d_SamplePSFyz[s*Samplepsf_Pixels], &d_SamplePSFxyz[s*Samplepsf_Pixels], SampledPSF_XYPixels, SampledPSF_XYPixels,
			d_PSFs, d_dPSFx, d_dPSFy, d_dPSFz, d_X, d_Y, d_Z, &d_transx[s*Nfit], &d_transy[s*Nfit],
			SampleSpacingXY, SampleSpacingZ, StartX, StartY, StartZ, N_int, PSFSizeOut, Nfit);
		CudaError("mLocalization::kernel_calcPSFPixel_wrapper");

		kernel_getdev_wrapper(dimGrid2, dimBlock2, &d_Data[s*Nfit*PSFSize], &d_GainR[s*Nfit*PSFSize], d_PSFs, d_dPSFx, d_dPSFy, d_dPSFz, 
			&d_I[s*Nfit], &d_Bg[s*Nfit], Nfit, PSFSize, &d_FirstDev[s*NP*Nfit], &d_SecondDev[s*NP*Nfit]);

		CudaError("mLocalization::kernel_getdev_wrapper");

	}

	//localization

	kernel_Localization_wrapper(dimGrid2, dimBlock2, d_ParamIn, d_ParamNext, d_Convergence, d_FirstDev, d_SecondDev,
		Nfit, N_int, PSFSizeOut, Lambda, SampleSpacingXY);

	cudaMemcpy(ParamNext, d_ParamNext, Nfit*NPL*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(Convergence, d_Convergence, Nfit*NPL*sizeof(float), cudaMemcpyDeviceToHost);

	CudaError("mLocalization::kernel_Localization_wrapper");

}

void Model_Localization::calPSF(){
	//PSF model generation
	int N_int = 4;		//There are N_int*N_int number of interpoloted points. 
	int BlockSize = BSZ;
	int NBlock = Nfit / BlockSize + 1;
	dim3 dimBlock(N_int, PSFSizeOut);
	dim3 dimGrid(PSFSizeOut, Nfit);
	dim3 dimBlock2(BlockSize);
	dim3 dimGrid2(NBlock);

	CudaError("calcPSF::setup Dims");
	//mexPrintf("%d %d %d \n", N_int, PSFSizeOut, NPSFs);

	//kernel calls
	int Samplepsf_Pixels = SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels;
	for (int s = 0; s < 4; s++)
	{
		 
		kernel_calcPSFPixel_wrapper(dimGrid, dimBlock, &d_SamplePSF[s*Samplepsf_Pixels], &d_SamplePSFx[s*Samplepsf_Pixels], &d_SamplePSFy[s*Samplepsf_Pixels], &d_SamplePSFz[s*Samplepsf_Pixels],
			&d_SamplePSFxy[s*Samplepsf_Pixels], &d_SamplePSFxz[s*Samplepsf_Pixels], &d_SamplePSFyz[s*Samplepsf_Pixels], &d_SamplePSFxyz[s*Samplepsf_Pixels], 
			SampledPSF_XYPixels, SampledPSF_XYPixels, d_PSFs, d_dPSFx, d_dPSFy, d_dPSFz, d_X, d_Y, d_Z, &d_transx[s*Nfit], &d_transy[s*Nfit],
			SampleSpacingXY, SampleSpacingZ, StartX, StartY, StartZ, N_int, PSFSizeOut, Nfit);

		kernel_calcPSFI_wrapper(dimGrid2, dimBlock2, d_PSFs, d_PSFI, &d_I[s*Nfit], &d_Bg[s*Nfit], Nfit, PSFSize);

		CudaError("calcPSF::kernel");
		cudaMemcpy(&PSFs[s*Nfit*PSFSizeOut*PSFSizeOut], d_PSFI, Nfit*PSFSizeOut*PSFSizeOut*sizeof(float), cudaMemcpyDeviceToHost);
		CudaError("calcPSF::memcpy");
	}
}

void Model_Localization::calCRLB(){
	//PSF model generation

	int N_int = 4;		//There are N_int*N_int number of interpoloted points. 
	int BlockSize = 256;
	int NBlock = Nfit / BlockSize + 1;
	dim3 dimBlock(N_int, PSFSizeOut);
	dim3 dimGrid(PSFSizeOut, Nfit);
	dim3 dimBlock2(BlockSize);
	dim3 dimGrid2(NBlock);

	CudaError("calCRLB::setup Dims");
	//mexPrintf("%d %d %d \n", N_int, PSFSizeOut, NPSFs);

	//kernel calls
	int Samplepsf_Pixels = SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels;
	for (int s = 0; s < 4; s++)
	{
		kernel_calcPSFPixel_wrapper(dimGrid, dimBlock, &d_SamplePSF[s*Samplepsf_Pixels], &d_SamplePSFx[s*Samplepsf_Pixels], &d_SamplePSFy[s*Samplepsf_Pixels], &d_SamplePSFz[s*Samplepsf_Pixels],
			&d_SamplePSFxy[s*Samplepsf_Pixels], &d_SamplePSFxz[s*Samplepsf_Pixels], &d_SamplePSFyz[s*Samplepsf_Pixels], &d_SamplePSFxyz[s*Samplepsf_Pixels], 
			SampledPSF_XYPixels, SampledPSF_XYPixels, d_PSFs, d_dPSFx, d_dPSFy, d_dPSFz, d_X, d_Y, d_Z, &d_transx[s*Nfit], &d_transy[s*Nfit],
			SampleSpacingXY, SampleSpacingZ, StartX, StartY, StartZ, N_int, PSFSizeOut, Nfit);

		CudaError("calCRLB::kernel_calcPSFPixel_wrapper");

		kernel_calFisherM_wrapper(dimGrid2, dimBlock2, d_PSFs, d_dPSFx, d_dPSFy, d_dPSFz, &d_I[s*Nfit], &d_Bg[s*Nfit], &d_GainR[s*Nfit*PSFSize], &d_ParamF[s*NPL*NPL*Nfit],
			s, Nfit, PSFSize);

		CudaError("calCRLB::kernel_calFisherM_wrapper");
	}

	//kernel calls
	kernel_calCRLB_wrapper(dimGrid2, dimBlock2, d_ParamF, d_ParamVar, Nfit);

	cudaMemcpy(ParamVar, d_ParamVar, Nfit*NPL*sizeof(float), cudaMemcpyDeviceToHost);
	CudaError("calCRLB::kernel_calCRLB_wrapper");

}

void Model_Localization::calErr(){
	//PSF model generation
	int N_int = 4;		//There are N_int*N_int number of interpoloted points. 
	dim3 dimBlock(N_int, PSFSizeOut);
	dim3 dimGrid(PSFSizeOut, Nfit);
	int BlockSize = BSZ;
	int NBlock = Nfit / BlockSize + 1;
	dim3 dimBlock2(BlockSize);
	dim3 dimGrid2(NBlock);


	CudaError("calcPSF::setup Dims");

	//kernel calls
	int Samplepsf_Pixels = SampledPSF_XYPixels*SampledPSF_XYPixels*SampledPSF_ZPixels;
	for (int s = 0; s < 4; s++)
	{
		kernel_calcPSFPixel_wrapper(dimGrid, dimBlock, &d_SamplePSF[s*Samplepsf_Pixels], &d_SamplePSFx[s*Samplepsf_Pixels], &d_SamplePSFy[s*Samplepsf_Pixels], &d_SamplePSFz[s*Samplepsf_Pixels],
			&d_SamplePSFxy[s*Samplepsf_Pixels], &d_SamplePSFxz[s*Samplepsf_Pixels], &d_SamplePSFyz[s*Samplepsf_Pixels], &d_SamplePSFxyz[s*Samplepsf_Pixels], 
			SampledPSF_XYPixels, SampledPSF_XYPixels, d_PSFs, d_dPSFx, d_dPSFy, d_dPSFz, d_X, d_Y, d_Z, &d_transx[s*Nfit], &d_transy[s*Nfit],
			SampleSpacingXY, SampleSpacingZ, StartX, StartY, StartZ, N_int, PSFSizeOut, Nfit);

		kernel_calcPSFI_wrapper(dimGrid2, dimBlock2, d_PSFs, d_PSFI, &d_I[s*Nfit], &d_Bg[s*Nfit], Nfit, PSFSize);

		kernel_Err_wrapper(dimGrid2, dimBlock2, d_PSFI, &d_Data[s*Nfit*PSFSize], &d_GainR[s*Nfit*PSFSize], d_Error, Nfit, PSFSize);
	}

	CudaError("calcPSF::kernel");

	//Max threads per block is 1024 = 32^2.  Since we might want bigger than 32*32, lets have each block do one row of pixels. 
	CudaError("kernel_Err_wrapper");
	cudaMemcpy(Error, d_Error, Nfit*2*sizeof(float), cudaMemcpyDeviceToHost);

}

void Model_Localization::CudaError(const char * InStr) {
	/*!
	*  \brief A simple function to dump Cuda errors and abort the run.
	*/
	cudaError_t errornum;
	const char *str = 0;
	if (errornum = cudaGetLastError()) {
		str = cudaGetErrorString(errornum);
		mexErrMsgIdAndTxt("CudaTemplate:CUDA", "%s: %s\nYou should clear this function in MATLAB for proper operation.\n", InStr, str);
	}
}