

// includes, system

#include <stdio.h>
// includes, project

//basic includes, others may be needed depending on application
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"
#include "cuda_runtime.h"

// Thread block size
#define BSZ 128
#define MEM 70
#define IMSZ 11
#define IMSZBIG 21
#define imMEM 4000
#define NK 256 //number of blocks to run in each kernel
#define pi 3.141592
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#define max(a,b)            (((a) > (b)) ? (a) : (b))

//kernel_MLEFit<<<dimGrid, dimBlock>>>(ii, sz, BlockSize, fitnum, d_xarray, d_yarray, d_Narray, d_barray, d_fishermatrix, BlockSize);


__global__ void kernel_guassiansampleblobs(int,int,int, float*,float*,float*, float*,float*,float*,float*,float*,float*);
__global__ void kernel_guassianintegrateblobs(int,int,int, float*,float*,float*, float*,float*,float*,float*,float*,float*);



//__device__ float PSF_xy(float x, int ii, float PSFSigma) {
//    norm=1.0/2.0/PSFSigma/PSFSigma;
//    return 1.0/2.0*(erf((ii-x+0.5)*sqrt(norm))-erf((ii-x-0.5)*sqrt(norm)));
//};
//
//__device__ float MODEL(float *x, float *y, float * Narray, int ii, int jj, float PSFSigma, float b, int N) {
//    float model;
//    model=b;
//    for (nn=0;nn<N;nn++)
//        model+=Narray[nn]*PSF_xy(xarray[nn], ii, PSFSigma)*PSF_xy(yarray[nn], jj, PSFSigma);
//    return model;
//
//};

void CUDAERRROR(const char *instr) {
	cudaError_t errornum;
	const char *str;
	if (errornum = cudaGetLastError()) {
		str = cudaGetErrorString(errornum);
		cudaThreadExit(); //release context so future cudaSetDevice calls work
		mexErrMsgIdAndTxt("CudaTemplate:CUDA", "%s: %s\nYou should clear this function in MATLAB for proper operation.\n", instr, str);
	}
}

void mexFunction(int nlhs, mxArray *plhs[],	int	nrhs, const	mxArray	*prhs[]) {
	int blockx;
	int threadx;
	int ii,iii,jj,kk,flag;
	int memblobsnum,ysz,xsz;
	float * xarray, * yarray, * Narray, *bg,*yt,*xl,*xsigma,*ysigma,*covariance,*im;
	float *d_xarray, *d_yarray, *d_Narray, *d_xsigma, *d_ysigma,*d_covariance,*d_im,*d_xl,*d_yt,*subim;
	const mwSize *datasize;
	int locr;
	mwSize imdim[2];



	if (nrhs<9)
		mexErrMsgTxt("xsize,ysize,x_array, y_array, N_array, sigmaX, sigmaY, covariance, UseIntegrated_FLAG\n");

	if (mxGetClassID(prhs[0])!=mxSINGLE_CLASS)
		mexErrMsgTxt("Data must be comprised of single floats!\n");
	if (mxGetClassID(prhs[1])!=mxSINGLE_CLASS)
		mexErrMsgTxt("Data must be comprised of single floats!\n");
	if (mxGetClassID(prhs[2])!=mxSINGLE_CLASS)
		mexErrMsgTxt("Data must be comprised of single floats!\n");
	if (mxGetClassID(prhs[3])!=mxSINGLE_CLASS)
		mexErrMsgTxt("Data must be comprised of single floats!\n");
	if (mxGetClassID(prhs[4])!=mxSINGLE_CLASS)
		mexErrMsgTxt("Data must be comprised of single floats!\n");
	if (mxGetClassID(prhs[5])!=mxSINGLE_CLASS)
		mexErrMsgTxt("Data must be comprised of single floats!\n");


	datasize=mxGetDimensions(prhs[2]);
	if (datasize[1]!=1)
		mexErrMsgTxt("xarray should be n X 1 array\n");

	datasize=mxGetDimensions(prhs[3]);

	if (datasize[1]!=1)
		mexErrMsgTxt("xarray should be n X 1 array\n");

	datasize=mxGetDimensions(prhs[4]);

	if (datasize[1]!=1)
		mexErrMsgTxt("xarray should be n X 1 array\n");

	datasize=mxGetDimensions(prhs[5]);
	if (datasize[1]!=1)
		mexErrMsgTxt("xarray should be n X 1 array\n");


	xsz =(float) mxGetScalar(prhs[0]);
	ysz =(float) mxGetScalar(prhs[1]);
	imdim[0]=xsz;
	imdim[1]=ysz;
	//PSFSigma=(float)mxGetScalar(prhs[1]); //matlab-dip_image convention
	xarray =(float *) mxGetData(prhs[2]);
	yarray =(float *) mxGetData(prhs[3]);
	Narray =(float *) mxGetData(prhs[4]);
	xsigma =(float *)mxGetData(prhs[5]);
	ysigma =(float *)mxGetData(prhs[6]);
	covariance =(float *)mxGetData(prhs[7]);
	flag =(float) mxGetScalar(prhs[8]);


	int blobn=datasize[0];
	float maxsigma=-1;
	float sigma;
	for(ii=0;ii<blobn;ii++){
		sigma=sqrt(pow(xsigma[ii],2)+pow(ysigma[ii],2));
		maxsigma=max(maxsigma,sigma);
	}

	int sz=(int) round(float(8*maxsigma));
	sz=min(sz,20);


	if ((flag!=1)&&(flag!=0))
		mexErrMsgTxt("flag can only be 0 or 1\n");

	// over allocate for additional thread reading error
	int BlockSize=min(ceil((float) 15000/4/sz/sz),64);
	memblobsnum=(int)ceil((float)datasize[0]/BlockSize)+128;

	//mexPrintf("Starting CUDA Malloc\n");
	
	CUDAERRROR("P1");

	cudaMalloc(&d_xarray, memblobsnum*BlockSize*sizeof(float));
	CUDAERRROR("M1");

	cudaMemset(d_xarray, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_xarray, xarray, datasize[0]*sizeof(float), cudaMemcpyHostToDevice);
	CUDAERRROR("S1");

	cudaMalloc((void**)&d_yarray, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_yarray, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_yarray, yarray,datasize[0]*sizeof(float), cudaMemcpyHostToDevice);
	CUDAERRROR("M2");

	cudaMalloc((void**)&d_Narray, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_Narray, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_Narray, Narray,datasize[0]*sizeof(float), cudaMemcpyHostToDevice);
	CUDAERRROR("M3");

	cudaMalloc((void**)&d_xsigma, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_xsigma, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_xsigma, xsigma,datasize[0]*sizeof(float), cudaMemcpyHostToDevice);
	CUDAERRROR("M4");

	cudaMalloc((void**)&d_ysigma, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_ysigma, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_ysigma, ysigma,datasize[0]*sizeof(float), cudaMemcpyHostToDevice);
	CUDAERRROR("M5");

	cudaMalloc((void**)&d_covariance, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_covariance, 0, memblobsnum*BlockSize*sizeof(float));
	cudaMemcpy(d_covariance, covariance,datasize[0]*sizeof(float), cudaMemcpyHostToDevice);
	CUDAERRROR("M6");


	cudaMalloc((void**)&d_im, sz*sz*memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_im, 0, sz*sz*memblobsnum*BlockSize*sizeof(float));

	cudaMalloc((void**)&d_xl, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_xl, 0, memblobsnum*BlockSize*sizeof(float));

	cudaMalloc((void**)&d_yt, memblobsnum*BlockSize*sizeof(float));
	cudaMemset(d_yt, 0, memblobsnum*BlockSize*sizeof(float));





	//only run NK blocks in each kernel
	int numK=(int)ceil((float)datasize[0]/BlockSize/NK);

	for (int ii=0;ii<numK;ii++) {

		blockx = min(ceil(((float)(((float)datasize[0])/BlockSize)-ii*NK)), NK);
		blockx = max(blockx,1);
		threadx= BlockSize;


		dim3 dimBlock(threadx);
		dim3 dimGrid(blockx);

		//printf("threadx: %d,blockx: %d\n", threadx, blockx);

		switch (flag)
		{
		case 0:
			kernel_guassiansampleblobs<<<dimGrid, dimBlock>>>(ii,BlockSize,sz, d_xarray,d_yarray,d_Narray, d_xsigma,d_ysigma,d_covariance,d_im,d_xl,d_yt);
			break;//15x15 images, 64 per block
		case 1:
			kernel_guassianintegrateblobs<<<dimGrid, dimBlock>>>(ii,BlockSize,sz, d_xarray,d_yarray,d_Narray, d_xsigma,d_ysigma,d_covariance,d_im,d_xl,d_yt);
			break;//15x15 images, 64 per block
		}

		CUDAERRROR("kernel");
		//mexEvalString("pause(0.001)");

	}

	subim= (float * )malloc(datasize[0]*sz*sz*sizeof(float));
	xl=(float * )malloc(datasize[0]*sizeof(float));
	yt=(float * )malloc(datasize[0]*sizeof(float));


	//reconstruct images
	plhs[0]=mxCreateNumericArray(2, imdim, mxSINGLE_CLASS, mxREAL);
	im=(float *)mxGetData(plhs[0]);

	cudaMemcpy(subim, d_im, datasize[0]*sz*sz*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(xl, d_xl, datasize[0]*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(yt, d_yt, datasize[0]*sizeof(float), cudaMemcpyDeviceToHost);


	for(kk=0;kk<blobn;kk++){
		for(jj=0;jj<sz;jj++){
			for(iii=0;iii<sz;iii++){
				if ((((int)xl[kk]+iii)<(xsz-1))&&(((int)yt[kk]+jj)<(ysz-1))){
					locr=((int)yt[kk]+jj)*xsz+(int)xl[kk]+iii;
					if((subim[kk*sz*sz+jj*sz+iii]>0)&&(subim[kk*sz*sz+jj*sz+iii]<100000)&&(locr>=0)&&(locr<=((xsz-1)*(ysz))))
					im[locr]+=subim[kk*sz*sz+jj*sz+iii];	
				}
			}
		}
	}



	free(subim);
	free(xl);
	free(yt);
	cudaFree(d_xarray);
	cudaFree(d_yarray);
	cudaFree(d_Narray);
	cudaFree(d_xsigma);
	cudaFree(d_ysigma);
	cudaFree(d_covariance);
	cudaFree(d_im);
	cudaFree(d_xl);
	cudaFree(d_yt);
	cudaDeviceReset();

}


//kernel_guassiansampleblobs<<<dimGrid, dimBlock>>>(ii,blockx,BlockSize,sz, d_xarray,d_yarray,d_Narray, d_xsigma,d_ysigma,d_covariance,d_im,d_xl,d_yt);   //15x15 images, 64 per block

__global__ void kernel_guassiansampleblobs(int iiK,int BlockSize, int sz, float *d_xarray,float *d_yarray,float *d_Narray, float *d_xsigma,float *d_ysigma,float *d_covariance,float *d_im,float *d_xl,float *d_yt  ) {
	int tx = threadIdx.x; //matrix number index
	int bx = blockIdx.x;
	float x,y,xsigma,ysigma,covariance,N;
	float xl;
	float yt;
	int ii,jj,pixelx,pixely;


	float model;//

	__shared__ float s_im[imMEM];


	bx=bx+iiK*NK;
	//import datas from device to shared memory

	x=d_xarray[bx*BlockSize+tx];
	y=d_yarray[bx*BlockSize+tx];
	N=d_Narray[bx*BlockSize+tx];
	xsigma=d_xsigma[bx*BlockSize+tx];
	ysigma=d_ysigma[bx*BlockSize+tx];
	covariance=d_covariance[bx*BlockSize+tx];
	xl=round(x)-round(float (sz/2-1));
	xl=max(xl,0);

	yt=round(y)-round(float (sz/2-1));
	yt=max(yt,0);


	for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {

		// generate model for pixel ii jj
		pixelx=ii;
		pixely=jj;
		s_im[tx*sz*sz+jj*sz+ii]=N/(2*pi*xsigma*ysigma*sqrt(1-pow(covariance,2)))*exp(-1/(2*(1-pow(covariance,2)))*(pow(x-xl-pixelx,2)/pow(xsigma,2)+pow(y-yt-pixely,2)/pow(ysigma,2)-2*covariance*(x-xl-pixelx)*(y-yt-pixely)/(xsigma*ysigma)));
	}



	for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++)
	{
		d_im[bx*BlockSize*sz*sz+tx*sz*sz+jj*sz+ii]=s_im[tx*sz*sz+jj*sz+ii];
		d_xl[bx*BlockSize+tx]=xl;
		d_yt[bx*BlockSize+tx]=yt;
	}

	return;



}



__global__ void kernel_guassianintegrateblobs(int iiK,int BlockSize, int sz, float *d_xarray,float *d_yarray,float *d_Narray, float *d_xsigma,float *d_ysigma,float *d_covariance,float *d_im,float *d_xl,float *d_yt  ) {
	int tx = threadIdx.x; //matrix number index
	int bx = blockIdx.x;
	float x,y,xsigma,ysigma,covariance,N;
	float xl;
	float yt;
	int ii,jj,pixelx,pixely;


	float model;//

	__shared__ float s_im[imMEM];


	bx=bx+iiK*NK;
	//import datas from device to shared memory

	x=d_xarray[bx*BlockSize+tx];
	y=d_yarray[bx*BlockSize+tx];
	N=d_Narray[bx*BlockSize+tx];
	xsigma=d_xsigma[bx*BlockSize+tx];
	ysigma=d_ysigma[bx*BlockSize+tx];
	covariance=d_covariance[bx*BlockSize+tx];

	xl=round(x)-round(float (sz/2-1));
	xl=max(xl,0);

	yt=round(y)-round(float (sz/2-1));
	yt=max(yt,0);

	for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {

		// generate model for pixel ii jj
		pixelx=ii;
		pixely=jj;
		s_im[tx*sz*sz+jj*sz+ii]=N/4*(erf((x-xl-pixelx-0.5)/sqrt(2*pow(xsigma,2)))-erf((x-xl-pixelx+0.5)/sqrt(2*pow(xsigma,2))))*(erf((y-yt-pixely-0.5)/sqrt(2*pow(ysigma,2)))-erf((y-yt-pixely+0.5)/sqrt(2*pow(ysigma,2))));  //exp(-1/(2*(1-pow(covariance,2)))*(pow(x-xl-pixelx,2)/pow(xsigma,2)+pow(y-yt-pixely,2)/pow(ysigma,2)-2*covariance*(x-xl-pixelx)*(y-yt-pixely)/(xsigma*ysigma)));
	}



	for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++)
	{
		d_im[bx*BlockSize*sz*sz+tx*sz*sz+jj*sz+ii]=s_im[tx*sz*sz+jj*sz+ii];
		d_xl[bx*BlockSize+tx]=xl;
		d_yt[bx*BlockSize+tx]=yt;
	}

	return;



}

//END OF KERNAL FUNCTION


