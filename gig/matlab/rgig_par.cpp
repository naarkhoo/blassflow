#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <stdlib.h>
#include "gig_link.h"
#include <iostream>

#define Calloc(n,type) (type *)calloc((size_t)(n),sizeof(type))

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	#define x_OUT           plhs[0]
	#define p prhs[0]
	#define a    prhs[1]
	#define b    prhs[2]
	#define seedMatlab	prhs[3]

	if(nrhs<4)
		mexErrMsgTxt("wrong number of input argumets.");
	if(nlhs != 1)
		mexErrMsgTxt("wrong number of output argumets.");



	double *pVEC, *aVEC, *bVEC, *x;
	int m;
	pVEC = mxGetPr(p);
	aVEC    = mxGetPr(a);
	bVEC    = mxGetPr(b);
	m = mxGetM(p);
	x_OUT = mxCreateDoubleMatrix(m, 1, mxREAL);
	x = mxGetPr(x_OUT);
	double *seed_double;	
	seed_double =mxGetPr(seedMatlab);
	unsigned long *seed = Calloc(6, unsigned long);
	for(int i = 0;i < 6 ; i++){
		seed[i] = (unsigned long) seed_double[i];
		}
	sample_gig(pVEC , aVEC, bVEC, x, m, seed);
	free(seed);
	return;
}
