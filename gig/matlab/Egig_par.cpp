#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <stdlib.h>
#include <iostream>
#include "gig_par.h"

// Calculate the Excptation of A GIG r.v EG , EG^{-1} , ELOGG
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	#define EG_OUT          plhs[0]
	#define EIG_OUT         plhs[1]
	#define ELOGG_OUT       plhs[2]
	#define p 		prhs[0]
	#define a    		prhs[1]
	#define b    		prhs[2]

	if(nrhs<3)
		mexErrMsgTxt("wrong number of input argumets.");
	if(nlhs != 3)
		mexErrMsgTxt("wrong number of output argumets.");
	double *pVEC, *aVEC, *bVEC, *EG, *EIG, *ELOGG;
	int m;
	pVEC = mxGetPr(p);
	aVEC    = mxGetPr(a);
	bVEC    = mxGetPr(b);
	m = mxGetM(p);
	EG_OUT = mxCreateDoubleMatrix(m, 1, mxREAL);
	EG = mxGetPr(EG_OUT);
	EIG_OUT = mxCreateDoubleMatrix(m, 1, mxREAL);
	EIG = mxGetPr(EIG_OUT);
	ELOGG_OUT = mxCreateDoubleMatrix(m, 1, mxREAL);
	ELOGG = mxGetPr(ELOGG_OUT);
	Egig(m, pVEC, aVEC, bVEC, EG ,EIG , ELOGG );
	return;
}
