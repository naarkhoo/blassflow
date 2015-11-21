#include "gig_link.h"
#include "gig_par.h"


#define Calloc(n,type) (type *)calloc((size_t)(n),sizeof(type))
void sample_gig(double* p,double* a,double* b,double* x,int n,unsigned long *seed)
{
	#ifdef _OPENMP
		int nP = omp_get_num_procs();
	#else
		int nP = 1;
	#endif

	


	RngStream * RngArray = new RngStream[nP];
	//for(int i = 0; i < nP ;i++)
	//		RngArray[i] = new RngStream();
	//gig * gigstream= Calloc(nP, gig); 
	gig *gigstream= new gig[nP];
	//U01s = new RngStream[nP];
	

	RngStream::SetPackageSeed(seed);
	for(int i = 0 ; i < nP; i++){
		gigstream[i].setRngStream(&(RngArray[i]));
	}
	#ifdef _OPENMP
		omp_set_num_threads(nP);
	#endif
	int thread = 0;
	int i;
	#pragma omp parallel private(i,thread) shared(x,p,a,b,n,gigstream)
	{
		#pragma omp for 	
		for(i = 0; i < n; i++){
			#ifdef _OPENMP
				thread = omp_get_thread_num();
			#endif
			// if using mexPrintf crashes things!!!
			//mexPrintf("thread, i = (%d,%d)\n",thread, i);
			gigstream[thread].rgig(x[i],
						&(p[i]), &(a[i]), &(b[i]));
		}
	}

	//free(gigstream);
	delete [] gigstream;	
	delete [] RngArray;
} 
