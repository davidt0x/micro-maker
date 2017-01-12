#include "mex.h"
#include <math.h>
#include <omp.h>

#include "ReconInfo.hpp"

// This function does all the work. It spawns all of the threads that do all the searching.
void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // Don't ever allow this mex file to be cleared. This prevents a hang
    // because this call utilizes the MIC processor.
    mexLock();

	if( nrhs < 2 ) {
		mexErrMsgTxt("result = SearchNNBMEX(SolidImage, ReconInfoStruct)");
	} else if( !( mxIsDouble(prhs[0])&&mxIsStruct(prhs[1])) ) {
		mexErrMsgTxt("SolidImage must be of type double. ReconInfoStruct must be a structure.");
	} 

    if( nlhs < 1 )
        mexErrMsgTxt("Too few output arguments specified.");

    // Get the reconstruction image
	double *S_start = mxGetPr(prhs[0]);

    // Get the number of dimensions
    mwSize NUM_DIMS = mxGetNumberOfDimensions(prhs[0]);
    const int *DIMS = mxGetDimensions(prhs[0]);

    Image4D S(S_start, DIMS[0], DIMS[1], DIMS[2], 1);

    if(NUM_DIMS < 2)
        mexErrMsgTxt("SolidImage must at least be 2D");

	int XSIZE = DIMS[0]; // Number of rows
	int YSIZE = DIMS[1]; // Number of columns
    int ZSIZE = (NUM_DIMS < 3) ? 1 : DIMS[2];

    int TOTAL_SIZE = XSIZE*YSIZE*ZSIZE*sizeof(double);

    // Extract all the fields from our structure.
    ReconInfo * Recon = new ReconInfo(prhs[1]);

    // Run the search, this function has side effects. It will
    // update the nearest neighborhoods table.
    unsigned int * nnb_table = Recon->SearchIt(S);

    // Copy the result to an output array.
    int NNB_DIMS[] = {XSIZE, YSIZE, ZSIZE, Recon->getNumExemplars()};
    plhs[0] = mxCreateNumericArray(4, NNB_DIMS, mxUINT32_CLASS, mxREAL);
    memcpy(mxGetData(plhs[0]), (void *)nnb_table, 
        XSIZE*YSIZE*ZSIZE*Recon->getNumExemplars()*sizeof(unsigned int));

    delete Recon;
    delete nnb_table;
    
}
