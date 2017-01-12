#include "mex.h"
#include "ReconInfo.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	if( nrhs < 2 ) {
		mexErrMsgTxt("result = OptimizeSolidMEX(InitialGuess, ReconInfoStruct)");
	} else if( !( mxIsDouble(prhs[0])&&mxIsStruct(prhs[1])) ) {
		mexErrMsgTxt("InitialGeuss must be of type double. ReconInfoStruct must be a structure.");
	} 

    if( nlhs < 1 )
        mexErrMsgTxt("Too few output arguments specified.");

    // Get the starting guess for the reconstruction
	double *S_start = mxGetPr(prhs[0]);

    // Get the number of dimensions
    mwSize NUM_DIMS = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *DIMS = mxGetDimensions(prhs[0]);

	int XSIZE = DIMS[0];
	int YSIZE = DIMS[1];
    int ZSIZE = (NUM_DIMS < 3) ? 1 : DIMS[2];
    int NUM_COMPONENTS = (NUM_DIMS < 4) ? 1 : DIMS[3];
    int TOTAL_SIZE = XSIZE*YSIZE*ZSIZE*NUM_COMPONENTS*sizeof(double);

    // Extract all the fields from our structure.
    ReconInfo * Recon = new ReconInfo(prhs[1]);

    // Check to make sure the starting guess the user passed in is the
    // right size.
    if(XSIZE != Recon->getSIZE_X() || 
       YSIZE != Recon->getSIZE_Y() ||
       ZSIZE != Recon->getSIZE_Z())
    {
        delete Recon;
        mexErrMsgTxt("Starting guess is wrong size!");
    }

    // Now we need to create a copy of the structure so we can optimize it.
    plhs[0] = mxCreateNumericArray(NUM_DIMS, DIMS, mxDOUBLE_CLASS, mxREAL);
    double *S = mxGetPr(plhs[0]);
    memcpy((void *)S, (void *)S_start, TOTAL_SIZE);

    // Wrap a 4D image class around this matrix
    Image4D img(S, XSIZE, YSIZE, ZSIZE, 1);

    mwSize cDIMS[2];
    cDIMS[0] = Recon->getNumExemplars();
    cDIMS[1] = 1;
    plhs[1] = mxCreateCellArray(2, cDIMS); 
    for(int i=0;i<cDIMS[0];i++)
        mxSetCell(plhs[1], i, Recon->getExemplarWeights()[i]); 

    // Run the optimization
    Recon->OptimizeIt(img);
    
    plhs[2] = Recon->getSourceTableArray();

    delete Recon;
    
}
