#include "mex.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <ctype.h>
#include <algorithm>
#include <set>

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{

	if( nrhs < 2 ) {
		mexErrMsgTxt("NBHoods = GetAllNHoodsMEX(Image, NB_SIZE)");
	} else if( !( mxIsDouble(prhs[0])&&mxIsDouble(prhs[1])) ) {
		mexErrMsgTxt("Image must be of type double and NB_SIZE must be double");
	} 

    if( nlhs < 1 )
        mexErrMsgTxt("Too few output arguments specified.");

    // Get the starting guess for the reconstruction
	double *I = mxGetPr(prhs[0]);
	int XSIZE = mxGetM(prhs[0]);
	int YSIZE = mxGetN(prhs[0]);
    int NB_SIZE = (int)mxGetScalar(prhs[1]);

    // Now we need to create a matrix big enough to store all the neighborhoods
    int NUM_NBS = (XSIZE-NB_SIZE+1)*(YSIZE-NB_SIZE+1);
    plhs[0] = mxCreateDoubleMatrix(NUM_NBS, NB_SIZE*NB_SIZE, mxREAL);
    double *Out = mxGetPr(plhs[0]);

    // If the user wants a lookup table saying where each neighborhood came from
    // then make that too.
    double *Lookup = NULL;
    if(nlhs > 1)
    {
        plhs[1] = mxCreateDoubleMatrix(NUM_NBS, 2, mxREAL);
        Lookup = mxGetPr(plhs[1]);
    }

    bool returnUniqueList = false;
    double *uniqOut = NULL;
    if(nlhs > 2)
    {
        returnUniqueList = true;
        plhs[2] = mxCreateDoubleMatrix(NUM_NBS, NB_SIZE*NB_SIZE, mxREAL);
        uniqOut = mxGetPr(plhs[2]);
    }

    int nbIndex = 0;

    for(int jj=0; jj<YSIZE-NB_SIZE+1; jj++)
    {
        for(int ii=0; ii<XSIZE-NB_SIZE+1; ii++)
        {
            int nbPixIndex = 0;
            set<double> uniqList;
            for (int ll=jj; ll<jj+NB_SIZE; ll++)
            {
                for (int kk=ii; kk<ii+NB_SIZE; kk++)
                {
                    Out[nbPixIndex*NUM_NBS+nbIndex] = I[ll*XSIZE+kk]; 

                    if(returnUniqueList)
                        uniqList.insert(I[ll*XSIZE+kk]);

                    nbPixIndex++;
                }
            }

            // If the user wants this, add it too the lookup table
            if(nlhs > 1)
            {
                Lookup[0*NUM_NBS+nbIndex] = ii+1;
                Lookup[1*NUM_NBS+nbIndex] = jj+1;
            }

            if(returnUniqueList)
            {
                int zz = 0;
                for(set<double>::iterator it=uniqList.begin(); it!=uniqList.end(); ++it, zz++)
                    uniqOut[nbIndex + zz*NUM_NBS] = *it;
            }
            
            nbIndex++;
        }
    }
    
}
