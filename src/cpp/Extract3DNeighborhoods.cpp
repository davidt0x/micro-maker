#include "mex.h"

#include <algorithm>
#include <set>

using namespace std;

inline int mod(int x, int m) {
//int __attribute__((target(mic))) mod(int x, int m) {
    return (x%m + m)%m;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	if( nrhs < 2 ) {
		mexErrMsgTxt("[nbX nbY nbZ nbXY nbXZ nbYZ] = Extract3DNeighborhoods(S, nbOffsets, NUM_CORES)");
	} else if( !( mxIsDouble(prhs[0]) && mxIsDouble(prhs[1])) )   {
		mexErrMsgTxt("S must be of type double and so must nbOffsets.");
	}

    // Get the 3D volume
    double *S = mxGetPr(prhs[0]);
 
    // Get the offsets for the neighborhood
    double *nbOffsets = mxGetPr(prhs[1]);

    if(mxGetNumberOfDimensions(prhs[1]) != 2)
        mexErrMsgTxt("nbOffsets must be a 2D array of voxel offsets for the neighborhood.");

    if(mxGetDimensions(prhs[1])[1] != 3)
        mexErrMsgTxt("nbOffsets must have NUM_NB_VOXELSx3 dimensions");

    // Get the neighborhood dimensions;
    int NUM_NB_VOXELS = mxGetDimensions(prhs[1])[0];
    int SIZE_OF_NB = NUM_NB_VOXELS*3;

    // Get the number of dimensions of the sample
    mwSize NUM_DIMS = mxGetNumberOfDimensions(prhs[0]);
    const mwSize *DIMS = mxGetDimensions(prhs[0]);

    if(NUM_DIMS < 3)
        mexErrMsgTxt("S must be a 3D matrix!");
    
    int XSIZE = DIMS[0];
    int YSIZE = DIMS[1];
    int ZSIZE = DIMS[2];

    // Each pixel has the same number of neighbrohoods. Lets figure out how many
    // pixels are in our image.
    int NUM_PIXELS = XSIZE*YSIZE*ZSIZE;
    int SIZE_OF_S = NUM_PIXELS;

    // Create the output arrays for each neighborhood
    mwSize NBS_DIMS[] = {NUM_PIXELS, NUM_NB_VOXELS};
    int SIZE_OF_OUT = NUM_PIXELS*NUM_NB_VOXELS;
    plhs[0] = mxCreateNumericArray(2, NBS_DIMS, mxDOUBLE_CLASS, mxREAL);
    double *output = mxGetPr(plhs[0]); 

	plhs[1] = mxCreateLogicalMatrix(NUM_PIXELS, 1);
    mxLogical *isPeriodic = mxGetLogicals(plhs[1]); 
	for(int i=0;i<NUM_PIXELS;i++)
		isPeriodic[i] = false;

    // The user has requested that we return list of unique values in each
    // neighborhood. This is useful in grain maps.
    bool returnUniqueList = false;
    double *uniqOutput = NULL;
    if(nlhs > 2)
    {
        returnUniqueList = true;
        plhs[2] = mxCreateNumericArray(2, NBS_DIMS, mxDOUBLE_CLASS, mxREAL);
        uniqOutput = mxGetPr(plhs[2]); 
    }


    // Get the number of cores to use if it is a parameter
    int NUM_CORES = 1;
    if(nrhs == 3)
        NUM_CORES = (int)mxGetScalar(prhs[2]);

    //#pragma offload target(mic:0) in(S[0:SIZE_OF_S]) in(nbOffsets[0:SIZE_OF_NB]) out(output[1:SIZE_OF_OUT])
    {
        #pragma omp parallel for num_threads(NUM_CORES)
        for(int i=0; i<XSIZE; i++)
        {
            for(int j=0; j<YSIZE; j++) 
            {
                for(int k=0; k<ZSIZE; k++)
                {
                    set<double> uniqList;

                    int pixi = k + j * (ZSIZE) + i * (YSIZE * ZSIZE);

                    for(int nb_vox=0; nb_vox < NUM_NB_VOXELS; nb_vox++)
                    {
                       
						int iu = i + nbOffsets[nb_vox + 0*NUM_NB_VOXELS];
						int ju = j + nbOffsets[nb_vox + 1*NUM_NB_VOXELS];
						int ku = k + nbOffsets[nb_vox + 2*NUM_NB_VOXELS];

						if(!isPeriodic[pixi])
							isPeriodic[pixi] = (iu < 0 || iu >= XSIZE || ju < 0 || ju >= YSIZE || ku < 0 || ku >= ZSIZE);
                        
						int nb_i = mod(iu, XSIZE);
                        int nb_j = mod(ju, YSIZE);
                        int nb_k = mod(ku, ZSIZE);
                        double val = S[nb_i + nb_j * (XSIZE) + nb_k *(XSIZE*YSIZE)];
                        output[pixi + nb_vox*NUM_PIXELS] = val;

                        if(returnUniqueList)
                            uniqList.insert(val);
                    }

                    if(returnUniqueList)
                    {
                        int zz = 0;
                        for(set<double>::iterator it=uniqList.begin(); it!=uniqList.end(); ++it, zz++)
                            uniqOutput[pixi + zz*NUM_PIXELS] = *it;
                    }
                 
                }
            }
        }
    }
}
