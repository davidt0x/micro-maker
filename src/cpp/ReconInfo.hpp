#ifndef _RECON_INFO_H
#define _RECON_INFO_H

#include "mex.h"
#include <math.h>
#include <omp.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <sstream>
#include <omp.h>

using namespace std;

#include "Histogram.hpp"
#include "Image4D.hpp"
#include "Exemplar.hpp"

int mod(int x, int m) {
    return (x%m + m)%m;
}

// Utility functions
// Extract a field from a matlab cell array. Make sure it exists too.
mxArray * getField(const mxArray * Recon, const char * field)
{
    mxArray * tmp = mxGetField(Recon, 0, field);
    if(tmp == NULL)
    {
		ostringstream sout;
		sout << field << " not found in ReconInfo" << std::endl;
        mexErrMsgTxt(sout.str().c_str());
    }
    return tmp;
}

int FindClosestValue(vector<double> vals, int num_candidates, double targetVal)
{
    double minDist = 1e200;
    int minI = 0;
    for(int zz=0; zz<num_candidates; zz++) {
        double dist = fabs(vals[zz] - targetVal); 
        if(dist < minDist) { 
            minI = zz;
            minDist = dist; 
        }
    }
             
    return minI;
}


inline int clamp_it(int x, double minx, double maxx)
{
    if(x < minx)
        x = minx;
    if(x > maxx)
        x = maxx; 

    return x;
}

// The reconstrution class. This class handles the reconstruction tasks.
class ReconInfo
{
    public:

        // Initialize the reconstruciton object with info sent from the matlab
        // environment. This allows the object to process matlab data and modify
        // it.
        ReconInfo(const mxArray * Recon)
        {
            mxArray *tmp;
            
            // Get the size of the neighboroods
            NB_SIZE = (int)mxGetScalar(getField(Recon, "NB_SIZE"));

            NUM_NB_VOXELS = NB_SIZE*NB_SIZE;

            // Get the number of cores to use
            NUM_CORES = (int)mxGetScalar(getField(Recon, "NUM_CORES"));

            // Get the size of the reconstruction
            double * reconsize = mxGetPr(getField(Recon, "RECON_SIZE"));
            SIZE_X = (int)reconsize[0]; 
            SIZE_Y = (int)reconsize[1];

            if(mxGetN(getField(Recon, "RECON_SIZE")) > 2)
                SIZE_Z = (int)reconsize[2];
            else
                SIZE_Z = 1;

            // Get the neighborhoods offsets array
            mxArray * nbOffsets_mxArray = getField(Recon, "nbOffsets");
            double *nbO = mxGetPr(nbOffsets_mxArray);
            const mwSize *NBO_DIMS= mxGetDimensions(nbOffsets_mxArray);
			int NUM_NB_VOXELS = NBO_DIMS[1];
            Image4D nbOff(nbO, NBO_DIMS[0], NBO_DIMS[1], NBO_DIMS[2], 1);

            mxArray *E_Edges = mxGetField(Recon, 0, "E_Edges");

            // If we have a edge map then this is grain map reconstruction
            isGrainMap =  E_Edges != NULL;

            // Load the angle map if it exists
            mxArray *mxAngMap = mxGetField(Recon, 0, "ANG_MAP");
            if(mxAngMap != NULL)
            {
                AngMap = mxGetPr(mxAngMap);
                ANGMAP_NUM_ORIENTS = mxGetM(mxAngMap);
                ANGMAP_N = mxGetN(mxAngMap);
                hasAngMap = true;
            }
            else
            {
                hasAngMap = false;
                ANGMAP_N = 1;
            }           

            // Get the exemplar images
            NumExemplars = mxGetNumberOfElements(getField(Recon, "EXEMPLARS"));
            ExemplarWeights.resize(NumExemplars);
            for(int i=0; i<NumExemplars; i++) {
                
                vector< vector<int> > nbOffsets;
                for(int voxi=0; voxi<NUM_NB_VOXELS; voxi++)
                {
                    vector<int> tmp(3);
                    tmp[0] = nbOff.at(i, voxi, 0, 0);
                    tmp[1] = nbOff.at(i, voxi, 1, 0);
                    tmp[2] = nbOff.at(i, voxi, 2, 0);
                    nbOffsets.push_back(tmp);
                }
                
                // Get the exemplar image and its dimensions
                mxArray *E = mxGetCell(getField(Recon, "EXEMPLARS"), i);
                double * data = mxGetPr(E);
                int szX = mxGetDimensions(E)[0];
                int szY = mxGetDimensions(E)[1];
                int numComponents = 1;
                if(mxGetNumberOfDimensions(E) > 2)
                    numComponents = mxGetDimensions(E)[2];

                // Copy the neighborhood lookup table to vector 
                mxArray *NBLook = 
                    mxGetCell(getField(Recon, "NB_ExemplarLookup"), i);
                double *look = mxGetPr(NBLook);
                vector< vector<int> > NBLookup;
                for(int p=0; p < mxGetM(NBLook); p++) {
                    vector<int> coords(2);
                    coords[0] = look[mxGetM(NBLook)*0 + p];
                    coords[1] = look[mxGetM(NBLook)*1 + p];
                    NBLookup.push_back(coords);
                }

                // Get the weight table for each texel.
                mwSize EX_DIMS[2];
                EX_DIMS[0] = szX;
                EX_DIMS[1] = szY;
                mxArray * txW = mxCreateNumericArray(2, EX_DIMS, mxDOUBLE_CLASS, mxREAL);
                ExemplarWeights[i] = txW;

                // If we have already calculated some weights for this exemplar
                // then copy them so we can use them again.
                mxArray * mxCellWeights = mxGetField(Recon, 0, "TexelWeights");
                if(mxCellWeights != NULL)
                {
                    mxArray * mxWeights = mxGetCell(mxCellWeights, i);
                    memcpy((void*)mxGetPr(txW),(void*)mxGetPr(mxWeights), 
                        sizeof(double)*szX*szY); 
                }
                else {
                    // If we have no weights, they should all start at 1.0
                    for(int w=0;w<szX*szY;w++)
                        mxGetPr(txW)[w] = 1.0;
                }

                // Construct and allocate this exemplar
                Exemplars.push_back(
                    new Exemplar(data, numComponents, NB_SIZE, NBLookup,
                        szX, szY, nbOffsets, mxGetPr(txW)));

            }   

            NUM_COMPONENTS = Exemplars[0]->getNUM_COMPONENTS();

            // Construct a 4D Image that is a copy of the nearest neighbor search
            // table passed from matlab.
            mwSize DIMS[4] = {SIZE_X, SIZE_Y, SIZE_Z, NumExemplars};
            mxArray *nnbCopy = mxCreateNumericArray(4, DIMS, mxDOUBLE_CLASS , mxREAL);
            NNB_Tables = Image4D(mxGetPr(nnbCopy), SIZE_X, SIZE_Y, SIZE_Z, NumExemplars);
            mxArray * tbl = getField(Recon, "NNB_Table");
            memcpy((void*)mxGetPr(nnbCopy),(void*)mxGetPr(tbl), sizeof(double)*SIZE_X*SIZE_Y*SIZE_Z*NumExemplars); 

            // For each exemplar, we want to keep a table the size of the
            // reconstruction that tells where the pixel came from in the 
            // exemplar.
            SourceTableArray = mxCreateNumericArray(4, DIMS, mxDOUBLE_CLASS, mxREAL);
            SourceTable = new Image4D(mxGetPr(SourceTableArray), SIZE_X, SIZE_Y, SIZE_Z, NumExemplars);
            
            // If the user passes a source table, copy it into the output table.
            mxArray * inSrcTbl = mxGetField(Recon, 0, "SourceTable");
            if(inSrcTbl != NULL)
                memcpy((void*)mxGetPr(SourceTableArray),(void*)mxGetPr(inSrcTbl), sizeof(double)*SIZE_X*SIZE_Y*SIZE_Z*NumExemplars); 

            NUM_PIXELS = SIZE_X*SIZE_Y*SIZE_Z;

            I = NULL;
            J = NULL;
            K = NULL;
        }

        // Destroy the reconstruction objects dynamically allocated memory
        ~ReconInfo()
        {
            // Delete all the exemplars
            for(int i=0;i<NumExemplars;i++)
                delete Exemplars[i];

            delete SourceTable;

        }

        void ExtractNeighborhoods(
                        int XSIZE, int YSIZE, int ZSIZE, 
                        int * nbOffsets, double * S, double * output)
        {
            #pragma omp parallel for num_threads(244)
            for(int i=0; i<XSIZE; i++)
                for(int j=0; j<YSIZE; j++) 
                    for(int k=0; k<ZSIZE; k++)
                        for(int nb_vox=0; nb_vox < NUM_NB_VOXELS; nb_vox++)
                        {
                            int nb_i = mod(i + nbOffsets[nb_vox + 0*NUM_NB_VOXELS], XSIZE);
                            int nb_j = mod(j + nbOffsets[nb_vox + 1*NUM_NB_VOXELS], YSIZE);
                            int nb_k = mod(k + nbOffsets[nb_vox + 2*NUM_NB_VOXELS], ZSIZE);
                            int pixi = k + j * (ZSIZE) + i * (YSIZE * ZSIZE);
                            output[pixi + nb_vox*(XSIZE*YSIZE*ZSIZE)] = S[nb_i + nb_j * (XSIZE) + nb_k *(XSIZE*YSIZE)];
                        }
        }
        
        void InitializeNNBSearch()
        {
             
        }

        // For every pixel in the structure, find the nearest neighborhoods
        // in the exemplar images.
        unsigned int * SearchIt(Image4D & S)
        {
            // The is the result of this function
            unsigned int * nnb_table = new unsigned int[S.getTotalSize()*NumExemplars];

            double * S_ptr = S.getData();
            int S_TOTAL_SIZE = S.getTotalSize();
            int HALF_NB_SIZE = (int)floor((double)NB_SIZE/2.0);

            // Convert the nbOffsets that are vectors in double * because STL
            // vectors are not working on Phi Architecture right now.
            int ** nbOffsets = new int*[NumExemplars];
            for(int exi=0;exi<NumExemplars;exi++)
            {
                nbOffsets[exi] = new int[NUM_NB_VOXELS*3];
                vector< vector<int> > nbO = Exemplars[exi]->getNBOffsets();
                for(int i=0;i<NB_SIZE*NB_SIZE;i++)
                {
                    nbOffsets[exi][i + 0*NUM_NB_VOXELS] = nbO[i][0];
                    nbOffsets[exi][i + 1*NUM_NB_VOXELS] = nbO[i][1];
                    nbOffsets[exi][i + 2*NUM_NB_VOXELS] = nbO[i][2];
                }
            }

            // Transfer S over to the MIC. This only needs to be done once.
            //#pragma offload_transfer target(mic:0) in(S_ptr[0:S_TOTAL_SIZE] : free_if(0)) 

            // For each exemplar, search every voxel in the solid to find 
            // the best matching neighborhood
            int XSIZE = S.getSizeX();
            int YSIZE = S.getSizeY();
            int ZSIZE = S.getSizeZ();
            {
                int NBS_SIZE = S_TOTAL_SIZE*NB_SIZE*NB_SIZE;
                for(int exi=0; exi < NumExemplars; exi++) {
                    int * nbO = nbOffsets[exi];
            //        #pragma offload target(mic:0) in(nbO[0:NUM_NB_VOXELS*3]) nocopy(S_ptr[0:S_TOTAL_SIZE]) 
                    {
                        vector< vector<double> > queries;
                        
                        double * nbs = new double[NBS_SIZE];
                        ExtractNeighborhoods(XSIZE, YSIZE, ZSIZE,
                                             nbO, S_ptr, nbs); 

                        delete nbs;
                    }

                } // End Exemplar For

            }
            
            for(int exi=0;exi<NumExemplars;exi++)
                delete nbOffsets[exi];
            delete nbOffsets;

            return nnb_table;
        }

        // Optimize the solid using the nearest neighbor table computed during
        // the search phase.
        void OptimizeIt(Image4D & S)
        {
            int HALF_NB_SIZE = (int)floor((double)NB_SIZE/2.0);

            // Generate a random ordering of pixels so we proceed through the
            // optimization in a random manner.
            GenerateRandomPixelSequence();

            // Get the histogram for each of the components
            vector<Histogram> SampleHist(NUM_COMPONENTS);
            vector<double> D(S.getSizeX()*S.getSizeY()*S.getSizeZ());
            for(int c=0;c<S.getNumComponents();c++)
            {
                int l = 0;
                for(int i=0;i<S.getSizeX();i++)
                    for(int j=0;j<S.getSizeY();j++)
                        for(int k=0;k<S.getSizeZ();k++,l++) 
                            D[l] = S.at(i,j,k,c);
                
                SampleHist[c] = Histogram(16);
                SampleHist[c].calcHistogram(D, S.getSizeX()*S.getSizeY()*S.getSizeZ() );
            }
            
            for(int exi=0; exi<NumExemplars; exi++)
            {
                for(int x=0;x<Exemplars[exi]->getSIZE_X(); x++)
                    for(int y=0;y<Exemplars[exi]->getSIZE_Y(); y++)
                    {
                        double val = Exemplars[exi]->at(x, y);
                        double weight = Exemplars[exi]->getTexelWeight(x, y) / (1 + 
                            std::max(0.0, SampleHist[0].getProbForValue(val) - Exemplars[exi]->getHistBinForValue(val, 0)));
                        Exemplars[exi]->setTexelWeight(x, y, weight); 
                    }
            }
            
            /*

            // For exemplar, update its texel weight table based on the current histograms
            vector< double * > posHists(NumExemplars);
            for(int exi=0; exi<NumExemplars; exi++)
            {
                posHists[exi] = new double[Exemplars[exi]->getSIZE_X()*Exemplars[exi]->getSIZE_Y()];
                for(int x=0;x<Exemplars[exi]->getSIZE_X(); x++)
                    for(int y=0;y<Exemplars[exi]->getSIZE_Y(); y++)
                        posHists[exi][x + y*Exemplars[exi]->getSIZE_X()] = 0;

                for(int i=0;i<S.getSizeX();i++)
                    for(int j=0;j<S.getSizeY();j++)
                        for(int k=0;k<S.getSizeZ();k++) 
                            posHists[exi][(int)SourceTable->at(i,j,k,exi)]++;
            }

            for(int exi=0; exi<NumExemplars; exi++)
            {
                for(int x=0;x<Exemplars[exi]->getSIZE_X(); x++)
                    for(int y=0;y<Exemplars[exi]->getSIZE_Y(); y++)
                    {
                        double NUM_EXEMPLAR_PIXELS = (double)(Exemplars[exi]->getSIZE_X()*Exemplars[exi]->getSIZE_Y());
                        double weight = Exemplars[exi]->getTexelWeight(x, y);
                        double sampleH = (double)posHists[exi][x + y*Exemplars[exi]->getSIZE_X()] / (double)NUM_PIXELS; 
                        double exH = (NUM_PIXELS / NUM_EXEMPLAR_PIXELS) / NUM_PIXELS;
                        weight = weight / (1 + std::max(0.0, sampleH - exH) );

                        Exemplars[exi]->setTexelWeight(x, y, weight); 
                    }
            }
            */
           
            // Main entry point of the optimization. Go through every pixel once
            #pragma omp parallel for num_threads(NUM_CORES)
            for(int pixi=0; pixi<NUM_PIXELS ; pixi++) 
            {
                int ii = I[pixi];
                int jj = J[pixi];
                int kk = K[pixi];
                
                vector<double> weighted_values(ANGMAP_N, 0.0);
                double sum_of_weights = 0;
           
                // At each voxel's optimization we need to look at what each neighborhood
                // it is involved with is suggesting it should be. We want to keep track
                // of the list of these values as well its location in the exemplar image.
                // This will allow us to do k-coherence search and discrete solving.
				int TotalCandidates = 0;
				vector<int> num_candidates(NumExemplars);
                vector< vector<double> > candidate_vals(NumExemplars, vector<double>(NB_SIZE*NB_SIZE));
                vector< vector<double> > candidate_locX(NumExemplars, vector<double>(NB_SIZE*NB_SIZE));
                vector< vector<double> > candidate_locY(NumExemplars, vector<double>(NB_SIZE*NB_SIZE));

                // For each exemplar, build the candidate list
                for(int exi=0; exi < NumExemplars; exi++)
                {
                    vector<double> VolFracsE = Exemplars[exi]->getVolumeFractions();

                    // Reset the candidate list for this voxel
                    num_candidates[exi] = 0;

                    // Get the coordinate offsets that describe this exemplars neighborhood.
                    vector< vector<int> > nbOffsets = Exemplars[exi]->getNBOffsets();
                    
                    // For each exemplar, we need to iterate over all voxels involved in the 
                    // neighborhood centered around the pixel we are optimizing.
                    for(int candi=0; candi<nbOffsets.size(); candi++) 
                    {
                        int pp = ii+nbOffsets[candi][0];
                        int ll = jj+nbOffsets[candi][1];
                        int qq = kk+nbOffsets[candi][2];

                        // tx, ty define the offset within the exemplar that the pixel
						// pp,ll,qq should overlap with. XXX Make General XXX
                        int	tx = candi % NB_SIZE;
						int ty = candi / NB_SIZE;

                        // If this pixel goes outside of the image lets skip it.
                        if(pp < 0 || pp > S.getSizeX()-1 || 
                           ll < 0 || ll > S.getSizeY()-1 || 
                           qq < 0 || qq > S.getSizeZ()-1)
                            continue;
                        
                        // At the voxel centered at  pp,ll,qq, what is the best matching neighborhood.
                        // in exemplar exi. Get the index. Need to subtract one because these are MATLAB
                        // indices and they are 1 based.
                        int nnb = NNB_Tables.at(pp, ll, qq, exi) - 1;

						// If the nearest neighbor index for this exemplar is negative then that means
						// we are near the edge and this neihgborhoods goes outside the bounds of the 
						// image. We need to skip it.
						if(nnb < 0)
							continue;

                        // Lookup the location of this neighborhood in the exemplar image
                        int nb_x = Exemplars[exi]->getXOfNB(nnb);
                        int nb_y = Exemplars[exi]->getYOfNB(nnb);

                        // Figure out this voxels location in the neighborhood
                        int vx = NB_SIZE - tx - 1;
                        int vy = NB_SIZE - ty - 1;
                        
                        // This is the value that the neighborhood centered 
                        // at pp,ll,qq says ii,jj,kk should be.
                        double val = Exemplars[exi]->at(nb_x+vx, nb_y+vy);
                        
                        // Get the weight for this exemplar texel
                        double weight = Exemplars[exi]->getTexelWeight(nb_x+vx, nb_y+vy);  

                        // If this isn't a grain map then just add the exemplar value to
                        // the candidate list
                        if(!isGrainMap)
                        {
                            candidate_vals[exi][num_candidates[exi]] = val;
                            weighted_values[0] += val*weight;
                        }
                        // If this neighborhood is a grain map. We need to figure out a grain id
                        // equivalency list for it placement in the current image. This involves
                        // scanning over the neighborhood and its overlapping region in the 
                        // reconstruction image.
                        else if(isGrainMap)
                        {

                            /*
                            vector<double> grainIdList;
							grainIdList.reserve(NB_SIZE*NB_SIZE);
                            for(int ci=0; ci<nbOffsets.size(); ci++)
                            {
                                int aa = pp+nbOffsets[ci][0];
                                int bb = ll+nbOffsets[ci][1];
                                int cc = qq+nbOffsets[ci][2];

                                // If this pixel goes outside of the image lets skip it.
                                if(aa < 0 || aa > S.getSizeX()-1 || 
                                   bb < 0 || bb > S.getSizeY()-1 || 
                                   cc < 0 || cc > S.getSizeZ()-1)
                                    continue;

								// tx, ty define the offset within the exemplar that the pixel
                                // pp,ll,qq should overlap with. XXX Make General XXX
                                int otx = ci % NB_SIZE;
                                int oty = ci / NB_SIZE;

                                // Figure out this voxels location in the neighborhood
                                //int ovx = NB_SIZE - otx - 1;
                                //int ovy = NB_SIZE - oty - 1;
                                
                                // This is the value that the neighborhood centered 
                                // at pp,ll,qq says aa,bb,cc should be.
                                double oval = Exemplars[exi]->at(nb_x+otx, nb_y+oty);
                                
                                if(oval == val)
                                    grainIdList.push_back(S.at(aa, bb, cc, 0));
								
                            }

							assert(grainIdList.size() > 0);
							
							sort(grainIdList.begin(), grainIdList.end());

                            val = getMode(grainIdList)[0];
                            */

                            int ang_index = val - 1;

                            candidate_vals[exi][num_candidates[exi]] = val;

                            for(int yy=0;yy<ANGMAP_N;yy++)
                                weighted_values[yy] += weight*AngMap[ang_index + yy*ANGMAP_NUM_ORIENTS];
                        }
                       

                        sum_of_weights = sum_of_weights + weight;

                        candidate_locX[exi][num_candidates[exi]] = nb_x+vx;
                        candidate_locY[exi][num_candidates[exi]] = nb_y+vy;
                        num_candidates[exi]++;
                        TotalCandidates++;
                    }

                }

                // Calculate the weighted least square solution
                for(int yy=0;yy<ANGMAP_N;yy++)
                    weighted_values[yy] = weighted_values[yy] / sum_of_weights;

				if(TotalCandidates == 0)
					continue;

                /*
                // Combined all the votes into one long vector
                vector<double> allVals;
                for(int exi=0; exi<NumExemplars; exi++)
					if(num_candidates[exi] > 0)
						allVals.insert(allVals.end(), candidate_vals[exi].begin(), candidate_vals[exi].begin()+num_candidates[exi]);
				
				sort(allVals.begin(), allVals.end());

                // Find the mode of all the votes
                vector<double> modes = getMode(allVals);

                double optimum_value = modes[0];
                */

                // For each exemplar, look at the portion of the candidate set that came
                // from it and see which value is closest to the weighted least square
                // solution.
                double optimum_value = 0;
                double minDist = 1e308;
                for(int exi=0; exi<NumExemplars; exi++)
                {
                    // Find the closest value for this exemplar
                    double tmpDist = 0.0;
                    //int tmpMinI = FindClosestCandidate(candidate_vals[exi], num_candidates[exi], weighted_values, tmpDist);    
                    int tmpMinI = FindClosestValue(candidate_vals[exi], num_candidates[exi], weighted_values[0]);

                    // We found the closest texel for this exemplar. Lets make this the source
                    // for this particular exemplar. 
                    int x = (int)candidate_locX[exi][tmpMinI];
                    int y = (int)candidate_locY[exi][tmpMinI];
                    SourceTable->at(ii, jj, kk, exi) = Exemplars[exi]->getSIZE_X() * y + x;
                    
                    if(tmpDist < minDist)
                    {
                        minDist = tmpDist; 
                        optimum_value = candidate_vals[exi][tmpMinI];
                    }

                }
                
                // Set the value of the pixel
                //S.at(ii, jj, kk, 0) = optimum_value;
                S.at(ii, jj, kk, 0) = optimum_value;

            } // For all pixels 

		//	for(int i=0; i<NumExemplars; i++)
		//		delete posHists[i];

        }


        int FindClosestCandidate(vector<double> vals, int num_candidates, vector<double> target, double & minDist)
        {
            minDist = 1e200;
            int minI = 0;
            for(int zz=0; zz<num_candidates; zz++) {
                double dist = 0.0;
                if(hasAngMap)
                    for(int yy=0;yy<target.size();yy++)
                        dist += fabs(AngMap[((int)vals[zz]-1) + yy*ANGMAP_NUM_ORIENTS] - target[yy]); 
                else
                    dist = fabs(vals[zz] - target[0]);
                    
                if(dist < minDist) { 
                    minI = zz;
                    minDist = dist; 
                }
            }
            return(minI);
        }

        int getNumExemplars() { return NumExemplars; }
        vector<Exemplar *> getExemplars() { return Exemplars; }
        vector<mxArray *> getExemplarWeights() { return ExemplarWeights; }

        // Return the dimensions of the reconstruction
        int getSIZE_X() { return SIZE_X; }
        int getSIZE_Y() { return SIZE_Y; }
        int getSIZE_Z() { return SIZE_Z; }

        // Get the number of cores has said we should run the code with.
        int getNumCores() { return NUM_CORES; }

        // Get the neighborhood size.
        int getNB_SIZE() { return NB_SIZE; }

        mxArray * getSourceTableArray() { return SourceTableArray; }

        private:

        // Round a number
        double round(double number)
        {
            return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
        }

        // Generate a random sequence of all pixels in the sample. We want to
        // optimize the sample in a random order because we don't want to 
        // introduce any directionality to the procedure.
        void GenerateRandomPixelSequence()
        {
            // Allocate space for the random sequence of pixels
            if(I == NULL) {
                I = new int[NUM_PIXELS];
                J = new int[NUM_PIXELS];
                K = new int[NUM_PIXELS];
            }

            // Make and array of all indices for the pixels in order
            int l=0;
            for (int i=0; i<SIZE_X; i++)
                for (int j=0; j<SIZE_Y; j++)
                    for (int k=0; k<SIZE_Z; k++)
                {
                    I[l] = i;
                    J[l] = j;
                    K[l] = k;
                    l++;
                }

            // Shuffle the ordering of the pixels to make it random
            shuffle(I, J, K, NUM_PIXELS);
        }

        /* Arrange the N elements of the ARRAYs in random order.
           effective if N is much smaller than RAND_MAX;
           this may not be the case, use a better random
           generator. */
        void shuffle(int *array1, int *array2, int *array3, size_t n)
        {
            if (n > 1) 
            {
                size_t i;
                for (i = 0; i < n - 1; i++) {
                    size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
                    int t = array1[j];
                    array1[j] = array1[i];
                    array1[i] = t;

                    t = array2[j];
                    array2[j] = array2[i];
                    array2[i] = t;
                    
                    t = array3[j];
                    array3[j] = array3[i];
                    array3[i] = t;
                }
            }
        }

        // Store the neighborhood size, it should be the same for all exemplars
        int NB_SIZE;

        // Number of pixels\voxels in the 2D neighborhood of NB_SIZE
        int NUM_NB_VOXELS;

        // The exemplar images we are using for reconstruction
        int NumExemplars;
        vector<Exemplar *> Exemplars;
        
        // The size of the reconstruction
        int SIZE_X, SIZE_Y, SIZE_Z;
        int NUM_PIXELS;

        // The number of components per pixel
        int NUM_COMPONENTS;

        // This is the nearest neighbor table. It tells us for
        // each pixel in the reconstruction, what is the best matching
        // neighborhood in the appropriate direction.
        Image4D NNB_Tables;

        // This is a matrix that contains which exemplar pixels where
        // copied to the output image. This will allow us to do k-coherence
        // search later.
        Image4D * SourceTable;
        mxArray * SourceTableArray;

        // All the pixels in the sample listed randomly. The list determines the
        // order we go through the sample during optimization.
        int *I;
        int *J;
        int *K;

        vector<mxArray*> ExemplarWeights;

        // The number of cores to use
        int NUM_CORES;

        // Is the image we a reconstructing a grain map
        bool isGrainMap;

        // If the exemplars are index images they may have a angle map
        double * AngMap;

        // The dimensions of the angle map
        int ANGMAP_NUM_ORIENTS;
        int ANGMAP_N;

        // Do the exemplars have and angle map
        bool hasAngMap;
};

#endif
