// This class stores all of the information for an exemplar image.
class Exemplar
{
    public:
        
        // Constructur that builds an Exemplar object 
        Exemplar(double *data, 
                 int NUM_COMPONENTS, 
                 int NB_SIZE,
                 vector< vector<int> > NBLookup,
                 int SIZE_X,
                 int SIZE_Y,
                 vector< vector<int> > nbOff,
                 double *texelWeights) 
            : data(data), NB_SIZE(NB_SIZE), NBLookup(NBLookup), nbOffsets(nbOff), texelWeights(texelWeights),
              NUM_COMPONENTS(NUM_COMPONENTS), SIZE_X(SIZE_X), SIZE_Y(SIZE_Y)
        {
            // Calculate the volume fraction of each component, this is unormalized.
            volFracs.resize(NUM_COMPONENTS);
            vector<double> compMax(NUM_COMPONENTS);
            vector<double> compMin(NUM_COMPONENTS);
            for(int i=0;i<NUM_COMPONENTS;i++)
            {
                volFracs[i] = 0.0;
                for(int j=0;j<SIZE_X*SIZE_Y;j++)
                    volFracs[i] = volFracs[i] + at(j);
                volFracs[i] /= SIZE_X*SIZE_Y;
            }

            // Calculate the histogram for each component
            CompHist.resize(NUM_COMPONENTS);
            vector<double> D(SIZE_X*SIZE_Y);
            for(int i=0;i<NUM_COMPONENTS;i++)
            {
                int j = 0;
                for(int q=0;q<SIZE_X;q++) {
                    for(int r=0;r<SIZE_Y;r++, j++) {
                        D[j] = at(q, r, i);
                    }
                }
                CompHist[i] = Histogram(16);
                CompHist[i].calcHistogram(D, SIZE_X*SIZE_Y);
            }

            // Create a lookup table to convert a linear index to X and Y coords
            // in the image
            ind2sub = new int[SIZE_X*SIZE_Y*2];
            int index = 0;
            for(int i=0; i<SIZE_X; i++) {
                for(int j=0; j<SIZE_Y; j++) {
                    ind2sub[index + 0*(SIZE_X*SIZE_Y)] = i;
                    ind2sub[index + 1*(SIZE_X*SIZE_Y)] = j;
                    index++;
                }
            }
            
        }

        ~Exemplar()
        {
            delete ind2sub;
        }


        // Get the array of values stored at location (row,col) in image
        inline double at(int row, int col, int component)
        {
            return data[row+col*SIZE_X+component*(SIZE_X*SIZE_Y)]; 
        }    

        // This is short hand if they only want the first pixel component
        // Useful when the image is MxNx1.
        inline double at(int row, int col)
        {
            return data[row+col*SIZE_X];
        }

        // This is if they want to use linear indices, zero based
        inline double at(int index)
        {
            return data[index];
        }

        double *getTexelWeights() { return texelWeights; }

        inline double getTexelWeight(int row, int col)
        {
            return texelWeights[row+col*SIZE_X];
        }

        inline void setTexelWeight(int row, int col, double value)
        {
            texelWeights[row+col*SIZE_X] = value;
        }

        // Getters for the private properties
        inline vector< vector<int> > & getNBOffsets() { return nbOffsets; }
        int getSIZE_X() { return SIZE_X; }
        int getSIZE_Y() { return SIZE_Y; }
        int getNUM_COMPONENTS() { return NUM_COMPONENTS; }
        int getNB_SIZE() { return NB_SIZE; }
        inline vector<double> & getVolumeFractions() { return volFracs; }

        int getNumNBs() { return NBLookup.size(); }
        int inline getXOfNB(int index) { return NBLookup[index][0]-1; }
        int inline getYOfNB(int index) { return NBLookup[index][1]-1; }

        vector<Histogram> getComponentHistograms() { return CompHist; }

        inline double getHistBinForValue(double value, int comp) { return CompHist[comp].getProbForValue(value); }

        int inline Index2SubscriptX(int index) { return ind2sub[index]; }
        int inline Index2SubscriptY(int index) { return ind2sub[index + (SIZE_X*SIZE_Y)]; }

    private:

        // The pixel data for this exemplar
        double *data;

        // The pixel dimensionality of this exemplar
        int NUM_COMPONENTS;
        
        // The size of the local neighborhood we are considering. 
        int NB_SIZE;

        // How many unique neighborhoods of size NB_SIZExNB_SIZE are in this
        // image.
        int NUM_NBs;

        // Number of rows in this exemplar image
        int SIZE_X;
    
        // Number of columns in this exemplar image
        int SIZE_Y;

        // This is a lookup table that maps indices of unique neighborhoods to
        // their actual X and Y locations in the exemplar image. 
        vector< vector<int> > NBLookup;

        // The neighborhood configuration of this exemplar within
        // the 3D reconstuction. This is basically the pattern of pixels
        // that define the orientation of this exemplar within the 3D sample.
        vector< vector<int> > nbOffsets;        

        // Volume fractions of each phase in the exemplar
        vector<double> volFracs;

        // A lookup table that converts linear index to 2D coords in the image.
        int * ind2sub;

        // A vector that contains each components histogram
        vector<Histogram> CompHist;

        double * texelWeights;

};


