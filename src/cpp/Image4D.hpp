// This structure makes it easier to access an image of 4 dimensions. That is,
// and image where each voxel can have multiple components.
class Image4D
{
    public:

        Image4D() : data(NULL), SIZE_X(0), SIZE_Y(0), SIZE_Z(0), NUM_COMPONENTS(0) {}

        // Construct and 4d image from pre-existing data. Do not copy
        Image4D(double *data, int szX, int szY, int szZ, int nComp)
        : data(data), SIZE_X(szX), SIZE_Y(szY), SIZE_Z(szZ), NUM_COMPONENTS(nComp), isAllocated(false) { }
        
        // Construct and 4D image, allocate it and handle memory management
        Image4D(int szX, int szY, int szZ, int nComp)
        : SIZE_X(szX), SIZE_Y(szY), SIZE_Z(szZ), NUM_COMPONENTS(nComp), isAllocated(true)
        { 
            data = new double[getTotalSize()];
        }

        ~Image4D()
        {
            if(isAllocated)
                delete data;
        }
        
        // Getters
        int getSizeX() { return SIZE_X; }
        int getSizeY() { return SIZE_Y; }
        int getSizeZ() { return SIZE_Z; }
        int getNumComponents() { return NUM_COMPONENTS; } // The fourth dimension
        int getTotalSize() { return SIZE_X*SIZE_Y*SIZE_Z*NUM_COMPONENTS; }
        double *getData() { return data; }

        // Get the data[row,col,level,component] value
        inline double & at(int row, int col, int level, int component) {
            return data[row + col*SIZE_X + level*SIZE_X*SIZE_Y + component*SIZE_X*SIZE_Y*SIZE_Z];
        }

        // When the image is 3D with one component we can use this.
        inline double & at(int row, int col, int level) {
            return data[row + col*SIZE_X + level*SIZE_X*SIZE_Y];
        }

        // When the Image is 2D with one component we can use this
        inline double & at(int row, int col) {
            return data[row + col*SIZE_X];
        }

    private:

        // Do we need to worry about memory management of this data
        bool isAllocated;

        // This is the real components of the data
        double *data;

        // The dimensions of the image
        int SIZE_X, SIZE_Y, SIZE_Z, NUM_COMPONENTS;
};


