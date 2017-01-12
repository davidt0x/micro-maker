#include <cassert>

typedef std::vector<double> vec_type ;

// PRE: v.empty() != true; v is sorted.
// POST: a non-empty vector of distinct modal values is returned.
vec_type getMode( const vec_type& v )
{
    vec_type modes ;
    vec_type::const_iterator it = v.begin() ;

    double runValue = *it++ ;
    unsigned runCount = 1 ;
    unsigned highestRunCount = runCount ;
    modes.push_back(runValue) ;

    while ( it != v.end() )
    {
        if ( runValue == *it )                      // run continuing.
        {
            if ( ++runCount > highestRunCount )     // current run is new high?
            {
                highestRunCount = runCount ;

                if ( modes.front() != runValue )    // then there should only be one mode
                {
                    modes.clear() ;
                    modes.push_back(runValue) ;
                }
            }
            else if ( runCount == highestRunCount ) // more than one mode, currently
                modes.push_back(runValue) ;
        }
        else                                        // new run beginning
            runValue = *it, runCount = 0 ;

        ++it ;
    }

    return modes ;
}


class Histogram
{

    public:

        Histogram() : NUM_BINS(16) {}

        Histogram(int num_bins) : NUM_BINS(num_bins) {}

        void calcHistogram(vector<double> & D, int lenD)
        {
            TOTAL_COUNTS = 0;

            // Initialize the histogram to all 0's
            hist.resize(NUM_BINS);
            for(int j=0;j<NUM_BINS;j++) hist[j] = 0;    

            for (int i=0;i<lenD; i++)
            {
                //assert(D[i] >= 0.0 && D[i] <= 1.0);
                valueMin = min(valueMin, D[i]);
                valueMax = max(valueMax, D[i]);
            }

            binWidth = 1.0 / (double)NUM_BINS;

            // Compute the histogram
            for (int i=0;i<lenD; i++)
                addCountForValue(D[i]);
        }

        int inline getCount(int i) { return hist[i]; }

        inline int getCountForValue(double value)
        {
            return hist[ getIndex(value) ];
        }

        inline double getProbForValue(double value)
        {
            return (double)hist[ getIndex(value) ] / (double)TOTAL_COUNTS;
        }

        inline void addCountForValue(double value)
        {
            hist[ getIndex(value) ]++;
            TOTAL_COUNTS++;
        }

        inline void removeCountForValue(double value)
        {
            hist[ getIndex(value) ]--;
            TOTAL_COUNTS--;
        }

        inline int getTotalCounts() { return TOTAL_COUNTS; }

    private:

        inline int getIndex(double value)
        {
            int index = (int)( value / binWidth );
            if(index >= NUM_BINS)
                index = NUM_BINS - 1;
            return index;
        }

        // The counts for the histogram
        vector<int> hist;

        // Number of bins in the histogram
        int NUM_BINS;

        int TOTAL_COUNTS;

        double valueMin;
        double valueMax;
        double binWidth;

};
