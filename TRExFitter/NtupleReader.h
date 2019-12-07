#ifndef NTUPLEREADER_H_
#define NTUPLEREADER_H_

class TRExFit;

/**
 * \class NtupleReader
 * \brief Class that processes ntuples and passes the information to TRExFit
 */

class NtupleReader {

    public:
        /**
          * The constructor
          * @param fitter a pointer to TRExFitter class
          */
        explicit NtupleReader(TRExFit* fitter);

        /**
          * The destructor doing nothing
          */
        ~NtupleReader();

        /**
          * Deleted constructors and assignment operators
          */
        NtupleReader(const NtupleReader& n) = delete;
        NtupleReader(NtupleReader&& n) = delete;
        NtupleReader& operator=(const NtupleReader& n) = delete;
        NtupleReader& operator=(NtupleReader&& n) = delete;

        /**
          * Method that reads the ntuples
          */
        void ReadNtuples();

    private: 
        /**
          * Pointer to the TRExFitClass
          */
        TRExFit* fFitter; 

        /**
          * A helper function to get a single variable
          * @param regIter Index of a region
          */   
        void DefineVariable(int regIter);
};

#endif
