#ifndef HISTOREADER_H_
#define HISTOREADER_H_

#include <memory>
#include <string>
#include <vector>

class Sample;
class SampleHist;
class Systematic;
class TH1D;
class TRExFit;

/**
  * \class HistoReader
  * \brief Class that reads histograms as an input
  */

class HistoReader {

    public:
        /**
          * The constructor
          * @param fitter A pointer to TRExFitter class
          */ 
        explicit HistoReader(TRExFit* fitter);

        /**
          * Destructor that does nothing
          */
        ~HistoReader();

        /**
          *
          * Deleted constructors and assignment operators
          */
        HistoReader() = delete;  
        HistoReader(const HistoReader& h) = delete;  
        HistoReader(HistoReader&& h) = delete;  
        HistoReader& operator=(const HistoReader& h) = delete;  
        HistoReader& operator=(HistoReader&& h) = delete;  

        /**
          * method that reads the histograms from user defined inputs
          */
        void ReadHistograms();
    
        /**
          * method that reads the histograms produced by TREx
          */
        void ReadTRExProducedHistograms();

    private:
        /// A pointer to the TRExFit class
        TRExFit* fFitter;
    
        /**
          *
          * A helper function to read singleHistogram from a file
          * @param vector of histo paths
          * @param pointer to systematic
          * @param index of channel
          * @param index of sample
          * @param flag if the systematic is up or down
          * @param flag if we are reading MC or Data
          * @return the read histogram
          */
         std::unique_ptr<TH1D> ReadSingleHistogram(const std::vector<std::string>& fullPaths,
                                                   Systematic* syst,
                                                   int i_ch,
                                                   int i_smp,
                                                   bool isUp,
                                                   bool isMC);

        /**
          * A helper function to read one region for data samples
          * @param index of the region
          * @param is data
          */  
        void ReadOneRegion(const int i_ch, const bool is_data);

        /**
          * A helper function to read shape and norm factors for samples
          * @param channel index
          * @param pointer to a sample
          */ 
        void ReadNormShape(SampleHist* sh,
                           const int i_ch,
                           const Sample* smp);

        /**
          * A helper function to read systeamtics
          * @param channel index
          * @param pointer to a systematic
          * @return true if systematics is overall
          */
        bool SetSystematics(const int i_ch,
                            Sample* ismp,
                            Systematic* syst);

}; 

#endif
