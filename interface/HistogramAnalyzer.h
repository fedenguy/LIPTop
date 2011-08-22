#ifndef histogramanalyzer_h
#define histogramanalyzer_h

#include <vector>
#include "TH1D.h" 
#include "TF1.h" 

class HistogramAnalyzer
{

 public:

  enum Measurements{kIntegral, kMPV, kMean, kRMS, kSkewness, kKurtosis, k10p, k25p, k75p, k90p };

  /**
     @short CTOR
   */
  HistogramAnalyzer() : fitFunc_ ( new TF1("mpvFitFunc","gaus",0,1000) ) { }
    
    /**
       @short returns the variables of interest (normalized to the MPV fit from gaussian)
     */
    std::vector<Double_t> analyzeHistogram(TH1F *h);

    /**
       @short DTOR
     */
  ~HistogramAnalyzer(){};

 private:

  TF1 *fitFunc_;

  std::vector<Double_t> histoMeasurements;
};


#endif
