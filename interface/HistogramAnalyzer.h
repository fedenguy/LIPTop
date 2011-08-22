#ifndef histogramanalyzer_h
#define histogramanalyzer_h

#include <vector>
#include <map>

#include "TH1D.h" 
#include "TF1.h" 
#include "TString.h"

class HistogramAnalyzer
{

 public:

  /**
     @short CTOR
   */
  HistogramAnalyzer() : fitFunc_ ( new TF1("mpvFitFunc","gaus",0,1000) ) { }
    
    /**
       @short returns the variables of interest (normalized to the MPV fit from gaussian)
     */
    std::map<TString,Double_t> analyzeHistogram(TH1F *h);
    
    std::map<TString,Double_t> getMeasurements() { return histoMeasurements; }

    /**
       @short DTOR
     */
  ~HistogramAnalyzer(){};

 private:

  TF1 *fitFunc_;

  std::map<TString,Double_t> histoMeasurements;
};


#endif
