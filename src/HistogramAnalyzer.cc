#include "LIP/Top/interface/HistogramAnalyzer.h"

using namespace std;

//
std::vector<Double_t> HistogramAnalyzer::analyzeHistogram(TH1F *h)
{
  histoMeasurements.clear();
  histoMeasurements.resize(9,0);
  if(h==0) return histoMeasurements;
  
  histoMeasurements[kIntegral] = h->Integral();
  if(histoMeasurements[kIntegral]<2) return histoMeasurements;

  //fit a gaussian near the most probable value
  Int_t iBin = h->GetMaximumBin();
  double mpv = h->GetXaxis()->GetBinCenter(iBin);
  fitFunc_->SetRange(mpv-25,mpv+25);
  fitFunc_->SetParLimits(1,mpv-10,mpv+10);
  h->Fit(fitFunc_,"LRQN");
  histoMeasurements[kMPV] = fitFunc_->GetParameter(1);

  histoMeasurements[kMean]     = h->GetMean()/histoMeasurements[kMPV];
  histoMeasurements[kRMS]      = h->GetRMS()/histoMeasurements[kMPV];
  histoMeasurements[kSkewness] = h->GetSkewness()/histoMeasurements[kMPV];
  histoMeasurements[kKurtosis] = h->GetKurtosis()/histoMeasurements[kMPV];

  Double_t xq[4];
  h->GetQuantiles(4,xq);
  histoMeasurements[k25p]=(xq[1]-xq[0])/histoMeasurements[kMPV];
  histoMeasurements[k50p]=(xq[2]-xq[1])/histoMeasurements[kMPV];
  histoMeasurements[k75p]=(xq[3]-xq[2])/histoMeasurements[kMPV];

  return histoMeasurements;
}
