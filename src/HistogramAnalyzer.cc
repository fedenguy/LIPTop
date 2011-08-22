#include "LIP/Top/interface/HistogramAnalyzer.h"

using namespace std;

//
std::vector<Double_t> HistogramAnalyzer::analyzeHistogram(TH1F *h)
{
  histoMeasurements.clear();
  histoMeasurements.resize(10,0);
  if(h==0) return histoMeasurements;

  histoMeasurements[kIntegral] = h->Integral();
  if(histoMeasurements[kIntegral]<2) return histoMeasurements;

  //quantiles for 10, 25, 50, 75, 90 %
  Double_t pq[5]={0.1,0.25,0.5,0.75,0.9};
  Double_t xq[5];
  h->GetQuantiles(5,xq,pq);
  Double_t median=xq[2];
  histoMeasurements[k10p]=(xq[0]-median)/median;
  histoMeasurements[k25p]=(xq[1]-median)/median;
  histoMeasurements[k75p]=(xq[3]-median)/median;
  histoMeasurements[k90p]=(xq[4]-median)/median;

  //fit a gaussian near the most probable value
  Int_t iBin = h->GetMaximumBin();
  double mpv = h->GetXaxis()->GetBinCenter(iBin);
  fitFunc_->SetRange(mpv-25,mpv+25);
  fitFunc_->SetParLimits(1,mpv-10,mpv+10);
  h->Fit(fitFunc_,"LRQN");
  histoMeasurements[kMPV]      = (fitFunc_->GetParameter(1)-median)/median;
  histoMeasurements[kMean]     = (h->GetMean()-median)/median;
  histoMeasurements[kRMS]      = h->GetRMS()/median;
  histoMeasurements[kSkewness] = h->GetSkewness()/median;
  histoMeasurements[kKurtosis] = h->GetKurtosis()/median;
  
  return histoMeasurements;
}
