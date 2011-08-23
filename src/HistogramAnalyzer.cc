#include "LIP/Top/interface/HistogramAnalyzer.h"

using namespace std;

//
std::map<TString, Double_t> HistogramAnalyzer::analyzeHistogram(TH1F *h)
{
  if(h==0) return histoMeasurements;

  histoMeasurements["kIntegral"] = h->Integral();
  if(histoMeasurements["kIntegral"]<2) return histoMeasurements;
  double mean=h->GetMean();

  //quantiles for 10, 25, 50, 75, 90 %
  Double_t pq[5]={0.1,0.25,0.5,0.75,0.9};
  Double_t xq[5];
  h->GetQuantiles(5,xq,pq);

  //width to the median
  histoMeasurements["k10p50p"]=(xq[2]-xq[0])/mean;
  histoMeasurements["k25p50p"]=(xq[2]-xq[1])/mean;
  histoMeasurements["k75p50p"]=(xq[3]-xq[2])/mean;
  histoMeasurements["k90p50p"]=(xq[4]-xq[2])/mean;

  //width asymmetries
  histoMeasurements["kAsymm10p"]=(histoMeasurements["k90p50p"]-histoMeasurements["k10p50p"])/mean;//(histoMeasurements["k90p50p"]+histoMeasurements["k10p50p"]);
  histoMeasurements["kAsymm25p"]=(histoMeasurements["k75p50p"]-histoMeasurements["k25p50p"])/mean;//(histoMeasurements["k75p50p"]+histoMeasurements["k25p50p"]);

  //width between different quantiles
  histoMeasurements["k90p10p"]=(xq[4]-xq[0])/mean;
  histoMeasurements["k75p25p"]=(xq[3]-xq[1])/mean;
  histoMeasurements["k75p10p"]=(xq[3]-xq[0])/mean;
  histoMeasurements["k90p25p"]=(xq[4]-xq[1])/mean;  
  histoMeasurements["k90p75p"]=(xq[4]-xq[3])/mean;

  //fit a gaussian near the most probable value
  Int_t iBin = h->GetMaximumBin();
  double mpv = h->GetXaxis()->GetBinCenter(iBin);
  fitFunc_->SetRange(mpv-25,mpv+25);
  fitFunc_->SetParLimits(1,mpv-10,mpv+10);
  h->Fit(fitFunc_,"LRQN");
  histoMeasurements["kMPV"]      = fitFunc_->GetParameter(1)/mean;
  histoMeasurements["kRMS"]      = h->GetRMS()/mean;
  histoMeasurements["kSkewness"] = h->GetSkewness()/mean;
  histoMeasurements["kKurtosis"] = h->GetKurtosis()/mean;
  
  return histoMeasurements;
}
