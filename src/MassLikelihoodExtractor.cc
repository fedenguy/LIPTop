#include "LIP/Top/interface/MassLikelihoodExtractor.h"

#include "TH1F.h"
#include "LIP/Top/interface/HistogramAnalyzer.h"

using namespace std;

MassLikelihoodExtractor::MassLikelihoodExtractor(bool debug = false):
  debug_(debug),
  hcomb1_(0), 
  hcomb2_(0)
{
  totalMassLikelihood_ = new TH1F("totalMassLikelihood","M_{t} likelihood from KINb solutions;M_{t} [GeV/c^{2}];Likelihood / (2.5 GeV/c^{2})",800,0,2000);
  totalMassLikelihood_->Reset();


  // hcomb1_ = new TH1F();
  // hcomb2_ = new TH1F();

}

MassLikelihoodExtractor::~MassLikelihoodExtractor()
{
  
}

void MassLikelihoodExtractor::processEvent(KinResultsHandler& kinHandler)
{
  HistogramAnalyzer histoAnalyzer;

  // Loading event: if not loaded, it is still not possible to use the same event already loaded outside of this
  //  t.GetEntry(event);

  // Eventually, acquire event info in order to label histograms and events and output files

  // Analyze first combination
  TH1F *h1=kinHandler.getHisto("mt",1);
  if(hcomb1_==0) { hcomb1_ = (TH1F *) h1->Clone("hcomb1"); hcomb1_->SetDirectory(0); }
  else           { hcomb1_->Add(h1); }
  std::map<TString,Double_t> res = histoAnalyzer.analyzeHistogram(hcomb1_);
  double mtop = kinHandler.getMPVEstimate(hcomb1_) [1];
  if(debug_)
    {
      cout << "Combination #1 has " << hcomb1_->Integral() << "solutions and mtop=" << mtop << endl;
      for(std::map<TString, Double_t>::iterator it = res.begin(); it != res.end(); it++) cout << "\t" << it->first << "=" << it->second << endl;
    }
  
  // Analyze second combination
  TH1F *h2=kinHandler.getHisto("mt",2);
  if(hcomb2_==0) { hcomb2_ = (TH1F *) h2->Clone("hcomb2"); hcomb2_->SetDirectory(0); }
  else           { hcomb2_->Add(h2); }
  res=histoAnalyzer.analyzeHistogram(hcomb2_);
  mtop = kinHandler.getMPVEstimate(hcomb2_) [1];
  if(debug_)
    {
      cout << "Combination #1 has " << hcomb2_->Integral() << "solutions and mtop=" << mtop << endl;
      for(std::map<TString, Double_t>::iterator it = res.begin(); it != res.end(); it++) cout << "\t" << it->first << "=" << it->second << endl;
    }

  // Baseline
  applyBaseline(*hcomb1_,1);
  applyBaseline(*hcomb2_,2);

  // Obtain event likelihood
  hcomb1_->Add(hcomb2_);
  
  evMassLikelihood_ = hcomb1_;

  for(Int_t i=0; i<=totalMassLikelihood_->GetNbinsX(); i++)
    totalMassLikelihood_->SetBinContent(i,TMath::Log(evMassLikelihood_->GetBinContent(i)));

  hcomb1_->Reset();
  hcomb2_->Reset();
}

TH1F MassLikelihoodExtractor::getEventLikelihood(/* event parameter */)
{
  return *evMassLikelihood_;
}

TH1F MassLikelihoodExtractor::getFinalLikelihood()
{
  return *totalMassLikelihood_;
}

void MassLikelihoodExtractor::applyBaseline(TH1F& h, int c)
{
  Double_t baseline = getBaselineFromHisto(h);

  if(debug_)
    {
      cout << "Baseline for combination #" << c << ": " << baseline << endl;
    }
  
  for(Int_t i=0; i<h.GetNbinsX(); i++)
    h.AddBinContent(i, baseline);
  
  return;
}

Double_t MassLikelihoodExtractor::getBaselineFromHisto(TH1F& h)
{
  Double_t min(0);
  bool minLock(false);
  Double_t max(0);

  for(Int_t i=0; i<h.GetNbinsX(); i++){
    if(h.GetBinContent(i) != 0){
      if(minLock == false){
	min = h.GetBinCenter(i);
	minLock = true;
      }
      max = h.GetBinCenter(i);
    }
  }

  Double_t result = 100 * (max - min) / fabs( 50000 - h.GetEntries() );
  
  if(debug_){
    cout << "Max: " << max << ", Min: " << min << " N: " << h.GetEntries() << endl;
    cout << "Baseline: " << result << endl;
  }
  
  return result;
}
