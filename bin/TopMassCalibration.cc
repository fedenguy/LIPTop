#include "TList.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TSystem.h"
#include "TProfile.h"

#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"
#include "LIP/Top/interface/MassMeasurement.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <vector>

using namespace std; 
using namespace RooFit ;

int maxPE=500;

std::map<TString,std::vector<Float_t> > evtYields;
typedef std::vector<TString> Proc_t;
std::map<TString,Proc_t> allSamples;
typedef std::vector<TH1D *> ProcHistos_t;
std::map<TString,ProcHistos_t> allHistos;
std::vector<TString> signalPts;

TString url,massParsUrl, syst;

TObjArray calibrate(TString mpoint);

//
int main(int argc, char* argv[])
{
  stringstream report;

  // load framework libraries                                                                                          
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  url=argv[1];
  gSystem->ExpandPathName(url);
  massParsUrl=argv[2];
  gSystem->ExpandPathName(massParsUrl);
  syst=argv[3];
  TString calibrationList="";
  if(argc>4)
    {
      calibrationList=argv[4];
      cout << "Calibrating for: " << calibrationList << endl;
    }
  else
    cout  << "Calibrating for all mass points" << endl;
    
  

  //
  // map the samples per type
  //
  std::vector<Float_t> procYields(4,0);    

  //signal samples first (all share the same yields)
  procYields[MassMeasurement::SF_EQ1BTAGS]=577.1;
  procYields[MassMeasurement::SF_GEQ2BTAGS]=1137.0;
  procYields[MassMeasurement::OF_EQ1BTAGS]=792.7;
  procYields[MassMeasurement::OF_GEQ2BTAGS]=1478.6;
  evtYields["Signal"]=procYields;
  if(syst=="baresignal" || syst=="effbup" || syst=="effbdown" || syst=="effqup" || syst=="effqdown")
    {
      allSamples["172.5"]=Proc_t(1,"TTJets_signal");
      allSamples["172.5"]=Proc_t(1,"TTJets");
    }
  if(syst=="matchingdown")
    {
      allSamples["172.5"]=Proc_t(1,"TTJets_matchingdown");
    }
  if(syst=="matchingup")
    {
      allSamples["172.5"]=Proc_t(1,"TTJets_matchingup");
    }
  if(syst=="scaledown")
    {
      allSamples["172.5"]=Proc_t(1,"TTJets_scaledown");
    }
  if(syst=="scaleup")
    {
      allSamples["172.5"]=Proc_t(1,"TTJets_scaleup");
    }
  if(syst=="std")
    {
      //signal samples
      allSamples["161.5"]=Proc_t(1,"TTJets_mass_161v5");
      allSamples["163.5"]=Proc_t(1,"TTJets_mass_163v5");
      allSamples["166.5"]=Proc_t(1,"TTJets_mass_166v5");
      allSamples["169.5"]=Proc_t(1,"TTJets_mass_169v5");
      allSamples["172.5"]=Proc_t(1,"TTJets_signal");
      allSamples["175.5"]=Proc_t(1,"TTJets_mass_175v5");
      allSamples["178.5"]=Proc_t(1,"TTJets_mass_178v5");
      allSamples["181.5"]=Proc_t(1,"TTJets_mass181v5");
      allSamples["184.5"]=Proc_t(1,"TTJets_mass184v5");
      
      //backgrounds
      Proc_t SingleTop;
      SingleTop.push_back("SingleTbar_tW");
      SingleTop.push_back("SingleT_tW");
      SingleTop.push_back("SingleTbar_t");
      SingleTop.push_back("SingleT_t");
      SingleTop.push_back("SingleTbar_s");
      SingleTop.push_back("SingleT_s");
      allSamples["SingleTop"]=SingleTop;
      procYields[MassMeasurement::SF_EQ1BTAGS]=38.5; 
      procYields[MassMeasurement::SF_GEQ2BTAGS]=33.2; 
      procYields[MassMeasurement::OF_EQ1BTAGS]=49.5; 
      procYields[MassMeasurement::OF_GEQ2BTAGS]=43.1;
      evtYields["SingleTop"]=procYields;
      
      Proc_t OtherTTbar;
      OtherTTbar.push_back("TTJets");
      allSamples["OtherTTbar"]=OtherTTbar;
      procYields[MassMeasurement::SF_EQ1BTAGS]=4.7+2.0; 
      procYields[MassMeasurement::SF_GEQ2BTAGS]=4.6+0.7; 
      procYields[MassMeasurement::OF_EQ1BTAGS]=6.3+4.7; 
      procYields[MassMeasurement::OF_GEQ2BTAGS]=8.4;
      evtYields["OtherTTbar"]=procYields;
      
      Proc_t DiBosons;
      DiBosons.push_back("WW");
      DiBosons.push_back("ZZ");
      DiBosons.push_back("WZ");
      allSamples["DiBosons"]=DiBosons;
      procYields[MassMeasurement::SF_EQ1BTAGS]=8.1;
      procYields[MassMeasurement::SF_GEQ2BTAGS]=1.7;
      procYields[MassMeasurement::OF_EQ1BTAGS]=2.0;
      procYields[MassMeasurement::OF_GEQ2BTAGS]=1.9;
      evtYields["DiBosons"]=procYields;
      
      Proc_t DY;
      DY.push_back("DYJetsToLL");
      allSamples["DY"]=DY;
      procYields[MassMeasurement::SF_EQ1BTAGS]=207.5;
      procYields[MassMeasurement::SF_GEQ2BTAGS]=56.8;
      procYields[MassMeasurement::OF_EQ1BTAGS]=29.3;
      procYields[MassMeasurement::OF_GEQ2BTAGS]=5.1;
      evtYields["DY"]=procYields;
    }

  //specific for b-tagging effiency systematics
  float newHeavyFlavorCut(0.207), newLightFlavorCut(0.233);
  if(syst=="effbup")   newHeavyFlavorCut-=0.09;
  if(syst=="effbdown") newHeavyFlavorCut+=0.08;
  if(syst=="effqup")   newLightFlavorCut-=0.02;
  if(syst=="effqdown") newLightFlavorCut+=0.02;
  
  //
  // build the base histograms to sample for the pseudo-experiments
  //
  TFile *f = TFile::Open(url);
  int ipt(1);
  cout << "[Filling templates]" << flush;
  for(std::map<TString,Proc_t>::iterator it = allSamples.begin(); it!= allSamples.end(); it++)
    {
      cout << "." << flush;
      TString sName("m"); sName += (ipt+1);      
      bool isSignal(it->first.IsFloat());
      if(isSignal && (calibrationList=="" || calibrationList.Contains(it->first)) ) signalPts.push_back(it->first);

      //instantiate the histograms
      ProcHistos_t ihistos;
      for(size_t icat=0; icat<4;icat++)
	{
	  TString hname(sName+"_"); hname += icat;
	  TH1D *h= new TH1D(hname,hname,80,100,500);      
	  h->SetDirectory(0);
	  TString tit(it->first); tit+= "(cat="; tit+=icat; tit +=")";
	  h->SetTitle(tit);
	  ihistos.push_back(h);
	}

      //now fill the histograms
      for(Proc_t::iterator pit = it->second.begin(); pit != it->second.end(); pit++)
	{
      	  TString tname=*pit + "/data";
	  TTree *t = (TTree *) f->Get(tname);
	  if(t==0) continue;

	  top::EventSummaryHandler evHandler;
	  if( !evHandler.attachToTree(t) ) continue;
	  for(Int_t i=0; i<evHandler.getEntries(); i++)
	    {
	      evHandler.getEntry(i);
	      top::EventSummary_t &ev = evHandler.getEvent();
	      
	      float mtop=ev.evmeasurements[0];
	      if(mtop==0) continue; 
	
	      int nbtags(ev.evmeasurements[5]);
	      if(syst=="effbup" || syst=="effbdown" || syst=="effqup" || syst=="effqdown")
		{
		  nbtags=0;
		  top::PhysicsEvent_t phys = getPhysicsEventFrom(ev);
		  for(unsigned int ijet=0; ijet<phys.jets.size(); ijet++)
		    {
		      if(phys.jets[ijet].pt()<30 || fabs(phys.jets[ijet].eta())>2.4) continue;
		      float btagDisc=phys.jets[ijet].btag7;
		      if(fabs(phys.jets[ijet].flavid)==5) nbtags += (btagDisc>newHeavyFlavorCut);
		      else                                nbtags += (btagDisc>newLightFlavorCut);
		    }
		}  
	      if(nbtags==0) continue;

	      bool isZcand= bool(ev.evmeasurements[2]);
	      if(isZcand) continue;

	      //check the category of the event
	      int catToFill(0);
	      bool isSF(ev.cat==EE || ev.cat==MUMU);
	      if(nbtags==1) catToFill=(isSF ? MassMeasurement::SF_EQ1BTAGS  : MassMeasurement::OF_EQ1BTAGS);
	      if(nbtags>1)  catToFill=(isSF ? MassMeasurement::SF_GEQ2BTAGS : MassMeasurement::OF_GEQ2BTAGS);

	      //fill the histograms
	      float evWeight=ev.weight*ev.xsecWeight;
	      ihistos[ catToFill ]->Fill(mtop,evWeight);
	    }
	}
      
      //save to collection
      allHistos[it->first]=ihistos;
    }    
  f->Close();
  cout  << endl;

  //
  // run calibration
  //
  cout << "[Running calibration fits] " << signalPts.size() << " tasks will be launched" << endl;

  //keep RooFit quiet
  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().getStream(0).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Fitting);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  
  for(size_t ipt=0; ipt<signalPts.size(); ipt++)
    {
      TString mpoint=signalPts[ipt];
      TObjArray results = calibrate(mpoint);
 
      //save summary of results to file
      TString outUrl=massParsUrl;
      outUrl.ReplaceAll(".txt","_calib.root");
      TFile *outF=TFile::Open(outUrl, "UPDATE");
      outF->rmdir(mpoint);
      TDirectory *outDir = outF->mkdir(mpoint);
      outDir->cd();
      for(Int_t ires=0; ires<results.GetEntriesFast(); ires++)
	{
	  TH1 *h=(TH1 *)results.At(ires);
	  h->SetDirectory(outDir);
	  h->Write();
	}
      outF->Close();
    }

  return -1;

}


// 
TObjArray calibrate(TString mpoint)
{
  float trueMass=mpoint.Atof();
  cout << "[Running calibration for " << trueMass << "]" << endl;

  //init a local fitter
  MassMeasurement mFitter(massParsUrl,mpoint);
  int nCategories=mFitter.getNumberOfCategories();
  TObjArray massFitH,biasH,pullH,statUncH;
  for(int i=0; i<=nCategories; i++)
    {
      TString postfix(""); 
      if(i) { postfix="_"; postfix+=i;}

      TH1D *h = new TH1D("massfit"+postfix,";m_{top}^{fit};Pseudo-experiments",50,trueMass*0.8,trueMass*1.2);  
      h->SetDirectory(0);
      h->SetMarkerStyle(20);
      massFitH.Add(h);

      h = new TH1D("bias"+postfix,";bias = (m_{top}^{fit}-m_{top}^{gen});Pseudo-experiments",50,-25.5,24.5);
      h->SetDirectory(0);
      h->SetMarkerStyle(20);
      biasH.Add(h);

      h = new TH1D("pull"+postfix,";pull = (m_{top}^{fit}-m_{top}^{gen}) / #sigma_{stat};Pseudo-experiments",25,-5.2,4.8);    
      h->SetDirectory(0);
      h->SetMarkerStyle(20);
      pullH.Add(h);

      h = new TH1D("statunc"+postfix,";#sigma_{stat};Pseudo-experiments",400,0,8);
      h->SetDirectory(0);
      h->SetMarkerStyle(20);
      statUncH.Add(h);
    }
        
  //generate n-ensembles and fit
  EnsembleMeasurement_t ensemble;
  for(Int_t ipe=0; ipe<maxPE; ipe++)
    {
      if(ipe%5==0) { printf("\r[\t%s\t] : %d / %d",mpoint.Data(),ipe,maxPE); cout << flush; }
      
      //fill the ensemble
      ensemble.nEvents=0;
      ensemble.status=false;
      ensemble.mass=0;
      ensemble.err=0;
      for(std::map<TString,ProcHistos_t>::iterator itt = allHistos.begin(); itt!=allHistos.end(); itt++)
	{
	  bool isSignal(itt->first.IsFloat());
	  if(isSignal && mpoint != itt->first) continue;
	  //cout << " " << itt->first << "(" << (isSignal ? "sig" : "bckg") << "): " << flush;
	  int procEvts(0);
	  for(size_t ih=0; ih<itt->second.size(); ih++)
	    {
	      //generate randomly from the template histograms
	      int nevExpected=evtYields[isSignal? "Signal" :itt->first][ih];
	      Int_t nevToGenerate=gRandom->Poisson(nevExpected);
	      procEvts+=nevToGenerate;
	      for(int iev=0; iev<nevToGenerate; iev++)
		{
		  ensemble.evMasses[ensemble.nEvents]=itt->second[ih]->GetRandom();
		  ensemble.evCategories[ensemble.nEvents]=Float_t(ih);
		  ensemble.nEvents++;
		}
	    }
	  //cout << procEvts << endl;
	}

      //now fit and monitor the bias, stat. uncertainty and pull
      MassFitResults_t result=mFitter.DoMassFit(ensemble,false);
      for(int i=0; i<=nCategories; i++)
	{
	  float mass    = (i==0 ? result.tMass    : result.iTmass[i-1]);
	  float statUnc = (i==0 ? result.tMassErr : result.iTmassErr[i-1]);
	  float bias    = mass-trueMass;
	  float pull    = bias/statUnc;
	  ((TH1*) massFitH.At(i))->Fill(mass);
	  ((TH1*) biasH.At(i))->Fill(bias);
	  ((TH1*) statUncH.At(i))->Fill(statUnc);
	  ((TH1*) pullH.At(i))->Fill(pull);
	}
    }
  printf("\n");

  //fit a gaussian and return the histograms
  TObjArray results; 
  for(int i=0; i<=nCategories; i++)
    {
      ((TH1*) massFitH.At(i))->Fit("gaus","LMEQ+");
      ((TH1*) biasH.At(i))->Fit("gaus","LMEQ+");
      ((TH1*) pullH.At(i))->Fit("gaus","LMEQ+");

      results.Add(massFitH.At(i));
      results.Add(biasH.At(i));
      results.Add(pullH.At(i));
      results.Add(statUncH.At(i));
    }

  return results;
}
