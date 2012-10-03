#include "TList.h"
#include "TRandom2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TSystem.h"
#include "TProfile.h"

#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/plotter.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"
#include "CMGTools/HtoZZ2l2nu/src/JSONWrapper.cc"
#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuSummaryHandler.h"
#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuPhysicsEvent.h"
#include "CMGTools/HtoZZ2l2nu/interface/SmartSelectionMonitor.h"
#include "CMGTools/HtoZZ2l2nu/interface/METUtils.h"

#include "LIP/Top/interface/MisassignmentMeasurement.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include <vector>

using namespace std;

stringstream report; 
float dataLumi         = 5041.0;
int maxJets            = 4;
int maxMixTries        = 100;
float sfMetCut         = 40;
TString jesUncFile     = "${CMSSW_BASE}/src/CMGTools/HtoZZ2l2nu/data/GR_R_42_V23_AK5PFchs_Uncertainty.txt";
JetCorrectionUncertainty *jetCorUnc=0;
TChain *dataChain, *signalChain;
std::map<TString, TChain *> bckgChain;
SmartSelectionMonitor controlHistos;
TRandom2 calibRndGen;

//
void printHelp();
void fillEventsChain(TString url,JSONWrapper::Object &jsonF);
void buildMljTemplates();
void fillMljTemplatesFrom(TChain *c,bool isData, bool isSignal, bool hasTop, bool isDY,bool isVV);

////CHECK ME
std::map<TString,Int_t> procYields;
std::map<TString,TTree *> procTrees,dataTrees;
std::map<std::pair<int,TString>, TGraphErrors*> calibCurves;
std::map<int, MisassignmentMeasurement *> exclusiveMisMeasurements;
std::map<std::pair<int, TString>, std::vector<double> > fCorrectTrue, fCorrectSigOnlyTrue, fCorrectData;
void runCalibration(TString url, JSONWrapper::Object &jsonF, Double_t lumi, int maxPE, Int_t jetbin, bool freezeResults);
void runPileupCheck(TString url);
////////


//                                                                                                                                                                                          
void printHelp()
{
  printf("--help    --> print this\n");
  printf("--in      --> input directory with the summary trees\n");
  printf("--json    --> json file with the description of the samples\n");
}

//
void fillEventsChain(TString url, JSONWrapper::Object &jsonF)
{
  //iterate over the processes
  std::vector<JSONWrapper::Object> Process = jsonF["proc"].daughters();
  for(unsigned int i=0;i<Process.size();i++)
    {
      TString procCtr(""); procCtr+=i;
      TString proc=(Process[i])["tag"].toString();
      bool isData(Process[i]["isdata"].toBool()); 
      bool isSignal(false);
      if(Process[i].isTag("issignal")) isSignal=Process[i]["issignal"].toBool();

      TChain *c=0;
      if(isData)
	{
	  if(dataChain==0) dataChain = new TChain;
	  c=dataChain;
	}
      else if(isSignal)
	{
	  if(signalChain==0) signalChain = new TChain;
	  c=signalChain;
	}
      
      //iterator over the sub-processes
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(size_t id=0; id<Samples.size(); id++)
	{
	  int split=1;
	  if(Samples[id].isTag("split"))split = Samples[id]["split"].toInt();
	  if(!isSignal && !isData)
	    {
	      c = new TChain;
	      bckgChain[ (Samples[id])["dtag"].toString() ] = c;
	    }
	  
	  //in case there was splitting of the analysis 
	  for(int isplit=0; isplit<split; isplit++)
	    {
	      TString tname((Samples[id])["dtag"].toString());
	      if(split>1) { tname += "_"; tname+=isplit; }
	      TString fname(url+tname+"_summary.root");
	      //logic is inverted in the ROOT API...
	      if(!gSystem->AccessPathName(fname))  c->AddFile(fname,TChain::kBigNumber,tname+"/data");
	    }
  
	}
    }

  //report what was found
  report << "[Samples used]" << endl;
  report << "\tData : ";     if(dataChain==0)   report << "n/a" << endl; else report << dataChain->GetEntries() << " events" << endl;
  report << "\tSignal : ";   if(signalChain==0) report << "n/a" << endl; else report << signalChain->GetEntries() << " events" << endl;
  for(std::map<TString,TChain *>::iterator cIt=bckgChain.begin(); cIt!=bckgChain.end(); cIt++)
    report << "\t" << cIt->first  << " : " << cIt->second->GetEntries() << " events" << endl;
  
}

//
void buildMljTemplates()
{
  TH1F *baseMlj= new TH1F("mlj",";Invariant Mass [GeV];Lepton-jet pairs",100,0,1000);
  TString systVars[] = {"","jesup","jesdown","jerup","jerdown"};
  TString mljVars[]  = {"","correct","misassigned","wrong","unmatched","rot","swap"};
  for(size_t iSystVar=0; iSystVar<sizeof(systVars)/sizeof(TString); iSystVar++)
    for(size_t imljVar=0; imljVar<sizeof(mljVars)/sizeof(TString); imljVar++)
      controlHistos.addHistogram( (TH1F *) baseMlj->Clone(mljVars[imljVar]+"mlj"+systVars[iSystVar]) );
  baseMlj->Delete();

  if(jetCorUnc==0)
    {
      gSystem->ExpandPathName(jesUncFile);  
      jetCorUnc = new JetCorrectionUncertainty(jesUncFile.Data());
    }
  
  fillMljTemplatesFrom(dataChain,  true,false,false,false,false);  cout << "." << flush;
  fillMljTemplatesFrom(signalChain,false,true,false,false,false);  cout << "." << flush;
  for(std::map<TString,TChain *>::iterator it = bckgChain.begin(); it != bckgChain.end(); it++)
    {
      bool hasTop(it->first.Contains("SingleT") || it->first.Contains("TTJets"));
      bool isDY(it->first.Contains("DYJets"));
      bool isVV(it->first.Contains("WW") || it->first.Contains("ZZ") || it->first.Contains("WZ"));
      fillMljTemplatesFrom(it->second,  false,false,hasTop,isDY,isVV);  cout << "." << flush;
    }
  cout << endl;
}

//
void fillMljTemplatesFrom(TChain *c,bool isData, bool isSignal, bool hasTop, bool isDY,bool isVV)
{
  //attach to tree
  if(c==0) return;
  const unsigned int iEntries=c->GetEntries();
  if(iEntries==0) return;
  ZZ2l2nuSummaryHandler evHandler;
  bool result=evHandler.attachToTree((TTree *)c);
  
  if(!result) return;
  TString chPrefix("");
  if(isData)                     chPrefix="data";  
  else if(isSignal)              chPrefix="signal";
  else if(isDY)                  chPrefix="dy"; 
  else if(isVV)                  chPrefix="vv"; 
  else if(hasTop && !isSignal)   chPrefix="stop";
  for(unsigned int i=0; i<iEntries; i++)
    {
      c->GetEntry(i,1);

      //decode the event
      ZZ2l2nuSummary_t &ev = evHandler.getEvent();
      float weight(isData ? 1.0 : ev.hptWeights[0]*ev.hptWeights[1]);  //for MC (xsec/NeventsGen)*puWeight
      
      TString ch("");
      if(ev.cat==MUMU)     ch="mumu"; 
      else if(ev.cat==EE)  ch="ee"; 
      else if(ev.cat==EMU) ch="emu";
      else                 continue;
      PhysicsEvent_t phys = getPhysicsEventFrom(ev);

      //leptons
      LorentzVector dil=phys.leptons[0]+phys.leptons[1];
      bool isZcand(ev.cat!=EMU && fabs(dil.mass()-91)<15);
      bool passMet(ev.cat==EMU || ((ev.cat==EE||ev.cat==MUMU) && phys.met[0].pt()>sfMetCut));
      PhysicsObjectLeptonCollection iRotLeptons = randomlyRotate(phys.leptons,phys.ajets,calibRndGen);

      //jets
      int nJPLtags(0);
      std::vector<PhysicsObjectJetCollection> jetsVar;
      LorentzVectorCollection metsVar;      
      METUtils::computeVariation(phys.ajets, phys.leptons, phys.met[0], jetsVar, metsVar, jetCorUnc);
      for(size_t iSystVar=0; iSystVar<jetsVar.size(); iSystVar++)
	{
	  TString systVar("");
	  if     (iSystVar==METUtils::JER_UP && !isData)    systVar="jerup";
	  else if(iSystVar==METUtils::JER_DOWN  && !isData) systVar="jerdown";
	  else if(iSystVar==METUtils::JES_UP && !isData)    systVar="jesup";
	  else if(iSystVar==METUtils::JES_DOWN && !isData)  systVar="jesdown";
	  else if(iSystVar!=METUtils::JER)                  continue;
	  
	  //analyze jets
	  int nGoodJets(0);
	  std::map<TString,std::vector<float> >mljs;
	  for(size_t ijet=0; ijet<jetsVar[iSystVar].size(); ijet++)
	    {
	      if(jetsVar[iSystVar][ijet].pt()<30 || fabs(jetsVar[iSystVar][ijet].eta())>2.5) continue;
	      nGoodJets++;
	      if(iSystVar==0) nJPLtags += (jetsVar[iSystVar][ijet].btag3>0.275);
	      int partonMatch=jetsVar[iSystVar][ijet].genid;
	      int flavorMatch=jetsVar[iSystVar][ijet].flavid;	      
	      
	      for(size_t ilep=0; ilep<2; ilep++)
		{
		  LorentzVector lj=phys.leptons[ilep]+jetsVar[iSystVar][ijet];
		  float imlj=lj.mass();

		  mljs[""].push_back( imlj );
		  
		  //rotated leptons model
		  if(iSystVar==0)
		    {
		      LorentzVector rotLJ=iRotLeptons[ilep]+jetsVar[iSystVar][ijet];
		      mljs["rot"].push_back(rotLJ.mass());
		    }
		  
		  //MC only
		  if(isData) continue;
		  
		  //check if the assignment is correct
		  int assignCode=(phys.leptons[ilep].genid*partonMatch);
		  bool isCorrect(assignCode<0  && fabs(flavorMatch)==5 && (isSignal || hasTop));
		  bool isFlipMatch(assignCode>0  && fabs(flavorMatch)==5 && (isSignal || hasTop));
		  if(isCorrect)       mljs["correct"].push_back(imlj);
		  else                mljs["misassigned"].push_back(imlj);
		  if(isFlipMatch)     mljs["wrong"].push_back(imlj);
		  else if(!isCorrect) mljs["unmatched"].push_back(imlj);
		}
	    }
	  if(nGoodJets<2 || nGoodJets>maxJets) continue;
	  
	  //do the mixing for signal and data only
	  if((isSignal || isData) && passMet && !isZcand)
	    {
	      int imixtry(1);
	      LorentzVectorCollection mixjets;
	      do{
		//unsigned int j=calibRndGen.Uniform(0,iEntries);
		unsigned int j=(i+imixtry)%iEntries;
		if(j==i) { imixtry++; continue; }
		else     imixtry++;
		if(imixtry>maxMixTries) {  cout << "Failed to mix 1 event" << endl; break; }
		c->GetEntry(j);
		
		ZZ2l2nuSummary_t &mixev = evHandler.getEvent();
		//if(evcat!= mixev.cat) continue;
		PhysicsEvent_t mixphys = getPhysicsEventFrom(mixev);

		//check that it is not a Z cand
		LorentzVector mixdil=mixphys.leptons[0]+mixphys.leptons[1];
		bool isZcand((mixev.cat==EE||mixev.cat==MUMU)&&fabs(mixdil.mass()-91)<15);   
		bool passMet(mixev.cat==EMU || ((mixev.cat==EE||mixev.cat==MUMU) && mixphys.met[0].pt()>sfMetCut));
		if(!passMet || isZcand) continue;

		//analyse the jets
		std::vector<PhysicsObjectJetCollection> mixjetsVar;
		LorentzVectorCollection mixmetsVar;
		METUtils::computeVariation(mixphys.ajets, mixphys.leptons, mixphys.met[0], mixjetsVar, mixmetsVar, jetCorUnc);
		for(size_t jjet=0; jjet< mixjetsVar[0].size(); jjet++)
		  {
		    if(mixjetsVar[0][jjet].pt()<30 || fabs(mixjetsVar[0][jjet].eta())>2.5) continue;
		    double dR1 = deltaR(phys.leptons[0].eta(),phys.leptons[0].phi(),mixjetsVar[0][jjet].eta(),mixjetsVar[0][jjet].phi());
		    double dR2 = deltaR(phys.leptons[1].eta(),phys.leptons[1].phi(),mixjetsVar[0][jjet].eta(),mixjetsVar[0][jjet].phi());
		    if(dR1<0.3 || dR2<0.3) continue; 
		    if(int(mixjets.size())<nGoodJets) mixjets.push_back( mixjetsVar[0][jjet] );
		    else                              break; 
		  }
		
		//continue until jet multiplicity is filled
		if(int(mixjets.size())<nGoodJets) continue;	
		break;
	      }while(1);
	
	      //compute the mljs for the mixed jets
	      for(size_t jjet=0; jjet<mixjets.size(); jjet++)
		{
		  for(size_t ilep=0; ilep<2; ilep++)
		    {
		      LorentzVector lj=phys.leptons[ilep]+mixjets[jjet];
		      float imlj=lj.mass();
		      mljs["swap"].push_back( imlj );
		    }
		}
	    }

	  //fill control histograms for what was computed
	  TString jetCat("eq"); jetCat+=nGoodJets; jetCat +="jets";
	  std::vector<TString> ctf;
	  if(passMet && !isZcand)
	    {
	      if(!isData)
		{
		  ctf.push_back("all");
		  ctf.push_back(ch);
		  ctf.push_back(ch+jetCat);
		}
	      ctf.push_back(chPrefix);
	      ctf.push_back(chPrefix+ch);
	      ctf.push_back(chPrefix+jetCat);
	      ctf.push_back(chPrefix+ch+jetCat);
	    }
	  else if(isZcand && nJPLtags==0)
	    {	   
	      ctf.push_back(chPrefix+ch+"z");
	      ctf.push_back(chPrefix+ch+"z"+jetCat);
	    }
	  for(std::map<TString,std::vector<float> >::iterator mIt=mljs.begin(); mIt!=mljs.end(); mIt++)
	    {
	      TString hname(mIt->first); hname+="mlj"; hname+=systVar;
	      for(size_t imlj=0; imlj<mIt->second.size(); imlj++)
		controlHistos.fillHisto(hname,ctf,(mIt->second)[imlj],weight);
	    }
	}
    }
}


//
int main(int argc, char* argv[])
{
  // load framework libraries                                                                                          
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  //configure
  TString url(""), json("");
  for(int i=1;i<argc;i++)
    {
      string arg(argv[i]);
      if(arg.find("--help")!=string::npos)              { printHelp();    return 0; }
      if(arg.find("--in")!=string::npos && i+1<argc)    { url=argv[i+1];   gSystem->ExpandPathName(url); i++;  printf("in      = %s\n", url.Data()); }
      if(arg.find("--json")!=string::npos && i+1<argc)  { json=argv[i+1];                                i++;  printf("json    = %s\n", json.Data()); }
      if(arg.find("--iLumi")!=string::npos && i+1<argc) { sscanf(argv[i+1],"%f",&dataLumi);              i++;  printf("lumi    = %f\n", dataLumi); }
    }
  if(url=="" || json=="") { printHelp(); return 0;}

  //readout the data
  JSONWrapper::Object jsonF(json.Data(), true);
  fillEventsChain(url, jsonF);

  //fill the templated distributions for Mlj
  buildMljTemplates();


  cout << report.str() << endl;

  //save everything
  TFile *outF=TFile::Open("MljAnalysisReport.root","recreate");
  outF->cd();
  controlHistos.Write();
  outF->Close();

  //calibrate method
  /*
  typedef std::pair<Double_t,Int_t> PseudoExperimentKey_t;
  std::vector<PseudoExperimentKey_t > lumiScenarios;
  //   lumiScenarios.push_back(PseudoExperimentKey_t(36,100));
  //   lumiScenarios.push_back(PseudoExperimentKey_t(50,100));
  //   lumiScenarios.push_back(PseudoExperimentKey_t(100,100));
  //   lumiScenarios.push_back(PseudoExperimentKey_t(200,100));
  //   lumiScenarios.push_back(PseudoExperimentKey_t(400,50));
  //   lumiScenarios.push_back(PseudoExperimentKey_t(600,25));
  //   lumiScenarios.push_back(PseudoExperimentKey_t(1000,10));
  lumiScenarios.push_back(PseudoExperimentKey_t(dataLumi,5));
  for(size_t ijet=0; ijet<njbins; ijet++)
    {
      //init the misassignment measurement handler from the mlj spectrum
      exclusiveMisMeasurements[jetBins[ijet]] = new MisassignmentMeasurement("${CMSSW_BASE}/src/CMGTools/HtoZZ2l2nu/data/GR_R_42_V23_AK5PFchs_Uncertainty.txt");
      for(int icat=0; icat<ncats; icat++)  misMeasurement->setBiasCorrections(cats[icat],0.0);

      for(size_t ipe=0; ipe<lumiScenarios.size(); ipe++)
	{
	  PseudoExperimentKey_t pe=lumiScenarios[ipe];
	  bool freezeResults=(pe.first==200);
	  runCalibration(url,pe.first,pe.second,jetBins[ijet],freezeResults );
	}
    }
  */

	  
// 	  TString iname(cats[icat]+jetBins[ijet]);
// 	  TGraphErrors *gr = new TGraphErrors;
// 	  gr->SetName(iname);
// 	  std::pair<int,TString> key(jetBins[ijet],cats[icat]);
// 	  calibCurves[key]=gr;
// 	}


}


/*


  typedef std::pair<Double_t,Int_t> PseudoExperimentKey_t;
  std::vector<PseudoExperimentKey_t > lumiScenarios;
  if(syst=="")// || syst=="metcut") 
    {
      lumiScenarios.push_back(PseudoExperimentKey_t(36,100));
      lumiScenarios.push_back(PseudoExperimentKey_t(50,100));
      lumiScenarios.push_back(PseudoExperimentKey_t(100,100));
      lumiScenarios.push_back(PseudoExperimentKey_t(200,100));
      lumiScenarios.push_back(PseudoExperimentKey_t(400,50));
      lumiScenarios.push_back(PseudoExperimentKey_t(600,25));
      lumiScenarios.push_back(PseudoExperimentKey_t(1000,10));
      lumiScenarios.push_back(PseudoExperimentKey_t(dataLumi,5));
    }
  else
    {
      lumiScenarios.push_back(PseudoExperimentKey_t(200,50));
    }



  //display the results
  setStyle();
  gStyle->SetOptFit(0);
  TCanvas *biasevolc = new TCanvas("biasevolc","biasevolc",true);
  biasevolc->SetCanvasSize(1200,1200);
  biasevolc->SetWindowSize(1200,1200);
  biasevolc->Divide(2,2);

  report << "Jet bin &" << flush;
  for(int icat=0; icat<4; icat++){
    report << catTits[icat] << flush;
    if(icat<3) report << " & " << flush; 
  };
  report << "\\\\\\hline" << endl;

  for(size_t ijet=0; ijet<3; ijet++)
    {
      report << "$" << jetBinTitles[ijet] << "$ &" << flush;
      for(int icat=0; icat<4; icat++)
	{
	  TPad *p=(TPad *) biasevolc->cd(icat+1);
	  std::pair<int,TString> key(jetBins[ijet],cats[icat]);
	  TGraphErrors *gr = calibCurves[key];
	  gr->SetTitle(jetBinTitles[ijet]);
	  gr->SetMarkerStyle(20+ijet*2);
	  gr->SetFillStyle(0);
	  gr->SetFillColor(0);
	  gr->Draw(ijet==0? "ap" : "p");
	  if(ijet==0)
	    {
	      gr->GetXaxis()->SetTitle("Integrated luminosity (pb^{-1})");
	      gr->GetYaxis()->SetTitle("Average bias");
	      gr->GetYaxis()->SetRangeUser(-0.25,0.25);
	      p->SetLogx(true);

	      TPaveText *pave = new TPaveText(0.20,0.8,0.35,0.9,"NDC");
	      pave->SetBorderSize(0);
	      pave->SetFillStyle(0);
	      pave->SetTextFont(42);
	      pave->AddText(catTits[icat]);
	      pave->Draw();
	    }
	  if(icat==0 && ijet==2)
	    {
	      TLegend *leg = p->BuildLegend();
	      formatForCmsPublic(p,leg,"CMS simulation",3);
	    }

	  //fit the bias using high lumi pseudo-experiments
	  gr->Fit("pol0","FMQR0+","",100,dataLumi);
	  Float_t avgBias=gr->GetFunction("pol0")->GetParameter(0);
	  Float_t avgBiasError=gr->GetFunction("pol0")->GetParError(0);
	  char buf[100];
	  sprintf(buf,"%3.3f $\\pm$ %3.3f", avgBias, avgBiasError);
	  report << buf;
	  if(icat<3) report << " & " << flush; 
	}

      report << "\\\\" << endl;
    }
  biasevolc->Modified();
  biasevolc->Update();
  biasevolc->SaveAs("MisassignmentBiasEvol.C");
  biasevolc->SaveAs("MisassignmentBiasEvol.png");
  biasevolc->SaveAs("MisassignmentBiasEvol.pdf");

  report << endl;
  for(size_t ijet=0; ijet<3; ijet++)
    {
      report << "$" << jetBinTitles[ijet] << "$ &" << flush;
      for(int icat=0; icat<4; icat++)
	{
	  report << catTits[icat] << flush;
	  if(icat<3) report << " & " << flush; 
	}
      report << "\\\\\\hline" << endl;


      report << "$f_{correct}^{MC}(signal) $ &" << flush;
      for(int icat=0; icat<4; icat++)
	{
	  std::pair<int,TString> key(jetBins[ijet],cats[icat]);
	  char buf[100];
	  sprintf(buf,"%3.3f $\\pm$ %3.3f", fCorrectSigOnlyTrue[key][0],fCorrectSigOnlyTrue[key][1]);
	  report << buf;
	  if(icat<3) report << " & " << flush; 	  
	}
      report << "\\\\" << endl;

      report << "$f_{correct}^{MC}$ &" << flush;
      for(int icat=0; icat<4; icat++)
	{
	  std::pair<int,TString> key(jetBins[ijet],cats[icat]);
	  char buf[100];
	  sprintf(buf,"%3.3f $\\pm$ %3.3f", fCorrectTrue[key][0],fCorrectTrue[key][1]);
	  report << buf;
	  if(icat<3) report << " & " << flush; 	  
	}
      report << "\\\\" << endl;

      report << "$f_{correct}^{data}$ &" << flush;
      for(int icat=0; icat<4; icat++)
	{
	  std::pair<int,TString> key(jetBins[ijet],cats[icat]);
	  char buf[100];
	  sprintf(buf,"%3.3f $\\pm$ %3.3f", fCorrectData[key][0],fCorrectData[key][1]);
	  report << buf;
	  if(icat<3) report << " & " << flush; 	  
	}
      report << "\\\\" << endl;
    }
  report << endl;
  cout << report.str() << endl;


  //runPileupCheck(url);

*/


	  
//
void runCalibration(TString url, JSONWrapper::Object &jsonF, Double_t lumi, int maxPE, Int_t jetbin, bool freezeResults)
{


  /*

  ZZ2l2nuSummaryHandler evHandler, ensembleHandler;
  // MisassignmentMeasurement *misMeasurement = exclusiveMisMeasurements[jetbin];
  // misMeasurement->resetHistograms(true);
  TString peTag("Jets: "); if(jetbin==JETINCLUSIVE) peTag += "inclusive"; else peTag += jetbin; peTag += " Lumi:"; peTag += int(lumi);
  for(int iPE=1; iPE<=maxPE; iPE++)
    {
      if(iPE%5==0) { printf("\r[ %s ]  completed %d/%d ",peTag.Data(),iPE,maxPE); cout << flush; }


	  

  //
  // map the samples per type
  //
  //signal
  
  TString syst("");
  if(syst=="matchingdown")                    procYields["TTJets_matchingdown"]=0;
  else if(syst=="matchingup")                 procYields["TTJets_matchingup"]=0;
  else if(syst=="scaledown")                  procYields["TTJets_scaledown"]=0;
  else if(syst=="scaleup")                    procYields["TTJets_scaleup"]=0;
  else if(syst=="flavor")                     procYields["TT2l_R0"]=0;
  else if(syst!="singletop" && syst!="notop") procYields["TTJets_signal"]=0;
  //other sources of t->Wb
  if(syst!="notop")
    {
      procYields["SingleTbar_tW"]=0;
      procYields["SingleT_tW"]=0;
      procYields["SingleTbar_t"]=0;
      procYields["SingleT_t"]=0;
      procYields["SingleTbar_s"]=0;
      procYields["SingleT_s"]=0;
      procYields["TTJets"]=0;
    }
  //other backgrounds
  procYields["WW"]=0;
  procYields["ZZ"]=0;
  procYields["WZ"]=0;
  procYields["DYJetsToLL"]=0;

  //
  // map the trees with the events for the pseudo-experiments
  //
  //  cout << "[Filling trees]" << endl;
  TFile *inF = TFile::Open(url);
  for(std::map<TString,Int_t>::iterator it = procYields.begin(); it!= procYields.end(); it++)
    {
      TString tname=it->first + "/data";
      TTree *t = (TTree *) inF->Get(tname);
      if(t==0) continue;
      Float_t xsecWeight;
      t->GetBranch("xsecWeight")->SetAddress(&xsecWeight);
      t->GetEntry(0);
      procTrees[it->first]  = t;
      double sf=1.0;
      if(it->first.Contains("DY") ) sf=(0.5*1.0+0.25*1.73+0.25*1.68);
      if(it->first.Contains("DY") && syst=="dyup")   sf*=1.15;
      if(it->first.Contains("DY") && syst=="dydown") sf*=0.85;
      procYields[it->first] = sf*lumi*xsecWeight*t->GetEntriesFast();
      //      cout << "\t" << it->first <<  " " <<  procYields[it->first] << " evts" << endl; 
    }
      bool initEnsemble(true);
      for(std::map<TString,TTree*>::iterator it = procTrees.begin(); it!= procTrees.end(); it++)
	{

	  bool hasSignal(it->first.Contains("TTJets") || it->first.Contains("SingleT") );
	  TTree *t=it->second;
	  bool result = evHandler.attachToTree(t);
	  if(!result) continue;
	  
	  //clone tree structure only and release to base directory
	  //NB the addresses of the branches are connected to the original tree
	  //   when filling the tree the values in the evHandler.evSummary_ structure will be copied
	  //   for details cf. http://root.cern.ch/root/html/TTree.html#TTree:CloneTree
	  if(initEnsemble)
	    {
	      ensembleHandler.initTree(t->CloneTree(0),false);
	      ensembleHandler.getTree()->SetDirectory(0);
	      initEnsemble=false;
	    }
	  
	  int nevts   = evHandler.getEntries();
	  if(nevts==0) continue;

	  //generate randomly according to expectations
	  float nevExpected = procYields[it->first];
	  Int_t nevToGenerate=gRandom->Poisson(nevExpected);
	  std::set<int> evtList;
	  for(int iev=0; iev<nevToGenerate; iev++)
	    {
	      int ientry=gRandom->Uniform(0,nevts);
	      if(evtList.find(ientry)!=evtList.end()) continue;
	      evtList.insert(ientry);
	      evHandler.getEntry(ientry);
	      evHandler.evSummary_.nmeasurements=1;
	      evHandler.evSummary_.evmeasurements[0]=hasSignal;
	      ensembleHandler.fillTree();
	    }	  
	}
      
      //take control of the filled ensemble and run the misassignment measurement
      ensembleHandler.attachToTree( ensembleHandler.getTree() );
      misMeasurement->measureMisassignments( ensembleHandler, 180, 40, false, jetbin, syst );

      //reset tree used ensemble
      ensembleHandler.getTree()->Delete("all");
    }
  

  //
  // apply a similar procedure to data
  //
  TString data[]={"DoubleMuMay10ReReco",  "DoubleElectronMay10ReReco",  "MuEGMay10ReReco",
		  "DoubleMuPromptRecov4", "DoubleElectronPromptRecov4", "MuEGPromptRecov4",
		  "DoubleMu05AugReReco",  "DoubleElectron05AugReReco",  "MuEG05AugReReco",
		  "DoubleMuPromptRecov6", "DoubleElectronPromptRecov6", "MuEGPromptRecov6"};
  const size_t ndata=sizeof(data)/sizeof(TString);
  bool initEnsemble(true);
  for(size_t idata=0; idata<ndata; idata++)
    {
      
      TTree *t=(TTree *) inF->Get(data[idata]+"/data");
      bool result = evHandler.attachToTree(t);
      if(!result) continue;
      if(initEnsemble)
	{
	  ensembleHandler.initTree(t->CloneTree(0),false);
	  ensembleHandler.getTree()->SetDirectory(0);
	  initEnsemble=false;
	}
      int nevts   = evHandler.getEntries();
      if(nevts==0) continue;
      for(int iev=0; iev<nevts; iev++)
	{
	  evHandler.getEntry(iev);
	  ensembleHandler.fillTree();
	}	  
    }
  ensembleHandler.attachToTree( ensembleHandler.getTree() );
  TString datamode(syst=="metcut"?syst:"");
  misMeasurement->measureMisassignments( ensembleHandler, 180, 40, true, jetbin, datamode );
  ensembleHandler.getTree()->Delete("all");
  misMeasurement->finishMonitoringHistograms();

  //save the bias (MC truth and data are retrieved if results are to be frozen)
  for(int icat=0; icat<4; icat++)
    {
      std::pair<int,TString> key(jetbin,cats[icat]);
      if(freezeResults)
	{
	  fCorrectTrue[key] = misMeasurement->getTrueCorrectPairsFraction(cats[icat]);
	  fCorrectData[key] = misMeasurement->getCorrectPairsFraction(cats[icat]);
	  fCorrectSigOnlyTrue[key] = misMeasurement->getTrueCorrectPairsFraction(cats[icat],true);
	}

      //if the number of pseudo-experiments is large get the bias estimate from the gaussian fit
      TH1 *biasH = misMeasurement->controlHistos.getHisto("bias",cats[icat]);
      Float_t fCorrBias=biasH->GetMean();
      Float_t fCorrBiasErr=biasH->GetMeanError();
      if(maxPE>20)
	{
	  TF1 *biasFit=biasH->GetFunction("gaus");
	  fCorrBias=biasFit->GetParameter(1);
	  fCorrBiasErr=biasFit->GetParError(1);
	}
      
      TGraphErrors *gr = calibCurves[key];
      Int_t ipt=gr->GetN();
      gr->SetPoint(ipt,lumi,fCorrBias);
      gr->SetPointError(ipt,0,fCorrBiasErr);



    }



  if(!freezeResults) return;
      
  //display the results
  char titBuf[250];
  sprintf(titBuf,"CMS preliminary, #sqrt{s}=7 TeV, #int L=%3.1f fb^{-1}",float(dataLumi/1000));
  setStyle();
  gStyle->SetOptFit(0);
  TCanvas *mljc = new TCanvas("mljc","mljc",true);
  mljc->SetCanvasSize(1200,1200);
  mljc->SetWindowSize(1200,1200);
  mljc->Divide(2,2);
  TCanvas *misc = new TCanvas("misc","misc",true);
  misc->SetCanvasSize(1200,1200);
  misc->SetWindowSize(1200,1200);
  misc->Divide(2,2);
  TCanvas *misresc = new TCanvas("misresc","misc",true);
  misresc->SetCanvasSize(1200,1200);
  misresc->SetWindowSize(1200,1200);
  misresc->Divide(2,2);
  TCanvas *biasc = new TCanvas("biasc","biasc",true);
  biasc->SetCanvasSize(1200,1200);
  biasc->SetWindowSize(1200,1200);
  biasc->Divide(2,2);
  TCanvas *pullc = new TCanvas("pullc","pullc",true);
  pullc->SetCanvasSize(1200,1200);
  pullc->SetWindowSize(1200,1200);
  pullc->Divide(2,2);

  for(int icat=0; icat<4; icat++)
    {
      //get the histograms
      TH1 *mcWrongH = misMeasurement->controlHistos.getHisto("avgwrongmlj",cats[icat]);
      formatPlot(mcWrongH,809,1,1,0,1001,true,false,1,809,809);
      mcWrongH->SetTitle("Wrong assignments");
      mcWrongH->Scale(dataLumi/lumi);

      TH1 *mcCorrectH= misMeasurement->controlHistos.getHisto("avgcorrectmlj",cats[icat]);
      formatPlot(mcCorrectH,614,1,1,0,1001,true,false,1,614,614);
      mcCorrectH->SetTitle("Correct assignments");
      mcCorrectH->Scale(dataLumi/lumi);

      TH1 *dataH=misMeasurement->controlHistos.getHisto("datainclusivemlj",cats[icat]);
      dataH->SetTitle("data");

      TH1 *mcmodelH=misMeasurement->controlHistos.getHisto("avgwrongmodelmlj",cats[icat]);
      mcmodelH->SetTitle("MC model");
      formatPlot(mcmodelH,8,9,2,1,0,true,false,8,8,8);
      mcmodelH->Scale(dataLumi/lumi);
      
      TH1 *datamodelH = misMeasurement->controlHistos.getHisto("datawrongmodelmlj",cats[icat]);
      datamodelH->SetTitle("data model");
      formatPlot(datamodelH,1,9,2,1,0,true,false,1,1,1);
      
      TH1 *mcSubtractedH = (TH1 *)mcCorrectH->Clone("mcres"+cats[icat]);
      mcSubtractedH->Add(mcWrongH);
      mcSubtractedH->Add(mcmodelH,-1);
      formatPlot(mcSubtractedH,8,9,2,1,0,true,false,8,8,8);
      mcSubtractedH->SetTitle("MC subtracted");

      TH1 *dataSubtractedH = (TH1 *)dataH->Clone("datares"+cats[icat]);
      dataSubtractedH->Add(datamodelH,-1);
      dataSubtractedH->SetTitle("data subtracted");

      TH1 *biasH = misMeasurement->controlHistos.getHisto("bias",cats[icat]);
      formatPlot(biasH,1,1,2,20,0,true,false,1,1,1);
      
      TH1 *pullH = misMeasurement->controlHistos.getHisto("pull",cats[icat]);
      formatPlot(pullH,1,1,2,20,0,true,false,1,1,1);
      
      TList nullList;
      TList dataList;   dataList.Add(dataH);
      TList mcList;     mcList.Add(mcWrongH);       mcList.Add(mcCorrectH);    
      TList modelList;  modelList.Add(mcmodelH);    modelList.Add(datamodelH);
      
      TPad *p=(TPad *) mljc->cd(icat+1);
      TLegend *leg = showPlotsAndMCtoDataComparison(p,mcList,nullList,dataList);     
      TPad *sp=(TPad *)p->cd(1);
      if(icat==0)  formatForCmsPublic(sp,leg,titBuf,3);
      else         leg->Delete();
      sp->SetLogx(true);
      sp->SetLogy(false);
      TPaveText *pave = new TPaveText(0.20,0.8,0.35,0.9,"NDC");
      pave->SetBorderSize(0);
      pave->SetFillStyle(0);
      pave->SetTextFont(42);
      pave->AddText(catTits[icat]);
      pave->Draw();
      sp=(TPad *)p->cd(2);
      sp->SetLogx(true);
      
      p=(TPad *)misc->cd(icat+1);
      leg = showPlotsAndMCtoDataComparison( p,mcList,modelList,dataList);     
      sp=(TPad *)p->cd(1);
      if(icat==0)  formatForCmsPublic(sp,leg,titBuf,4);
      else         leg->Delete();
      sp->SetLogx(true);
      sp->SetLogy(false);
      pave->DrawClone();
      sp=(TPad *)p->cd(2);
      sp->SetLogx(true);

      p=(TPad *)misresc->cd(icat+1);
      TList mcCorList;    mcCorList.Add(mcCorrectH); 
      TList mcresList;    mcresList.Add(mcSubtractedH);
      TList dataresList;  dataresList.Add(dataSubtractedH);
      leg = showPlotsAndMCtoDataComparison(p,mcCorList,mcresList,dataresList);     
      sp=(TPad *)p->cd(1);
      if(icat==0)  formatForCmsPublic(sp,leg,titBuf,2);
      else         leg->Delete();
      sp->SetLogx(true);
      sp->SetLogy(false);
      pave->DrawClone();
      sp=(TPad *)p->cd(2);
      sp->SetLogx(true);

      p=(TPad *)biasc->cd(icat+1);
      biasH->Draw("e1");
      leg=p->BuildLegend();
      if(icat==0)  formatForCmsPublic(p,leg,"CMS simulation",2);
      leg->Delete();
      pave->DrawClone();

      p=(TPad *)pullc->cd(icat+1);
      pullH->Draw("e1");
      leg=p->BuildLegend();
      if(icat==0)  formatForCmsPublic(p,leg,"CMS simulation",2);
      leg->Delete();
      pave->DrawClone();
    }

  TString postfix("");
  if(jetbin>0) { postfix+="_"; postfix += jetbin; }

  mljc->Modified();
  mljc->Update();
  mljc->SaveAs("Mlj"+postfix+".C");
  mljc->SaveAs("Mlj"+postfix+".png");
  mljc->SaveAs("Mlj"+postfix+".pdf");

  misc->Modified();
  misc->Update();
  misc->SaveAs("MljWithModel"+postfix+".C");
  misc->SaveAs("MljWithModel"+postfix+".png");
  misc->SaveAs("MljWithModel"+postfix+".pdf");

  misresc->Modified();
  misresc->Update();
  misresc->SaveAs("MljSubtracted"+postfix+".C");
  misresc->SaveAs("MljSubtracted"+postfix+".png");
  misresc->SaveAs("MljSubtracted"+postfix+".pdf");

  biasc->Modified();
  biasc->Update();
  biasc->SaveAs("MisassignmentBias"+postfix+".C");
  biasc->SaveAs("MisassignmentBias"+postfix+".png");
  biasc->SaveAs("MisassignmentBias"+postfix+".pdf");

  pullc->Modified();
  pullc->Update();
  pullc->SaveAs("MisassignmentPull"+postfix+".C");
  pullc->SaveAs("MisassignmentPull"+postfix+".png");
  pullc->SaveAs("MisassignmentPull"+postfix+".pdf");

  inF->Close();    
	      */
}



//
void runPileupCheck(TString url)
{
  /*
  SelectionMonitor controlHistos;
  
  //generate fake data with different pileup average
  double massAxis[]={10,20,30,40,50,60,70,80,90,100,115,130,145,160,185,200,250,300,400,500,1000,2000};
  int nMassBins=sizeof(massAxis)/sizeof(double)-1;
  for(int i=0; i<12; i++)
    {
      TString name("gen"); name+=i;
      controlHistos.addHistogram( new TH1D(name+"pu",name+" PU",36,-0.5,35.5) );
      TH1 *h=controlHistos.getHisto(name+"pu");
      for(int ibin=1; ibin<=h->GetXaxis()->GetNbins(); ibin++)
	{
	  Float_t binCenter=h->GetXaxis()->GetBinCenter(ibin);
	  h->SetBinContent(ibin,TMath::Poisson(binCenter,ibin));
	}
      controlHistos.addHistogram( new TH1D(name+"mlj","Lepton-jet spectrum;Invariant Mass(l,j) [GeV/c^{2}];Lepton-jet pairs",nMassBins,massAxis) );
    }

  //now build the mlj spectrum
  TFile *inF = TFile::Open(url);
  for(std::map<TString,Int_t>::iterator it = procYields.begin(); it!= procYields.end(); it++)
    {
      TString tname=it->first + "/data";
      TTree *t = (TTree *) inF->Get(tname);
      if(t==0) continue;

      top::EventSummaryHandler evSummaryHandler;
      if( !evSummaryHandler.attachToTree( t) ) continue;
      const Int_t totalEntries=evSummaryHandler.getEntries();

      //1. fill the PU in MC
      controlHistos.addHistogram( new TH1D(it->first + "pumc",it->first , 36,-0.5,35.5) );
      TH1 *MCH=controlHistos.getHisto(it->first+"pumc");
      float procWeight=dataLumi*totalEntries;
      for(int ientry=0; ientry<totalEntries; ientry++)
	{
	  evSummaryHandler.getEntry(ientry);
	  top::EventSummary_t &ev = evSummaryHandler.getEvent();
	  procWeight *= ev.xsecWeight;
	  MCH->Fill(ev.ngenpu);
	}
      
      //2. instantiate 12 reweighters
      edm::LumiReWeighting *LumiWeights[12];
      for(int i=0; i<12; i++)
	{
	  std::vector< float > MC_distr,Lumi_distr;
	  TString name("gen"); name+=i;
	  TH1 *LumiH=controlHistos.getHisto(name+"pu");
	  for(int ibin=1; ibin<=MCH->GetXaxis()->GetNbins(); ibin++)
	    {
	      MC_distr.push_back( MCH->GetBinContent(i) );
	      Lumi_distr.push_back( LumiH->GetBinContent(i) );
	    }
	  LumiWeights[i] = new edm::LumiReWeighting(MC_distr,Lumi_distr);
	}


      //3. construct the weighted spectra
      std::vector<float> weight(12,0);
      for(int ientry=0; ientry<totalEntries; ientry++)
	{
	  evSummaryHandler.getEntry(ientry);
	  top::EventSummary_t &ev = evSummaryHandler.getEvent();
	  top::PhysicsEvent_t phys = getPhysicsEventFrom(ev);
	  
	  //get the weights
	  for(int i=0; i<12; i++)  weight[i]=LumiWeights[i]->weight( ev.ngenpu )*procWeight;

	  //inclusive mlj spectrum
	  for(top::PhysicsObjectJetCollection::iterator it=phys.jets.begin(); it!=phys.jets.end(); it++)
	    {
              if(it->pt()<=30 || fabs(it->eta())>=2.5) continue;
             
	      for(int i=0; i<12; i++) 
		{
		  TString name("gen"); name+=i;
		  for(int ilep=0; ilep<2; ilep++)
		    {
		      LorentzVector lj=phys.leptons[ilep]+*it;
		      float mlj=lj.mass();
		      controlHistos.fillHisto(name+"mlj","all",mlj,weight[i],true);
		    }
		}
	    }
	  
	}

      //delete the lumi reweighters
      for(int i=0; i<12; i++) delete LumiWeights[i];
    }
  
  //now display the results
  TCanvas *mljc = new TCanvas("mljc","mljc",true);
  mljc->SetCanvasSize(600,600);
  mljc->SetWindowSize(600,600);
  TGraphErrors *mljavg = new TGraphErrors; mljavg->SetMarkerStyle(20);
  TGraphErrors *mljrms = new TGraphErrors; mljrms->SetMarkerStyle(20);
  for(int i=0; i<12; i++) 
    {
      TString name("gen"); name+=i;
      TH1 *h=controlHistos.getHisto(name+"mlj","all");
      h->SetLineColor(kGreen-9+i);
      h->SetMarkerColor(kGreen-9+i);
      h->SetFillStyle(0);
      h->SetMarkerStyle(1);
      h->Draw(i==0?"hist":"histsame");

      mljavg->SetPoint(i,i,h->GetMean());
      mljavg->SetPointError(i,i,h->GetMeanError());

      mljrms->SetPoint(i,i,h->GetRMS());
      mljrms->SetPointError(i,i,h->GetRMSError());
    }
  TLegend *leg=mljc->BuildLegend();
  formatForCmsPublic(mljc,leg,"CMS simulation",10);
  mljc->Modified();
  mljc->Update();
  mljc->SaveAs("MljPU.C");
  mljc->SaveAs("MljPU.png");
  mljc->SaveAs("MljPU.pdf");

  TCanvas *mljmomc = new TCanvas("mljmomc","mljmomc",true);
  mljmomc->SetCanvasSize(600,600);
  mljmomc->SetWindowSize(600,600);
  mljmomc->Divide(1,2);
  mljmomc->cd(1);
  mljavg->Draw("ap");
  mljavg->GetXaxis()->SetTitle("Average pileup");
  mljavg->GetYaxis()->SetTitle("Average M_{lj}");
  mljmomc->cd(2);
  mljrms->Draw("ap");
  mljrms->GetXaxis()->SetTitle("Average pileup");
  mljrms->GetYaxis()->SetTitle("RMS M_{lj}");
  mljmomc->Modified();
  mljmomc->Update();
  mljmomc->SaveAs("MljPUMomenta.C");
  mljmomc->SaveAs("MljPUMomenta.png");
  mljmomc->SaveAs("MljPUMomenta.pdf");

  */
}
