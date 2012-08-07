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
#include "CMGTools/HtoZZ2l2nu/interface/plotter.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/MisassignmentMeasurement.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include <vector>

using namespace std;
 
Double_t dataLumi=2165.0;
std::map<TString,Int_t> procYields;
std::map<TString,TTree *> procTrees,dataTrees;
std::map<std::pair<int,TString>, TGraphErrors*> calibCurves;
std::map<int, MisassignmentMeasurement *> exclusiveMisMeasurements;
TString cats[]    = {"ee", "mumu",   "emu",  "all"};
TString catTits[] = {"ee", "#mu#mu", "e#mu", "inclusive"};
enum JetBinCodes       {JETINCLUSIVE=0, JETS2=2,   JETS3=3};
TString jetBinTitles[]={"inclusive",    "=2 jets", "=3 jets"};
stringstream report;
std::map<std::pair<int, TString>, std::vector<double> > fCorrectTrue, fCorrectSigOnlyTrue, fCorrectData;

//
void runCalibration(TString url,Double_t lumi=50, int maxPE=100,Int_t jetbin=JETINCLUSIVE, TString syst="", bool freezeResults=false);
void runPileupCheck(TString url);


//                                                                                                                                                                                          
void printHelp()
{
  printf("--help    --> print this\n");
  printf("--in      --> input file with the summary trees\n");
  printf("--syst    --> systematic to be evaluated (dyup/down, jesup/down, jer, metcut, scaleup/down, matchingup/down\n");
}

//
int main(int argc, char* argv[])
{
  // load framework libraries                                                                                          
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  TString url(""), syst("");
  for(int i=1;i<argc;i++)
    {
      string arg(argv[i]);
      if(arg.find("--help")!=string::npos)             { printHelp();    return 0; }
      if(arg.find("--in")!=string::npos && i+1<argc)   { url=argv[i+1];         gSystem->ExpandPathName(url);         i++;  printf("in      = %s\n", url.Data()); }
      if(arg.find("--syst")!=string::npos && i+1<argc) { syst=argv[i+1];                                            i++;  printf("syst    = %s\n", syst.Data()); }
    }
  if(url=="") { printHelp(); return 0;}
  
  int jetBins[]={JETINCLUSIVE,JETS2,JETS3};
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


  for(size_t ijet=0; ijet<3; ijet++)
    {
      //init the misassignment measurement handler from the mlj spectrum
      MisassignmentMeasurement *misMeasurement = new MisassignmentMeasurement;
      misMeasurement->initJetUncertainties("${CMSSW_BASE}/src/CMGTools/HtoZZ2l2nu/data/GR_R_42_V20_AK5PFchs_Uncertainty.txt",
					   "${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_PtResolution_AK5PF.txt",
					   "${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_EtaResolution_AK5PF.txt",
					   "${CMSSW_RELEASE_BASE}/src/CondFormats/JetMETObjects/data/Spring10_PhiResolution_AK5PF.txt");
      for(int icat=0; icat<4; icat++) 
	{
	  misMeasurement->setBiasCorrections(cats[icat],0.0);
	  TGraphErrors *gr = new TGraphErrors;
	  TString iname(cats[icat]); iname+=jetBins[ijet];
	  gr->SetName(iname);
	  std::pair<int,TString> key(jetBins[ijet],cats[icat]);
	  calibCurves[key]=gr;
	}
      exclusiveMisMeasurements[jetBins[ijet]]=misMeasurement;

      for(size_t ipe=0; ipe<lumiScenarios.size(); ipe++)
	{
	  PseudoExperimentKey_t pe=lumiScenarios[ipe];
	  bool freezeResults=(pe.first==200);//ipe==lumiScenarios.size()-1);
	  runCalibration(url,pe.first,pe.second,jetBins[ijet],syst, freezeResults );
	}
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
}
	  
//
void runCalibration(TString url,Double_t lumi, int maxPE,Int_t jetbin, TString syst, bool freezeResults)
{
  //
  // map the samples per type
  //
  //signal
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

  //
  // run calibration
  //
  top::EventSummaryHandler evHandler, ensembleHandler;
  MisassignmentMeasurement *misMeasurement = exclusiveMisMeasurements[jetbin];
  misMeasurement->resetHistograms(true);
  TString peTag("Jets: "); if(jetbin==JETINCLUSIVE) peTag += "inclusive"; else peTag += jetbin; peTag += " Lumi:"; peTag += int(lumi);
  for(int iPE=1; iPE<=maxPE; iPE++)
    {
      if(iPE%5==0) { printf("\r[ %s ]  completed %d/%d ",peTag.Data(),iPE,maxPE); cout << flush; }
      
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
}



//
void runPileupCheck(TString url)
{

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
}
