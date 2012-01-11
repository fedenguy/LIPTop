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
#include "LIP/Top/interface/HFCMeasurement.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"

#include <vector>

using namespace std;
using namespace RooFit;

typedef std::vector<TString> Proc_t;
typedef std::vector<Float_t> ProcYields_t;
typedef std::vector<TH1D *> ProcHistos_t;


stringstream report;
 
//
void runCalibration(TString url, int fitType, TString btagWP, std::map<TString,Double_t> &fitPars, int maxPE, Double_t dataLumi, TString syst,  bool freezeResults);

//
void printHelp()
{
  printf("--help    --> print this\n");
  printf("--iLumi   --> integrated luminosity to be used /pb  (default:2165)\n");
  printf("--in      --> input file with the summary trees\n");
  printf("--par     --> fitter configuration file\n");
  printf("--fit     --> fit type code                         (default: 0 = R\n");
  printf("--btag    --> btag algorithm                        (default: TCHEL)\n");
  printf("--syst    --> systematic uncertainty to evalute     (default:none)\n");
  printf("--npe     --> number of pseudo-experiments to throw (default:1)\n");
  printf("command line example: HeavyFlavorFitCalibration --iLumi 2165 --in EventSummary.root --par hfcFitter_cfg.txt --npe==500\n");
}

//
std::map<TString,Double_t> parseParametersFrom(TString parfileURL)
{
  std::map<TString,Double_t> fitPars;
  ifstream in;
  in.open(parfileURL);
  TString line;
  while (1) {
    in >> line;
    if (!in.good()) break;
    TObjArray *tokens = line.Tokenize(":");
    if(tokens->GetEntriesFast()<2) continue;
    TString key = ((TObjString *)tokens->At(0))->GetString();
    TString val = ((TObjString *)tokens->At(1))->GetString();
    fitPars[key]=val.Atof();
  }
  in.close();
  return fitPars;
}


//
int main(int argc, char* argv[])
{
  // load framework libraries                                                                                          
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  TString url("");
  TString fitParFile("");
  Int_t fitType(HFCMeasurement::FIT_R);
  TString btagWP("TCHEL");
  TString syst("");
  int maxPE(1);
  Double_t dataLumi(2165);
  for(int i=1;i<argc;i++)
    {
      string arg(argv[i]);
      if(arg.find("--help")!=string::npos)              { printHelp();	  return 0; }
      if(arg.find("--iLumi")!=string::npos && i+1<argc) { sscanf(argv[i+1],"%lf",&dataLumi);                          i++;  printf("lumi    = %f\n", dataLumi);          }
      if(arg.find("--in")!=string::npos && i+1<argc)    { url=argv[i+1];         gSystem->ExpandPathName(url);        i++;  printf("in      = %s\n", url.Data());        }
      if(arg.find("--par")!=string::npos && i+1<argc)   { fitParFile=argv[i+1];  gSystem->ExpandPathName(fitParFile); i++;  printf("parFile = %s\n", fitParFile.Data()); }
      if(arg.find("--syst")!=string::npos && i+1<argc)  { syst=argv[i+1];                                             i++;  printf("syst    = %s\n", syst.Data());       }
      if(arg.find("--npe")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&maxPE);                              i++;  printf("N_{PE}  = %d\n", maxPE);             }
      if(arg.find("--fit")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&fitType);                            i++;  printf("fitType  = %d\n", fitType);          }
      if(arg.find("--btag")!=string::npos && i+1<argc)  { btagWP=argv[i+1];                                           i++;  printf("btag WP  = %s\n", btagWP.Data());    }
    } 
  if(url=="" || fitParFile=="") { printHelp(); return 0; }
  

  bool freezeResults(true);
  std::map<TString,Double_t> fitPars=parseParametersFrom(fitParFile);
  runCalibration(url, fitType, btagWP, fitPars, maxPE, dataLumi, syst, freezeResults );
  
  cout << report.str() << endl;
}
	  
//
void runCalibration(TString url, int fitType, TString btagWP, std::map<TString,Double_t> &fitPars, int maxPE, Double_t dataLumi, TString syst,  bool freezeResults)
{

  TH1D *fitResH = new TH1D("fitres",";Fit result;Pseudo-experiments",100,0,2);             fitResH->SetDirectory(0);
  TH1D *biasH = new TH1D("bias",";bias;Pseudo-experiments",25,-0.52,0.48);                 biasH->SetDirectory(0);
  TH1D *pullH = new TH1D("pull",";pull = bias / #sigma;Pseudo-experiments",25,-5.2,4.8);   pullH->SetDirectory(0);
  TH1D *statUncH = new TH1D("staterr",";#sigma;Pseudo-experiments",50,-0.05,0.05);         statUncH->SetDirectory(0);

  //
  // map the samples per type
  //
  std::map<TString, Proc_t>        descForProc;
  std::map<TString, ProcYields_t > yieldsForProc;
  std::map<TString, ProcHistos_t > histosForProc;
  
  reweight::PoissonMeanShifter *puShifter=0;
  if(syst=="puup") puShifter = new reweight::PoissonMeanShifter(+0.6);
  if(syst=="pudown") puShifter = new reweight::PoissonMeanShifter(-0.6);

  //signal
  Double_t rToGen(1.0);
  ProcYields_t procYields(6,0);
  procYields[HFCMeasurement::EE_2JETS]=546.8;
  procYields[HFCMeasurement::EE_3JETS]=273.4;
  procYields[HFCMeasurement::MUMU_2JETS]=698.3;
  procYields[HFCMeasurement::MUMU_3JETS]=335.8;
  procYields[HFCMeasurement::EMU_2JETS]=1882.3;
  procYields[HFCMeasurement::EMU_3JETS]=915.4;
  yieldsForProc["Signal"]=procYields;

  if(syst=="matchingdown")        descForProc["Signal"] = Proc_t(1,"TTJets_matchingdown");
  else if(syst=="matchingup")     descForProc["Signal"] = Proc_t(1,"TTJets_matchingup");
  else if(syst=="scaledown")      descForProc["Signal"] = Proc_t(1,"TTJets_scaledown");
  else if(syst=="scaleup")        descForProc["Signal"] = Proc_t(1,"TTJets_scaleup");
  else if(syst.Contains("flavor")) 
    {
      sscanf(syst.Data(),"flavor%lf",&rToGen);

      descForProc["Signal1"]     = Proc_t(1,"TT2l_R1");     
      yieldsForProc["Signal1"]   = procYields;
      descForProc["Signal05a"]   = Proc_t(1,"TT2l_R05a");
      yieldsForProc["Signal05a"] = procYields;
      descForProc["Signal05b"]   = Proc_t(1,"TT2l_R05b");
      yieldsForProc["Signal05b"] = procYields;
      descForProc["Signal0"]     = Proc_t(1,"TT2l_R0");
      yieldsForProc["Signal0"]   = procYields;

      for(size_t icat=0; icat<6; icat++)
	{
	  yieldsForProc["Signal1"][icat]   *= pow(rToGen,2);
	  yieldsForProc["Signal05a"][icat] *= rToGen*(1-rToGen);
	  yieldsForProc["Signal05b"][icat] *= rToGen*(1-rToGen);
	  yieldsForProc["Signal0"][icat]   *= pow(1-rToGen,2);
	}
    }
  else                     descForProc["Signal"] = Proc_t(1,"TTJets_signal");

  //backgrounds
  Proc_t SingleTop;
  SingleTop.push_back("SingleTbar_tW");
  SingleTop.push_back("SingleT_tW");
  SingleTop.push_back("SingleTbar_t");
  SingleTop.push_back("SingleT_t");
  SingleTop.push_back("SingleTbar_s");
  SingleTop.push_back("SingleT_s");
  procYields[HFCMeasurement::EE_2JETS]=33.5;
  procYields[HFCMeasurement::EE_3JETS]=10.7;
  procYields[HFCMeasurement::MUMU_2JETS]=43.2;
  procYields[HFCMeasurement::MUMU_3JETS]=12.6;
  procYields[HFCMeasurement::EMU_2JETS]=111.3;
  procYields[HFCMeasurement::EMU_3JETS]=34.7;
  descForProc["SingleTop"]=SingleTop;
  yieldsForProc["SingleTop"]=procYields;

  Proc_t OtherTTbar;
  OtherTTbar.push_back("TTJets");
  procYields[HFCMeasurement::EE_2JETS]=4.8;
  procYields[HFCMeasurement::EE_3JETS]=3.2;
  procYields[HFCMeasurement::MUMU_2JETS]=0.9;
  procYields[HFCMeasurement::MUMU_3JETS]=2.1;
  procYields[HFCMeasurement::EMU_2JETS]=7.8;
  procYields[HFCMeasurement::EMU_3JETS]=8.2;
  descForProc["OtherTTbar"]=OtherTTbar;
  yieldsForProc["OtherTTbar"]=procYields;
  
  Proc_t DiBosons;
  DiBosons.push_back("WW");
  DiBosons.push_back("ZZ");
  DiBosons.push_back("WZ");
  procYields[HFCMeasurement::EE_2JETS]=11.9;
  procYields[HFCMeasurement::EE_3JETS]=2.4;
  procYields[HFCMeasurement::MUMU_2JETS]=13.6;
  procYields[HFCMeasurement::MUMU_3JETS]=2.8;
  procYields[HFCMeasurement::EMU_2JETS]=38.5;
  procYields[HFCMeasurement::EMU_3JETS]=6.6;
  descForProc["DiBosons"]=DiBosons;
  yieldsForProc["DiBosons"]=procYields;
  
  Proc_t DY;
  DY.push_back("DYJetsToLL");
  DY.push_back("DYJetsToMuMu_M20to50");
  DY.push_back("DYJetsToHFCMeasurement::EEM20to50");
  DY.push_back("DYJetsToTauTau_M20to50");

  descForProc["DY"]=DY;
  float uncScaleFactor(1.0);
  if(syst=="dyup")   uncScaleFactor *= 1.3;
  if(syst=="dydown") uncScaleFactor *= 0.7;
  procYields[HFCMeasurement::EE_2JETS]=298.0*uncScaleFactor;
  procYields[HFCMeasurement::EE_3JETS]=84.8*uncScaleFactor;
  procYields[HFCMeasurement::MUMU_2JETS]=451.4*uncScaleFactor;
  procYields[HFCMeasurement::MUMU_3JETS]=121.7*uncScaleFactor;
  procYields[HFCMeasurement::EMU_2JETS]=157.3*uncScaleFactor;
  procYields[HFCMeasurement::EMU_3JETS]=33.3*uncScaleFactor;
  yieldsForProc["DY"]=procYields;

  //
  // fill the histograms to be used for the pseudo-experiments
  //
  int maxJets=int(fitPars["maxJets"]);
  Double_t btagWPcut(fitPars[btagWP+"_cut"]);
  TFile *inF = TFile::Open(url);
  top::EventSummaryHandler evHandler;
  cout << "[Filling histograms]" << flush;
  Double_t avgeff[2]     = {0.0, 0.0};
  Double_t avgeffNorm[2] = {0.0, 0.0};
  for(std::map<TString, Proc_t>::iterator pIt=descForProc.begin(); pIt!=descForProc.end(); pIt++)
    {
      bool isSignalRequired(pIt->first.Contains("Signal"));

      //instantiate the histograms
      ProcHistos_t procHistos;
      for(size_t icat=0; icat<procYields.size(); icat++)
	{
	  TString name(pIt->first); name+=icat;
	  TH1D *bmh=new TH1D(name,";b-tag multiplicity;Events",maxJets+1,0,maxJets+1);
	  bmh->SetDirectory(0);
	  procHistos.push_back(bmh);
	}
      
      //fill histograms from the weighted sum of each contribution to the process
      cout << endl << "\t" << pIt->first << " has " << pIt->second.size() << " sources" << endl;
      for(Proc_t::iterator it = pIt->second.begin(); it != pIt->second.end(); it++)
	{
	  TTree *t=(TTree *)inF->Get(*it+"/data");
	  if(t==0) continue;
	  bool result = evHandler.attachToTree(t);
	  if(!result) continue;
	  
	  cout << endl;
	  const int nentries=evHandler.getEntries();
	  for(int ientry=0; ientry<nentries; ientry++)
	    {
	      printf("\r\t\t [ %3.0f/100 ] completed",100*float(ientry)/float(nentries));
	      evHandler.getEntry(ientry);
	      top::EventSummary_t &ev = evHandler.getEvent();
	      if(isSignalRequired && !ev.isSignal) continue;
	      if(!isSignalRequired && ev.isSignal) continue;

	      //event weight
	      Double_t weight =ev.weight*ev.xsecWeight*dataLumi;
	      if(puShifter) weight *= puShifter->ShiftWeight( ev.ngenpu );
	      
	      top::PhysicsEvent_t phys = getPhysicsEventFrom(ev);
	      
	      //sanity checks for the selection
	      LorentzVector dilepton = phys.leptons[0]+phys.leptons[1];
	      float dilcharge=(phys.leptons[0].id/fabs(phys.leptons[0].id)) *(phys.leptons[1].id/fabs(phys.leptons[1].id));
	      if((ev.cat==EE || ev.cat==MUMU) && (fabs(dilepton.mass()-91)<15 || phys.met.pt()<30)) continue;
	      if(dilcharge>0) continue;

	      //count b-tags
	      int njets=0;
	      int nbtags=0;
	      for(unsigned int ijet=0; ijet<phys.jets.size(); ijet++)
		{
		  if(phys.jets[ijet].pt()<30 || fabs(phys.jets[ijet].eta())>2.4) continue;

		  njets++;
		  double btag(-9999.);
		  if(btagWP.Contains("TCHE") )        btag = phys.jets[ijet].btag1;
		  else if(btagWP.Contains("TCHP") )   btag = phys.jets[ijet].btag2;
		  else if(btagWP.Contains("SSVHE") )  btag = phys.jets[ijet].btag3;
		  else if(btagWP.Contains("JBP") )    btag = phys.jets[ijet].btag4;
		  else if(btagWP.Contains("JP") )     btag = phys.jets[ijet].btag5;
		  else if(btagWP.Contains("SSVHP") )  btag = phys.jets[ijet].btag6;
		  else if(btagWP.Contains("CSV") )    btag = phys.jets[ijet].btag7;
		  bool isBtagged(btag>btagWPcut);
		  bool isMatchedToB( fabs(phys.jets[ijet].flavid)==5 );
		  avgeff[isMatchedToB]     += isBtagged*weight;
		  avgeffNorm[isMatchedToB] += weight;
		  nbtags += isBtagged;
		}
	      if(njets>maxJets || njets<2) continue;
	      

	      //check event category
	      int icat=-1;
	      if(ev.cat==EE &&   njets==2)  icat=HFCMeasurement::EE_2JETS;
	      if(ev.cat==EE &&   njets==3)  icat=HFCMeasurement::EE_3JETS;
	      if(ev.cat==MUMU && njets==2)  icat=HFCMeasurement::MUMU_2JETS;
	      if(ev.cat==MUMU && njets==3)  icat=HFCMeasurement::MUMU_3JETS;
	      if(ev.cat==EMU &&  njets==2)  icat=HFCMeasurement::EMU_2JETS;
	      if(ev.cat==EMU &&  njets==3)  icat=HFCMeasurement::EMU_3JETS;
	      if(icat<0) continue;

	      //fill histogram
	      procHistos[icat]->Fill(nbtags,weight);
	    }
	}
	 
      //save histograms for sampling
      histosForProc[pIt->first]=procHistos;
    }
  inF->Close();    
  cout << endl;
  
  //
  //create the histos for the ensemble to be fit
  //
  ProcHistos_t ensembleHistos;
  for(size_t icat=0; icat<procYields.size(); icat++)
    {
      TString name("ensemble"); name+=icat;
      TH1D *bmh=new TH1D(name,";b-tag multiplicity;Events",maxJets+1,0,maxJets+1);
      bmh->SetDirectory(0);
      ensembleHistos.push_back(bmh);
    }

  //
  // configure the fitter
  //
  HFCMeasurement hfcFitter(fitType,maxJets,rToGen);
  hfcFitter.configureBtagAlgo   (btagWP, btagWPcut);
  double mceb=fitPars[btagWP+"_effb"];
  double mceq=fitPars[btagWP+"_effq"];
  if(fitType==HFCMeasurement::FIT_R)
    {
      for(size_t i=0; i<2; i++) avgeff[i] /= avgeffNorm[i];
      mceq=avgeff[0];
      mceb=avgeff[1];
    }
  hfcFitter.setBtagEfficiency   (mceb, fitPars[btagWP+"_sfb"], fitPars[btagWP+"_sfbunc"]);
  hfcFitter.setMistagEfficiency (mceq, fitPars[btagWP+"_sfq"], fitPars[btagWP+"_sfqunc"]);

  TString channels[]={"ee","mumu","emu"};
  for(int ich=0; ich<3; ich++)
    {
      for(int ijets=2; ijets<=maxJets; ijets++)
	{
	  TString key(channels[ich]+"_"); key += ijets; key += "jets";
	  hfcFitter.setSelectionFractions( fitPars[key+"_fcorrect"],   fitPars[key+"_fcorrectunc"],
					   fitPars[key+"_fttbar"],     fitPars[key+"_fttbarunc"],
					   fitPars[key+"_fsingletop"], fitPars[key+"_fsingletop"],
					   ijets, channels[ich]);
	}
    }
  hfcFitter.printConfiguration(report);

  //keep RooFit quiet                                                                                                                                                          
  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().getStream(0).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Fitting);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);

 
  //
  // run the calibration
  //
  cout << "[Running calibration]" << endl;
  for(int iPE=1; iPE<=maxPE; iPE++)
    {
      if(iPE%5==0) { printf("\rCompleted %d/%d ",iPE,maxPE); cout << flush; }
      
      //fill randomly the ensemble histos
      for(size_t ih=0; ih<ensembleHistos.size(); ih++) ensembleHistos[ih]->Reset("ICE");
      for(std::map<TString,ProcHistos_t>::iterator itt = histosForProc.begin(); itt!=histosForProc.end(); itt++)
        {
	  
	  int procEvts(0);
	  for(size_t ih=0; ih<itt->second.size(); ih++)
	    {
	      //generate randomly from the template histograms                                                                                                                            
	      int nevExpected=yieldsForProc[itt->first][ih];
	      Int_t nevToGenerate=gRandom->Poisson(nevExpected);
	      procEvts+=nevToGenerate;
	      for(int iev=0; iev<nevToGenerate; iev++) ensembleHistos[ih]->Fill( itt->second[ih]->GetRandom() );
	    }
	}

      //fit
      hfcFitter.fitHFCtoMeasurement(ensembleHistos,false);
    
      double effb(hfcFitter.model.abseb->getVal()) , effq(hfcFitter.model.abseq->getVal());
      report << "[Ensemble #" << iPE << "]" << endl; 
      report << "\t R   = " << hfcFitter.model.r->getVal() << " +" << hfcFitter.model.r->getAsymErrorHi() << " " << hfcFitter.model.r->getAsymErrorLo() << endl;
      for(int icat=0; icat<6; icat++) report << "\t\t R_{" << icat << "}  = " << hfcFitter.model.rFit[icat] << " +" << hfcFitter.model.rFitAsymmErrHi[icat] << " " << hfcFitter.model.rFitAsymmErrLo[icat] << endl;
      report << "\t e_b = " << effb*hfcFitter.model.sfeb->getVal() << " +" << effb*hfcFitter.model.sfeb->getAsymErrorHi() << " " << effb*hfcFitter.model.sfeb->getAsymErrorLo() << endl;
      for(int icat=0; icat<6; icat++) report << "\t\t eb_{" << icat << "}  = " << hfcFitter.model.ebFit[icat] << " +" << hfcFitter.model.ebFitAsymmErrHi[icat] << " " << hfcFitter.model.ebFitAsymmErrLo[icat] << endl;
      report << "\t e_q = " << effq*hfcFitter.model.sfeq->getVal() << " +" << effq*hfcFitter.model.sfeq->getAsymErrorHi() << " " << effq*hfcFitter.model.sfeq->getAsymErrorLo() << endl;
      for(int icat=0; icat<6; icat++) report << "\t\t eq_{" << icat << "}  = " << hfcFitter.model.eqFit[icat] << " +" << hfcFitter.model.eqFitAsymmErrHi[icat] << " " << hfcFitter.model.eqFitAsymmErrLo[icat] << endl;

      float valFit   ( fitType==HFCMeasurement::FIT_R ? hfcFitter.model.r->getVal()   : effb*hfcFitter.model.sfeb->getVal() );
      float valFitErr( fitType==HFCMeasurement::FIT_R ? hfcFitter.model.r->getError() : effb*hfcFitter.model.sfeb->getError() );
      float trueVal  ( fitType==HFCMeasurement::FIT_R ? rToGen                        : mceb );
      fitResH->Fill(valFit);
      statUncH->Fill(valFitErr);
      float bias=valFit-trueVal;
      biasH->Fill(bias);
      float statUnc=hfcFitter.model.r->getError();
      if(statUnc>0) pullH->Fill(bias/statUnc);
    }
  

  //save the templates and results to a file for posterity
  TFile *outF = TFile::Open("HeavyFlavorFitCalibration.root","RECREATE");
  outF->cd();
  fitResH->Write();
  statUncH->Write();
  statUncH->Write();
  biasH->Write();
  pullH->Write();
  for(std::map<TString, ProcHistos_t >::iterator it = histosForProc.begin();  it != histosForProc.end(); it++)
    {
      TDirectory *procDir = outF->mkdir(it->first);
      procDir->cd();
      for(ProcHistos_t::iterator hit=it->second.begin(); hit!=it->second.end(); hit++) (*hit)->Write(); 
    }
  outF->Close();


  /*      
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
  */


}



