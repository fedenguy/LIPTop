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

#include "LIP/Top/interface/HFCMeasurement.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <vector>

using namespace std;
using namespace RooFit;

typedef std::vector<TString> Proc_t;
typedef std::vector<Float_t> ProcYields_t;
typedef std::vector<TH1D *> ProcHistos_t;


stringstream report;
 
//
void runCalibration(TString url, int fitType, int nuisanceType, int runMode, TString btagWP, 
		    std::map<TString,Double_t> &fitPars, int maxPE, Double_t dataLumi, TString syst,  bool freezeResults);

//
void printHelp()
{
  printf("--help    --> print this\n");
  printf("--iLumi   --> integrated luminosity to be used /pb  (default:2165)\n");
  printf("--in      --> input file with the summary trees\n");
  printf("--par     --> fitter configuration file\n");
  printf("--fit     --> fit type code                         (default: 0 = R)\n");
  printf("--nui     --> nuisance for fcorrect                 (default: 0 = Gaussian)\n");
  printf("--run     --> run mode                              (default: 0 = inclusive only)\n");
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
  Int_t nuisanceType(HFCMeasurement::GAUSSIAN);
  int runMode(HFCMeasurement::FITINCLUSIVEONLY);
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
      if(arg.find("--syst")!=string::npos)              { syst=argv[i+1];                                                   printf("syst    = %s\n", syst.Data());       }
      if(arg.find("--npe")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&maxPE);                              i++;  printf("N_{PE}  = %d\n", maxPE);             }
      if(arg.find("--fit")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&fitType);                            i++;  printf("fitType  = %d\n", fitType);          }
      if(arg.find("--nui")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&nuisanceType);                       i++;  printf("nuisanceType = %d\n", nuisanceType); }
      if(arg.find("--run")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&runMode);                            i++;  printf("runMode = %d\n", runMode);           }
      if(arg.find("--btag")!=string::npos && i+1<argc)  { btagWP=argv[i+1];                                           i++;  printf("btag WP  = %s\n", btagWP.Data());    }
    } 
  //  if(url=="" || fitParFile=="") { printHelp(); return 0; }
  if(fitParFile=="") { printHelp(); return 0; }
  

  HFCMeasurement fitter(HFCMeasurement::FIT_R, fitParFile);

  /*  
      bool freezeResults(true);
      std::map<TString,Double_t> fitPars=parseParametersFrom(fitParFile);
      runCalibration(url, fitType, nuisanceType, runMode, btagWP, fitPars, maxPE, dataLumi, syst, freezeResults );
      
      cout << report.str() << endl;
  */
}

/*
	  
//
void runCalibration(TString url, int fitType, int nuisanceType, int runMode, TString btagWP, 
		    std::map<TString,Double_t> &fitPars, int maxPE, Double_t dataLumi, TString syst,  bool freezeResults)
{

  TH1D *fitResH = new TH1D("fitres",";Fit result;Pseudo-experiments",100,0,2);             fitResH->SetDirectory(0);
  TH1D *fitDiffH[6];
  for(size_t icat=0; icat<6; icat++) 
    {
      TString catStr(""); catStr += icat;
      fitDiffH[icat] = new TH1D("fitdiff"+catStr,";R_{i}-R;Pseudo-experiments",40,-0.195,0.205);
      fitDiffH[icat]->SetDirectory(0);
    }
  TH1D *biasH = new TH1D("bias",";bias;Pseudo-experiments",50,-0.052,0.048);               biasH->SetDirectory(0);
  TH1D *pullH = new TH1D("pull",";pull = bias / #sigma;Pseudo-experiments",50,-5.2,4.8);   pullH->SetDirectory(0);
  TH1D *statUncH = new TH1D("staterr",";#sigma;Pseudo-experiments",50,-0.05,0.05);         statUncH->SetDirectory(0);
  TH1D *lowerPointH = new TH1D("lowerpoint",";bias;Pseudo-experiments",100,0.5,1.);       lowerPointH->SetDirectory(0);
  TH1D *upperPointH = new TH1D("uppperpoint",";bias;Pseudo-experiments",100,0.5,1.);      upperPointH->SetDirectory(0);

  //
  // map the samples per type
  //
  std::map<TString, Proc_t>        descForProc;
  std::map<TString, ProcYields_t > yieldsForProc;
  std::map<TString, ProcHistos_t > histosForProc;
  
  reweight::PoissonMeanShifter *puShifter=0;
  if(syst=="puup")   puShifter = new reweight::PoissonMeanShifter(+0.6);
  if(syst=="pudown") puShifter = new reweight::PoissonMeanShifter(-0.6);

  //signal
  Double_t rToGen(1.0);
  ProcYields_t procYields(6,0);
  procYields[HFCMeasurement::EE_2JETS]=548;
  procYields[HFCMeasurement::EE_3JETS]=274;
  procYields[HFCMeasurement::MUMU_2JETS]=702;
  procYields[HFCMeasurement::MUMU_3JETS]=352;
  procYields[HFCMeasurement::EMU_2JETS]=1884;
  procYields[HFCMeasurement::EMU_3JETS]=1046;
  yieldsForProc["Signal"]=procYields;

  if(syst=="matchingdown")        descForProc["Signal"] = Proc_t(1,"TTJets_matchingdown");
  else if(syst=="matchingup")     descForProc["Signal"] = Proc_t(1,"TTJets_matchingup");
  else if(syst=="scaledown")      descForProc["Signal"] = Proc_t(1,"TTJets_scaledown");
  else if(syst=="scaleup")        descForProc["Signal"] = Proc_t(1,"TTJets_scaleup");
  else if(syst.Contains("flavor")) 
    {
      sscanf(syst.Data(),"flavor%lf",&rToGen);
      cout << "*********************" << rToGen << " ********************* " << endl;

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
  procYields[HFCMeasurement::EE_2JETS]=33.7;
  procYields[HFCMeasurement::EE_3JETS]=10.8;
  procYields[HFCMeasurement::MUMU_2JETS]=43;
  procYields[HFCMeasurement::MUMU_3JETS]=13.5;
  procYields[HFCMeasurement::EMU_2JETS]=111.6;
  procYields[HFCMeasurement::EMU_3JETS]=36.5;
  if(!syst.Contains("flavor"))
    {
      descForProc["SingleTop"]=SingleTop;
      yieldsForProc["SingleTop"]=procYields;
    }

  Proc_t OtherTTbar;
  OtherTTbar.push_back("TTJets");
  procYields[HFCMeasurement::EE_2JETS]=4.8;
  procYields[HFCMeasurement::EE_3JETS]=3.2;
  procYields[HFCMeasurement::MUMU_2JETS]=0.9;
  procYields[HFCMeasurement::MUMU_3JETS]=2.2;
  procYields[HFCMeasurement::EMU_2JETS]=7.6;
  procYields[HFCMeasurement::EMU_3JETS]=8.2;
  if(!syst.Contains("flavor"))
    {
      descForProc["OtherTTbar"]=OtherTTbar;
      yieldsForProc["OtherTTbar"]=procYields;
    }

  Proc_t DiBosons;
  DiBosons.push_back("WW");
  DiBosons.push_back("ZZ");
  DiBosons.push_back("WZ");
  procYields[HFCMeasurement::EE_2JETS]=11.8;
  procYields[HFCMeasurement::EE_3JETS]=2.0;
  procYields[HFCMeasurement::MUMU_2JETS]=14;
  procYields[HFCMeasurement::MUMU_3JETS]=2.6;
  procYields[HFCMeasurement::EMU_2JETS]=38;
  procYields[HFCMeasurement::EMU_3JETS]=6.7;
  if(!syst.Contains("flavor"))
    {
      descForProc["DiBosons"]=DiBosons;
      yieldsForProc["DiBosons"]=procYields;
    }

  Proc_t DY;
  DY.push_back("DYJetsToLL");
  DY.push_back("DYJetsToLL_M10to50");
  float uncScaleFactor(1.0);
  if(syst=="dyup")   uncScaleFactor *= 1.15;
  if(syst=="dydown") uncScaleFactor *= 0.85;
  procYields[HFCMeasurement::EE_2JETS]=291*uncScaleFactor;
  procYields[HFCMeasurement::EE_3JETS]=96*uncScaleFactor;
  procYields[HFCMeasurement::MUMU_2JETS]=424*uncScaleFactor;
  procYields[HFCMeasurement::MUMU_3JETS]=139*uncScaleFactor;
  procYields[HFCMeasurement::EMU_2JETS]=159.6*uncScaleFactor;
  procYields[HFCMeasurement::EMU_3JETS]=30.3*uncScaleFactor;
  if(!syst.Contains("flavor"))
    {
      descForProc["DY"]=DY;
      yieldsForProc["DY"]=procYields;
    }


  //MC truth for the fit parameters
  Double_t fcorrect[6], fcorrectNorm[6], fttbar[6], fttbarNorm[6], fsingletop[6], fsingletopNorm[6];
  Double_t avgbeffPerCat[6], avgbeffPerCatNorm[6], avgleffPerCat[6], avgleffPerCatNorm[6];
  Double_t avgbeff(0.),avgbeffNorm(0.),avgleff(0.), avgleffNorm(0.);
  for(size_t i=0; i<6; i++) 
    { 
      fcorrect[i]=0;      fcorrectNorm[i]=0; 
      fttbar[i]=0;        fttbarNorm[i]=0; 
      fsingletop[i]=0;    fsingletopNorm[i]=0; 
      avgbeffPerCat[i]=0; avgbeffPerCatNorm[i]=0;
      avgleffPerCat[i]=0; avgleffPerCatNorm[i]=0;
    }


  //
  // fill the histograms to be used for the pseudo-experiments
  //
  int maxJets=int(fitPars["maxJets"]);
  Double_t btagWPcut(fitPars[btagWP+"_cut"]);
  TFile *inF = TFile::Open(url);
  top::EventSummaryHandler evHandler;
  cout << "[Filling histograms]" << flush;
  for(std::map<TString, Proc_t>::iterator pIt=descForProc.begin(); pIt!=descForProc.end(); pIt++)
    {
      bool isSignalRequired(pIt->first.Contains("Signal"));
      bool isSingleTop(pIt->first.Contains("SingleT"));
      bool isTTbar(pIt->first.Contains("TT") || isSignalRequired);

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
	      int nbtags(0);
	      int ncorrect(0);
	      float nbstagged(0), nlstagged(0), nbs(0), nls(0);
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
		  nbtags += isBtagged;
		  nbstagged += isMatchedToB*isBtagged;     nbs += isMatchedToB;
		  nlstagged += (!isMatchedToB)*isBtagged;  nls += !isMatchedToB;
		  ncorrect += ((isTTbar || isSingleTop) && fabs(phys.jets[ijet].genid)==5 && isMatchedToB);
		}
	      if(njets>maxJets || njets<2) continue;
	      if(ncorrect>2) ncorrect=2;
	      
	      
	      //check event category
	      int icat=-1;
	      if(ev.cat==EE &&   njets==2)  icat=HFCMeasurement::EE_2JETS;
	      if(ev.cat==EE &&   njets==3)  icat=HFCMeasurement::EE_3JETS;
	      if(ev.cat==MUMU && njets==2)  icat=HFCMeasurement::MUMU_2JETS;
	      if(ev.cat==MUMU && njets==3)  icat=HFCMeasurement::MUMU_3JETS;
	      if(ev.cat==EMU &&  njets==2)  icat=HFCMeasurement::EMU_2JETS;
	      if(ev.cat==EMU &&  njets==3)  icat=HFCMeasurement::EMU_3JETS;
	      if(icat<0) continue;

	      fttbarNorm[icat]+=weight;
	      if(isTTbar)     { fttbar[icat] += weight;    fsingletopNorm[icat] += weight; }
	      if(isSingleTop) { fsingletop[icat] += weight; }
	      //	      if(it->Contains("TT2l_R1") || it->Contains("TTJets") || isSingleTop)
	      fcorrect[icat]     += ncorrect*weight;
	      fcorrectNorm[icat] += 2*njets*weight;
	      avgleff             += nlstagged*weight;   avgleffNorm             += nls*weight; 
	      avgleffPerCat[icat] += nlstagged*weight;   avgleffPerCatNorm[icat] += nls*weight;
	      avgbeff             += nbstagged*weight;   avgbeffNorm             += nbs*weight; 
	      avgbeffPerCat[icat] += nbstagged*weight;   avgbeffPerCatNorm[icat] += nbs*weight;

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
  HFCMeasurement hfcFitter(fitType, nuisanceType,maxJets,rToGen);
  hfcFitter.configureBtagAlgo   (btagWP, btagWPcut);
  double mceb=fitPars[btagWP+"_effb"];
  double mceq=fitPars[btagWP+"_effq"];
  if(fitType==HFCMeasurement::FIT_R)
    {
      avgbeff /= avgbeffNorm;
      mceb=avgbeff;
      avgleff /= avgleffNorm;
      mceq=avgleff;
      report << "MC truth for fit parameters" << endl
	     << "[Inclusive]"
	     << " effb=" << avgbeff 
	     << " effq=" << avgleff << endl;
      
      for(size_t i=0; i<6; i++)
	{
	  report << "\t cat #" << i << ")"
		 << " effb=" << avgbeffPerCat[i]/avgbeffPerCatNorm[i]
		 << " effq=" << avgleffPerCat[i]/avgleffPerCatNorm[i]
		 << " fcorrect=" << fcorrect[i]/fcorrectNorm[i] 
		 << " fttbar=" << fttbar[i]/fttbarNorm[i] 
		 << " fsingletop=" << fsingletop[i]/fsingletopNorm[i] << endl;
	}
    }
  hfcFitter.setBtagEfficiency   (mceb, fitPars[btagWP+"_sfb"], fitPars[btagWP+"_sfbunc"]);
  hfcFitter.setMistagEfficiency (mceq, fitPars[btagWP+"_sfq"], fitPars[btagWP+"_sfqunc"]);

  TString channels[]={"ee","mumu","emu"};
  for(int ich=0; ich<3; ich++)
    {
      for(int ijets=2; ijets<=maxJets; ijets++)
	{
	  TString key(channels[ich]+"_"); key += ijets; key += "jets";
	  hfcFitter.setParametersForCategory( fitPars[key+"_fcorrect"],   fitPars[key+"_fcorrectunc"],
					      fitPars[key+"_fttbar"],     fitPars[key+"_fttbarunc"],
					      fitPars[key+"_fsingletop"], fitPars[key+"_fsingletopunc"],
					      fitPars[key+"_effb"]/mceb,  fitPars[key+"_effq"]/mceq, 
					      ijets, channels[ich]);

	}
    }
  //  hfcFitter.printConfiguration(report);

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

      if(syst=="vareffs" && iPE>1)
	{
	  Double_t varsfeb=fitPars[btagWP+"_sfb"]+gRandom->Gaus(0,fitPars[btagWP+"_sfbunc"]);
	  hfcFitter.model.sfeb->setVal(varsfeb);
	  Double_t varsfeq=fitPars[btagWP+"_sfq"]+gRandom->Gaus(0,fitPars[btagWP+"_sfqunc"]);
	  hfcFitter.model.sfeq->setVal(varsfeq);
	  cout << fitPars[btagWP+"_sfb"] << " ---> " << varsfeb << endl
	       << fitPars[btagWP+"_sfq"] << " ---> " << varsfeq << endl;
	}
      
      //fit
      hfcFitter.fitHFCtoMeasurement(ensembleHistos,runMode,false);

      if(syst=="vareffs")
	{
	  hfcFitter.model.sfeb->setVal(fitPars[btagWP+"_sfb"]);
	  hfcFitter.model.sfeq->setVal(fitPars[btagWP+"_sfq"]);
	}

      double effb(hfcFitter.model.abseb->getVal()) , effq(hfcFitter.model.abseq->getVal());
      report << "[Ensemble #" << iPE << "]" << endl; 
      report << "\t R   = " << hfcFitter.model.r->getVal() << " +" << hfcFitter.model.r->getAsymErrorHi() << " " << hfcFitter.model.r->getAsymErrorLo() << endl;
      //for(int icat=0; icat<6; icat++) report << "\t\t R_{" << icat << "}  = " << hfcFitter.model.rFit[icat] << " +" << hfcFitter.model.rFitAsymmErrHi[icat] << " " << hfcFitter.model.rFitAsymmErrLo[icat] << endl;
      // report << "\t e_b = " << effb*hfcFitter.model.sfeb->getVal() << " +" << effb*hfcFitter.model.sfeb->getAsymErrorHi() << " " << effb*hfcFitter.model.sfeb->getAsymErrorLo() << endl;
      // for(int icat=0; icat<6; icat++) report << "\t\t eb_{" << icat << "}  = " << hfcFitter.model.ebFit[icat] << " +" << hfcFitter.model.ebFitAsymmErrHi[icat] << " " << hfcFitter.model.ebFitAsymmErrLo[icat] << endl;
      //    report << "\t e_q = " << effq*hfcFitter.model.sfeq->getVal() << " +" << effq*hfcFitter.model.sfeq->getAsymErrorHi() << " " << effq*hfcFitter.model.sfeq->getAsymErrorLo() << endl;
      // for(int icat=0; icat<6; icat++) report << "\t\t eq_{" << icat << "}  = " << hfcFitter.model.eqFit[icat] << " +" << hfcFitter.model.eqFitAsymmErrHi[icat] << " " << hfcFitter.model.eqFitAsymmErrLo[icat] << endl;

      float valFit   ( fitType==HFCMeasurement::FIT_R ? hfcFitter.model.rFitResult           : effb*hfcFitter.model.sfeb->getVal() );
      float valFitErr( fitType==HFCMeasurement::FIT_R ? hfcFitter.model.rFitResultAsymmErrHi : effb*hfcFitter.model.sfeb->getError() );
      float trueVal  ( fitType==HFCMeasurement::FIT_R ? rToGen                               : mceb );
      fitResH->Fill(valFit);
      if(runMode==HFCMeasurement::FITEXCLUSIVECATEGORIES)
	for(int icat=0; icat<6; icat++) fitDiffH[icat]->Fill(hfcFitter.model.rFit[icat]-valFit);
      statUncH->Fill(valFitErr);
      float bias=valFit-trueVal;
      biasH->Fill(bias);
      if(valFitErr>0) pullH->Fill(bias/valFitErr);
      if(fitType==HFCMeasurement::FIT_R_CONSTRAINED)
	{
	  lowerPointH->Fill(hfcFitter.model.rFitLowerLimit);
	  upperPointH->Fill(hfcFitter.model.rFitUpperLimit);
	}

    }
  

  //save the templates and results to a file for posterity
  TFile *outF = TFile::Open("HeavyFlavorFitCalibration.root","UPDATE");
  outF->cd();
  TString baseDir=syst;
  if(baseDir=="") baseDir="std";
  if(nuisanceType!=HFCMeasurement::GAUSSIAN) { baseDir += "_"; baseDir += nuisanceType; }
  outF->rmdir(baseDir);   //remove if already existing
  TDirectory *outDir = outF->mkdir(baseDir);
  outDir->cd();
  fitResH->Write();
  for(int icat=0; icat<6; icat++)
    {
      fitDiffH[icat]->Fit("gaus","LMEQ+");
      fitDiffH[icat]->Write();
    }
  statUncH->Write();
  statUncH->Write();
  biasH->Write();
  pullH->Write();
  lowerPointH->Write();
  upperPointH->Write();
  for(std::map<TString, ProcHistos_t >::iterator it = histosForProc.begin();  it != histosForProc.end(); it++)
    {
      TDirectory *procDir = outDir->mkdir(it->first);
      procDir->cd();
      for(ProcHistos_t::iterator hit=it->second.begin(); hit!=it->second.end(); hit++) (*hit)->Write(); 
    }
  outF->Close();

}


*/
