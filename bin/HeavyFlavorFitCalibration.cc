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
  //
  // map the samples per type
  //
  std::map<TString,Int_t> procYields;
  std::map<TString,TTree *> procTrees,dataTrees;

  //signal
  float rToGen(1);
  if(syst=="matchingdown")                    procYields["TTJets_matchingdown"]=0;
  else if(syst=="matchingup")                 procYields["TTJets_matchingup"]=0;
  else if(syst=="scaledown")                  procYields["TTJets_scaledown"]=0;
  else if(syst=="scaleup")                    procYields["TTJets_scaleup"]=0;
  else if(syst.Contains("flavor")) 
    {
      procYields["TT2l_R1"]=0;
      procYields["TT2l_R05a"]=0;
      procYields["TT2l_R05b"]=0;
      procYields["TT2l_R0"]=0;
      sscanf(syst.Data(),"flavor%f",&rToGen);
    }
  else                                        procYields["TTJets_signal"]=0;
  procYields["SingleTbar_tW"]=0;
  procYields["SingleT_tW"]=0;
  procYields["SingleTbar_t"]=0;
  procYields["SingleT_t"]=0;
  procYields["SingleTbar_s"]=0;
  procYields["SingleT_s"]=0;
  procYields["TTJets"]=0;
  
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
      if(it->first.Contains("DY") ) sf=1.5;
      if(it->first.Contains("DY") && syst=="dyup")   sf*=1.3;
      if(it->first.Contains("DY") && syst=="dydown") sf*=0.7;
      if(it->first=="TT2l_R1"   && syst.Contains("flavor") ) sf *= pow(rToGen,2);
      if(it->first=="TT2l_R05a" && syst.Contains("flavor") ) sf *= rToGen*(1-rToGen);
      if(it->first=="TT2l_R05b" && syst.Contains("flavor") ) sf *= rToGen*(1-rToGen);
      if(it->first=="TT2l_R0"   && syst.Contains("flavor") ) sf *= pow(1-rToGen,2);

      procYields[it->first] = sf*dataLumi*xsecWeight*t->GetEntriesFast();
      //  cout << "\t" << it->first <<  " " <<  procYields[it->first] << " evts" << endl; 
    }

  //
  // configure the fitter
  //
  int maxJets=int(fitPars["maxJets"]);
  HFCMeasurement hfcFitter(fitType,maxJets,int(fitPars["smR"]));
  hfcFitter.configureBtagAlgo   (btagWP, fitPars[btagWP+"_cut"]);
  hfcFitter.setBtagEfficiency   (fitPars[btagWP+"_effb"], fitPars[btagWP+"_sfb"], fitPars[btagWP+"_sfbunc"]);
  hfcFitter.setMistagEfficiency (fitPars[btagWP+"_effq"], fitPars[btagWP+"_sfq"], fitPars[btagWP+"_sfqunc"]);
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

  return;
  //
  // run the calibration
  //
  top::EventSummaryHandler evHandler, ensembleHandler;
  for(int iPE=1; iPE<=maxPE; iPE++)
    {
      if(iPE%5==0) { printf("\rCompleted %d/%d ",iPE,maxPE); cout << flush; }
      
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
      

      //reset tree used ensemble
      ensembleHandler.getTree()->Delete("all");
    }
  

  if(!freezeResults) return;
  /*      
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

  inF->Close();    
}



