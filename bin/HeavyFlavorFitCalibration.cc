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
 
TH1 *getObservedBtagFrom(TString url, TString wp, TString sampleType);

//check me
void runCalibration(TString url, int fitType, int nuisanceType, int runMode, TString btagWP, 
		    std::map<TString,Double_t> &fitPars, int maxPE, Double_t dataLumi, TString syst,  bool freezeResults);

//
void printHelp()
{
  printf("--help    --> print this\n");
  printf("--in      --> input file with the summary trees\n");
  printf("--ws      --> input file with the workspace\n");
  printf("--par     --> fitter configuration file\n");
  printf("--btag    --> tag efficiencies configuration file\n");
  printf("--fit     --> fit type code                         (default: R=0,eb=1,ebvsR=2,vtb=3)\n");
  printf("--npe     --> number of pseudo-experiments to throw (default:1)\n");
  printf("command line examples: HeavyFlavorFitCalibration --in plotter.root --par hfcParams_cfg.json --btag wpParams_cfg.json\n");
  printf("                       HeavyFlavorFitCalibration --ws WS.root --npe 100\n");
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

  TString url(""),wsurl("");
  TString fitParFile(""),btagParFile("");
  Int_t fitType(HFCMeasurement::FIT_R);
  TString syst("");
  int maxPE(1);
  for(int i=1;i<argc;i++)
    {
      string arg(argv[i]);
      if(arg.find("--help")!=string::npos)              { printHelp();	  return 0; }
      if(arg.find("--in")!=string::npos && i+1<argc)    { url=argv[i+1];         gSystem->ExpandPathName(url);         i++;  printf("in      = %s\n", url.Data());             }
      if(arg.find("--ws")!=string::npos && i+1<argc)    { wsurl=argv[i+1];       gSystem->ExpandPathName(wsurl);       i++;  printf("ws      = %s\n", wsurl.Data());           }
      if(arg.find("--par")!=string::npos && i+1<argc)   { fitParFile=argv[i+1];  gSystem->ExpandPathName(fitParFile);  i++;  printf("parFile = %s\n", fitParFile.Data());      }
      if(arg.find("--btag")!=string::npos && i+1<argc)  { btagParFile=argv[i+1]; gSystem->ExpandPathName(btagParFile); i++;  printf("btagFile  = %s\n", btagParFile.Data());   }
      if(arg.find("--npe")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&maxPE);                               i++;  printf("N_{PE}  = %d\n", maxPE);                  }
      if(arg.find("--fit")!=string::npos && i+1<argc)   { sscanf(argv[i+1],"%d",&fitType);                             i++;  printf("fitType  = %d\n", fitType);               }
    } 
  if(url=="" && wsurl=="")      { printHelp(); return 0; }
  if(url!="" && fitParFile=="") { printHelp(); return 0; }

  setStyle();
  gStyle->SetOptFit(1111);

  //keep RooFit quiet                                                                                                                                                          
  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().getStream(0).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Minimization);
  RooMsgService::instance().getStream(1).removeTopic(Fitting);
  RooMsgService::instance().getStream(0).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Eval);
  RooMsgService::instance().getStream(1).removeTopic(Integration);
  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(ObjectHandling);


  //
  // STANDARD FIT OR MC CLOSURE
  //
  if(wsurl=="")
    {
      //run pseudo-experiments
      HFCMeasurement *fitter=new HFCMeasurement(fitType, fitParFile,btagParFile);
      TString sampleType = fitter->getSampleType(); 
      TString wp         = fitter->getWP();
      
      //get observed/expected b-tag multiplicity from plotter
      TH1 *btagObs=getObservedBtagFrom(url,wp,sampleType);
      if(btagObs==0) {
	cout << "Unable to find " << sampleType << " b-tag distribution for " << wp << " in " << url << endl;
	return 0;
      }
      
      if(sampleType!="data")
	{
	  //fitter.fitHFCfrom(btagObs,true);
	  
	  std::map<std::string, TH1F *> closureTests;
	  TH1 *btagObsSample=0;
	  for(int ipe=1; ipe<=maxPE; ipe++)
	    {
	      
	      //instantiate a new fitter
	      fitter =new HFCMeasurement(HFCMeasurement::FIT_R, fitParFile,btagParFile);
	      if(btagObsSample==0){
		btagObsSample=(TH1 *)btagObs->Clone("btagobssample");
		btagObsSample->SetDirectory(0);
	      }
	      
	      //sample distribution
	      btagObsSample->Reset("ICE");
	      btagObsSample->FillRandom(btagObs,btagObs->Integral());
	      
	      //fit it
	      fitter->fitHFCfrom(btagObsSample,false);
	      
	      //fill closure test histograms
	      std::map<std::string,HFCMeasurement::FitResult_t> &res=fitter->getResults();
	      for(std::map<std::string,HFCMeasurement::FitResult_t>::iterator rit = res.begin(); rit!= res.end(); rit++)
		{
		  if(rit->first.find("jets")==string::npos) continue;
		  if(rit->second.status==false) continue;
		  TH1F *biasH=0; 
		  if(closureTests.find(rit->first)==closureTests.end())
		    {
		      biasH = new TH1F(("bias_"+rit->first).c_str(),";bias;Pseudo-experiments",50,-0.052,0.048);               
		      biasH->SetDirectory(0);
		      closureTests[rit->first]=biasH;
		    }
		  else
		    biasH=closureTests[rit->first];
		  if(biasH) biasH->Fill(rit->second.poiFit-1.0);
		}
	      
	      delete fitter;
	    }
	  
	  //show results
	  TCanvas *c = new TCanvas("c","c",1200,1200);
	  c->Divide(3,3);
	  int icat=0;
	  for(std::map<std::string, TH1F *>::iterator rit=closureTests.begin(); rit!=closureTests.end(); rit++,icat++)
	    {
	      string tag=rit->first;
	      c->cd(icat+1);
	      rit->second->Draw("e1");
	      rit->second->GetXaxis()->SetNdivisions(10);
	      rit->second->GetYaxis()->SetTitleOffset(1.2);
	      rit->second->Fit("gaus");
	      
	      //category
	      TPaveText *pt = new TPaveText(0.2,0.85,0.45,0.94,"brNDC");
	      pt->SetBorderSize(0);
	      pt->SetFillColor(0);
	      pt->SetFillStyle(0);
	      pt->SetTextFont(42);
	      TString caption=tag.c_str();
	      caption=caption.ReplaceAll("mu","#mu");
	      caption=caption.ReplaceAll("2",",=2 ");
	      caption=caption.ReplaceAll("3",",=3 ");
	      caption=caption.ReplaceAll("4",",=4 ");
	      pt->AddText(caption);
	      pt->Draw();
	      
	      //overall header
	      if(icat==0)
		{
		  pt = new TPaveText(0.1,0.96,0.9,1.0,"brNDC");
		  pt->SetBorderSize(0);
		  pt->SetFillColor(0);
		  pt->SetFillStyle(0);
		  pt->AddText("CMS simulation");
		  pt->Draw();
		}
	    }

	  c->SaveAs("HFCClosureTest.png");
	  c->SaveAs("HFCClosureTest.pdf");
	  c->SaveAs("HFCClosureTest.C");
	} 
      else
	{
	  fitter->fitHFCfrom(btagObs,true);
	  
	  std::map<std::string,HFCMeasurement::FitResult_t> &res=fitter->getResults();
	  
	  //display results it a tex table
	  string sample[]={"inclusive","ee","mumu","emu"};
	  string sampleTitle[]={"inclusive","$ee$","$\\mu\\mu$","$e\\mu$"};
	  
	  //table header
	  cout << endl
	       << "\\begin{table}" << endl
	       << "\\begin{center}" << endl
	       << "\\caption{}" << endl
	       << "\\label{tab:hfcsummarytable}" << endl
	       << "\\begin{tabular}{lcccc}\\hline\\hline" << endl
	       << "Sample" << flush;
	  for(int isample=0; isample<4; isample++) cout << " & " << sampleTitle[isample];
	  cout << "\\\\\\hline" << endl;
	  
	  //result
	  cout << "$\\mathcal{" << fitter->title() << "}$" << flush; 
	  for(int isample=0; isample<4; isample++) cout << " & " << res[ sample[isample] ].poiFit << "$\\pm$" << res[ sample[isample] ].poiErr;
	  cout << "\\\\\\hline" << endl;
	  
	  //uncertainties
	  /*
	    std::map<std::string,Double_t> &nuiList=res[sample[0]].uncBreakup;
	    for(std::map<std::string,Double_t>::iterator it=nuiList.begin(); it!=nuiList.end(); it++)
	    {
	    if(it->first=="stat")      cout << "Stat. unc." << flush;
	    else if(it->first=="syst") cout << "Syst. unc." << flush;
	    else                       cout << "~~~~" << it->first << flush; 
	    for(int isample=0; isample<4; isample++) 
	    {
	    if(res[sample[isample]].uncBreakup.find(it->first)==res[sample[isample]].uncBreakup.end())
	    cout << " & ";
	    else
	    cout << " & " << res[ sample[isample] ].uncBreakup[it->first];
	    }
	    
	    if(it->first=="stat" || it->first=="syst")  cout << "\\\\\\hline" << endl;
	    else                                        cout << "\\\\" << endl;
	    }
	  */
	  cout << "\\hline\\end{tabular}" << endl
	       << "\\end{center}" << endl
	       << "\\end{table}" << endl << endl;
	}
    }
  //
  // UNCERTAINTY BREAKUP
  //
  else
    {
      //linearity check
      TCanvas *c=new TCanvas("c","c",1200,400);
      c->Divide(5,1);
      for(int ip=0; ip<5; ip++)
	{
	  float rgen=ip*1./5;

	  TH1F *templateH=0, *ensembleH=0;
	  for(int ipe=0; ipe<=maxPE; ipe++)
	    {
	      TFile *inF = TFile::Open(wsurl);
	      RooWorkspace *ws=(RooWorkspace *) inF->Get("w");
	      
	      HFCMeasurement *fitter=new HFCMeasurement(ws,4,HFCMeasurement::FIT_R);
	      RooRealVar* firstPOI = (RooRealVar*) ws->set("poi")->first();
	      ws->loadSnapshot("default");

	      TIterator *nuis_params_itr = ws->set("nuisances")->createIterator();
	      TObject *nuis_params_obj;
	      while((nuis_params_obj=nuis_params_itr->Next())){
		RooRealVar *nuiVar=(RooRealVar *)nuis_params_obj;
		//TString name(nuiVar->GetName());
		//if(name.Contains("fcor_") || name.Contains("tt_")) 
		// {
		nuiVar->setConstant(kTRUE);
		// }
	      }
	      firstPOI->setVal(rgen);
	      if(templateH==0) { 
		templateH=fitter->generateBtagObs();
		ensembleH=(TH1F *)templateH->Clone("ensemble");
		ensembleH->SetDirectory(0);
		delete fitter;
		continue;
	      }
	      
	      ensembleH->Reset("ICE");
	      ensembleH->FillRandom(templateH,templateH->Integral());
	      fitter->plrFit(ensembleH);
	      if(ipe==1)
		{
		  c->cd(ip+1);
		  templateH->DrawClone("hist");
		  TH1F *fitH=fitter->generateBtagObs();
		  fitH->SetLineColor(kBlue);
		  fitH->DrawClone("histsame");
		}
	      delete fitter;
	      
	      inF->Close();
	    }
	  
	  delete templateH;
	  delete ensembleH;
	}

      c->SaveAs("tmp.png");
    }
}

//
TH1 *getObservedBtagFrom(TString url, TString wp, TString sampleType)
{
  
  TH1 *btagObs=0;
  TFile *fIn = TFile::Open(url);
  if(fIn==0) return btagObs;

  TList *dirs=fIn->GetListOfKeys();
  for(int iproc=0; iproc<dirs->GetEntries(); iproc++)
    {
      TString iDir(dirs->At(iproc)->GetName());
      if(sampleType!="data" && iDir=="data") continue;
      if(sampleType=="data" && iDir!="data") continue;
      
      TH1 *h    = (TH1 *)fIn->Get(iDir+"/"+wp+"btagsextended");
      if(h==0) continue;
      if(btagObs==0) { btagObs = (TH1 *)h->Clone("btagobs"); btagObs->SetDirectory(0); }
      else           { btagObs->Add(h); }
    }
  fIn->Close();
  
  return btagObs;
}  
