#include <iostream>
#include <boost/shared_ptr.hpp>
#include <fstream>

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/KinResultsHandler.h"
#include "LIP/Top/interface/KinAnalysis.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/ObjectFilters.h"
#include "CMGTools/HtoZZ2l2nu/interface/SelectionMonitor.h"

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TSystem.h"
#include "TFile.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TCanvas.h"
#include "TString.h"
#include "TDirectory.h"

#include "CMGTools/HtoZZ2l2nu/interface/TMVAUtils.h"


using namespace std;

//
int main(int argc, char* argv[])
{
  SelectionMonitor controlHistos; //plot storage

  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  //check arguments
  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
    return 0;
  }
  
  //configure                                                                                                                                                                                                                 
  const edm::ParameterSet &runProcess = edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("runProcess");
  TString evurl               = runProcess.getParameter<std::string>("input");
  TString outUrl              = runProcess.getParameter<std::string>("outdir");
  TString kindir              = runProcess.getParameter<std::string>("kindir");
  bool isMC                   = runProcess.getParameter<bool>("isMC");
  int mcTruthMode             = runProcess.getParameter<int>("mctruthmode");
  TString dirname             = runProcess.getParameter<std::string>("dirName");
  bool saveSummaryTree        = runProcess.getParameter<bool>("saveSummaryTree");
  bool useMVA                 = runProcess.getParameter<bool>("useMVA");
  if(mcTruthMode!=2) useMVA=false;
  edm::ParameterSet tmvaInput = runProcess.getParameter<edm::ParameterSet>("tmvaInput");
  TString studyTag            = runProcess.getParameter<std::string>("studyTag");
  TString weightsDir          = runProcess.getParameter<std::string>("weightsDir");
  std::vector<std::string> methodList = runProcess.getParameter<std::vector<std::string> >("methodList");
  std::vector<std::string> varsList    = runProcess.getParameter<std::vector<std::string> >("varsList");
  bool trainMVA                        = runProcess.getParameter<bool>("doTrain");

  //book histos
  controlHistos.addHistogram( new TH1F ("njets", ";Jets;Events", 6, 0.,6.) );
  controlHistos.addHistogram( new TH1F ("btags", ";b-tag multiplicity;Events", 6, 0.,6.) );
  for(int ibin=1; ibin<=controlHistos.getHisto("njets","all")->GetXaxis()->GetNbins(); ibin++)
    {
      TString label(""); label += ibin-1; 
      if(controlHistos.getHisto("njets","all")->GetXaxis()->GetNbins()==ibin) label ="#geq" + label;
      controlHistos.getHisto("njets","all")->GetXaxis()->SetBinLabel(ibin,label + " jets");
      controlHistos.getHisto("btags","all")->GetXaxis()->SetBinLabel(ibin,label + " b-tags");
    }
  
  controlHistos.addHistogram( new TH1F ("leadjet", "; Leading jet p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("subleadjet", "; Sub-leading jet p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("leadlepton", "; Leading lepton p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("subleadlepton", "; Sub-leading lepton p_{T} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("met", "; #slash{E}_{T} [GeV/c]; Events / (10 GeV/c)", 40, 0.,400.) );
  controlHistos.addHistogram( new TH1F ("ht", "; #sum_{jets} [GeV/c]; Events / (20 GeV/c)",40, 0.,800.) );
  controlHistos.addHistogram( new TH1F ("st", "; #sum_{leptons,E_{T}^{miss}} p_{T} [GeV/c]; Events / (20 GeV/c)",40, 0.,800.) );
  controlHistos.addHistogram( new TH1F ("sumpt", "; #sum_{leptons} p_{T} [GeV/c]; Events / (20 GeV/c)",25, 0.,500.) );
  controlHistos.addHistogram( new TH1F ("htlep", "; #sum_{jets,leptons,E_{T}^{miss}} [GeV/c]; Events / (20 GeV/c)",70, 0.,1400.) );
  controlHistos.addHistogram( new TH1F ("mtop", "; m_{Top} [GeV/c^{2}]; Events / (15 GeV/c^{2})", 40, 0.,450.) );
  controlHistos.addHistogram( new TH1F ("ptttbar", "; p_{t#bar{t}} [GeV/c]; Events / (10 GeV/c)", 25, 0.,250.) );
  controlHistos.addHistogram( new TH1F ("afb", "; #Delta #eta(t,#bar{t}) ); Events / (0.1)", 100, -5.,5.) );
  controlHistos.addHistogram( new TH1F ("mttbar", "; Mass(t,#bar{t}) [GeV/c^{2}]; Events / (50 GeV/c^{2})", 40, 0.,2000.) );
  controlHistos.addHistogram( new TH2F ("mtopvsdilmass", "; m_{Top} [GeV/c^{2}]; Mass(l,l') [GeV/c^{2}]; Events", 100, 0.,500.,100, 0.,500.) );
  controlHistos.addHistogram( new TH2F ("mtopvsmlj", "; m_{Top} [GeV/c^{2}]; Mass(l,j) [GeV/c^{2}]; Events", 100, 0.,500.,100, 0.,500.) );
  controlHistos.addHistogram( new TH2F ("mtopvsmet", "; m_{Top} [GeV/c^{2}]; #slash{E}_{T} [GeV/c]; Events", 100, 0.,500.,50, 0.,250.) );
  controlHistos.addHistogram( new TH2F ("mtopvsmttbar", "; m_{Top} [GeV/c^{2}]; Mass(t,#bar{t}) [GeV/c^{2}]; Events", 100, 0.,500.,100, 0.,2000.) );
  controlHistos.addHistogram( new TH2F ("mtopvsafb", "; m_{Top} [GeV/c^{2}]; #Delta #eta(t,#bar{t}) ); Events", 100, 0.,500.,100, -5.,5.) );
  controlHistos.addHistogram( new TH2F ("mttbarvsafb", "; Mass(t,#bar{t}) [GeV/c^{2}];#Delta #eta(t,#bar{t}); Events", 100, 0.,2000.,100,-5.,5.) );
  
  TString cats[]={"all","ee","mumu","emu","etau","mutau"};
  size_t ncats=sizeof(cats)/sizeof(TString);
  TString subcats[]={"","eq0btags","eq1btags","geq2btags"};
  size_t nsubcats=sizeof(subcats)/sizeof(TString);
  for(size_t icat=0; icat<ncats; icat++)
    for(size_t jcat=0; jcat<nsubcats; jcat++)
      controlHistos.initMonitorForStep(cats[icat]+subcats[jcat]);
  

  gSystem->Exec("mkdir -p " + outUrl);
  outUrl += "/";
  outUrl += gSystem->BaseName(evurl);
  TFile *file=TFile::Open(outUrl, "recreate");

  //
  // TMVA
  //  
  TMVA::Factory *tmvaFactory=0;
  TMVA::Reader *tmvaReader=0;
  const unsigned int nVariables = varsList.size()+1;
  std::vector<Double_t> tmvaVarsD( nVariables,0 );
  std::vector<Float_t> tmvaVarsF( nVariables,0 );
  if(useMVA)
    {
      TMVA::Tools::Instance();
      
      if(trainMVA)
	{
	  //create the factory object. 
	  tmvaFactory = new TMVA::Factory( studyTag.Data(), file,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
	 
	  //variables to study
	  TString varsString("");
	  for(std::vector<std::string>::iterator it = varsList.begin(); it != varsList.end(); it++) 
	    {
	      if(it != varsList.begin()) varsString += ":";
	      varsString += *it;
	      tmvaFactory->AddVariable( *it, *it, "", 'F');
	    }
	  tmvaFactory->AddSpectator( "eventCategory" );
	  
	  cout << "==> Start TMVAClassification with " << methodList.size() << " methods and " << nVariables-1 << " variables" << endl;
	}
      else
	{
	  //reader for the methods
	  tmvaReader = new TMVA::Reader( "!Color:!Silent" );
	  
	  //variables to use
	  for(size_t ivar=0; ivar<varsList.size(); ivar++)   tmvaReader->AddVariable( varsList[ivar], &tmvaVarsF[ivar] );
	  tmvaReader->AddSpectator("eventCategory", &tmvaVarsF[varsList.size()]);

	  //read the methods already trained
	  for(size_t imet=0; imet<methodList.size(); imet++)
	    {
	      //open the file with the method description                                                                                                                                                                                                 
	      TString weightFile = weightsDir + "/" + studyTag + ".weights.xml";
	      gSystem->ExpandPathName(weightFile);
	      tmvaReader->BookMVA(methodList[imet], weightFile);
	      TH1 *h=tmva::getHistogramForDiscriminator( methodList[imet] );
	      controlHistos.addHistogram( h );
	    }
	}
    }
  

  //fix entries flag
  ofstream *outf=0;
  if(!isMC) outf=new ofstream("highmassevents.txt",ios::app);

  //process events file
  gSystem->ExpandPathName(evurl);
  TFile *evfile = TFile::Open(evurl);
  if(evfile==0) return -1;
  if(evfile->IsZombie()) return -1;
  EventSummaryHandler evSummaryHandler;
  if( !evSummaryHandler.attachToTree( (TTree *)evfile->Get(dirname) ) ) 
    {
      evfile->Close();
      return -1;
    }  
  TTree *evTree=evSummaryHandler.getTree();

  //init event spy
  EventSummaryHandler *spyEvents=0;
  TFile *spyFile=0;
  TDirectory *spyDir=0;  
  if(saveSummaryTree)
    {
      spyEvents = new EventSummaryHandler;
      spyFile = TFile::Open("EventSummaries.root","UPDATE");
      TString evtag=gSystem->BaseName(evurl);
      evtag.ReplaceAll(".root","");
      spyFile->rmdir(evtag);
      spyDir = spyFile->mkdir(evtag);
      TTree *outT = new TTree("data","Event summary");
      spyEvents->initTree(outT);
    }

  //process kin file
  TString kinUrl(evurl);
  kinUrl.ReplaceAll(".root","/"+kindir);
  gSystem->ExpandPathName(kinUrl);
  cout << "Kin results from " << kinUrl << " to be processed with summary from " << evurl << endl;
  KinResultsHandler kinHandler;
  kinHandler.init(kinUrl,false);
  TChain *t=kinHandler.getResultsChain();
 
  //loop over events
  std::map<TString,int> selEvents;
  for (int inum=0; inum < evTree->GetEntriesFast(); ++inum)
  {
    evTree->GetEvent(inum);
    EventSummary_t &ev = evSummaryHandler.getEvent();
    TString key(""); key+= ev.run; key+="-"; key += ev.lumi; key+="-"; key += ev.event;
    selEvents[key]=inum;
  }

  //loop over kin results
  int nresults(0),neventsused(0);
  int nsigtrain(0), nsigtest(0), nbkgtrain(0), nbkgtest(0);
  for (int inum=0; inum < t->GetEntries(); ++inum){
    t->GetEvent(inum);
 
    //get original event
    Int_t irun,ievent,ilumi;
    kinHandler.getEventInfo(irun,ievent,ilumi);

    TString key("");  key+= irun; key+="-"; key += ilumi;  key+="-"; key += ievent;

    if(selEvents.find(key)==selEvents.end()) continue;
    nresults++;
    evTree->GetEntry( selEvents[key] );

    //get event summary
    EventSummary_t &ev = evSummaryHandler.getEvent();
    if(isMC)
      {
	if(mcTruthMode==1 && !ev.isSignal) continue;
	if(mcTruthMode==2 && ev.isSignal) continue;
      }
    
    std::vector<TString> categs;
    categs.push_back("all");
    if(ev.cat==dilepton::MUMU)  categs.push_back("mumu");
    if(ev.cat==dilepton::EE)  categs.push_back("ee");
    if(ev.cat==dilepton::EMU) categs.push_back("emu");
    if(ev.cat==dilepton::ETAU) categs.push_back("etau");
    if(ev.cat==dilepton::MUTAU) categs.push_back("mutau");
    
    //get particles from the event
    int njets(0),nbtags(0);
    KinCandidateCollection_t leptons, jets, mets,vtx;
    for(Int_t ipart=0; ipart<ev.nparticles; ipart++)
      {
	TLorentzVector p4(ev.px[ipart],ev.py[ipart],ev.pz[ipart],ev.en[ipart]);
	if(isnan(p4.Pt()) || isinf(p4.Pt())) continue;
	switch( ev.id[ipart] )
	  {
	  case 0:
	    mets.push_back( KinCandidate_t(p4,p4.Pt()) );
            break;
          case 1:
            jets.push_back( KinCandidate_t(p4, ev.info1[ipart]) );
	    njets++;
	    if(ev.info1[ipart]>1.7) nbtags++;
            break;
	  case 500:
	    vtx.push_back( KinCandidate_t(p4,p4.Pt()) );
	    break;
          default:
            leptons.push_back( KinCandidate_t(p4,ev.id[ipart]) );
            break;
	  }
      }
    sort(leptons.begin(),leptons.end(),KinAnalysis::sortKinCandidates);
    sort(jets.begin(),jets.end(),KinAnalysis::sortKinCandidates);
    sort(mets.begin(),mets.end(),KinAnalysis::sortKinCandidates);

    TString subcat="eq0btags";
    if(nbtags==1) subcat="eq1btags";
    if(nbtags>=2) subcat="geq2btags";

    //get the combination preferred by KIN
    TH1F *h1=kinHandler.getHisto("mt",1), *h2=kinHandler.getHisto("mt",2);
    h1->Rebin(2); h2->Rebin(2);
    Int_t icomb=(h1->Integral()< h2->Integral())+1;
    TH1F *mpref=kinHandler.getHisto("mt",icomb);
    double mtop = kinHandler.getMPVEstimate(mpref) [1];
    TH1F *mttbarpref=kinHandler.getHisto("mttbar",icomb);
    double mttbar = kinHandler.getMPVEstimate(mttbarpref)[1];
    TH1F *afbpref=kinHandler.getHisto("afb",icomb);
    double afb = kinHandler.getMPVEstimate(afbpref)[1];

    //
    // analyze histos with MVA
    //
    if(useMVA)
      {
	std::map<TH1 *,bool> assignedResults;
	std::map<TH1 *,std::map<int,double> > discriResults;
	//	    assignedResults[h1] = (goodAssignment==1);
	//	    assignedResults[h2] = (goodAssignment==2);
	for(std::map<TH1 *,bool>::iterator resIt = assignedResults.begin();
	    resIt != assignedResults.end();
	    resIt++)
	  {
	    
	    //std::vector<float> histoMomenta = analyzeHistogram(resIt->first);
	    int varCounter(0);
	    for(std::vector<std::string>::iterator it = varsList.begin(); it != varsList.end(); it++)
	      {
		//if(*it=="myvar")        { tmvaVarsF[varCounter++]=histoMomenta[kMyVar]; }
	      }
	    tmvaVarsF[varCounter++] = ev.cat;
	    for(size_t ivar=0; ivar<tmvaVarsF.size(); ivar++)  tmvaVarsD[ivar]=tmvaVarsF[ivar];
	    
	    if(trainMVA)
	      {
		if(resIt->second)
		  {
		    if ( inum%2 == 0 ){ tmvaFactory->AddSignalTrainingEvent( tmvaVarsD,1. ); nsigtrain++; }
		    else              { tmvaFactory->AddSignalTestEvent    ( tmvaVarsD, 1. ); nsigtest++; }
		  }
		else
		  {
		    if ( inum%2 == 0 ){ tmvaFactory->AddBackgroundTrainingEvent( tmvaVarsD, 1. ); nbkgtrain++; }
		    else              { tmvaFactory->AddBackgroundTestEvent    ( tmvaVarsD, 1. ); nbkgtest++; }
		  }
	      }
	    else
	      {
		std::map<int,double> hDiscriResults;
		for(size_t imet=0; imet<methodList.size(); imet++) 
		  hDiscriResults[imet]=tmvaReader->EvaluateMVA( methodList[imet] );
		discriResults[resIt->first] = hDiscriResults;
	      }
	  }
      }
    
    //Compute dilepton/dijet invariant mass
    TLorentzVector dil = leptons[0].first+leptons[1].first;
    float dilmass = dil.M();
    if(fabs(dilmass-91)<15 && (ev.cat==dilepton::EE || ev.cat==dilepton::MUMU))continue;
    double ptlep1(max(leptons[0].first.Pt(),leptons[1].first.Pt())), ptlep2(min(leptons[0].first.Pt(),leptons[1].first.Pt()));    
    TLorentzVector dij = jets[0].first+jets[1].first;
    float mjj=dij.M();
    double ptjet1(max(jets[0].first.Pt(),jets[1].first.Pt())), ptjet2(min(jets[0].first.Pt(),jets[1].first.Pt()));

    TLorentzVector ptttbar=leptons[0].first+leptons[1].first+jets[0].first+jets[1].first+mets[0].first;
    
    //get the lepton-jet pairs
    TLorentzVector lj1=leptons[0].first+jets[icomb==1?0:1].first;
    TLorentzVector lj2=leptons[1].first+jets[icomb==1?1:0].first;

    //ht
    double ht(0);
    for(size_t ijet=0; ijet<jets.size(); ijet++) ht += jets[ijet].first.Pt();
    double sumptlep(leptons[0].first.Pt()+leptons[1].first.Pt()+mets[0].first.Pt());
    double st(sumptlep+mets[0].first.Pt());
    double htlep(st+ht);

    //fill histos
    float weight = ev.weight;

    for(std::vector<TString>::iterator cIt = categs.begin(); cIt != categs.end(); cIt++)
      {
	if(mtop>0)
	  {
// 	    results[*cIt+"_njets"]->Fill(njets,weight);
// 	    results[*cIt+"_btags"]->Fill(nbtags,weight);
// 	    results[*cIt+"_leadjet"]->Fill(ptjet1,weight);
// 	    results[*cIt+"_subleadjet"]->Fill(ptjet2,weight);
// 	    results[*cIt+"_leadlepton"]->Fill(ptlep1,weight);
// 	    results[*cIt+"_subleadlepton"]->Fill(ptlep2,weight);
// 	    results[*cIt+"_met"]->Fill(mets[0].first.Pt(),weight);
// 	    results[*cIt+"_ht"]->Fill(ht,weight);
// 	    results[*cIt+"_st"]->Fill(st,weight);
// 	    results[*cIt+"_sumpt"]->Fill(sumptlep,weight);
// 	    results[*cIt+"_htlep"]->Fill(htlep,weight);
// 	    results[*cIt+"_ptttbar"]->Fill(ptttbar.Pt(),weight);
	    
// 	    results[*cIt+subcat+"_mtop"]->Fill(mtop,weight);
// 	    results[*cIt+"_mtop"]->Fill(mtop,weight);
// 	    results[*cIt+subcat+"_dilmass"]->Fill(dilmass,weight);
// 	    ((TH2F *)results[*cIt+"_mtopvsdilmass"])->Fill(mtop,dilmass,weight);
// 	    ((TH2F *)results[*cIt+"_mtopvsmlj"])->Fill(mtop,lj1.M(),weight);
// 	    ((TH2F *)results[*cIt+"_mtopvsmlj"])->Fill(mtop,lj2.M(),weight);
// 	    ((TH2F *)results[*cIt+"_mtopvsmet"])->Fill(mtop,mets[0].first.Pt(),weight);
// 	    ((TH2F *)results[*cIt+"_mtopvsmttbar"])->Fill(mtop,mttbar,weight);
// 	    ((TH2F *)results[*cIt+"_mtopvsafb"])->Fill(mtop,afb,weight);
// 	    ((TH2F *)results[*cIt+"_mttbarvsafb"])->Fill(mttbar,afb,weight);
// 	    results[*cIt+"_afb"]->Fill(afb);
// 	    results[*cIt+"_mttbar"]->Fill(mttbar);
	  }
      }
    neventsused++;

    //save for further study
    //if(mtop>0 && spyEvents && ev.normWeight==1)
    if(mtop>0 && spyEvents && ev.weight>1.7)
      {
	std::vector<float> measurements;
	measurements.push_back(mtop);
	measurements.push_back(mttbar);
	measurements.push_back(afb);
	measurements.push_back(ptttbar.Pt());
	measurements.push_back(nbtags);
	measurements.push_back(njets);
	measurements.push_back(htlep);
	spyEvents->fillTreeWithEvent( ev, measurements );
      }

    //for data only
    if (!isMC) 
      *outf << "| " << irun << ":" << ilumi << ":" << ievent 
	    << " | " << categs[1] 
	    << " | " << mtop 
	    << " | " << mttbar 
	    << " | " << ptlep1 << ";" << ptlep2  << " | " << dilmass
	    << " | " << ptjet1 << ";" << ptjet2  << " | " << mjj
	    << " | " << mets[0].first.Pt() << " | " << htlep << endl; 
  }
  kinHandler.end();

  if(useMVA)
    {
      if(trainMVA)
	{
	  cout << "\t Good assignments split in " << nsigtrain << " training events + " << nsigtest  << " test events" << endl
	       << "\t Wrong assignments split in " << nbkgtrain << " training events + " << nbkgtest  << " test events" << endl;

	  //prepare the training
	  tmvaFactory->PrepareTrainingAndTestTree("","","nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
	  for(size_t i=0; i<methodList.size(); i++)  tmva::bookMethod(tmvaFactory,methodList[i]);

	  // now train, test and evaluate                 
	  std::cout << "\t Start training" << std::endl;
	  tmvaFactory->TrainAllMethods();
	  std::cout << "\t Start test" << std::endl;
	  tmvaFactory->TestAllMethods();
	  std::cout << "\t Start evaluate" << std::endl;
	  tmvaFactory->EvaluateAllMethods();

	  std::cout << "\t Results available in: " << file->GetName() << std::endl;
	  std::cout << "\t runMVAStudy is done!" << std::endl;
	  
	  delete tmvaFactory;
	}
    }

  if(!isMC) 
    {
      outf->close(); 
      delete outf;
    }


  //if MC: rescale to number of selected events and to units of pb
  cout << "From " << selEvents.size() << "original events found " << nresults << " kin results - used " << neventsused << endl; 

  float cnorm=1;
  if(isMC && nresults)
    {
      double scaleFactor=double(selEvents.size())/double(nresults);
      TH1F *cutflowH = (TH1F *) evfile->Get("evAnalyzer/top/cutflow");
      if(cutflowH)
	{
	  cnorm=cutflowH->GetBinContent(1);
	  if(cnorm>0) scaleFactor/=cnorm;
	}
      cout << selEvents.size() << " " << nresults << " " << scaleFactor << endl;
      //      for(std::map<TString,TH1 *>::iterator hIt = results.begin(); hIt != results.end(); hIt++) hIt->second->Scale(scaleFactor);
    }


  //save to file
  TDirectory *baseOutDir=file->mkdir("localAnalysis");
  SelectionMonitor::StepMonitor_t &mons=controlHistos.getAllMonitors();
  std::map<TString, TDirectory *> outDirs;
  outDirs["all"]=baseOutDir->mkdir("all");
  outDirs["ee"]=baseOutDir->mkdir("ee");
  outDirs["emu"]=baseOutDir->mkdir("emu");
  outDirs["mumu"]=baseOutDir->mkdir("mumu");
  for(SelectionMonitor::StepMonitor_t::iterator it =mons.begin(); it!= mons.end(); it++)
    {
      TString icat("all");
      if(it->first.BeginsWith("ee")) icat="ee";
      if(it->first.BeginsWith("emu")) icat="emu";
      if(it->first.BeginsWith("mumu")) icat="mumu";
      outDirs[icat]->cd();
      for(SelectionMonitor::Monitor_t::iterator hit=it->second.begin(); hit!= it->second.end(); hit++)
        {
	  if(isMC && cnorm>0) hit->second->Scale(1./cnorm);
	  if( !((TClass*)hit->second->IsA())->InheritsFrom("TH2")
	      && !((TClass*)hit->second->IsA())->InheritsFrom("TGraph") )
	    fixExtremities(hit->second,true,true);
	  hit->second->Write();

        }
    }
  file->Close(); 


  if(spyEvents)
    {
      spyDir->cd();
      spyEvents->getTree()->Write();
      spyFile->Write();
      spyFile->Close();
    }
}  
