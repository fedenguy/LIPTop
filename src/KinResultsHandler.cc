#include "LIP/Top/interface/KinResultsHandler.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "TSystem.h"

using namespace std;
using namespace top;

//
KinResultsHandler::KinResultsHandler()
  : kinFile_(0),  kinTree_(0), kinChain_(0)
{
}

//
void KinResultsHandler::init(TString outpath,bool doWrite, int maxJetMult) 
{
  doWrite_=doWrite;

  typedef std::pair<TString, int> KinHistoKey;

  //allocate memory
  int ncombs=maxJetMult*(maxJetMult-1);

  //open file
  if(!doWrite_) 
    {
      //kinFile_ = TFile::Open(outpath);
      //kinTree_ = (TTree *) kinFile_->Get("kin");

      kinChain_ = new TChain("kin");
      if(!outpath.Contains(".root"))
	{
	  std::map<int,TString> files;
	  void *dirp=gSystem->OpenDirectory(outpath);
	  while(dirp){
	    TString entry(gSystem->GetDirEntry(dirp));
	    if(entry=="") break;
	    if(!entry.Contains(".root")) continue;
	    int ifile;
	    sscanf(entry.Data(),"KinAnalysis_%d.root",&ifile);
	    files[ifile]=outpath+"/"+entry;
	  }
	  for(std::map<int,TString>::iterator it = files.begin(); it != files.end(); it++) kinChain_->AddFile(it->second);
	}
      else  kinChain_->AddFile(outpath);
            
      TObjArray *branches = kinChain_->GetListOfBranches();
      for(int ibranch=0; ibranch<branches->GetEntriesFast(); ibranch++)
	{
	  TBranch *br = (TBranch *) branches->At(ibranch);
	  TString name(br->GetName());
	  
	  if( name=="run" )        kinChain_->SetBranchAddress(name,&iRun_);
	  else if( name=="event" ) kinChain_->SetBranchAddress(name,&iEvent_);
	  else if( name=="lumi" )  kinChain_->SetBranchAddress(name,&iLumi_);
	  else
	    {
	      TObjArray *tkns=name.Tokenize("_");
	      KinHistoKey key( ((TObjString *)tkns->At(0))->GetString(),
			       ((TObjString *)tkns->At(1))->GetString().Atoi() );
	      kinHistos_[ key ] = 0;
	      kinChain_->SetBranchAddress(name, &kinHistos_[key] );
	    }
	}     
    } 
  else
    {
      kinFile_ = TFile::Open(outpath, doWrite ? "RECREATE" : "");
      kinFile_->SetCompressionLevel( 9 );
      kinTree_ = new TTree( "kin","Kinematics analysis of top dilepton events" );
      kinTree_->SetAutoSave();
      kinTree_->Branch("run",  &iRun_, "run/I");
      kinTree_->Branch("lumi", &iLumi_, "lumi/I");
      kinTree_->Branch("event", &iEvent_, "event/I");
      for(int icomb=1; icomb<=ncombs; icomb++)
	{
	  TString cat(""); cat += icomb;
	  TString title("Combination #"); title+=cat;
	  
	  kinHistos_[KinHistoKey("mt",icomb)] = new TH1F("mt_"+cat,title+";M_{t} [GeV/c^{2}];N_{solutions} / (5 GeV/c^{2})",400,0,2000);
	  kinHistos_[KinHistoKey("mttbar",icomb)]  = new TH1F("mttbar_"+cat,title+";M_{t#bar{t}} [GeV/c^{2}];N_{solutions} / (20 GeV/c^{2})", 250,0,5000);   
	  kinHistos_[KinHistoKey("mt2",icomb)] = new TH1F("mt2_"+cat,title+";M_{T2} [GeV/c^{2}];N_{solutions} / (5 GeV/c^{2})",100,0,500);	 
	  kinHistos_[KinHistoKey("afb",icomb)] = new TH1F("afb_"+cat,title+";A_{fb};N_{solutions} / (0.05)",100,-4.95,5.05);
	}

      for(std::map<KinHistoKey,TH1F *>::iterator it = kinHistos_.begin(); it != kinHistos_.end(); it++)
	{
	  it->second->Sumw2();
	  kinTree_->Branch(it->second->GetName(),"TH1F",&it->second,32000,0);
	}
    }
  //500
  fitFunc_ =  new TF1("fitFunc","gaus",0,1000);
}

//
void KinResultsHandler::addResults(EventSummary_t &ev)
{
  iRun_=ev.run;
  iLumi_=ev.lumi;
  iEvent_=ev.event;
  kinTree_->Fill();
}

//
void KinResultsHandler::end()
{
  if(kinTree_==0 || kinFile_==0) return;

  if(doWrite_)
    {
      kinFile_->cd();
      kinTree_->Print();
      kinFile_->Write();
    }
  kinFile_->Close();
}

//    
std::vector<double> KinResultsHandler::getMPVEstimate(TH1 *h)
{
  std::vector<double> res(4,0);
  if(h==0) return res;
  if(h->Integral()==0) return res;

  //fit a gaussian near the most probable value
  Int_t iBin = h->GetMaximumBin();	  
  double mpv = h->GetXaxis()->GetBinCenter(iBin);
  fitFunc_->SetRange(mpv-25,mpv+25);
  fitFunc_->SetParLimits(1,mpv-10,mpv+10);
  h->Fit(fitFunc_,"LRQN");
  for(size_t iparam=0; iparam<3; iparam++) res[iparam]= fitFunc_->GetParameter(iparam);
  if(fitFunc_->GetNDF()>0) res[3] = fitFunc_->GetChisquare()/fitFunc_->GetNDF();
  return res;
}

//
void KinResultsHandler::resetHistos()
{
  for(std::map<std::pair<TString, int>,TH1F *>::iterator it = kinHistos_.begin();
      it != kinHistos_.end(); it++)
    it->second->Reset("ICE");
}

//
TH1F *KinResultsHandler::getHisto(TString var, int nComb)
{
  std::pair<TString, int> key(var,nComb);
  if(kinHistos_.find(key)==kinHistos_.end()) return 0;
  return kinHistos_[key];
}


//
KinResultsHandler::~KinResultsHandler()
{
//   for(std::map<std::pair<TString, int>,TH1F *>::iterator it = kinHistos_.begin();
//       it != kinHistos_.end(); it++)
//     it->second->Delete();
}
