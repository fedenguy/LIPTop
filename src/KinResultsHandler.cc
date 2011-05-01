#include "LIP/Top/interface/KinResultsHandler.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"

using namespace std;

//
KinResultsHandler::KinResultsHandler()
  : kinFile_(0),  kinTree_(0)
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
  kinFile_ = TFile::Open(outpath, doWrite ? "RECREATE" : "");

  if(!doWrite_) 
    {
      kinTree_ = (TTree *) kinFile_->Get("kin"); 
      TObjArray *branches = kinTree_->GetListOfBranches();
      for(int ibranch=0; ibranch<branches->GetEntriesFast(); ibranch++)
	{
	  TBranch *br = (TBranch *) branches->At(ibranch);
	  TString name(br->GetName());
	  if( name=="run" )   br->SetAddress(&iRun_);
	  else if( name=="event" ) br->SetAddress(&iEvent_);
	  else if( name=="lumi" )  br->SetAddress(&iLumi_);
	  else
	    {
	      TObjArray *tkns=name.Tokenize("_");
	      KinHistoKey key( ((TObjString *)tkns->At(0))->GetString(),
			       ((TObjString *)tkns->At(1))->GetString().Atoi() );
	      kinHistos_[ key ] = 0;
	      br->SetAddress( &kinHistos_[key] );
	    }
	}     
    } 
  else
    {
      kinFile_->SetCompressionLevel( 9 );
      kinTree_ = new TTree( "kin","Kinematics analysis of top dilepton events" );
      kinTree_->SetAutoSave();
      kinTree_->Branch("run",  &iRun_, "run/I");
      kinTree_->Branch("lumi", &iEvent_, "lumi/I");
      kinTree_->Branch("event", &iLumi_, "event/I");
      for(int icomb=1; icomb<=ncombs; icomb++)
	{
	  TString cat(""); cat += icomb;
	  TString title("Combination #"); title+=cat;
	  
	  kinHistos_[KinHistoKey("mt",icomb)] = new TH1F("mt_"+cat,title+";M_{t} [GeV/c^{2}];N_{solutions} / (2 GeV/c^{2})",200,100,500);
	  kinHistos_[KinHistoKey("mttbar",icomb)]  = new TH1F("mttbar_"+cat,title+";M_{t#bar{t}} [GeV/c^{2}];N_{solutions} / (10 GeV/c^{2})", 200,0,2000);   
	  kinHistos_[KinHistoKey("mt2",icomb)] = new TH1F("mt2_"+cat,title+";M_{T2} [GeV/c^{2}];N_{solutions} / (2 GeV/c^{2})",200,0,400);	 
	  kinHistos_[KinHistoKey("afb",icomb)] = new TH1F("afb_"+cat,title+";A_{fb};N_{solutions} / (0.05)",100,-4.95,5.05);
	}

      for(std::map<KinHistoKey,TH1F *>::iterator it = kinHistos_.begin(); it != kinHistos_.end(); it++)
	{
	  it->second->Sumw2();
	  kinTree_->Branch(it->second->GetName(),"TH1F",&it->second,32000,0);
	}
    }
  
  fitFunc_ =  new TF1("fitFunc","gaus",0,500);
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

  //close file and delete allocated memory
  if(doWrite_)
    {
      kinFile_->cd();
      //kinTree_->SetDirectory(kinFile_);
      kinTree_->Print();
      kinFile_->Write();
      //kinTree_->Write();
    }
  kinFile_->Close();
  //kinFile_->Delete();
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
