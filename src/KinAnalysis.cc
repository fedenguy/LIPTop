#include "LIP/Top/interface/KinAnalysis.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"

using namespace std;

//
KinAnalysis::KinAnalysis(TString &scheme,int maxTries, int maxJetMult,float mw, float mb)
  : scheme_(scheme),
    maxTries_(maxTries),
    maxJetMult_(maxJetMult),
    mw_(mw),
    mb_(mb)
{
  //seed
  TTimeStamp timeStamp;
  UInt_t theSeed = timeStamp.GetNanoSec() % timeStamp.GetSec() + timeStamp.GetDate();
  rndGen_.SetSeed(theSeed);
  gRandom->SetSeed(theSeed & 0xfd1c);

  //pz parameterization
  deltaPzFunc_ = new TF1("dpzfunc","gaus(0)+gaus(3)",-2000,2000);
  deltaPzFunc_->FixParameter(0,2.79);  deltaPzFunc_->FixParameter(3,1.53);
  deltaPzFunc_->FixParameter(1,0);     deltaPzFunc_->FixParameter(4,0);
  deltaPzFunc_->FixParameter(2,323);   deltaPzFunc_->FixParameter(5,690); 

  //check variations
  if( scheme_.Contains("mdpz") )   { deltaPzFunc_->FixParameter(2,323*0.5);   deltaPzFunc_->FixParameter(5,690*0.5); }
  if( scheme_.Contains("pdpz") )   { deltaPzFunc_->FixParameter(2,323*2);   deltaPzFunc_->FixParameter(5,690*2); }
  if( scheme_.Contains("pmw") )    mw_ += 0.025;
  if( scheme_.Contains("mmw") )    mw_ -= 0.025;
}

//
void KinAnalysis::runOn(EventSummary_t &ev, JetResolution *ptResol, JetResolution *etaResol, JetResolution *phiResol, JetCorrectionUncertainty *jecUnc)
{
  try{
    
    KinCandidateCollection_t leptons, jets, mets;
    for(Int_t ipart=0; ipart<ev.nparticles; ipart++)
      {
	TLorentzVector p4(ev.px[ipart],ev.py[ipart],ev.pz[ipart],ev.en[ipart]);
	switch( ev.id[ipart] )
	  {
	  case 0:
	    mets.push_back( KinCandidate_t(p4,p4.Pt()) );
	    break;
	  case 1:
	    jets.push_back( KinCandidate_t(p4, ev.info1[ipart]) );
	    break;
	  default:
	    leptons.push_back( KinCandidate_t(p4,p4.Pt()) );
	    break;
	  }
      }
    sort(leptons.begin(),leptons.end(),KinAnalysis::sortKinCandidates);
    sort(jets.begin(),jets.end(),KinAnalysis::sortKinCandidates);
    sort(mets.begin(),mets.end(),KinAnalysis::sortKinCandidates);


    //set to the base values
    int nComb=0;
    resetHistos();
 
    //jet energy scale
    double jet1Scale(1.0), jet2Scale(1.0);
    for(int ijet=0; ijet<maxJetMult_; ijet++) 
      {
	for(int jjet=0; jjet<maxJetMult_; jjet++)
	  {
	    if(ijet==jjet) continue;
	    nComb++;
		  
	    cout << "\t - Comb #" << nComb << " (" << leptons[0].first.Pt() << ";" << jets[ijet].first.Pt() << ":" << jets[ijet].second << ")"
		 << " (" << leptons[1].first.Pt() << ";" << jets[jjet].first.Pt() << ":" << jets[jjet].second << ")" << endl;

	    for(int ivar=0; ivar<3; ivar++)
	      {
		if(ivar>0)
		  {
		    jecUnc->setJetEta(jets[ijet].first.Eta());
		    jecUnc->setJetPt(jets[ijet].first.Pt());
		    jet1Scale = 1.0 + (jet1Scale<1?-1:+1)*jecUnc->getUncertainty(true);
		    
		    jecUnc->setJetEta(jets[jjet].first.Eta());
		    jecUnc->setJetPt(jets[jjet].first.Pt());
		    jet2Scale = 1.0 + (jet2Scale<1?-1:+1)*jecUnc->getUncertainty(true);
		  }
				  
		for(int itry=1; itry<maxTries_; itry++)
		  {		  
		    //leptons
		    TLorentzVector pl1 = leptons[0].first;
		    TLorentzVector pl2 = leptons[1].first;

		    //jets;
		    double deltaPz = deltaPzFunc_->GetRandom();
		    float ptScaleRes = (ptResol->resolutionEtaPt(jets[ijet].first.Eta(),jets[ijet].first.Pt())->GetRandom()-1.0);
		    float etaRes = etaResol->resolutionEtaPt(jets[ijet].first.Eta(),jets[ijet].first.Pt())->GetRandom();
		    float phiRes = phiResol->resolutionEtaPt(jets[ijet].first.Eta(),jets[ijet].first.Pt())->GetRandom();
		    float newpt = (1.0+ptScaleRes)*jets[ijet].first.Pt();
		    float neweta = etaRes+jets[ijet].first.Eta();
		    float newphi = phiRes+jets[ijet].first.Phi();
		    TLorentzVector pb1;
		    pb1.SetPtEtaPhiM(newpt,neweta,newphi,mb_);
		    pb1 *= jet1Scale;

		    ptScaleRes = (ptResol->resolutionEtaPt(jets[jjet].first.Eta(),jets[jjet].first.Pt())->GetRandom()-1.0);
		    etaRes = etaResol->resolutionEtaPt(jets[jjet].first.Eta(),jets[jjet].first.Pt())->GetRandom();
		    phiRes = phiResol->resolutionEtaPt(jets[jjet].first.Eta(),jets[jjet].first.Pt())->GetRandom();
		    newpt = (1.0+ptScaleRes)*jets[jjet].first.Pt();
		    neweta = etaRes+jets[jjet].first.Eta();
		    newphi = phiRes+jets[jjet].first.Phi();
		    TLorentzVector pb2;
		    pb2.SetPtEtaPhiM(newpt,neweta,newphi,mb_);
		    pb2 *= jet2Scale;
		    
		    //MET (must correct for jet energy scale/resolution smearing)
		    TLorentzVector met(mets[0].first);
		    float metResol= 1.0+rndGen_.Gaus(0,0.1);
		    double dPhiMET = rndGen_.Gaus(0,0.1);
		    float metx = metResol*( met.Px()-(pb1.Px()-jets[ijet].first.Px())-(pb2.Px()-jets[jjet].first.Px())-(pl1.Px()-leptons[0].first.Px())-(pl2.Px()-leptons[1].first.Px()) );
		    float mety = metResol*( met.Py()-(pb1.Py()-jets[ijet].first.Py())-(pb2.Py()-jets[jjet].first.Py())-(pl1.Py()-leptons[0].first.Py())-(pl2.Py()-leptons[1].first.Py()) );
		    TVector3 metConstraint( metx, mety, deltaPz-pb1.Pz()-pl1.Pz()-pb2.Pz()-pl2.Pz());
		    metConstraint.RotateZ(dPhiMET);
		      
		    //prevent strange values
		    if(pl1.Pt()<1 || pb1.Pt()<1 || pl2.Pt()<1 || pb2.Pt()<1 || metConstraint.Pt()<1) continue;
		      
		    TTbarSolutionCollection_t sols = kin_.findSolutions(pl1,pb1,pl2,pb2,metConstraint,mw_);    
		    if(sols.size()==0) continue;
		    TTbarSolution_t *sol = &(sols.back());
		    TLorentzVector ttbar = sol->pt1+sol->pt2;
		    float avgMtop = (sol->pt1.M()+sol->pt2.M())*0.5;
		    float mttbar = ttbar.M();
		    std::vector<double> mt2 = getMT2( *sol );
		    float afb = sol->pt1.Eta()-sol->pt2.Eta();
			  
		    getHisto("mt", nComb)->Fill( avgMtop );
		    getHisto("mttbar",nComb)->Fill(mttbar);
		    getHisto("minmt2",nComb)->Fill(mt2[1]);
		    getHisto("maxmt2",nComb)->Fill(mt2[0]);
		    getHisto("afb",nComb)->Fill(afb);
		  }
	      }
	  }
      }
  }
  catch(std::exception &e){
    cout << e.what() << endl;
  }
}


//
void KinAnalysis::resetHistos()
{
  //  for(std::map<KinHistoKey,TH1 *>::iterator it = kinHistos_.begin();
  //      it != kinHistos_.end(); it++)
  //    it->second->Reset("ICE");
}

//
TH1 *KinAnalysis::getHisto(TString var, int nComb)
{
  //  KinHistoKey key(var,nComb);
  //  if(kinHistos_.find(key)==kinHistos_.end()) return 0;
  //  return kinHistos_[key];
  return 0;
}

//
void KinAnalysis::bookHistos()
{
  int ncombs=maxJetMult_*(maxJetMult_-1);
  for(int icomb=1; icomb<=ncombs; icomb++)
    {
      TString cat(""); cat += icomb;
      TString title("Combination #"); title+=cat;
      /*
      kinHistos_[KinHistoKey("mt",icomb)]=  new TH1D("mt_"+cat,title+";M_{t} [GeV/c^{2}];N_{solutions} / (2 GeV/c^{2})",200,100,500);
      kinHistos_[KinHistoKey("mttbar",icomb)] = new TH1D("mttbar_"+cat,title+";M_{t#bar{t}} [GeV/c^{2}];N_{solutions} / (10 GeV/c^{2})", 200,0,2000);
      kinHistos_[KinHistoKey("minmt2",icomb)] = new TH1D("minmt2_"+cat,title+";min M_{T2} [GeV/c^{2}];N_{solutions} / (2 GeV/c^{2})",200,0,400);
      kinHistos_[KinHistoKey("maxmt2",icomb)] = new TH1D("maxmt2_"+cat,title+";max M_{T2} [GeV/c^{2}];N_{solutions} / (2 GeV/c^{2})",200,0,400);
      kinHistos_[KinHistoKey("afb",icomb)] = new TH1D("afb_"+cat,title+";A_{fb};N_{solutions} / (0.05)",100,-4.95,5.05);
      */
    }

  //  for(std::map<KinHistoKey,TH1 *>::iterator it = kinHistos_.begin();
  //      it != kinHistos_.end(); it++)  
  //    formatPlot( it->second, 1,1,1,20,0,true,true,1,1,1);
}
