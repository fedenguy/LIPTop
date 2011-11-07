#include "LIP/Top/interface/KinAnalysis.h"
#include "TRandom.h"

using namespace std;

//
TLorentzVectorCollection randomlyRotate(TLorentzVectorCollection &leptons, TLorentzVectorCollection &jets)
{
  //get rotated leptons
  TLorentzVectorCollection rotLeptons;

  //create a rotated copy of all leptons
  for(TLorentzVectorCollection::iterator lit = leptons.begin(); lit != leptons.end(); lit++)
    {
      do
	{
	  //rotate lepton
	  double en    = lit->E();
	  double pabs  = lit->P();
	  double phi   = gRandom->Uniform(0,2*TMath::Pi());
	  double theta = TMath::ACos( gRandom->Uniform(-1,1) );
	  TLorentzVector rotLep =  TLorentzVector( pabs*TMath::Cos(phi)*TMath::Sin(theta),
						   pabs*TMath::Sin(phi)*TMath::Sin(theta),
						   pabs*TMath::Cos(theta),
						   en);

	  //require selectable kinematics
	  if( TMath::Abs(rotLep.Eta())>2.4 || rotLep.Pt()<20 ) continue;  

	  //require object separation
	  double minDR(1000);
	  for(TLorentzVectorCollection::iterator jit = jets.begin(); jit != jets.end(); jit++)
	    {
	      double dR = jit->DeltaR(rotLep);
	      if(dR>minDR) continue;
	      minDR=dR;
	    }
	  if(minDR<0.4) continue;
	  
	  //save lepton
	  rotLeptons.push_back(rotLep);
	  break;
	} while( 1 );
    }

  //all done
  return rotLeptons;
} 

//
KinAnalysis::KinAnalysis(TString &scheme,int maxTries, int maxJetMult,float mw, float mb, TString outpath, bool doWrite)
  : scheme_(scheme),
    maxTries_(maxTries),
    maxJetMult_(maxJetMult),
    mw_(mw),
    mb_(mb)
{
  //init the handler
  resHandler_.init(outpath,doWrite,maxJetMult);

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
void KinAnalysis::runOn(top::EventSummary_t &ev, JetResolution *ptResol, JetResolution *etaResol, JetResolution *phiResol, JetCorrectionUncertainty *jecUnc)
{
  try{
    
    KinCandidateCollection_t leptons, jets, mets;
    TLorentzVectorCollection leptonsp4,jetsp4;
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
	    jetsp4.push_back(p4);
	    break;
	  default:
	    leptons.push_back( KinCandidate_t(p4,ev.id[ipart]) );
	    leptonsp4.push_back(p4);
	    break;
	  }
      }

    //random rotation of leptons
    if(scheme_=="randrot") 
      {
	TLorentzVector deltaLep(0,0,0,0);
	leptonsp4 = randomlyRotate(leptonsp4,jetsp4);
	for(size_t ilep=0; ilep<leptons.size(); ilep++) 
	  {
	    deltaLep += leptons[ilep].first-leptonsp4[ilep];
	    leptons[ilep].first = leptonsp4[ilep];
	  }
	mets[0].first += deltaLep;
      }
    
    //order collections
    sort(leptons.begin(),leptons.end(),KinAnalysis::sortKinCandidates);
    sort(jets.begin(),jets.end(),KinAnalysis::sortKinCandidates);
    sort(mets.begin(),mets.end(),KinAnalysis::sortKinCandidates);
    if(leptons.size()<2 || jets.size()<2 || mets.size()<1) return;


    //debug
    cout << "[KinAnalysis][runOn] " << ev.run << " : " << ev.lumi << " : " << ev.event << endl
	 << "Scheme is: " << scheme_ << endl
	 << "Leptons #1 : (" << leptons[0].first.Pt() << ";" << leptons[0].first.Eta() << ";" << leptons[0].first.Phi() << ") q:" << leptons[0].second << endl  
	 << "        #2 : (" << leptons[1].first.Pt() << ";" << leptons[1].first.Eta() << ";" << leptons[1].first.Phi() << ") q:" << leptons[1].second << endl  
	 << "Jets    #1 : (" << jets[0].first.Pt() << ";" << jets[0].first.Eta() << ";" << jets[0].first.Phi() << ") btag:" << jets[0].second  << endl  
	 << "        #2 : (" << jets[1].first.Pt() << ";" << jets[1].first.Eta() << ";" << jets[1].first.Phi() << ") btag:" << jets[1].second << endl  
	 << "MET        : (" << mets[0].first.Pt() << ";" << mets[0].first.Phi() << ")" << endl;
    
    //set to the base values
    int nComb=0;
    resHandler_.resetHistos();
 
    //jet energy scale
    double jet1Scale(1.0), jet2Scale(1.0);
    for(int ijet=0; ijet<maxJetMult_; ijet++) 
      {
	for(int jjet=0; jjet<maxJetMult_; jjet++)
	  {
	    if(ijet==jjet) continue;
	    nComb++;
		  
// 	    for(int ivar=0; ivar<3; ivar++)
// 	      {
// 		if(ivar>0)
// 		  {
// 		    jecUnc->setJetEta(jets[ijet].first.Eta());
// 		    jecUnc->setJetPt(jets[ijet].first.Pt());
// 		    jet1Scale = 1.0 + (jet1Scale<1?-1:+1)*jecUnc->getUncertainty(true);
		    
// 		    jecUnc->setJetEta(jets[jjet].first.Eta());
// 		    jecUnc->setJetPt(jets[jjet].first.Pt());
// 		    jet2Scale = 1.0 + (jet2Scale<1?-1:+1)*jecUnc->getUncertainty(true);
// 		  }
				  
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

		    //compute the full kinematics obtained
		    TTbarSolution_t *sol = &(sols.back());
		    TLorentzVector ttbar = sol->pt1+sol->pt2;
		    float avgMtop = (sol->pt1.M()+sol->pt2.M())*0.5;
		    float mttbar = ttbar.M();
		    std::vector<double> mt2 = getMT2( *sol );
		    float afb = sol->pt1.Eta()-sol->pt2.Eta();

		    //fill histos
		    resHandler_.getHisto("mt", nComb)->Fill( avgMtop );
		    resHandler_.getHisto("mttbar",nComb)->Fill(mttbar);
		    resHandler_.getHisto("mt2",nComb)->Fill(mt2[0]);
		    resHandler_.getHisto("afb",nComb)->Fill(afb);
		  }
// 	      }
 	  }
      }
    
    //save resuls
    resHandler_.addResults( ev );
  }
  catch(std::exception &e){
    cout << e.what() << endl;
  }
}



//
KinAnalysis::~KinAnalysis()
{
}
