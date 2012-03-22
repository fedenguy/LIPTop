#include "LIP/Top/interface/KinAnalysis.h"
#include "CMGTools/HtoZZ2l2nu/interface/MacroUtils.h"
#include "TRandom.h"

using namespace std;

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
void KinAnalysis::runOn(top::EventSummary_t &ev, JetResolution *ptResol, JetResolution *etaResol, JetResolution *phiResol, JetCorrectionUncertainty *jecUnc, bool isMC)
{
  try{

    top::PhysicsEvent_t phys = getPhysicsEventFrom(ev);
    
    //dilepton
    KinCandidateCollection_t leptons;
    LorentzVectorCollection finalLeptonsP4;
    LorentzVector deltaLep(0,0,0,0);
    for(size_t ilep=0; ilep<phys.leptons.size(); ilep++)
      {
	double les(1.0);
	LorentzVector p4(phys.leptons[ilep]);
	if( (scheme_=="lesup" || scheme_=="lesdown") && fabs(phys.leptons[ilep].id)==11)
	  {
	    bool isEB( fabs(p4.Eta()) < 1.442 );
	    if(scheme_=="lesup")   les = isEB ? 1.005 : 1.025;
	    if(scheme_=="lesdown") les = isEB ? 0.995 : 0.975;
	  }
	p4 *= les; 
	deltaLep += (p4-phys.leptons[ilep]);
	finalLeptonsP4.push_back(p4);
	
	if(p4.pt()<20 || fabs(p4.eta())>2.5) continue;
	leptons.push_back( KinCandidate_t(p4,phys.leptons[ilep].id) );
      }
    if(leptons.size()<2) return;

    //jets
    std::vector<LorentzVector> origJetsP4; for(size_t ijet=0; ijet<phys.jets.size(); ijet++) origJetsP4.push_back(phys.jets[ijet]);
    std::vector<LorentzVector> finalJetsP4(origJetsP4);
    if(scheme_=="jesup" || scheme_=="jesdown")
      {
	std::vector<LorentzVectorCollection> jetsVar;
	LorentzVectorCollection metsVar;
	jet::computeVariation(origJetsP4,phys.met,jetsVar,metsVar,ptResol,etaResol,phiResol,jecUnc);
	if(scheme_=="jesup")   { for(size_t ijet=0; ijet<origJetsP4.size(); ijet++) finalJetsP4[ijet]=jetsVar[jet::JES_UP][ijet];   }
	if(scheme_=="jesdown") { for(size_t ijet=0; ijet<origJetsP4.size(); ijet++) finalJetsP4[ijet]=jetsVar[jet::JES_DOWN][ijet]; }
      }
    else if(scheme_=="jesWup" || scheme_=="jesWdown" || scheme_=="jesW")
      {
	//cf. http://cms-physics.web.cern.ch/cms-physics/public/TOP-11-015-pas.pdf
	int shiftSign=0;
	if(scheme_=="jesWup") shiftSign=1;
	if(scheme_=="jesWdown") shiftSign=-1;
	float baseJetScale(isMC ? 1.0 :  1.004);
	float jetScale=baseJetScale+shiftSign*(sqrt(pow(0.005,2)+pow(0.012,2)));
	for(size_t ijet=0; ijet<origJetsP4.size(); ijet++) finalJetsP4[ijet]=jetScale*origJetsP4[ijet];
      }
    LorentzVector deltaJets(0,0,0,0);
    KinCandidateCollection_t jets;
    for(size_t ijet=0; ijet<phys.jets.size(); ijet++)
      {
	deltaJets += finalJetsP4[ijet]-origJetsP4[ijet];
	if(finalJetsP4[ijet].pt()<30 || fabs(finalJetsP4[ijet].eta())>2.5) continue;
	jets.push_back(KinCandidate_t(finalJetsP4[ijet],phys.jets[ijet].btag7));
      }
    if(jets.size()<2) return;
    

    //met
    LorentzVector finalMet=phys.met-deltaLep-deltaJets;
    LorentzVector clusteredP4(0,0,0,0);
    for(size_t ilep=0; ilep<finalLeptonsP4.size(); ilep++) clusteredP4 += finalLeptonsP4[ilep];
    for(size_t ijet=0; ijet<finalJetsP4.size(); ijet++)    clusteredP4 += finalJetsP4[ijet];
    LorentzVector unclusteredP4(finalMet+clusteredP4);
    if(scheme_=="unclustup" || scheme_=="unclustdown")
      {
	float sf(scheme_=="unclustdown" ? 0.9 : 1.1);
	unclusteredP4 *= sf;
      }
    finalMet =unclusteredP4-clusteredP4;
    KinCandidateCollection_t mets;  
    mets.push_back( KinCandidate_t(finalMet,finalMet.pt()));    
    LorentzVector deltaMet(finalMet-phys.met);
    if(mets[0].first.pt()<30) return;
    
    //order collections
    sort(leptons.begin(),leptons.end(),KinAnalysis::sortKinCandidates);
    sort(jets.begin(),jets.end(),KinAnalysis::sortKinCandidates);
    sort(mets.begin(),mets.end(),KinAnalysis::sortKinCandidates);

    //debug
    cout << "[KinAnalysis][runOn] " << ev.run << " : " << ev.lumi << " : " << ev.event << endl
	 << "Scheme is: " << scheme_ << endl
	 << "Leptons #1 : (" << leptons[0].first.pt() << ";" << leptons[0].first.eta() << ";" << leptons[0].first.phi() << ") q:" << leptons[0].second << endl  
	 << "        #2 : (" << leptons[1].first.pt() << ";" << leptons[1].first.eta() << ";" << leptons[1].first.phi() << ") q:" << leptons[1].second << endl  
	 << "Jets    #1 : (" << jets[0].first.pt() << ";" << jets[0].first.eta() << ";" << jets[0].first.phi() << ") btag:" << jets[0].second  << endl  
	 << "        #2 : (" << jets[1].first.pt() << ";" << jets[1].first.eta() << ";" << jets[1].first.phi() << ") btag:" << jets[1].second << endl  
	 << "MET        : (" << mets[0].first.pt() << ";" << mets[0].first.phi() << ")" << endl;
    cout << "Deltas: Lep: " << deltaLep.pt() << " Jet:" << deltaJets.pt() << " MET: " << deltaMet.pt() << endl;
    
    //set to the base values
    int nComb=0;
    resHandler_.resetHistos();

    //try all combinatorics
    for(int ijet=0; ijet<maxJetMult_; ijet++) 
      {
	for(int jjet=0; jjet<maxJetMult_; jjet++)
	  {
	    if(ijet==jjet) continue;
	    nComb++;
	    
	    for(int itry=1; itry<maxTries_; itry++)
	      {		  
		//get a new pzttbar hypothesis
		double deltaPz = deltaPzFunc_->GetRandom();		
		
		//leptons
		LorentzVector pl1 = leptons[0].first;
		LorentzVector pl2 = leptons[1].first;
	
		//smear jets and MET using the resolution functions
		std::vector<LorentzVector> origJetsP4; origJetsP4.push_back(jets[0].first); origJetsP4.push_back(jets[1].first);
		std::vector<LorentzVectorCollection> jetsVar;
		LorentzVectorCollection metsVar;
		jet::computeVariation(origJetsP4,mets[0].first,jetsVar,metsVar,ptResol,etaResol,phiResol,jecUnc);
		LorentzVector pb1 = jetsVar[jet::JER][0];
		LorentzVector pb2 = jetsVar[jet::JER][1];
		LorentzVector met(metsVar[jet::JER]);
		
		//MET (must correct for jet energy scale/resolution smearing)
		float metResol= 1.0;//+rndGen_.Gaus(0,0.1);
		//double dPhiMET = rndGen_.Gaus(0,0.1);
		float metx = metResol*met.px();
		float mety = metResol*met.py();
		TVector3 metConstraint( metx, mety, deltaPz-pb1.Pz()-pl1.Pz()-pb2.Pz()-pl2.Pz());
		//metConstraint.RotateZ(dPhiMET);

		//prevent strange values
		if(pl1.Pt()<1 || pb1.Pt()<1 || pl2.Pt()<1 || pb2.Pt()<1 || metConstraint.Pt()<1) continue;
		
		TLorentzVector Tpl1(pl1.px(),pl1.py(),pl1.pz(),pl1.energy());
		TLorentzVector Tpl2(pl2.px(),pl2.py(),pl2.pz(),pl2.energy());
		TLorentzVector Tpb1(pb1.px(),pb1.py(),pb1.pz(),pb1.energy());
		TLorentzVector Tpb2(pb2.px(),pb2.py(),pb2.pz(),pb2.energy());
		TTbarSolutionCollection_t sols = kin_.findSolutions(Tpl1,Tpb1,Tpl2,Tpb2,metConstraint,mw_);    
		if(sols.size()==0) continue;
		
		//compute the full kinematics obtained
		TTbarSolution_t *sol = &(sols.back());
		TLorentzVector ttbar = sol->pt1+sol->pt2;
		float avgMtop = (sol->pt1.M()+sol->pt2.M())*0.5;
		float mttbar = ttbar.M();
		//std::vector<double> mt2 = getMT2( *sol );
		//float afb = sol->pt1.Eta()-sol->pt2.Eta();
		
		//fill histos
		resHandler_.getHisto("mt", nComb)->Fill( avgMtop );
		resHandler_.getHisto("mttbar",nComb)->Fill(mttbar);
		//resHandler_.getHisto("mt2",nComb)->Fill(mt2[0]);
		//resHandler_.getHisto("afb",nComb)->Fill(afb);
	      }
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
