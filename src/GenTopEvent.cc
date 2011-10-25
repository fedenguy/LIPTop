#include "LIP/Top/interface/GenTopEvent.h"
#include "DataFormats/Math/interface/deltaR.h"


using namespace std;

namespace gen
{
  namespace top
  {
     
    //
    const reco::Candidate * Event::findFirstMotherOf(const reco::Candidate *p)
    {
      if(p==0) return 0;
      const reco::Candidate *nextP = p->mother();
      if(nextP==0) return 0;
      int pdgId = p->pdgId();
      int motherPdgId = nextP->pdgId();
      int iTry(0);
      while(pdgId==motherPdgId)
	{
	  iTry++;
	  nextP = nextP->mother();
	  if(nextP==0) return 0;
	  motherPdgId = nextP->pdgId();
	  if(iTry>10) return 0;
	}
      return nextP;
    }
      
    //
    const reco::Candidate * Event::getFinalStateFor(const reco::Candidate * p)
    {
      if(p==0) return 0;
	
      const reco::Candidate * prevState=p;
      do{	
	const reco::Candidate * nextState=0;
	int nDaughters = prevState->numberOfDaughters();
	for(int iDaughter=0; iDaughter<nDaughters; iDaughter++)
	  {
	    const reco::Candidate * dau = prevState->daughter(iDaughter);
	    if(dau==0) continue;
	    if(dau->pdgId()!= p->pdgId()) continue;
	    nextState=dau;	   
	    break;
	  }
	if(nextState==0) break;
	if(nextState==prevState) break;
	prevState=nextState;
      }while(1);
	
      return prevState;
    }
      
      
    //
    const reco::Candidate * Event::getCandidateFor(const reco::Candidate * recoCandidate, int id, double matchCone)
    {
      if(recoCandidate==0) return 0;
      std::list<reco::CandidatePtr > &candList = (abs(id)>6 ? leptons : quarks );
	
      //find the closest candidate (compare also with the final state)
      double minDeltaR(100.);
      const reco::Candidate * matchCandidate=0;
      for(std::list<reco::CandidatePtr >::iterator it = candList.begin();
	  it != candList.end();
	  it++)
	{
	  if(abs(id)>6 && TMath::Abs((*it)->pdgId()) != TMath::Abs(id)) continue;
	  const reco::Candidate *  fs = getFinalStateFor( it->get() );
	  double dR = deltaR( **it, *recoCandidate );
	  double dR2 = deltaR( *fs, *recoCandidate );
	  if(dR>minDeltaR && dR2>minDeltaR) continue;
	  matchCandidate = it->get();
	  minDeltaR=(dR<dR2 ? dR : dR2);
	}
	
      //return non null if found inside a reasonable match cone
      if(minDeltaR>matchCone) matchCandidate=0;
	
      return matchCandidate;
    }


    //
    const reco::Candidate * Event::getNearestCandidateFor(double eta, double phi, int id, double matchCone)
    {
      std::list<reco::CandidatePtr > &candList = (abs(id)>6 ? leptons : quarks );
	
      //find the closest candidate (compare also with the final state)
      double minDeltaR(100.);
      const reco::Candidate * matchCandidate=0;
      for(std::list<reco::CandidatePtr >::iterator it = candList.begin();
	  it != candList.end();
	  it++)
	{
	  if(abs(id)>6 && TMath::Abs((*it)->pdgId()) != TMath::Abs(id)) continue;
	  const reco::Candidate *  fs = getFinalStateFor( it->get() );
	  double dR = deltaR( it->get()->eta(), it->get()->phi(), eta, phi );
	  double dR2 = deltaR( fs->eta(), fs->phi(), eta, phi );
	  if(dR>minDeltaR && dR2>minDeltaR) continue;
	  matchCandidate = it->get();
	  minDeltaR=(dR<dR2 ? dR : dR2);
	}
      
      //return non null if found inside a reasonable match cone
      if(minDeltaR>matchCone) matchCandidate=0;
	
      return matchCandidate;
    }
      
    //
    int Event::assignTTEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup)
    {
      //reset all
      int quarkCounter(0),chLeptonCounter(0),electronCounter(0),muonCounter(0),tauCounter(0);  
      tops.clear(); leptons.clear();  neutrinos.clear(); quarks.clear();         
	
      //retrieve the generator level particle collection      
      edm::Handle<edm::View<reco::Candidate> > genParticles;
      iEvent.getByLabel(genLabel_, genParticles);
      
      //iterate over the collection and store t->Wb decays
      std::vector<const reco::Candidate *> W_Sons, theQuarks;
      for(size_t i = 0; i < genParticles.product()->size(); ++ i) {
	reco::CandidatePtr p = genParticles->ptrAt(i);
	int id_p   = p->pdgId();   
	  
	//select tops
	if(abs(id_p) != 6) continue;   
	if( p->status()==3) tops.push_back( p );
	  
	for(size_t b = 0; b < p->numberOfDaughters(); ++ b){
	  const reco::Candidate *top_daughter = p->daughter(b);
	    
	  //store W daughters
	  if(abs(top_daughter->pdgId()) == 24)
	    {                                 
	      for(size_t c = 0; c < top_daughter->numberOfDaughters(); ++c)
		{
		  const reco::Candidate *W_daughter = top_daughter->daughter(c);
		  if(abs(W_daughter->pdgId()) != 24) W_Sons.push_back(W_daughter);
		}
	    }
	  else if(abs(top_daughter->pdgId())<6 && abs(top_daughter->pdgId())!=0 
		  && abs(top_daughter->pdgId())!=22 ) {
	    theQuarks.push_back( top_daughter );
	  }
	}
      }

      //check W decay chain
      std::vector<const reco::Candidate *> theLeptons, theNeutrinos;
      for(std::vector<const reco::Candidate *>::const_iterator Wson_i = W_Sons.begin(); Wson_i != W_Sons.end(); Wson_i++)
	{
	  int id_Wson = (*Wson_i)->pdgId();
	  
	  //taus must be traced further
	  if(abs((*Wson_i)->pdgId()) == 15)
	    {	    	         
	      int tempLepton = 0;
	      for(size_t c = 0; c < (*Wson_i)->numberOfDaughters(); ++ c){
		const reco::Candidate *TauSon = (*Wson_i)->daughter(c);
		if(abs(TauSon->pdgId()) == 15)
		  {
		    for(size_t d = 0; d < TauSon->numberOfDaughters(); ++ d){
		      const reco::Candidate *TauSon2 = TauSon->daughter(d);
		      int tausonid = abs(TauSon2->pdgId());
		      if     (tausonid == 11){tempLepton = 11; theLeptons.push_back ( TauSon2 ); }
		      else if(tausonid == 13){tempLepton = 13; theLeptons.push_back ( TauSon2 ); }
		      else if(tausonid == 12 || tausonid==14 || tausonid==16) theNeutrinos.push_back(TauSon2);
		    }
		  }
		else if(abs(TauSon->pdgId()) ==11 ) {tempLepton = 11; theLeptons.push_back ( TauSon ); }
		else if(abs(TauSon->pdgId()) ==13 ) {tempLepton = 13; theLeptons.push_back ( TauSon ); }
		else if(abs(TauSon->pdgId())==12 || abs(TauSon->pdgId())==14 || abs(TauSon->pdgId())==16) theNeutrinos.push_back(TauSon);
	      }
	      
	      if     (tempLepton == 11) { chLeptonCounter++; electronCounter++;  }
	      else if(tempLepton == 13) { chLeptonCounter++; muonCounter++;      }
	      else                      { chLeptonCounter++; tauCounter++;       theLeptons.push_back( *Wson_i );  }	    
	    }
	  else{	  
	    switch(abs(id_Wson)){
	    case 1:  quarkCounter++;  break;
	    case 2:  quarkCounter++;  break;
	    case 3:  quarkCounter++;  break;
	    case 4:  quarkCounter++;  break;
	    case 5:  quarkCounter++;  break;
	    case 6:  quarkCounter++;  break;
	    case 11: chLeptonCounter++; electronCounter++; theLeptons.push_back( *Wson_i); break;
	    case 13: chLeptonCounter++; muonCounter++;     theLeptons.push_back( *Wson_i); break;
	    case 12: case 14: case 16: theNeutrinos.push_back( *Wson_i ); break;
	    }
	  }
	}


      //now get the candidate ptrs... (should exist a more efficient way of doing this)
      for(size_t i = 0; i < genParticles.product()->size(); ++ i) 
	{
	  reco::CandidatePtr p = genParticles->ptrAt(i);
	  for(std::vector<const reco::Candidate *>::iterator it = theQuarks.begin(); it != theQuarks.end(); it++)
	    if( *it == p.get() ) quarks.push_back(p);
	  for(std::vector<const reco::Candidate *>::iterator it = theLeptons.begin(); it != theLeptons.end(); it++)
	    if( *it == p.get() ) leptons.push_back(p);
	  for(std::vector<const reco::Candidate *>::iterator it = theNeutrinos.begin(); it != theNeutrinos.end(); it++)
	    if( *it == p.get() ) neutrinos.push_back(p);
      }
	    
      //assign the ttbar decay channel
      int ttChannel = UNKNOWN;
      if     (quarkCounter == 4) {ttChannel = ALLJETS;}
      else if(quarkCounter == 2 && electronCounter == 1) {ttChannel = EJETS;}
      else if(quarkCounter == 2 && muonCounter == 1) {ttChannel = MUJETS;}
      else if(quarkCounter == 2 && tauCounter == 1) {ttChannel = TAUJETS;}
      else if(quarkCounter == 0 && electronCounter == 2) {ttChannel = EE;}
      else if(quarkCounter == 0 && electronCounter == 1 && muonCounter == 1){ttChannel = EMU;}
      else if(quarkCounter == 0 && electronCounter == 1 && tauCounter  == 1){ttChannel = ETAU;}
      else if(quarkCounter == 0 && muonCounter == 2) {ttChannel = MUMU;}
      else if(quarkCounter == 0 && muonCounter == 1 && tauCounter == 1) {ttChannel = MUTAU;}
      else if(quarkCounter == 0 && tauCounter == 2) {ttChannel = TAUTAU;}
	
      tops.sort();      tops.unique();
      leptons.sort();   leptons.unique();
      quarks.sort();    quarks.unique();
      neutrinos.sort(); neutrinos.unique();
	
      return ttChannel;
    }
    

    //
    float Event::getDYMass(const edm::Event& iEvent, const edm::EventSetup& iSetup)
    {
      //retrieve the generator level particle collection     
      edm::Handle<edm::View<reco::GenParticle> > genParticles;
      iEvent.getByLabel(genLabel_, genParticles);

      float mass(0.);
      for(size_t i = 0; i < genParticles->size(); ++ i) 
	{
	  const reco::GenParticle & p = (*genParticles)[i];
	  int id_p   = p.pdgId();   
	  
	  //select Z
	  if(abs(id_p) != 23 && abs(id_p)!=22) continue;   
	  if( p.status()!=3) continue;
	  
	  mass = p.mass();
	  break;
	}
      return mass;
    }
  }
}

