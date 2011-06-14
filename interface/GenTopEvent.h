#ifndef _gen_top_event_hh_
#define _gen_top_event_hh_

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TVector3.h"

namespace gen
{
  namespace top
  {
    class Event
    {
    public:
      Event() { };
      ~Event() { };
      
      enum TTChannel { UNKNOWN=0, ALLJETS, EJETS, MUJETS, TAUJETS, EE, EMU, ETAU, MUMU, MUTAU, TAUTAU };

      std::list<const reco::Candidate *> tops;
      std::list<const reco::Candidate *> neutrinos;
      std::list<const reco::Candidate *> leptons,dyDaughters;
      std::list<const reco::Candidate *> quarks;
      edm::InputTag genLabel_;

      /**
	 @short filters according to the pthat
	 @return true if minpthat<pthat<maxpthat or if pthat info is not found
       */
      inline bool filterDYBin(double mindymass=0,double maxdymass=14000)
      {
	if(dyDaughters.size()!=2) return true;
	math::XYZTLorentzVector totalMom(0,0,0,0);
	for(std::list<const reco::Candidate *>::iterator it = dyDaughters.begin();
	    it != dyDaughters.end();
	    it++)
	  {
	    const math::XYZTLorentzVector &imom = (*it)->p4();
	    totalMom += imom;
	  }
	double mass = totalMom.mass();
	if(mass<mindymass || mass> maxdymass) return false;
	return true;
      }

      /**
	 @short iterates back in the particles list until the mother is found
      */
      const reco::Candidate *findFirstMotherOf(const reco::Candidate *p);

      /**
	 @short iterates further down in the particles list until the last state is found
       */
      const reco::Candidate *getFinalStateFor(const reco::Candidate *p);

      /**
	 @short finds a match candidate in the generator level particle lists
      */
      const reco::Candidate *getCandidateFor(const reco::Candidate *recoCandidate, int id, double matchCone=0.1);

      /**
	 @short finds a match candidate in the generator level particle lists
      */
      const reco::Candidate *getNearestCandidateFor(double eta, double phi, int id, double matchCone=0.1);
      
      /**
	 @short finds TTchannel
      */
      int assignTTEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup);
    };
  }
}

#endif
