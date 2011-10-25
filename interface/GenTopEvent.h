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

      std::list<reco::CandidatePtr > tops;
      std::list<reco::CandidatePtr > neutrinos;
      std::list<reco::CandidatePtr > leptons;
      std::list<reco::CandidatePtr > quarks;
      edm::InputTag genLabel_;

      /**
	 @short iterates back in the particles list until the mother is found
      */
      const reco::Candidate * findFirstMotherOf(const reco::Candidate * p);

      /**
	 @short iterates further down in the particles list until the last state is found
       */
      const reco::Candidate * getFinalStateFor(const reco::Candidate * p);

      /**
	 @short finds a match candidate in the generator level particle lists
      */
      const reco::Candidate * getCandidateFor(const reco::Candidate * recoCandidate, int id, double matchCone=0.1);

      /**
	 @short finds a match candidate in the generator level particle lists
      */
      const reco::Candidate * getNearestCandidateFor(double eta, double phi, int id, double matchCone=0.1);
      
      /**
	 @short finds TTchannel
      */
      int assignTTEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup);


      /**
	 @short get DY->ll mass
       */
      float getDYMass(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      
    };
  }
}

#endif
