#ifndef misassignmentmeasurement_hh
#define misassignmentmeasurement_hh

#include <iostream>
#include <string>
#include <vector>

#include "LIP/Top/interface/EventSummaryHandler.h"

#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/SelectionMonitor.h"

#include "TFile.h"
#include "TRandom2.h"
#include "TH1D.h"


/**
   @short 
   takes care of measuring the missassignments from the M_{lj} spectrum
   while bookeeping the statistics for the pseudo-experiments
 */
class MisassignmentMeasurement
{
 public:

  MisassignmentMeasurement() : nMeasurements(0)  { bookMonitoringHistograms();  }

  ~MisassignmentMeasurement()  { }

  TRandom2 rndGen;
  SelectionMonitor controlHistos;

  /**
     @short books the histograms
  */
  void bookMonitoringHistograms();

  //
  void saveMonitoringHistograms();

  /**
     @short returns a collection of randomly rotated leptons
  */
  PhysicsObjectLeptonCollection randomlyRotate( PhysicsObjectLeptonCollection &leptons, PhysicsObjectJetCollection &jets);

  /**
     @short runs the measurement
  */
  void measureMisassignments(EventSummaryHandler &evHandler, double mcut=190, double minMlj=40, bool isData=false);
  
  /**
     @short getters
  */
  std::vector<double> getCorrectPairsFraction() 
    {
      std::vector<double> res(2,0);
      res[0]=fCorrectPairsEst;
      res[1]=fCorrectPairsEstErr;
      return res;
    }
  std::vector<double> getAlpha() 
    {
      std::vector<double> res(2,0);
      res[0]=alphaEst;
      res[1]=alphaEstErr;
      return res;
    }
  double getNorm() { return kNorm; }

 private:
  /**
   */
  void resetHistograms();

  int nMeasurements;

  double kNorm;
  double fCorrectPairsEst, fCorrectPairsEstErr;
  double alphaEst,alphaEstErr;
};


#endif
