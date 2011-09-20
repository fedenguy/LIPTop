#ifndef liptop_toolslinkdef_h
#define liptop_toolslinkdef_h

#include "LIP/Top/interface/EventSummaryHandler.h"
#include "LIP/Top/interface/HeavyFlavorPDF.h"
#include "LIP/Top/interface/TopKinSolver.h"
#include "LIP/Top/interface/KinAnalysis.h"
#include "LIP/Top/interface/HFCMeasurement.h"
#include "LIP/Top/interface/MisassignmentMeasurement.h"
#include "LIP/Top/interface/MassMeasurement.h"
#include "CMGTools/HtoZZ2l2nu/interface/setStyle.h"
#include "CMGTools/HtoZZ2l2nu/interface/plotter.h"
#include "CMGTools/HtoZZ2l2nu/interface/SelectionMonitor.h"

#ifdef __CINT__

#pragma link off all class; 
#pragma link off all function; 
#pragma link off all global; 
#pragma link off all typedef;

#pragma link C++ function setStyle;
#pragma link C++ function getNewCanvas;
#pragma link C++ function formatPlot;
#pragma link C++ function fixExtremities;
#pragma link C++ function formatForCmsPublic;
#pragma link C++ function showPlots;
#pragma link C++ function showStackPlot;
#pragma link C++ function showSimplePlot;
#pragma link C++ function showMCtoDataComparison;
#pragma link C++ function getPlotAsTable;
#pragma link C++ function getProjections;
#pragma link C++ function showPlotsAndMCtoDataComparison;
#pragma link C++ class SelectionMonitor;

#pragma link C++ struct EventSummary_t;
#pragma link C++ class EventSummaryHandler;
#pragma link C++ typedef LorentzVector;
#pragma link C++ typedef LorentzVectorCollection;
#pragma link C++ PhysicsObject;
#pragma link C++ PhysicsObject_Lepton;
#pragma link C++ PhysicsObject_Jet;
#pragma link C++ typedef PhysicsObjectCollection;
#pragma link C++ typedef PhysicsObjectLeptonCollection;
#pragma link C++ typedef PhysicsObjectJetCollection;
#pragma link C++ struct PhysicsEvent_t;
#pragma link C++ function getPhysicsEventFrom;
#pragma link C++ class HeavyFlavorPDF;
#pragma link C++ struct TTbarSolution_t;
#pragma link C++ typedef TTbarSolutionCollection_t;
#pragma link C++ function mTTbarOrder;
#pragma link C++ function getMT;
#pragma link C++ function getMT2;
#pragma link C++ class TopKinSolver;
#pragma link C++ class KinAnalysis;
#pragma link C++ function randomlyRotate;
#pragma link C++ struct CombinedHFCModel_t;
#pragma link C++ class HFCMeasurement;
#pragma link C++ class MisassignmentMeasurement;
#pragma link C++ struct EnsembleMeasurement_t;
#pragma link C++ struct MassFitResults_t;
#pragma link C++ class MassMeasurement;

#endif

#endif

// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 8
// End:

