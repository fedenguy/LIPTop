#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"


#include "LIP/Top/interface/DileptonEventCleaner.h"
#include "PhysicsTools/UtilAlgos/interface/EDAnalyzerWrapper.h"
typedef edm::AnalyzerWrapper<DileptonEventCleaner> WrappedEDDileptonEventCleaner;
DEFINE_FWK_MODULE(DileptonEventCleaner);  
