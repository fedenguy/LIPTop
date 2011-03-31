#include "LIP/Top/interface/DileptonEventCleaner.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"
#include "PhysicsTools/UtilAlgos/interface/FWLiteAnalyzerWrapper.h"

typedef fwlite::AnalyzerWrapper<DileptonEventCleaner> WrappedFWLiteCleanEventAnalyzer;

int main(int argc, char* argv[]) 
{
  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();
  
  // only allow one argument for this simple example which should be the python cfg file
  if ( argc < 2 ) {
    std::cout << "Usage : " << argv[0] << " [parameters.py]" << std::endl;
    return 0;
  }
  if( !edm::readPSetsFrom(argv[1])->existsAs<edm::ParameterSet>("process") ){
    std::cout << " Error : 'process' should be ParametersSet to be run as fwlite in your configuration file" << std::endl; 
    exit(0);
  }
  
  // get the python configuration                                                                                                      
  WrappedFWLiteCleanEventAnalyzer ana(edm::readPSetsFrom(argv[1])->getParameter<edm::ParameterSet>("process"), std::string("evAnalyzer"), std::string("evAnalyzer"));
  ana.beginJob();
  ana.analyze();
  ana.endJob();
  return 0;
}
