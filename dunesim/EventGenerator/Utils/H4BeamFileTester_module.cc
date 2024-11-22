////////////////////////////////////////////////////////////////////////
// Class:       H4BeamFileTester
// Plugin Type: analyzer (Unknown Unknown)
// File:        H4BeamFileTester_module.cc
//
// Generated at Fri Mar  8 15:24:47 2024 by Jacob Calcutt using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "H4BeamFileService.h"

namespace dune {
  class H4BeamFileTester;
}


class dune::H4BeamFileTester : public art::EDAnalyzer {
public:
  explicit H4BeamFileTester(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  H4BeamFileTester(H4BeamFileTester const&) = delete;
  H4BeamFileTester(H4BeamFileTester&&) = delete;
  H4BeamFileTester& operator=(H4BeamFileTester const&) = delete;
  H4BeamFileTester& operator=(H4BeamFileTester&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

private:

  // Declare member data here.

};


dune::H4BeamFileTester::H4BeamFileTester(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void dune::H4BeamFileTester::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  art::ServiceHandle<dune::H4BeamFileService> beamFileService;
  std::cout << "Art event: " << e.run() << " " << e.subRun() << " " <<
               e.id().event() << std::endl;
  std::cout << "BFS: " << beamFileService->GetCurrentEvent() << std::endl;
}

DEFINE_ART_MODULE(dune::H4BeamFileTester)
