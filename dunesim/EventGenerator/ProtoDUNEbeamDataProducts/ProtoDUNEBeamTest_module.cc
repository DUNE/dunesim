////////////////////////////////////////////////////////////////////////
// Class:       ProtoDUNEBeamTest
// Plugin Type: analyzer (art v2_07_03)
// File:        ProtoDUNEBeamTest_module.cc
//
// Generated at Mon Sep  4 06:55:33 2017 by Leigh Whitehead using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "dune/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEbeamsim.h"

namespace sim {
  class ProtoDUNEBeamTest;
}


class sim::ProtoDUNEBeamTest : public art::EDAnalyzer {
public:

  explicit ProtoDUNEBeamTest(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ProtoDUNEBeamTest(ProtoDUNEBeamTest const &) = delete;
  ProtoDUNEBeamTest(ProtoDUNEBeamTest &&) = delete;
  ProtoDUNEBeamTest & operator = (ProtoDUNEBeamTest const &) = delete;
  ProtoDUNEBeamTest & operator = (ProtoDUNEBeamTest &&) = delete;

  virtual void beginJob() override;
  virtual void endJob() override;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

};


sim::ProtoDUNEBeamTest::ProtoDUNEBeamTest(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
{

}

void sim::ProtoDUNEBeamTest::beginJob()
{

}

void sim::ProtoDUNEBeamTest::analyze(art::Event const & evt)
{

  // Get the reconstructed tracks
  auto beamsim = evt.getValidHandle<std::vector<sim::ProtoDUNEbeamsim> >("generator");
  const sim::ProtoDUNEbeamsim temp = (*beamsim)[1];
  unsigned short nInst = temp.NInstruments();
  std::cout << "Number of beam instruments read from the file = " << nInst << std::endl;

}

void sim::ProtoDUNEBeamTest::endJob()
{

}

DEFINE_ART_MODULE(sim::ProtoDUNEBeamTest)

