////////////////////////////////////////////////////////////////////////
// Class:       ProtoDUNEBeamTest
// Plugin Type: analyzer (art v2_07_03)
// File:        ProtoDUNEBeamTest_module.cc
//
// Generated at Mon Sep  4 06:55:33 2017 by Leigh Whitehead using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////
//// Modified by Pablo and Leigh H. Howard, 
//Smear important variables of beam monitors, mimic Cherenkov monitors response
//and store some histograms through art utilities
//// pablo.fer@cern.ch
//// July 2018
///////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
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
//Histogram utilities
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include <fstream>
#include "TH1F.h"


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
TH1F *fT0FHist_zoom;
TH1F *fT0FsHist_zoom;
TH1F *fUxHist_zoom;
TH1F *fUxsHist_zoom;
TH1F *fUyHist_zoom;
TH1F *fUysHist_zoom;
TH1F *fUzHist_zoom;
TH1F *fUzsHist_zoom;
TH1F *fT0FHist;
TH1F *fT0FsHist;
TH1F *fUxHist;
TH1F *fUxsHist;
TH1F *fUyHist;
TH1F *fUysHist;
TH1F *fUzHist;
TH1F *fUzsHist;
};


sim::ProtoDUNEBeamTest::ProtoDUNEBeamTest(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
{

}

void sim::ProtoDUNEBeamTest::beginJob()
{
// Declaring histograms
  art::ServiceHandle<art::TFileService> tfs;
  fT0FHist_zoom = tfs->make<TH1F>("TOF_TRUE_zoom","T0F (ns)",100,95,96);
  fT0FsHist_zoom = tfs->make<TH1F>("TOF_SMEARED_zoom","T0F (ns)",100,95,96);
  fUxHist_zoom = tfs->make<TH1F>("X_DIRECTION_TRUE_zoom","Ux",100,-0.05,0.05);
  fUyHist_zoom = tfs->make<TH1F>("Y_DIRECTION_TRUE_zoom","Uy",100,-0.05,0.05);
  fUzHist_zoom = tfs->make<TH1F>("Z_DIRECTION_TRUE_zoom","Uz",100,-1,-0.95);
  fUxsHist_zoom = tfs->make<TH1F>("X_DIRECTION_SMEARED_zoom","Ux",100,-0.05,0.05);
  fUysHist_zoom = tfs->make<TH1F>("Y_DIRECTION_SMEARED_zoom","Uy",100,-0.05,0.05);
  fUzsHist_zoom = tfs->make<TH1F>("Z_DIRECTION_SMEARED_zoom","Uz",100,-1,-0.95);
  fT0FHist = tfs->make<TH1F>("TOF_TRUE","T0F (ns)",100,90,100);
  fT0FsHist = tfs->make<TH1F>("TOF_SMEARED","T0F (ns)",100,90,100);
  fUxHist = tfs->make<TH1F>("X_DIRECTION_TRUE","Ux",100,-1,1);
  fUyHist = tfs->make<TH1F>("Y_DIRECTION_TRUE","Uy",100,-1,1);
  fUzHist = tfs->make<TH1F>("Z_DIRECTION_TRUE","Uz",100,-1,1);
  fUxsHist = tfs->make<TH1F>("X_DIRECTION_SMEARED","Ux",100,-1,1);
  fUysHist = tfs->make<TH1F>("Y_DIRECTION_SMEARED","Uy",100,-1,1);
  fUzsHist = tfs->make<TH1F>("Z_DIRECTION_SMEARED","Uz",100,-1,1);
}

void sim::ProtoDUNEBeamTest::analyze(art::Event const & evt)
{

  // Get the reconstructed tracks
  auto beamsim = evt.getValidHandle<std::vector<sim::ProtoDUNEbeamsim> >("generator");
  const sim::ProtoDUNEbeamsim temp = (*beamsim)[0];
//  unsigned short nInst = temp.NInstruments();
// Computing true and smeared TOF (between TOF1 and TRIG2)
  sim::ProtoDUNEBeamInstrument tof1 = temp.GetInstrument("TOF1");
  sim::ProtoDUNEBeamInstrument trig2 = temp.GetInstrument("TRIG2");
  float tof = trig2.GetT() - tof1.GetT();
  float tof_smeared = trig2.GetSmearedVar1() -tof1.GetSmearedVar1();
//  std::cout << "TOF " << tof << std::endl;
//  std::cout << "TOFS " << tof_smeared << std::endl;

  sim::ProtoDUNEBeamInstrument bprofext = temp.GetInstrument("BPROFEXT");
  sim::ProtoDUNEBeamInstrument bprof4 = temp.GetInstrument("BPROF4");
  double dir[3],dir_smeared[3],dummy;
  dummy = pow(pow(bprofext.GetX() -bprof4.GetX(),2)+pow(bprofext.GetY() -bprof4.GetY(),2)+pow(bprofext.GetZ() -bprof4.GetZ(),2),0.5);
  dir[0] = (bprofext.GetX() -bprof4.GetX())/dummy;
  dir[1] = (bprofext.GetY() -bprof4.GetY())/dummy;
  dir[2] = (bprofext.GetZ() -bprof4.GetZ())/dummy;
  dummy = pow(pow(bprofext.GetSmearedVar1() -bprof4.GetSmearedVar1(),2)+pow(bprofext.GetSmearedVar2() -bprof4.GetSmearedVar2(),2)+pow(bprofext.GetZ() -bprof4.GetZ(),2),0.5);
  dir_smeared[0] = (bprofext.GetSmearedVar1() -bprof4.GetSmearedVar1())/dummy;
  dir_smeared[1] = (bprofext.GetSmearedVar2() -bprof4.GetSmearedVar2())/dummy;
  dir_smeared[2] = (bprofext.GetZ() -bprof4.GetZ())/dummy;
//  std::cout << "cos" << dir[0]*dir_smeared[0]+dir[1]*dir_smeared[1]+dir[2]*dir_smeared[2] << std::endl;
//  std::cout << bprofext.GetX() << "," << bprofext.GetSmearedVar1() << std::endl;
//  std::cout << dir[0] << "," << dir_smeared[0] << std::endl;

//Filling hisotgrams
  fT0FHist_zoom->Fill(tof);
  fT0FsHist_zoom->Fill(tof_smeared);
  fUxHist_zoom->Fill(dir[0]);
  fUxsHist_zoom->Fill(dir_smeared[0]);
  fUyHist_zoom->Fill(dir[1]);
  fUysHist_zoom->Fill(dir_smeared[1]);
  fUzHist_zoom->Fill(dir[2]);
  fUzsHist_zoom->Fill(dir_smeared[2]);
  fT0FHist->Fill(tof);
  fT0FsHist->Fill(tof_smeared);
  fUxHist->Fill(dir[0]);
  fUxsHist->Fill(dir_smeared[0]);
  fUyHist->Fill(dir[1]);
  fUysHist->Fill(dir_smeared[1]);
  fUzHist->Fill(dir[2]);
  fUzsHist->Fill(dir_smeared[2]);
}

void sim::ProtoDUNEBeamTest::endJob()
{

}

DEFINE_ART_MODULE(sim::ProtoDUNEBeamTest)

