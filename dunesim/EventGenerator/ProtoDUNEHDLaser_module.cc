////////////////////////////////////////////////////////////////////////
// Class:       ProtoDUNEHDLaser
// Plugin Type: producer (Unknown Unknown)
// File:        ProtoDUNEHDLaser_module.cc
//
// Generated at Mon Apr 18 11:20:27 2022 by Jacob Calcutt using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

#include <memory>

namespace evgen {
  class ProtoDUNEHDLaser;
}


class evgen::ProtoDUNEHDLaser : public art::EDProducer {
public:
  explicit ProtoDUNEHDLaser(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ProtoDUNEHDLaser(ProtoDUNEHDLaser const&) = delete;
  ProtoDUNEHDLaser(ProtoDUNEHDLaser&&) = delete;
  ProtoDUNEHDLaser& operator=(ProtoDUNEHDLaser const&) = delete;
  ProtoDUNEHDLaser& operator=(ProtoDUNEHDLaser&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  void GenerateLaserPulse(simb::MCTruth & truth);

  //double fJouleToeV = 6.241509e18;

  // Declare member data here.
  TRandom3 fRNG;

  double fCenterX, fCenterY, fCenterZ; //Central Beam Position
  double fWidthX, fWidthY; //Beam Profile Shape
  //Laser settings 
  size_t fIntensity; //N photons
  double fWavelength; //wavelength in nm
  std::string fInputFile;
  int fCurrentEvent;

  TTree * fInputTree;

  double fX, fY, fZ, fTheta, fPhi;
  double fToRad = TMath::Pi()/180.;
};


evgen::ProtoDUNEHDLaser::ProtoDUNEHDLaser(fhicl::ParameterSet const& p)
  : EDProducer{p}, fRNG(p.get<int>("Seed", 0)),
    fCenterX(p.get<double>("CenterX")), fCenterY(p.get<double>("CenterY")),
    fCenterZ(p.get<double>("CenterZ")), fWidthX(p.get<double>("WidthX")),
    fWidthY(p.get<double>("WidthY")), fIntensity(p.get<size_t>("Intensity")),
    fWavelength(p.get<double>("Wavelength")),
    fInputFile(p.get<std::string>("InputFile")),
    fCurrentEvent(p.get<int>("StartEvent")) {

  //MCTruth to hold set of photons
  produces<std::vector<simb::MCTruth>>();
  auto * input_file = TFile::Open(fInputFile.c_str());
  fInputTree = (TTree*)input_file->Get("tree");
  //fInputTree->SetDirectory(0);

  fInputTree->SetBranchAddress("x", &fX);
  fInputTree->SetBranchAddress("y", &fY);
  fInputTree->SetBranchAddress("z", &fZ);
  fInputTree->SetBranchAddress("theta", &fTheta);
  fInputTree->SetBranchAddress("phi", &fPhi);
}

void evgen::ProtoDUNEHDLaser::produce(art::Event& e) {
  // Define the truth collection for this event.
  auto truthcol = std::make_unique<std::vector<simb::MCTruth>>();
  simb::MCTruth truth;
  truth.SetOrigin(simb::kSingleParticle);

  GenerateLaserPulse(truth);

  // Add the MCTruth to the vector
  truthcol->push_back(truth);

  // Finally, add the MCTruth to the event
  e.put(std::move(truthcol));
}

void evgen::ProtoDUNEHDLaser::GenerateLaserPulse(simb::MCTruth & truth) {

  fInputTree->GetEntry(fCurrentEvent); 
  
  double momentum = 1.0;//GeV/c?
  simb::MCParticle laser(0, 66613, "primary");
  double energy = sqrt(laser.Mass()*laser.Mass() + momentum*momentum);

  double px = momentum*sin(fTheta*fToRad)*cos(fPhi*fToRad);
  double py = momentum*sin(fTheta*fToRad)*sin(fPhi*fToRad);
  double pz = momentum*cos(fTheta*fToRad);

  std::cout << "(X, Y, Z)" << fX << " " << fY << " " << fZ << std::endl;
  std::cout << "(PX, PY, PZ)" << px << " " << py << " " << pz << std::endl;
  const TLorentzVector mom(px, py, pz, energy);
  const TLorentzVector pos(fX, fY, fZ, 0.);
  laser.AddTrajectoryPoint(pos, mom);
  truth.Add(laser);
 
  ++fCurrentEvent;
}

DEFINE_ART_MODULE(evgen::ProtoDUNEHDLaser)
