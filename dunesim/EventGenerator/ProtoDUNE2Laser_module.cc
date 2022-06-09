////////////////////////////////////////////////////////////////////////
// Class:       ProtoDUNE2Laser
// Plugin Type: producer (Unknown Unknown)
// File:        ProtoDUNE2Laser_module.cc
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
#include "TRandom3.h"

#include <memory>

namespace evgen {
  class ProtoDUNE2Laser;
}


class evgen::ProtoDUNE2Laser : public art::EDProducer {
public:
  explicit ProtoDUNE2Laser(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ProtoDUNE2Laser(ProtoDUNE2Laser const&) = delete;
  ProtoDUNE2Laser(ProtoDUNE2Laser&&) = delete;
  ProtoDUNE2Laser& operator=(ProtoDUNE2Laser const&) = delete;
  ProtoDUNE2Laser& operator=(ProtoDUNE2Laser&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  void GenerateLaserPulse(simb::MCTruth & truth);

  double fJouleToeV = 6.241509e18;

  // Declare member data here.
  TRandom3 fRNG;

  double fCenterX, fCenterY, fCenterZ; //Central Beam Position
  double fWidthX, fWidthY; //Beam Profile Shape
  //Laser settings 
  size_t fIntensity; //N photons
  double fWavelength; //wavelength in nm

};


evgen::ProtoDUNE2Laser::ProtoDUNE2Laser(fhicl::ParameterSet const& p)
  : EDProducer{p}, fRNG(p.get<int>("Seed", 0)),
    fCenterX(p.get<double>("CenterX")), fCenterY(p.get<double>("CenterY")),
    fCenterZ(p.get<double>("CenterZ")), fWidthX(p.get<double>("WidthX")),
    fWidthY(p.get<double>("WidthY")), fIntensity(p.get<size_t>("Intensity")),
    fWavelength(p.get<double>("Wavelength")) {

  //MCTruth to hold set of photons
  produces<std::vector<simb::MCTruth>>();
}

void evgen::ProtoDUNE2Laser::produce(art::Event& e) {
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

void evgen::ProtoDUNE2Laser::GenerateLaserPulse(simb::MCTruth & truth) {

  
  double energy = 1.0;//GeV/c?
  double x = fRNG.Gaus(fCenterX, fWidthX);
  double y = fRNG.Gaus(fCenterY, fWidthY);
  const TLorentzVector pos(x, y, fCenterZ, 0.);

  simb::MCParticle muon(0, -13, "primary");
  const TLorentzVector mom(0., 0., sqrt(energy*energy - muon.Mass()*muon.Mass()), energy);
  muon.AddTrajectoryPoint(pos, mom);
  truth.Add(muon);
 

  /*for (size_t i = 0; i < fIntensity; ++i) {
    simb::MCParticle photon(i, 22, "primary", -1, 0., 1);
    double energy = 1.e6*TMath::HC()*fJouleToeV/fWavelength;
    const TLorentzVector mom(0., 0., energy, energy);

    double x = fRNG.Gaus(fCenterX, fWidthX);
    double y = fRNG.Gaus(fCenterY, fWidthY);
    const TLorentzVector pos(x, y, fCenterZ, 0.);

    std::cout << "Generating photon at (" << x << ", " << y << ", " <<
                 fCenterZ << ") with energy " << energy << std::endl;
    photon.AddTrajectoryPoint(pos, mom);
    truth.Add(photon);
  }*/
}

DEFINE_ART_MODULE(evgen::ProtoDUNE2Laser)
