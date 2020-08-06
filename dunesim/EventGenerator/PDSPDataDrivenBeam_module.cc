////////////////////////////////////////////////////////////////////////
// Class:       PDSPDataDrivenBeam
// Plugin Type: producer (art v3_05_01)
// File:        PDSPDataDrivenBeam_module.cc
//
// Generated at Thu Aug  6 08:22:36 2020 by Jacob Calcutt using cetskelgen
// from cetlib version v3_10_00.
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

#include <memory>

class PDSPDataDrivenBeam;


class PDSPDataDrivenBeam : public art::EDProducer {
public:
  explicit PDSPDataDrivenBeam(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDSPDataDrivenBeam(PDSPDataDrivenBeam const&) = delete;
  PDSPDataDrivenBeam(PDSPDataDrivenBeam&&) = delete;
  PDSDataDrivenBeam&
      operator=(PDSPDataDrivenBeam const&) = delete;
  PDSPDataDrivenBeam& operator=(PDSPDataDrivenBeam&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;
  void beginJob() override;

private:
  
  TFile* fInputFile = 0x0;
  std::string fFileName;
  std::vector<std::string> fParticleTypes;

};


PDSPDataDrivenBeam::PDSPDataDrivenBeam(
    fhicl::ParameterSet const& p)
  : EDProducer{p},
    fFileName(p.get<std::string>("FileName"),
    fParticleTypes(p.get<std::vector<std::string>>("ParticleTypes") {

}

void PDSPDataDrivenBeam::produce(art::Event& e) {

}

void PDSPDataDrivenBeam::beginJob() {
  //Open input file
  fInputFile = new TFile(fFileName.c_str(), "READ");

  //Get all the PDFs  
  //These will be 1D distributions in momentum for each
  //4D monitor bin (x,y for up/downstream)
  //Stored in directories according to their particle type
  for (size_t i = 0; i < fParticleTypes.size(); ++i) {
    std::cout << "Fetching directory for " << fParticleTypes[i] << std::endl;
    TDirectory * current_dir = fInputFile->cd(fParticleTypes[i].c_str());
    TList * hists = current_dir->GetListOfKeys();
    for (size_t j = 0; j < hists->GetSize(); ++j) {
      std::cout << hists->At(j)->GetName() << std::endl;
    }
  }

}

DEFINE_ART_MODULE(PDSPDataDrivenBeam)
