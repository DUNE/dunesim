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

#include "TFile.h"

class PDSPDataDrivenBeam;

/*
struct MonitorRegion {
  int upstream_low;
  int upstream_high;
  int downstream_low;
  int downstream_high;

  MonitorRegion(int up_low, int up_high, int down_low, int down_high)
    : upstream_low(up_low), upstream_high(up_high),
      downtream_low(down_low), downstream_high(down_high) {};

  check(int up, int down) {
    return ((upstream_low <= up && up <= upstream_high) &&
            (downstream_low <= down && down <= downstream_high));
  };
};
*/

class PDSPDataDrivenBeam : public art::EDProducer {
public:
  explicit PDSPDataDrivenBeam(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDSPDataDrivenBeam(PDSPDataDrivenBeam const&) = delete;
  PDSPDataDrivenBeam(PDSPDataDrivenBeam&&) = delete;
  PDSPDataDrivenBeam&
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
    fFileName(p.get<std::string>("FileName")) {

  fParticleTypes = p.get<std::vector<std::string>>("ParticleTypes");

}

void PDSPDataDrivenBeam::produce(art::Event& e) {

}

void PDSPDataDrivenBeam::beginJob() {
  //Open input file
  fInputFile = new TFile(fFileName.c_str(), "READ");

  //Get all the PDFs  
  //These will be a 1D distributions in momentum for each
  //4D monitor bin (x,y for up/downstream)
  //Stored in directories according to their particle type
  for (size_t i = 0; i < fParticleTypes.size(); ++i) {
    std::cout << "Fetching directory for " << fParticleTypes[i] << std::endl;
    TDirectory * current_dir = (TDirectory*)fInputFile->Get(fParticleTypes[i].c_str());
    TList * hists = (TList*)current_dir->GetListOfKeys();
    for (int j = 0; j < hists->GetSize(); ++j) {
      std::cout << hists->At(j)->GetName() << std::endl;
      std::string title = hists->At(j)->GetTitle();
      std::string h_bins = title.substr(0, title.find("_"));
      std::string v_bins = title.substr(title.find("_")+1, title.size()-1);
      std::cout << title << " " << h_bins << " " << v_bins << std::endl;
    }
  }

}

DEFINE_ART_MODULE(PDSPDataDrivenBeam)
