////////////////////////////////////////////////////////////////////////

// Class:       BackgroundValidator

// Plugin Type: analyzer (art v3_06_03)

// File:        BackgroundValidator_module.cc

//

// Generated at Mon Jul 12 05:19:35 2021 by Pierre Lasorak using cetskelgen

// from cetlib version v3_11_01.

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

#include "lardataobj/Simulation/sim.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "canvas/Persistency/Common/FindMany.h"

#include "canvas/Persistency/Common/FindManyP.h"

#include "canvas/Persistency/Common/FindOneP.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "nusimdata/SimulationBase/MCTruth.h"

#include "art_root_io/TFileService.h"

#include "TTree.h"

class BackgroundValidator;

class BackgroundValidator : public art::EDAnalyzer {

public:

  explicit BackgroundValidator(fhicl::ParameterSet const& p);

  // The compiler-generated destructor is fine for non-base

  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.

  BackgroundValidator(BackgroundValidator const&) = delete;

  BackgroundValidator(BackgroundValidator&&) = delete;

  BackgroundValidator& operator=(BackgroundValidator const&) = delete;

  BackgroundValidator& operator=(BackgroundValidator&&) = delete;

  void beginJob();

  void endJob() {}

  // Required functions.

  void analyze(art::Event const& e) override;

private:

  const std::string input_label;

  TTree* tree;

  int    Run     ;

  int    SubRun  ;

  int    Event   ;

  int    PDG     ;

  double Energy  ;

  double Momentum;

  double DirX    ;

  double DirY    ;

  double DirZ    ;

  double StartX  ;

  double StartY  ;

  double StartZ  ;

  double EndX    ;

  double EndY    ;

  double EndZ    ;

};

BackgroundValidator::BackgroundValidator(fhicl::ParameterSet const& p)

  : EDAnalyzer{p},

  input_label{p.get<std::string>("InputLabel")},

  tree{nullptr} {

  

  }

void BackgroundValidator::beginJob() {

  art::ServiceHandle<art::TFileService> tfs;

  tree = tfs->make<TTree>("BackgroundValidation","Tree to validate background");

  tree->Branch("Run"     , &Run     );

  tree->Branch("SubRun"  , &SubRun  );

  tree->Branch("Event"   , &Event   );

  tree->Branch("PDG"     , &PDG     );

  tree->Branch("Energy"  , &Energy  );

  tree->Branch("Momentum", &Momentum);

  tree->Branch("DirX"    , &DirX    );

  tree->Branch("DirY"    , &DirY    );

  tree->Branch("DirZ"    , &DirZ    );

  tree->Branch("StartX"  , &StartX  );

  tree->Branch("StartY"  , &StartY  );

  tree->Branch("StartZ"  , &StartZ  );

  tree->Branch("EndX"    , &EndX    );

  tree->Branch("EndY"    , &EndY    );

  tree->Branch("EndZ"    , &EndZ    );

}

void BackgroundValidator::analyze(art::Event const& e) {

  Run    = e.run();

  SubRun = e.subRun();

  Event  = e.event();

  

  art::Handle<std::vector<simb::MCTruth>> truths;

  e.getByLabel(input_label, truths);

  for (size_t i=0; i<truths->size(); ++i) {

    const simb::MCTruth truth = truths->at(i);

    for (int it=0; it<truth.NParticles(); ++it) {

      simb::MCParticle p = truth.GetParticle(it);

      PDG      = p.PdgCode();

      Energy   = p.E();

      Momentum = p.P();

      DirX     = p.Px() / p.P();

      DirY     = p.Py() / p.P();

      DirZ     = p.Pz() / p.P();

      StartX   = p.Vx();

      StartY   = p.Vy();

      StartZ   = p.Vz();

      EndX     = p.EndX();

      EndY     = p.EndY();

      EndZ     = p.EndZ();

      tree->Fill();

    }

    

  }

}

DEFINE_ART_MODULE(BackgroundValidator)
