#ifndef H4BeamInput_h
#define H4BeamInput_h

#include "H4BeamFileService.h"
#include "art/Framework/Core/InputSourceMacros.h" 
#include "art/Framework/IO/Sources/Source.h" 
#include "art/Framework/IO/Sources/SourceTraits.h"
#include "art/Framework/Core/fwd.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Provenance/FileFormatVersion.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"



namespace dune {
//Forward declare the class
class H4BeamInputDetail;
}

class dune::H4BeamInputDetail {
 public:
  H4BeamInputDetail(fhicl::ParameterSet const & ps,
                    art::ProductRegistryHelper & rh,
                    art::SourceHelper const & sh);

  void readFile(std::string const & filename, art::FileBlock*& fb);

  bool readNext(art::RunPrincipal const* const inR,
                art::SubRunPrincipal const* const inSR,
                art::RunPrincipal*& outR,
                art::SubRunPrincipal*& outSR,
                art::EventPrincipal*& outE);

  void closeCurrentFile() {
    art::ServiceHandle<dune::H4BeamFileService> beamFileService;
    beamFileService->Reset();
  };

 private:
  art::SourceHelper const& pmaker;
  //int fLogLevel;
  int fSkipEvents = 0;
  size_t fNEventsAvailable = 0;
  bool fAccessedFirstEvent = false;
 };
#endif
