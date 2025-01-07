#include "H4BeamInput.h"
#include "fhiclcpp/ParameterSetRegistry.h"

dune::H4BeamInputDetail::H4BeamInputDetail(
    fhicl::ParameterSet const & ps,
    art::ProductRegistryHelper & rh,
    art::SourceHelper const & sh) 
  : pmaker(sh),
    fSkipEvents(ps.get<int>("skipEvents", 0)) {

}

void dune::H4BeamInputDetail::readFile(
    std::string const & filename, art::FileBlock*& fb) {
  MF_LOG_INFO("H4Beam")
      << "Opening Beam Input File " <<
         filename;

  art::ServiceHandle<dune::H4BeamFileService> beamFileService;
  beamFileService->OpenFile(filename);
  beamFileService->IncrementEvent(fSkipEvents);
  fNEventsAvailable = beamFileService->GetEventIDs()->size();
  MF_LOG_INFO("H4Beam")
      << "Has events: " << fNEventsAvailable;

  fb = new art::FileBlock(art::FileFormatVersion(1, "RawEvent2011"),
                          filename);
}

bool dune::H4BeamInputDetail::readNext(art::RunPrincipal const* const inR,
                                       art::SubRunPrincipal const* const inSR,
                                       art::RunPrincipal*& outR,
                                       art::SubRunPrincipal*& outSR,
                                       art::EventPrincipal*& outE) {

  // Establish default 'results'
  outR = 0;
  outSR = 0;
  outE = 0;

  art::ServiceHandle<dune::H4BeamFileService> beamFileService;
  if (fAccessedFirstEvent) {
    //This class can do this because it's a friend
    beamFileService->IncrementEvent();
  }
  if (beamFileService->GetCurrentEvent() >= fNEventsAvailable) return false;
  size_t run_id = beamFileService->GetRun();
  size_t subrun_id = beamFileService->GetSubRun();

  // make new run if inR is 0 or if the run has changed
  if (inR == 0 || inR->run() != run_id) {
    outR = pmaker.makeRunPrincipal(run_id, 0);
  }
  // make new subrun if inSR is 0 or if the subrun has changed
  art::SubRunID subrun_check(run_id, subrun_id);
  if (inSR == 0 || subrun_check != inSR->subRunID()) {
    outSR = pmaker.makeSubRunPrincipal(run_id, subrun_id, 0);
  }

  size_t event = beamFileService->GetCurrentEvent();
  outE = pmaker.makeEventPrincipal(run_id, subrun_id, event, 0, false);

  //Increment here
  //This class can do this because it's a friend
  //beamFileService->IncrementEvent();

  fAccessedFirstEvent = true;
  return true;
}

//typedef for shorthand
namespace dune {
  using H4BeamInputSource = art::Source<H4BeamInputDetail>;
}

DEFINE_ART_INPUT_SOURCE(dune::H4BeamInputSource)
