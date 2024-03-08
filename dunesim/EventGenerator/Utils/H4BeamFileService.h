#ifndef H4BeamFileService_H
#define H4BeamFileService_H

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "TFile.h"
#include "dunesim/EventGenerator/Utils/ProtoDUNETriggeredBeamUtils.h"

namespace dune
{

  class H4BeamFileService {
  public:
    explicit H4BeamFileService(fhicl::ParameterSet const& p, art::ActivityRegistry& areg);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // sets the pointer and assumes ownership of it

    //void SetPtr(std::unique_ptr<TFile> fileptr);
    void OpenFile(const std::string & filename);

    // gets a non-owning copy of the file pointer

    //TFile * GetPtr();
    std::map<int, evgen::BeamEvent> * GetBeamEvents() {
      return fAllBeamEvents.get(); 
    };
    std::vector<int> * GetEventIDs() {
      return fFinalTriggeredEventIDs.get();
    };

    void Reset();
    size_t GetRun() {
      return fRun;
    };

    size_t GetCurrentEvent() {
      return fCurrentEvent;
    };

    void IncrementEvent(int n=1) {fCurrentEvent += n;}

  private:

    size_t fCurrentEvent = 0;
    std::unique_ptr<std::map<int, evgen::BeamEvent>> fAllBeamEvents;
    std::unique_ptr<std::vector<int>> fFinalTriggeredEventIDs;
    evgen::ProtoDUNETriggeredBeamUtils fBeamUtils;
    size_t fRun;

    std::string fTRIG1TreeName, fTRIG2TreeName, fTOF1TreeName, fBPROF1TreeName,
                fBPROF2TreeName, fBPROF3TreeName, fBPROF4TreeName,
                fBPROFEXTTreeName, fNP04frontTreeName;
    bool fIsNP02;

  };

}

DECLARE_ART_SERVICE(dune::H4BeamFileService, LEGACY)


#endif
