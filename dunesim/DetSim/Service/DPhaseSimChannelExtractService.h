// DPhaseSimChannelExtractService.h
//
// Vyacheslav Galymov
// February 2016
//
// Simulate signals from dual-phase detector:
//      - charge amplification in gas
//      - charge collection in two views
//

#ifndef _DPhaseSimChannelExtractService_H_
#define _DPhaseSimChannelExtractService_H_

#include <vector>
#include <string>
#include "dune/DuneInterface/Service/SimChannelExtractService.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/Utilities/LArFFT.h"
#include "dune/Utilities/SignalShapingServiceDUNEDPhase.h"
#include "dune/Utilities/CrpGainService.h"
//#include "dune/DetSim/Tool/CrpGainSimTool.h"

namespace sim {
  class SimChannel;
}

//class CrpGainSimulator;

class DPhaseSimChannelExtractService : public SimChannelExtractService {

public:

  DPhaseSimChannelExtractService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  int extract(detinfo::DetectorClocksData const& clockData,
              const sim::SimChannel* psc, AdcSignalVector& sig) const;

  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

  //float GainPerView(){ return fDPGainPerView; };

private:

  // standard larsoft FFT service
  art::ServiceHandle<util::LArFFT> m_pfft;

  // dual-phase signal response service
  art::ServiceHandle<util::SignalShapingServiceDUNEDPhase> m_psss;

  //
  art::ServiceHandle<util::CrpGainService> m_crpgain;
  // tool to simulate charge gain in CRPs
  //std::string m_CrpGainToolName;
  //std::unique_ptr<CrpGainSimTool> m_CrpGainTool;

  unsigned int m_ntick;
  //float fDPGainPerView; // gain in dual-phase
};

DECLARE_ART_SERVICE_INTERFACE_IMPL(DPhaseSimChannelExtractService, SimChannelExtractService, LEGACY)

#endif
