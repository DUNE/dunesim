// StuckBitAdcDistortionService.h
//
// David Adams
// December 2015
//
// Service that adds stuck bits to ADCs.
//
// FCL parameters:
//    StuckBitsProbabilitiesFname - Input file name
//    StuckBitsOverflowProbHistoName  - Name of histogram with overflow probablilities
//    StuckBitsUnderflowProbHistoName - Name of histogram with underflow probablilities
//    RandomSeed - Overrides NuRandomService if set nonzero.
//    LogLevel - (0=none, 1=init only, ...)

#ifndef StuckBitAdcDistortionService_H
#define StuckBitAdcDistortionService_H

#include "dune/DuneInterface/AdcDistortionService.h"
#include <string>
#include <iostream>
#include "dune/DuneInterface/AdcTypes.h"

namespace fhicl {
class ParameterSet;
}
namespace art {
class ActivityRegistry;
}
namespace CLHEP {
class HepRandomEngine;
}

class StuckBitAdcDistortionService : public AdcDistortionService {

public:

  typedef unsigned int Channel;

  // Ctor.
  StuckBitAdcDistortionService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Dtor.
  ~StuckBitAdcDistortionService();

  // Modify an input ADC vector.
  int modify(Channel chan, AdcCountVector& adcvec) const;

  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="  ") const;

private:

  // Configuration.
  std::string fStuckBitsProbabilitiesFname;     ///< file holding ADC stuck code probabilities 
  std::string fStuckBitsOverflowProbHistoName;  ///< Name of hist with ADC stuck code overflow probs 
  std::string fStuckBitsUnderflowProbHistoName; ///< Name of hist with ADC stuck code underflow probs 
  int m_RandomSeed;
  int m_LogLevel;

  double      fOverflowProbs[64];     ///< array of probs for LSF bits getting stuck at 000000
  double      fUnderflowProbs[64];    ///< array of probs for LSF bits getting stuck at 111111

  CLHEP::HepRandomEngine* m_pran;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(StuckBitAdcDistortionService, AdcDistortionService, LEGACY)

#endif
