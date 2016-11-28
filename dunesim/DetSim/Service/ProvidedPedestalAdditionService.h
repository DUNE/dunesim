// ProvidedPedestalAdditionService.h

// David Adams
// December 2015
//
// Pedestal additon service that obtains pedestals and their
// RMS values from a pedestal service.
// Pedestal noise is assigned using independent random Gaussian values
// for each tick with sigma = NoiseScale*RMS
// Parameters:
//    NoiseScale - The above noise scale. (<= 0 for no pedestal noise).
//    RandomSeed - Overrides NuRandomService if set nonzero.
//    LogLevel - (0=none, 1=init only, ...)

#ifndef ProvidedPedestalAdditionService_H
#define ProvidedPedestalAdditionService_H

#include "dune/DuneInterface/PedestalAdditionService.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

class TH1;
namespace CLHEP {
class HepRandomEngine;
}
namespace lariov {
class DetPedestalProvider;
}

class ProvidedPedestalAdditionService : public PedestalAdditionService {

public:

  ProvidedPedestalAdditionService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  ~ProvidedPedestalAdditionService();

  // Add pedestals and pedestal noise to a signal array.
  int addPedestal(Channel chan, AdcSignalVector& sigs, float& ped, float& pedrms) const;

  // Print parameters.
  virtual std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  // Configuration.
  float m_NoiseScale;
  int m_RandomSeed;
  int m_LogLevel;

  TH1* m_PedNoiseHist;   ///< Histogram of pedestal noise counts

  CLHEP::HepRandomEngine* m_pran;

  const lariov::DetPedestalProvider& m_PedestalProvider;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(ProvidedPedestalAdditionService, PedestalAdditionService, LEGACY)

#endif

