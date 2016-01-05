// WhiteChannelNoiseService

// Implementation of TPC channel noise model with white noise.
// Same as the nose model 2 in SimWireDUNE35t, e.g. from dunetpc v04_29_01.

#ifndef WhiteChannelNoiseService_H
#define WhiteChannelNoiseService_H

#include "dune/DuneInterface/ChannelNoiseService.h"
#include <vector>
#include <iostream>

class TH1;
namespace CLHEP {
class HepRandomEngine;
}

class WhiteChannelNoiseService : public ChannelNoiseService {

public:

  // Ctor.
  WhiteChannelNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Add noise to a signal array.
  int addNoise(Channel chan, AdcSignalVector& sigs) const;

  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:
 
  // Histograms.
  TH1* fNoiseHist;      ///< distribution of noise counts

  CLHEP::HepRandomEngine* m_pran;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(WhiteChannelNoiseService, ChannelNoiseService, LEGACY)

#endif
