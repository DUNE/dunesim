// WhiteChannelNoiseService
//
// David Adams
// April 2016
// Implementation of TPC channel noise model with white noise.
// Same as the nose model 2 in SimWireDUNE35t, e.g. from dunetpc v04_29_01.
// FCL parameters:
//    RandomSeed - Overrides NuRandomService if set nonzero.
//    LogLevel - (0=none, 1=init only, ...)

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

  // Dtor.
  ~WhiteChannelNoiseService();

  // Add noise to a signal array.
  int addNoise(Channel chan, AdcSignalVector& sigs) const;

  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:
 
  // Configuration.
  int m_RandomSeed;
  int m_LogLevel;

  // Histograms.
  TH1* fNoiseHist;      ///< distribution of noise counts

  CLHEP::HepRandomEngine* m_pran;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(WhiteChannelNoiseService, ChannelNoiseService, LEGACY)

#endif
