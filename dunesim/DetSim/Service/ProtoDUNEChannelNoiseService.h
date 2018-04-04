// ProtoDUNEChannelNoiseService

// David Adams
// January 2016
//
// Implementation of TPC channel noise model with an exponential
// shape in frequency.
// Same as the nose model 1 in SimWireDUNE35t, e.g. from dunetpc v04_29_01.
//
// DLA Feb 2016: Change normalization so RMS does not vary with FFT size.
// See https://cdcvs.fnal.gov/redmine/issues/11470.

#ifndef ProtoDUNEChannelNoiseService_H
#define ProtoDUNEChannelNoiseService_H

#include "dune/DuneInterface/ChannelNoiseService.h"
#include <vector>
#include <iostream>

class TH1;
namespace CLHEP {
class HepRandomEngine;
}

class ProtoDUNEChannelNoiseService : public ChannelNoiseService {

public:

  // Ctor.
  ProtoDUNEChannelNoiseService(fhicl::ParameterSet const& pset);

  // Ctor.
  ProtoDUNEChannelNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Dtor.
  ~ProtoDUNEChannelNoiseService();

  // Add noise to a signal array.
  int addNoise(Channel chan, AdcSignalVector& sigs) const;

  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

  // Fill a noise vector.
  // Input vector contents are lost.
  // The size of the vector is obtained from the FFT service.
  void generateNoise(float wirelength, float ENOB, float aLowCutoff,
                     AdcSignalVector& noise, TH1* aNoiseHist) const;

private:
 
  // Fill the noise vectors.
  void generateNoise();

  // Parameters.
  float        fLowCutoffZ;        ///< low frequency filter cutoff (kHz) for Z (collection) plane
  float        fLowCutoffU;        ///< low frequency filter cutoff (kHz)  for U plane
  float        fLowCutoffV;        ///< low frequency filter cutoff (kHz)  for V plane
  unsigned int fNoiseArrayPoints;  ///< number of points in randomly generated noise array
  bool         fOldNoiseIndex;     ///< Use old selection of noise array index
  float        fWhiteNoiseZ;       ///< Level (per freq bin) for white noise for Z.
  float        fWhiteNoiseU;       ///< Level (per freq bin) for white noise for U.
  float        fWhiteNoiseV;       ///< Level (per freq bin) for white noise for V.
  int          fRandomSeed;        ///< Seed for random number service. If absent or zero, use SeedSvc.
  int          fLogLevel;          ///< Log message level: 0=quiet, 1=init only, 2+=every event
  float        fWirelengthZ;
  float        fWirelengthU;
  float        fWirelengthV;
  float        fENOB;

  // Noise arrays.
  AdcSignalVectorVector fNoiseZ;  ///< noise on each channel for each time for Z (collection) plane
  AdcSignalVectorVector fNoiseU;  ///< noise on each channel for each time for U plane
  AdcSignalVectorVector fNoiseV;  ///< noise on each channel for each time for V plane

  // Histograms.
  //TH1* fNoiseHist;      ///< distribution of noise counts // unused
  TH1* fNoiseHistZ;     ///< distribution of noise counts for Z
  TH1* fNoiseHistU;     ///< distribution of noise counts for U
  TH1* fNoiseHistV;     ///< distribution of noise counts for V
  TH1* fNoiseChanHist;  ///< distribution of accessed noise samples

  CLHEP::HepRandomEngine* m_pran;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(ProtoDUNEChannelNoiseService, ChannelNoiseService, LEGACY)

#endif
