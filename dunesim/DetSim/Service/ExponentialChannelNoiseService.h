// ExponentialChannelNoiseService

// Implementation of TPC channel noise model with an exponential
// shape in frequency.
// Same as the nose model 1 in SimWireDUNE35t, e.g. from dunetpc v04_29_01.

#ifndef ExponentialChannelNoiseService_H
#define ExponentialChannelNoiseService_H

#include "dune/DuneInterface/ChannelNoiseService.h"
#include <vector>
#include <iostream>

class TH1;
namespace CLHEP {
class HepRandomEngine;
}

class ExponentialChannelNoiseService : public ChannelNoiseService {

public:

  // Ctor.
  ExponentialChannelNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Add noise to a signal array.
  int addNoise(Channel chan, AdcSignalVector& sigs) const;

  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

  // Fill a noise vector.
  // Input vector contents are lost.
  // The size of the vector is obtained from the FFT service.
  void generateNoise(float aNoiseFact, float aNoiseWidth, float aLowCutoff, AdcSignalVector& noise) const;

private:
 
  // Fill the noise vectors.
  void generateNoise();

  // Parameters.
  float        fNoiseFactZ;        ///< noise scale factor for Z (collection) plane
  float        fNoiseWidthZ;       ///< exponential noise width (kHz)  for Z (collection) plane
  float        fLowCutoffZ;        ///< low frequency filter cutoff (kHz) for Z (collection) plane
  float        fNoiseFactU;        ///< noise scale factor  for U plane
  float        fNoiseWidthU;       ///< exponential noise width (kHz)   for U plane
  float        fLowCutoffU;        ///< low frequency filter cutoff (kHz)  for U plane
  float        fNoiseFactV;        ///< noise scale factor   for V plane
  float        fNoiseWidthV;       ///< exponential noise width (kHz)   for V plane
  float        fLowCutoffV;        ///< low frequency filter cutoff (kHz)  for V plane
  unsigned int fNoiseArrayPoints;  ///< number of points in randomly generated noise array
  bool         fOldNoiseIndex;     ///< Use old selection of noise array index

  // Noise arrays.
  AdcSignalVectorVector fNoiseZ;  ///< noise on each channel for each time for Z (collection) plane
  AdcSignalVectorVector fNoiseU;  ///< noise on each channel for each time for U plane
  AdcSignalVectorVector fNoiseV;  ///< noise on each channel for each time for V plane

  // Histograms.
  TH1* fNoiseHist;      ///< distribution of noise counts
  TH1* fNoiseChanHist;  ///< distribution of accessed noise samples

  CLHEP::HepRandomEngine* m_pran;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(ExponentialChannelNoiseService, ChannelNoiseService, LEGACY)

#endif
