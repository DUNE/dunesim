// Implementation of TPC channel noise model with an exponential
// plus constant shape in frequency.
// 29 Sep 2016

// m.thiesse@sheffield.ac.uk


#ifndef ExpoWhiteChannelNoiseService_H
#define ExpoWhiteChannelNoiseService_H

#include "dune/DuneInterface/ChannelNoiseService.h"
#include <vector>
#include <iostream>

class TH1;
namespace CLHEP {
class HepRandomEngine;
}

class ExpoWhiteChannelNoiseService : public ChannelNoiseService {

public:

  // Ctor.
  ExpoWhiteChannelNoiseService(fhicl::ParameterSet const& pset);

  // Ctor.
  ExpoWhiteChannelNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Dtor.
  ~ExpoWhiteChannelNoiseService();

  // Add noise to a signal array.
  int addNoise(Channel chan, AdcSignalVector& sigs) const;

  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

  // Fill a noise vector.
  // Input vector contents are lost.
  // The size of the vector is obtained from the FFT service.
  void generateNoise(float aNoiseNorm, float aNoiseWidth, float aLowCutoff,
                     AdcSignalVector& noise, TH1* aNoiseHist) const;

private:
 
  // Fill the noise vectors.
  void generateNoise();

  // Parameters.
  float        fNoiseNormZ;        ///< noise scale factor for Z (collection) plane
  float        fNoiseWidthZ;       ///< exponential noise width (kHz)  for Z (collection) plane
  float        fLowCutoffZ;        ///< low frequency filter cutoff (kHz) for Z (collection) plane
  float        fNoiseNormU;        ///< noise scale factor  for U plane
  float        fNoiseWidthU;       ///< exponential noise width (kHz)   for U plane
  float        fLowCutoffU;        ///< low frequency filter cutoff (kHz)  for U plane
  float        fNoiseNormV;        ///< noise scale factor   for V plane
  float        fNoiseWidthV;       ///< exponential noise width (kHz)   for V plane
  float        fLowCutoffV;        ///< low frequency filter cutoff (kHz)  for V plane
  unsigned int fNoiseArrayPoints;  ///< number of points in randomly generated noise array
  bool         fOldNoiseIndex;     ///< Use old selection of noise array index
  int          fRandomSeed;        ///< Seed for random number service. If absent or zero, use SeedSvc.
  int          fLogLevel;          ///< Log message level: 0=quiet, 1=init only, 2+=every event

  // Noise arrays.
  AdcSignalVectorVector fNoiseZ;  ///< noise on each channel for each time for Z (collection) plane
  AdcSignalVectorVector fNoiseU;  ///< noise on each channel for each time for U plane
  AdcSignalVectorVector fNoiseV;  ///< noise on each channel for each time for V plane

  // Histograms.
  TH1* fNoiseHist;      ///< distribution of noise counts
  TH1* fNoiseHistZ;     ///< distribution of noise counts for Z
  TH1* fNoiseHistU;     ///< distribution of noise counts for U
  TH1* fNoiseHistV;     ///< distribution of noise counts for V
  TH1* fNoiseChanHist;  ///< distribution of accessed noise samples

  CLHEP::HepRandomEngine* m_pran;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(ExpoWhiteChannelNoiseService, ChannelNoiseService, LEGACY)

#endif
