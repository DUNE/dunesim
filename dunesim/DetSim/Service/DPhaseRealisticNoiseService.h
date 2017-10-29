////////////////////////////////////////////////////////////////////////////////
// DPhaseRealisticNoiseService
//
// Andrea Scarpelli (andrea.scarpelli@cern.ch)
// Ocotober 2017
//
// Implementation of realistic noise from frequency spectrum in 3x1x1 data
// using same logic of ExponentialChannelNoiseService
///////////////////////////////////////////////////////////////////////////////

#ifndef DPhaseRealisticNoiseService_H
#define DPhaseRealisticNoiseService_H

#include "dune/DuneInterface/ChannelNoiseService.h"
#include <vector>
#include <iostream>

class TH1;
namespace CLHEP {
class HepRandomEngine;
}

class DPhaseRealisticNoiseService : public ChannelNoiseService {

public:

  // Ctor.
  DPhaseRealisticNoiseService(fhicl::ParameterSet const& pset);

  // Ctor.
  DPhaseRealisticNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Dtor.
  ~DPhaseRealisticNoiseService();

  // Add noise to a signal array.
  int addNoise(Channel chan, AdcSignalVector& sigs) const;

  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

  //import a model of realistic noise fft from a root file
  //Store it in an array for each view 
  importNoiseModel(std::string noiseModel, std::vector<double> & frequencyArrayX,
                                    std::vector<double> & frequencyArrayY) const

  // Fill a noise vector.
  // Input vector contents are lost.
  // The size of the vector is obtained from the FFT service.
  void generateNoise(float aNoiseNorm, float aNoiseWidth, float aLowCutoff,
                     AdcSignalVector& noise, TH1* aNoiseHist) const;

private:

  // Fill the noise vectors.
  void generateNoise();

  // Parameters.
  std::string  fNoiseModel;        ///< noise model root file
  unsigned int fNoiseArrayPoints;  ///< number of points in randomly generated noise array
  double       fRandomizeX;        ///< randomization of the average frequency spectrum  (on kX or kZ)
  double       fRandomizeY;        ///<< randomization of the average frequency spectrum (on kY)
  bool         fOldNoiseIndex;     ///< Use old selection of noise array index
  float        fWhiteNoiseX;       ///< Level (per freq bin) for white noise for X.
  float        fWhiteNoiseY;       ///< Level (per freq bin) for white noise for Y.
  int          fRandomSeed;        ///< Seed for random number service. If absent or zero, use SeedSvc.
  int          fLogLevel;          ///< Log message level: 0=quiet, 1=init only, 2+=every event

  //frequency arrays imported
  std::vector<double> fNoiseModelFrequenciesX   ///< Array storing the frequencies imported from the model in kHz for plane kX (kZ)
  std::vector<double> fNoiseModelFrequenciesY   ///< Array storing the frequencies imported from the model in kHz for plane kY (kZ)

  // Noise arrays.
  AdcSignalVectorVector fNoiseX;  ///< noise on each channel for each time for X plane
  AdcSignalVectorVector fNoiseY;  ///< noise on each channel for each time for Y plane

  // Histograms.
  TH1* fNoiseHist;      ///< distribution of noise counts
  TH1* fNoiseHistX;     ///< distribution of noise counts for X
  TH1* fNoiseHistY;     ///< distribution of noise counts for Y
  TH1* fNoiseChanHist;  ///< distribution of accessed noise samples

  CLHEP::HepRandomEngine* m_pran;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(DPhaseRealisticNoiseService, ChannelNoiseService, LEGACY)

#endif
