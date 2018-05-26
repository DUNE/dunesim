////////////////////////////////////////////////////////////////////////////////
// DPhaseRealisticNoiseService
//
// Andrea Scarpelli (andrea.scarpelli@cern.ch)
// October 2017
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

  // Fill a noise vector.
  // Input vector contents are lost.
  // The size of the vector is obtained from the FFT service.
  void generateNoise(std::vector<double> frequencyVector, AdcSignalVector& noise,
                                    TH1* aNoiseHist, double customRandom);

private:

  //tool getters for use the model size everywhere in the code
  unsigned int GetModelSize() const;

  //tool setter for use the model size everywhere in the code
  void SetModelSize(unsigned int size);

  // Tool function used  mirrorWaveform();
  double GetShift(AdcSignalVector time_vector, int window_length) const;

  // Tool producing associations between channels and semi-random phase;
  void Chan2Phase(std::map<Channel, double> &PhaseChannelMap) const;

  // Fill the noise vectors.
  void generateNoise();

  //Mirror the output waveform to match geometries dofferent from 3x1x1 geo
  void mirrorWaveform(AdcSignalVector& noise, int TimeSamples) const;

  //import a model of realistic noise fft from a root file
  //Store it in an array for each view
  void importNoiseModel(std::string noiseModel, std::vector<double> & frequencyArrayX,
                                    std::vector<double> & frequencyArrayY) const;

  // Parameters.
  std::string  fNoiseModel;        ///< noise model root file
  unsigned int fNoiseArrayPoints;  ///< number of points in randomly generated noise array
  double       fRandomizeX;        ///< randomization of the average frequency spectrum  (on kX or kZ)
  double       fRandomizeY;        ///<< randomization of the average frequency spectrum (on kY)
  double       fSmooth;
  bool         fSetFirst0;         ///<< set the first bin of the frequency array to 0
  bool         fSetBaseline;       ///<< Sum baseline model to the data
  bool         fOldNoiseIndex;     ///< Use old selection of noise array index
  int          fRandomSeed;        ///< Seed for random number service. If absent or zero, use SeedSvc.
  int          fLogLevel;          ///< Log message level: 0=quiet, 1=init only, 2+=every event

  //frequency arrays imported
  std::vector<double> fNoiseModelFrequenciesX;   ///< Array storing the frequencies imported from the model in kHz for plane kX (kZ)
  std::vector<double> fNoiseModelFrequenciesY;   ///< Array storing the frequencies imported from the model in kHz for plane kY (kZ)

  // Noise arrays.
  AdcSignalVectorVector fNoiseX;  ///< noise on each channel for each time for X plane
  AdcSignalVectorVector fNoiseY;  ///< noise on each channel for each time for Y plane

  // Histograms.
  //TH1* fNoiseHist;      ///< distribution of noise counts // unused
  TH1* fNoiseHistX;     ///< distribution of noise counts for X
  TH1* fNoiseHistY;     ///< distribution of noise counts for Y
  TH1* fNoiseChanHist;  ///< distribution of accessed noise samples

  //more
  //double fPhase; // unused
  unsigned int fModelsize;

  CLHEP::HepRandomEngine* m_pran;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(DPhaseRealisticNoiseService, ChannelNoiseService, LEGACY)

#endif
