///////////////////////////////////////////////////////////////////////////////
// DPhaseCoherentNoiseService
//
// Andrea Scarpelli (mailto:andrea.scarpelli@cern.ch)
// August 2018
//
// Implementation of realistic noise from frequency spectrum in 3x1x1 data
// Phase and amplitude randomization following the correlation patterns seen in
// the 3x1x1 detector
//
///////////////////////////////////////////////////////////////////////////////

#ifndef DPhaseCoherentNoiseService_H
#define DPhaseCoherentNoiseService_H

#include "dune/DuneInterface/ChannelNoiseService.h"
#include <vector>
#include <iostream>

class TH1;
namespace CLHEP {
class HepRandomEngine;
}

class DPhaseCoherentNoiseService : public ChannelNoiseService {

public:

  //rename the most used std::map definition into Map
  typedef std::map<int, std::vector< float >>  Map;
  typedef std::map<Channel, AdcSignalVector>  Mask;

  //FFT frequency and amplitude maps
  Map fChFrequencyMap;   ///< Map storing the frequency vector for each channel
  Map fChAmplitudeMap;   ///< Map holding the amplitude of the frequency for each channel

  // Ctor.
  DPhaseCoherentNoiseService(fhicl::ParameterSet const& pset);

  // Ctor.
  DPhaseCoherentNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Dtor.
  ~DPhaseCoherentNoiseService();

  // Add noise to a signal array.
  int addNoise(Channel chan, AdcSignalVector& sigs) const;

  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:

  //
  void getAmplitudeArray( Channel chan, std::vector<float> & array ) const;

  //
  void getFrequencyArray( Channel chan, std::vector<float> & array ) const;

  //generate the array with phases
  void makePhaseMap( Map & phaseMap, int size, float minShift, float maxShift);

  //get the correct map number (one new per event)
  int getNumber( Channel chan ) const;

  //fill noise
  void getNoiseArray( std::vector< float > & noiseArray,
     std::vector< float > ampArray,  std::vector< float > freqArray, std::vector< float > phaseArray, float randAmp ) const;

  void getPhases();

  //import a model of realistic noise fft from a root file
  //Store it in an array for each view
  void importNoiseModel(std::string noiseModel,
                  Map & chFrequencyMap, Map & chAmplitudeMap, double cut, int &normalization ) const;

  // Parameters.
  std::string                 fNoiseModel;         ///<< noise model root file
  double                      fAmplitudeCut;       ///<< only frequencies with amplitude above this cut will be considered
  int                         fNormalization;      ///<< Normalization factor ( similar to the one for the InFFT )
  std::vector< float >        fRandomize;          ///<< randomization of the amplitude
  std::vector< float >        fPhaseShift;          ///<< Phase shift for each group of 30 channels
  std::vector<int>            fChannelGroup;       ///<< Channels in the same group get the same phase
  std::vector< float >        fInchoerentNoise;    ///<< Mean and std of the incoherent noise
  int                         fNumberOfPhases;     ///<< Number of pregenerated phase shift maps
  int                         fRandomSeed;         ///<< Seed for random number service. If absent or zero, use SeedSvc.
  int                         fLogLevel;           ///<< Log message level: 0=quiet, 1=init only, 2+=every event


  int fMaxFrequencySize;
  std::vector< Map > fPhaseMap;                   ///< Pregenerated phase shift maps
  mutable int fNum=0;                             ///< Hold the correct event number


  CLHEP::HepRandomEngine* m_pran;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(DPhaseCoherentNoiseService, ChannelNoiseService, LEGACY)

#endif
