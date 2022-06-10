// SPhaseUniqChannelNoiseService

// Jingbo Wang (jiowang@ucdavis.edu)
// March 2019
//
// Implementation of a general TPC channel noise model with:
// (1) white noise
// (2) Inherent Gaussian noise in frequency
// (3) MicroBooNE noise in frequency
// (4) Coherent noise (exponential + Gaussian) in frequency 
//     (Note a: phase at each frequency bin is randamized at the moment. Will be updated soon
//      Note b: Currently, consecutive offline channels (configurable) are grouped together and 
//              the same coherent noise waveform is assigned to channels within the same group. )
//
// The default parameters are obtained from the ProtoDUNE-SP data (run 4096)
// fcl file: dunetpc/fcl/protodune/detsim/protoDUNE_detsim_data_driven_noise.fcl
//

#ifndef SPhaseUniqChannelNoiseService_H
#define SPhaseUniqChannelNoiseService_H

#include "dunecore/DuneInterface/Service/UniqChannelNoiseService.h"
#include <vector>
#include <iostream>

class TH1;
class TF1;
namespace CLHEP {
class HepRandomEngine;
}

class SPhaseUniqChannelNoiseService : public UniqChannelNoiseService {

public:

  // Ctor.
  SPhaseUniqChannelNoiseService(fhicl::ParameterSet const& pset);

  // Ctor.
  SPhaseUniqChannelNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&);

  // Dtor.
  ~SPhaseUniqChannelNoiseService();

  // Add noise to a signal array.
  int addNoise(detinfo::DetectorClocksData const&,
               detinfo::DetectorPropertiesData const&,
               Channel chan, AdcSignalVector& sigs) const;

  void newEvent();
  unsigned int getGroupNumberFromOfflineChannel(unsigned int offlinechan) const;
  void setNumNoiseArrayPoints(unsigned int n_np, unsigned int n_cnp);
  void setIdxNoiseArrayPoints(unsigned int i_n, unsigned int i_cn);

  // Print the configuration.
  std::ostream& print(std::ostream& out =std::cout, std::string prefix ="") const;

private:
 
  // Fill the noise vectors.
  void generateNoise(detinfo::DetectorClocksData const& clockData);
  
  // Fill a noise vector.
  // Input vector contents are lost.
  // The size of the vector is obtained from the FFT service.
  void generateMicroBooNoise(detinfo::DetectorClocksData const& clockData,
                             float wirelength, float ENOB,
                     AdcSignalVector& noise, TH1* aNoiseHist) const;
  void generateGaussianNoise(detinfo::DetectorClocksData const& clockData,
                             AdcSignalVector& noise, std::vector<float> gausNorm,
	                    std::vector<float> gausMean, std::vector<float> gausSigma,
	                    TH1* aNoiseHist) const;
  void generateCoherentNoise(detinfo::DetectorClocksData const& clockData,
                             AdcSignalVector& noise, std::vector<float> gausNorm,
	                    std::vector<float> gausMean, std::vector<float> gausSigma,
	                    float cohExpNorm, float cohExpWidth, float cohExpOffset, 
	                    TH1* aNoiseHist) const;
  
  // Make coherent groups
  void makeCoherentGroupsByOfflineChannel(unsigned int nchpergroup);
  std::vector<unsigned int> fChannelGroupMap;   ///< assign each channel a group number
  std::vector<int> fGroupCoherentNoiseMap; ///< assign each group a noise 
  //unsigned int getGroupNumberFromOfflineChannel(unsigned int offlinechan) const;
  unsigned int getCohNoiseChanFromGroup(unsigned int cohgroup) const;
  
  // General parameters
  bool         fGenUniqNoise;      ///< generate unique noise per channel per event
  bool         fSpecifyNoisePoint; ///< generate specified number of unique noise points 
  unsigned int fNoiseArrayPoints;  ///< number of points in randomly generated noise array
  unsigned int fIdxNoiseArrayPoint;
  int          fRandomSeed;        ///< Seed for random number service. If absent or zero, use SeedSvc.
  int          fLogLevel;          ///< Log message level: 0=quiet, 1=init only, 2+=every event
  
  // Inherent White noise parameters
  bool         fEnableWhiteNoise;
  float        fWhiteNoiseZ;       ///< Level (per freq bin) for white noise for Z.
  float        fWhiteNoiseU;       ///< Level (per freq bin) for white noise for U.
  float        fWhiteNoiseV;       ///< Level (per freq bin) for white noise for V.  
  
  // Inherent Gaussian noise
  bool         fEnableGaussianNoise;
  std::vector<float>	fGausNormU;		///<  noise scale factor for the gaussian component in coherent noise
  std::vector<float>	fGausMeanU;		///<  mean of the gaussian component in coherent noise
  std::vector<float>	fGausSigmaU;	///<  sigma of the gaussian component in coherent noise
  std::vector<float>	fGausNormV;		///<  noise scale factor for the gaussian component in coherent noise
  std::vector<float>	fGausMeanV;		///<  mean of the gaussian component in coherent noise
  std::vector<float>	fGausSigmaV;	///<  sigma of the gaussian component in coherent noise
  std::vector<float>	fGausNormZ;		///<  noise scale factor for the gaussian component in coherent noise
  std::vector<float>	fGausMeanZ;		///<  mean of the gaussian component in coherent noise
  std::vector<float>	fGausSigmaZ;	///<  sigma of the gaussian component in coherent noise
  
  // Inherent MicroBoone noise parameters
  bool         fEnableMicroBooNoise;  ///< enable MicroBooNE noise model
  float        fENOB;                 ///< Effective number of bits
  float        fWirelengthZ;
  float        fWirelengthU;
  float        fWirelengthV;
  std::vector<float>  fNoiseFunctionParameters;  ///< Parameters in the MicroBooNE noise model
  TF1* fUBNoisePoisFnc;
  TF1* fUBNoiseWrLDFnc;
  TF1* fUBNoiseGainFnc;

  // Coherent Noise parameters
  bool         fEnableCoherentNoise;
  unsigned int fNChannelsPerCoherentGroup;
  unsigned int fNCoherentGroups;
  unsigned int fExpNoiseArrayPoints;  ///< number of points in randomly generated noise array
  unsigned int fCohNoiseArrayPoints;  ///< number of points in randomly generated noise array
  unsigned int fIdxCohNoiseArrayPoint;
  float        fCohExpNorm;           ///< noise scale factor for the exponential component component in coherent noise
  float        fCohExpWidth;          ///< width of the exponential component in coherent noise
  float        fCohExpOffset;         ///< Amplitude offset of the exponential background component in coherent noise
  std::vector<float>	fCohGausNorm;		///<  noise scale factor for the gaussian component in coherent noise
  std::vector<float>	fCohGausMean;		///<  mean of the gaussian component in coherent noise
  std::vector<float>	fCohGausSigma;	///<  sigma of the gaussian component in coherent noise 
  TF1 *fCohNoiseFnc;
  TF1* fCohPoisFnc;
  
  // Inherent Gausian noise arrays
  AdcSignalVectorVector fGausNoiseZ;
  AdcSignalVectorVector fGausNoiseU;
  AdcSignalVectorVector fGausNoiseV;
  
  // Inherent MicroBoo noise arrays
  AdcSignalVectorVector fMicroBooNoiseZ;  
  AdcSignalVectorVector fMicroBooNoiseU; 
  AdcSignalVectorVector fMicroBooNoiseV;
  
  // Coherent Noise array.
  AdcSignalVectorVector fCohNoise;  ///< noise on each channel for each time for all planes
  

  // Histograms.
  
  TH1* fGausNoiseHistZ;     ///< distribution of noise counts for Z
  TH1* fGausNoiseHistU;     ///< distribution of noise counts for U
  TH1* fGausNoiseHistV;     ///< distribution of noise counts for V
  TH1* fGausNoiseChanHist;  ///< distribution of accessed noise samples

  TH1* fMicroBooNoiseHistZ;     ///< distribution of noise counts for Z
  TH1* fMicroBooNoiseHistU;     ///< distribution of noise counts for U
  TH1* fMicroBooNoiseHistV;     ///< distribution of noise counts for V
  TH1* fMicroBooNoiseChanHist;  ///< distribution of accessed noise samples
  
  TH1* fCohNoiseHist;      ///< distribution of noise counts
  TH1* fCohNoiseChanHist;  ///< distribution of accessed noise samples

  CLHEP::HepRandomEngine* m_pran;

};

DECLARE_ART_SERVICE_INTERFACE_IMPL(SPhaseUniqChannelNoiseService, UniqChannelNoiseService, LEGACY)

#endif
