////////////////////////////////////////////////////////////////////////////////
//
// Implement a service including realistic noise using 3x1x1 data
//
// A realistic nose model can be read from .root file, containing the FFT transform
// for each channel. In the fft array is stored average frequency
// on the chosen plane.
//
// Custom randomization of the power spectrum can be configured
//
// NB: Hardcoded rotated geometry! when rotated geo is committed,
// line ** has to be replaced in the following way:     if ( view==geo::kY ) {
//
// TODO: Average frequency spectrum of random ch frequency spectrum?
// TODO: Effect of padding and normalisation
// TODO: Randomisation
//
////////////////////////////////////////////////////////////////////////////////

#include "dune/DetSim/Service/DPhaseRealisticNoiseService.h"
#include <sstream>
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/LArFFT.h"
#include "larcore/Geometry/Geometry.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "TH1F.h"
#include "TFile.h"
#include "TKey.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TRandom3.h"

using std::cout;
using std::ostream;
using std::endl;
using std::string;
using std::ostringstream;
using rndm::NuRandomService;
using CLHEP::HepJamesRandom;

//**********************************************************************

DPhaseRealisticNoiseService::
DPhaseRealisticNoiseService(fhicl::ParameterSet const& pset)
: fRandomSeed(0), fLogLevel(1),
  fNoiseHistX(nullptr), fNoiseHistY(nullptr),
  fNoiseChanHist(nullptr),
  m_pran(nullptr) {
  const string myname = "DPhaseRealisticNoiseService::ctor: ";
  fNoiseModel        = pset.get<string>("NoiseModel");
  fRandomizeX        = pset.get<double>("RandomizeX");
  fRandomizeY        = pset.get<double>("RandomizeY");
  fWhiteNoiseX       = pset.get<double>("WhiteNoiseX");
  fWhiteNoiseY       = pset.get<double>("WhiteNoiseY");
  fNoiseArrayPoints  = pset.get<unsigned int>("NoiseArrayPoints");
  fOldNoiseIndex     = pset.get<bool>("OldNoiseIndex");
  bool haveSeed = pset.get_if_present<int>("RandomSeed", fRandomSeed);
  if ( fRandomSeed == 0 ) haveSeed = false;
  pset.get_if_present<int>("LogLevel", fLogLevel);
  fNoiseX.resize(fNoiseArrayPoints);
  fNoiseY.resize(fNoiseArrayPoints);
  int seed = fRandomSeed;
  art::ServiceHandle<art::TFileService> tfs;
  fNoiseHistX = tfs->make<TH1F>("unoise", ";X Noise [ADC counts];", 1000,   -10., 10.);
  fNoiseHistY = tfs->make<TH1F>("vnoise", ";Y Noise [ADC counts];", 1000,   -10., 10.);
  fNoiseChanHist = tfs->make<TH1F>("NoiseChan", ";Noise channel;", fNoiseArrayPoints, 0, fNoiseArrayPoints);
  // Assign a unique name for the random number engine ExponentialChannelNoiseServiceVIII
  // III = for each instance of this class.
  string rname = "DPhaseRealisticNoiseService";
  if ( haveSeed ) {
    if ( fLogLevel > 0 ) cout << myname << "WARNING: Using hardwired seed." << endl;
    m_pran = new HepJamesRandom(seed);
  } else {
    if ( fLogLevel > 0 ) cout << myname << "Using NuRandomService." << endl;
    art::ServiceHandle<NuRandomService> seedSvc;
    m_pran = new HepJamesRandom;
    if ( fLogLevel > 0 ) cout << myname << "    Initial seed: " << m_pran->getSeed() << endl;
    seedSvc->registerEngine(NuRandomService::CLHEPengineSeeder(m_pran), rname);
  }
  if ( fLogLevel > 0 ) cout << myname << "  Registered seed: " << m_pran->getSeed() << endl;

  importNoiseModel(fNoiseModel, fNoiseModelFrequenciesX, fNoiseModelFrequenciesY);

  //pregenerated noise model waveforms based on the realistic noise frequency pattern
  for ( unsigned int isam=0; isam<fNoiseArrayPoints; ++isam ) {
    generateNoise(fNoiseModelFrequenciesX, fRandomizeX, fNoiseX[isam], fNoiseHistX);
    generateNoise(fNoiseModelFrequenciesY, fRandomizeY, fNoiseY[isam], fNoiseHistY);
  }
  if ( fLogLevel > 1 ) print() << endl;
} //m_pran(nullptr)

//**********************************************************************

DPhaseRealisticNoiseService::
DPhaseRealisticNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: DPhaseRealisticNoiseService(pset) { }

//**********************************************************************

DPhaseRealisticNoiseService::~DPhaseRealisticNoiseService() {
  const string myname = "DPhaseRealisticNoiseService::dtor: ";
  if ( fLogLevel > 0 ) {
    cout << myname << "Deleting random engine with seed " << m_pran->getSeed() << endl;
  }
  delete m_pran;
}

//**********************************************************************

int DPhaseRealisticNoiseService::addNoise(Channel chan, AdcSignalVector& sigs) const {
  //TODO: can we introduce here more realistic simulation of noise for the noisy channels

  CLHEP::RandFlat flat(*m_pran);
  CLHEP::RandGauss gaus(*m_pran);
  unsigned int noisechan = 0;
  if ( fOldNoiseIndex ) {
    // Keep this strange way of choosing noise channel to be consistent with old results.
    // The relative weights of the first and last channels are 0.5 and 0.6.
    noisechan = nearbyint(flat.fire()*(1.*(fNoiseArrayPoints-1)+0.1));
  } else {
    noisechan = flat.fire()*fNoiseArrayPoints;
    if ( noisechan == fNoiseArrayPoints ) --noisechan;
  }
  fNoiseChanHist->Fill(noisechan);
  art::ServiceHandle<geo::Geometry> geo;
  const geo::View_t view = geo->View(chan);
  for ( unsigned int itck=0; itck<sigs.size(); ++itck ) {
    double tnoise = 0.0;
    double wnoise = 0.0;
    if ( view==geo::kZ ) {
      tnoise = fNoiseX[noisechan][itck];
      wnoise = fWhiteNoiseX;
    } else if ( view==geo::kY ) {
      tnoise = fNoiseY[noisechan][itck];
      wnoise = fWhiteNoiseY;
    }
    if ( wnoise != 0.0 ) tnoise += wnoise*gaus.fire();
    sigs[itck] += tnoise;
  }
  return 0;
}

//**********************************************************************

ostream& DPhaseRealisticNoiseService::print(ostream& out, string prefix) const {
  out << prefix << "ExponentialChannelNoiseService: " << endl;

  out << prefix << "        NoiseModel:     " << fNoiseNormZ   << endl;
  out << prefix << "        RandomizationX: " << fRandomizeX   << endl;
  out << prefix << "        RandomizationY: " << fRandomizeY   << endl;
  out << prefix << "        WhiteNoiseU:    " << fWhiteNoiseX  << endl;
  out << prefix << "        WhiteNoiseV:    " << fWhiteNoiseY  << endl;
  out << prefix << "        RandomSeed:     " << fRandomSeed   << endl;
  out << prefix << "        LogLevel:       " << fLogLevel     << endl;
  out << prefix << "    Actual random seed: " << m_pran->getSeed();

  return out;
}

//**********************************************************************

void DPhaseRealisticNoiseService::DPhaseRealisticNoiseService::
importNoiseModel(std::string noiseModel, std::vector<double> & frequencyArrayX,
                                  std::vector<double> & frequencyArrayY) const {
  //Read TFile using the string. Not sure what is the best way to do it within larsoft
  TFile *fin = TFile::Open(noiseModel.c_str(), "READ");
  if(!fin->IsOpen()){
      cout << "ERROR!: Can't open file: " << noiseModel << endl;
      return 0;
  }
  else{
    if(fLogLevel > 1)
      cout << "File: " << noiseModel << " successfully opened!" << endl;
  }
  //get the histogram
  TIter next(fin->GetListOfKeys());
  TKey *key;
  while(key = (TKey*)next()){
    if( key->GetClassName() == "TH2D"){
      TH2D *inputHist = (TH1D*)key->Get(key->GetName());
    }else if (key->GetClassName() == "TProfile2D" ){
      TH2D *inputHist = (TH1D*)key->Get(key->GetName());
    }else{
      cout << "ERROR! Object: " << key->GetName() << " in file " << noiseModel
                                        << "has not the right format!" << endl;
    }
  }

  //TODO: Add Crosscheck to verify the histogram is in the correct form
  //Only histograms containing 311 fft map with format ch vs frequencies in MHz
  //can be read at the moment

  //histogram based on 311 contains 1280 bins in x (num of channels for the 311)
  //and 1025 bins along Y (compatible with )
  //sample freqency is 2.5 MHz

  art::ServiceHandle<geo::Geometry> geo;
  const geo::View_t view;

  frequencyArrayX.resize(h->GetNbinsY());
  frequencyArrayY.resize(h->GetNbinsX());

  double averagedFreq=0.;

  for(size_t f =0; f < h->GetNbinsY(); f++){
    double sum =0.;
    for(size_t ch =0; f< h->GetNbinsX(); ch++){
      sum+=h->GetBinContent(ch, f);
    }
    averagedFreq = sum/h->GetNbinsX();
    view = geo->View(ch);
    if(view == geo::kZ){
      frequencyArrayX.at(f) = avergatedFreq;
    }else if(view == geo::kY){
      frequencyArrayX.at(f) = avergatedFreq;
    }else{
      cout << "ERROR: invalid view " << view << endl;
    }
  }


  fin->Close();
  return 0;
}

//**********************************************************************

void DPhaseRealisticNoiseService::
generateNoise(std::vector<double> frequencyVector, AdcSignalVector& noise,
                                  TH1* aNoiseHist, double customRandom) const {
  const string myname = "ExponentialChannelNoiseService::generateNoise: ";
  if ( fLogLevel > 1 ) {
    cout << myname << "Generating noise." << endl;
  }
  // Fetch sampling rate.
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  float sampleRate = detprop->SamplingRate();
  // Fetch FFT service and # ticks.
  art::ServiceHandle<util::LArFFT> pfft;
  unsigned int ntick = pfft->FFTSize();
  CLHEP::RandFlat flat(*m_pran);
  // Manipulate realistic noise spectrum.
  unsigned nbin = ntick/2 + 1;
  std::vector<TComplex> noiseFrequency(nbin, 0.);
  double pval = 0.;
  double lofilter = 0.;
  double phase = 0.;
  double rnd[2] = {0.};

  //Padding the frequency vector to handle a time domain larger than the input
  //frequency domain
  //TODO: study the effect of padding
  //TODO: study the effect of normalisation
  frequencyVector.resize(ntick/2+1, 0.);

  for ( unsigned int i=0; i<ntick/2+1; ++i ) {
    //read noise spectrum
    pval = frequencyVector[i];
    // randomize (%) //<< TODO: understand correcty oh this could be done
    flat.fireArray(2, rnd, 0, 1);
    pval *= (1 + customRandom*rnd[0]); //<<TODO now customRandom is the max possible fluctuaton from the average. UNREALISTIC!!
    // random phase angle
    phase = rnd[1]*2.*TMath::Pi();
    TComplex tc(pval*cos(phase),pval*sin(phase));
    noiseFrequency[i] += tc;
  }
  // Obtain time spectrum from frequency spectrum.
  noise.clear();
  noise.resize(ntick,0.0);
  std::vector<double> tmpnoise(noise.size());
  pfft->DoInvFFT(noiseFrequency, tmpnoise);
  noiseFrequency.clear();
  // Multiply each noise value by ntick as the InvFFT
  // divides each bin by ntick assuming that a forward FFT
  // has already been done.
  // As in the ExponentialChannelNoiseService, the normalisation is done using the sqrt of the ticks
  for ( unsigned int itck=0; itck<noise.size(); ++itck ) {
    noise[itck] = sqrt(ntick)*tmpnoise[itck];
  }
  for ( unsigned int itck=0; itck<noise.size(); ++itck ) {
    aNoiseHist->Fill(noise[itck]);
  }
}

//**********************************************************************

void DPhaseRealisticNoiseService::generateNoise() {
  fNoiseX.resize(fNoiseArrayPoints);
  fNoiseY.resize(fNoiseArrayPoints);
  for ( unsigned int inch=0; inch<fNoiseArrayPoints; ++inch ) {
    generateNoise(fNoiseModelFrequencies, fNoiseX[inch], fNoiseHistX);
    generateNoise(fNoiseModelFrequencies, fNoiseY[inch], fNoiseHistY);
  }
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(DPhaseRealisticNoiseService, ChannelNoiseService)

//**********************************************************************
