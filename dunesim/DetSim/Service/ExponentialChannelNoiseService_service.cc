// ExponentialChannelNoiseService.cxx

#include "dune/DetSim/Service/ExponentialChannelNoiseService.h"
#include <sstream>
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArFFT.h"
#include "Geometry/Geometry.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "artextensions/SeedService/SeedService.hh"
#include "CLHEP/Random/RandFlat.h"
#include "TH1F.h"

using std::cout;
using std::ostream;
using std::endl;
using std::string;
using std::ostringstream;

#undef UseSeedService

//**********************************************************************

ExponentialChannelNoiseService::
ExponentialChannelNoiseService(fhicl::ParameterSet const& pset)
: fRandomSeed(-1), fLogLevel(1),
  fNoiseHistZ(nullptr), fNoiseHistU(nullptr), fNoiseHistV(nullptr),
  fNoiseChanHist(nullptr),
  m_pran(nullptr) {
  const string myname = "ExponentialChannelNoiseService::ctor: ";
  fNoiseNormZ        = pset.get<double>("NoiseNormZ");
  fNoiseWidthZ       = pset.get<double>("NoiseWidthZ");
  fLowCutoffZ        = pset.get<double>("LowCutoffZ");
  fNoiseNormU        = pset.get<double>("NoiseNormU");
  fNoiseWidthU       = pset.get<double>("NoiseWidthU");
  fLowCutoffU        = pset.get<double>("LowCutoffU");
  fNoiseNormV        = pset.get<double>("NoiseNormV");
  fNoiseWidthV       = pset.get<double>("NoiseWidthV");
  fLowCutoffV        = pset.get<double>("LowCutoffV");
  fNoiseArrayPoints  = pset.get<unsigned int>("NoiseArrayPoints");
  fOldNoiseIndex     = pset.get<bool>("OldNoiseIndex");
  bool haveSeed = pset.get_if_present<int>("RandomSeed", fRandomSeed);
  pset.get_if_present<int>("LogLevel", fLogLevel);
  fNoiseZ.resize(fNoiseArrayPoints);
  fNoiseU.resize(fNoiseArrayPoints);
  fNoiseV.resize(fNoiseArrayPoints);
  int seed = fRandomSeed;
#ifdef UseSeedService
  art::ServiceHandle<artext::SeedService> seedSvc;
  int seed = seedSvc->getSeed("ExponentialChannelNoiseService");
#else
  if ( ! haveSeed ) {
    cout << myname << "WARNING: Using hardwired seed." << endl;
    seed = 1005;
  } else {
    cout << myname << "WARNING: Using seed from FCL." << endl;
  }
#endif
  art::ServiceHandle<art::TFileService> tfs;
  fNoiseHistZ = tfs->make<TH1F>("znoise", ";Z Noise [ADC counts];", 1000,   -10., 10.);
  fNoiseHistU = tfs->make<TH1F>("unoise", ";U Noise [ADC counts];", 1000,   -10., 10.);
  fNoiseHistV = tfs->make<TH1F>("vnoise", ";V Noise [ADC counts];", 1000,   -10., 10.);
  fNoiseChanHist = tfs->make<TH1F>("NoiseChan", ";Noise channel;", fNoiseArrayPoints, 0, fNoiseArrayPoints);
  art::EngineCreator ecr;
  // Assign a unique name for the random number engine for each instance of this class.
  string basename = "ExponentialChannelNoiseService";
  for ( unsigned int itry=0; itry<1000; ++itry ) {
    try {
      ostringstream ssnam;
      ssnam << basename;
      if ( itry ) {
        ssnam << "V";
        if ( itry < 10 ) ssnam << 0;
        if ( itry < 100 ) ssnam << 0;
        if ( itry < 1000 ) ssnam << 0;
        ssnam << itry;
      }
      m_pran = &ecr.createEngine(seed, "HepJamesRandom", ssnam.str());
      break;
    } catch (...) {
      if ( ++itry >= 1000 ) {
        cout << myname << "Too many random number engines." << endl;
        m_pran = &ecr.createEngine(seed, "HepJamesRandom", basename);
        abort();  // Preceding should raise an exception.
      }
    }
  }
  for ( unsigned int isam=0; isam<fNoiseArrayPoints; ++isam ) {
    generateNoise(fNoiseNormZ, fNoiseWidthZ, fLowCutoffZ, fNoiseZ[isam], fNoiseHistZ);
    generateNoise(fNoiseNormU, fNoiseWidthU, fLowCutoffU, fNoiseU[isam], fNoiseHistU);
    generateNoise(fNoiseNormV, fNoiseWidthV, fLowCutoffV, fNoiseV[isam], fNoiseHistV);
  }
  if ( fLogLevel > 0 ) print();
}

//**********************************************************************

ExponentialChannelNoiseService::
ExponentialChannelNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: ExponentialChannelNoiseService(pset) { }

//**********************************************************************

int ExponentialChannelNoiseService::addNoise(Channel chan, AdcSignalVector& sigs) const {
  CLHEP::HepRandomEngine& engine = *m_pran;
  CLHEP::RandFlat flat(engine);
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
    if      ( view==geo::kU ) tnoise = fNoiseU[noisechan][itck];
    else if ( view==geo::kV ) tnoise = fNoiseV[noisechan][itck];
    else                      tnoise = fNoiseZ[noisechan][itck];
    sigs[itck] += tnoise;
  }
  return 0;
}

//**********************************************************************

ostream& ExponentialChannelNoiseService::print(ostream& out, string prefix) const {
  string myprefix = prefix + "  ";
  out << myprefix << "       NoiseNormZ: " << fNoiseNormZ   << endl;
  out << myprefix << "      NoiseWidthZ: " << fNoiseWidthZ  << endl;
  out << myprefix << "       LowCutoffZ: " << fLowCutoffZ << endl;
  out << myprefix << "       NoiseNormU: " << fNoiseNormU   << endl;
  out << myprefix << "      NoiseWidthU: " << fNoiseWidthU  << endl;
  out << myprefix << "       LowCutoffU: " << fLowCutoffU << endl;
  out << myprefix << "       NoiseNormV: " << fNoiseNormV   << endl;
  out << myprefix << "      NoiseWidthV: " << fNoiseWidthV  << endl;
  out << myprefix << "       LowCutoffV: " << fLowCutoffV << endl;
  out << myprefix << " NoiseArrayPoints: " << fNoiseArrayPoints << endl;
  out << myprefix << "    OldNoiseIndex: " << fOldNoiseIndex << endl;
  out << myprefix << "       RandomSeed: " <<  fRandomSeed << endl;
  out << myprefix << "         LogLevel: " <<  fLogLevel << endl;
  out << myprefix << " Actual random seed: " << m_pran->getSeed() << endl;
  return out;
}

//**********************************************************************

void ExponentialChannelNoiseService::
generateNoise(float aNoiseNorm, float aNoiseWidth, float aLowCutoff,
              AdcSignalVector& noise, TH1* aNoiseHist) const {
  // Fetch sampling rate.
  art::ServiceHandle<util::DetectorProperties> detprop;
  float sampleRate = detprop->SamplingRate();
  // Fetch FFT service and # ticks.
  art::ServiceHandle<util::LArFFT> pfft;
  unsigned int ntick = pfft->FFTSize();
  // Fetch random number engine.
#ifdef UseSeedService
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine& engine = rng->getEngine("ExponentialChannelNoiseService");
#else
  CLHEP::HepRandomEngine& engine = *m_pran;
#endif
  CLHEP::RandFlat flat(engine);
  // Create noise spectrum in frequency.
  unsigned nbin = ntick/2 + 1;
  std::vector<TComplex> noiseFrequency(nbin, 0.);
  double pval = 0.;
  double lofilter = 0.;
  double phase = 0.;
  double rnd[2] = {0.};
  // width of frequencyBin in kHz
  double binWidth = 1.0/(ntick*sampleRate*1.0e-6);
  for ( unsigned int i=0; i<ntick/2+1; ++i ) {
    // exponential noise spectrum 
    pval = aNoiseNorm*exp(-(double)i*binWidth/aNoiseWidth);
    // low frequency cutoff     
    lofilter = 1.0/(1.0+exp(-(i-aLowCutoff/binWidth)/0.5));
    // randomize 10%
    flat.fireArray(2, rnd, 0, 1);
    pval *= lofilter*(0.9 + 0.2*rnd[0]);
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
  // DLA Feb 2016: Change factor from ntick --> sqrt(ntick) so that the RMS 
  // does not depend on ntick (FFT size).
  for ( unsigned int itck=0; itck<noise.size(); ++itck ) {
    noise[itck] = sqrt(ntick)*tmpnoise[itck];
  }
  for ( unsigned int itck=0; itck<noise.size(); ++itck ) {
    aNoiseHist->Fill(noise[itck]);
  }
}

//**********************************************************************

void ExponentialChannelNoiseService::generateNoise() {
  fNoiseZ.resize(fNoiseArrayPoints);
  fNoiseU.resize(fNoiseArrayPoints);
  fNoiseV.resize(fNoiseArrayPoints);
  for ( unsigned int inch=0; inch<fNoiseArrayPoints; ++inch ) {
    generateNoise(fNoiseNormZ, fNoiseWidthZ, fLowCutoffZ, fNoiseZ[inch], fNoiseHistZ);
    generateNoise(fNoiseNormU, fNoiseWidthU, fLowCutoffU, fNoiseZ[inch], fNoiseHistU);
    generateNoise(fNoiseNormV, fNoiseWidthV, fLowCutoffV, fNoiseZ[inch], fNoiseHistV);
  }
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(ExponentialChannelNoiseService, ChannelNoiseService)

//**********************************************************************
