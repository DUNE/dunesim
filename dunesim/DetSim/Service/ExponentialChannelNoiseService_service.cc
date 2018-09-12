// ExponentialChannelNoiseService.cxx

#include "dune/DetSim/Service/ExponentialChannelNoiseService.h"
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
#include "TRandom3.h"

using std::cout;
using std::ostream;
using std::endl;
using std::string;
using std::ostringstream;
using rndm::NuRandomService;
using CLHEP::HepJamesRandom;

//**********************************************************************

ExponentialChannelNoiseService::
ExponentialChannelNoiseService(fhicl::ParameterSet const& pset)
: fRandomSeed(0), fLogLevel(1),
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
  fWhiteNoiseZ       = pset.get<double>("WhiteNoiseZ");
  fWhiteNoiseU       = pset.get<double>("WhiteNoiseU");
  fWhiteNoiseV       = pset.get<double>("WhiteNoiseV");
  fNoiseArrayPoints  = pset.get<unsigned int>("NoiseArrayPoints");
  fOldNoiseIndex     = pset.get<bool>("OldNoiseIndex");
  bool haveSeed = pset.get_if_present<int>("RandomSeed", fRandomSeed);
  if ( fRandomSeed == 0 ) haveSeed = false;
  pset.get_if_present<int>("LogLevel", fLogLevel);
  fNoiseZ.resize(fNoiseArrayPoints);
  fNoiseU.resize(fNoiseArrayPoints);
  fNoiseV.resize(fNoiseArrayPoints);
  int seed = fRandomSeed;
  art::ServiceHandle<art::TFileService> tfs;
  fNoiseHistZ = tfs->make<TH1F>("znoise", ";Z Noise [ADC counts];", 1000,   -10., 10.);
  fNoiseHistU = tfs->make<TH1F>("unoise", ";U Noise [ADC counts];", 1000,   -10., 10.);
  fNoiseHistV = tfs->make<TH1F>("vnoise", ";V Noise [ADC counts];", 1000,   -10., 10.);
  fNoiseChanHist = tfs->make<TH1F>("NoiseChan", ";Noise channel;", fNoiseArrayPoints, 0, fNoiseArrayPoints);
  // Assign a unique name for the random number engine ExponentialChannelNoiseServiceVIII
  // III = for each instance of this class.
  string rname = "ExponentialChannelNoiseService";
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
  for ( unsigned int isam=0; isam<fNoiseArrayPoints; ++isam ) {
    generateNoise(fNoiseNormZ, fNoiseWidthZ, fLowCutoffZ, fNoiseZ[isam], fNoiseHistZ);
    generateNoise(fNoiseNormU, fNoiseWidthU, fLowCutoffU, fNoiseU[isam], fNoiseHistU);
    generateNoise(fNoiseNormV, fNoiseWidthV, fLowCutoffV, fNoiseV[isam], fNoiseHistV);
  }
  if ( fLogLevel > 1 ) print() << endl;
}

//**********************************************************************

ExponentialChannelNoiseService::
ExponentialChannelNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: ExponentialChannelNoiseService(pset) { }

//**********************************************************************

ExponentialChannelNoiseService::~ExponentialChannelNoiseService() {
  const string myname = "ExponentialChannelNoiseService::dtor: ";
  if ( fLogLevel > 0 ) {
    cout << myname << "Deleting random engine with seed " << m_pran->getSeed() << endl;
  }
  delete m_pran;
}

//**********************************************************************

int ExponentialChannelNoiseService::addNoise(Channel chan, AdcSignalVector& sigs) const {
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
    if ( view==geo::kU ) {
      tnoise = fNoiseU[noisechan][itck];
      wnoise = fWhiteNoiseU;
    } else if ( view==geo::kV ) {
      tnoise = fNoiseV[noisechan][itck];
      wnoise = fWhiteNoiseV;
    } else {
      tnoise = fNoiseZ[noisechan][itck];
      wnoise = fWhiteNoiseZ;
    }
    if ( wnoise != 0.0 ) tnoise += wnoise*gaus.fire();
    sigs[itck] += tnoise;
  }
  return 0;
}

//**********************************************************************

ostream& ExponentialChannelNoiseService::print(ostream& out, string prefix) const {
  out << prefix << "ExponentialChannelNoiseService: " << endl;
  out << prefix << "        NoiseNormZ: " << fNoiseNormZ   << endl;
  out << prefix << "       NoiseWidthZ: " << fNoiseWidthZ  << endl;
  out << prefix << "        LowCutoffZ: " << fLowCutoffZ << endl;
  out << prefix << "        NoiseNormU: " << fNoiseNormU   << endl;
  out << prefix << "       NoiseWidthU: " << fNoiseWidthU  << endl;
  out << prefix << "        LowCutoffU: " << fLowCutoffU << endl;
  out << prefix << "        NoiseNormV: " << fNoiseNormV   << endl;
  out << prefix << "       NoiseWidthV: " << fNoiseWidthV  << endl;
  out << prefix << "        LowCutoffV: " << fLowCutoffV << endl;
  out << prefix << "  NoiseArrayPoints: " << fNoiseArrayPoints << endl;
  out << prefix << "     OldNoiseIndex: " << fOldNoiseIndex << endl;
  out << prefix << "       WhiteNoiseZ: " << fWhiteNoiseZ << endl;
  out << prefix << "       WhiteNoiseU: " << fWhiteNoiseU << endl;
  out << prefix << "       WhiteNoiseV: " << fWhiteNoiseV << endl;
  out << prefix << "        RandomSeed: " <<  fRandomSeed << endl;
  out << prefix << "          LogLevel: " <<  fLogLevel << endl;
  out << prefix << "  Actual random seed: " << m_pran->getSeed();
  return out;
}

//**********************************************************************

void ExponentialChannelNoiseService::
generateNoise(float aNoiseNorm, float aNoiseWidth, float aLowCutoff,
              AdcSignalVector& noise, TH1* aNoiseHist) const {
  const string myname = "ExponentialChannelNoiseService::generateNoise: ";
  if ( fLogLevel > 1 ) {
    cout << myname << "Generating noise." << endl;
    if ( fLogLevel > 2 ) {
      cout << myname << "    Norm: " << aNoiseNorm << endl;
      cout << myname << "   Width: " << aNoiseWidth << endl;
      cout << myname << "  Cutoff: " << aLowCutoff << endl;
      cout << myname << "    Seed: " << m_pran->getSeed() << endl;
    }
  }
  // Fetch sampling rate.
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  float sampleRate = detprop->SamplingRate();
  // Fetch FFT service and # ticks.
  art::ServiceHandle<util::LArFFT> pfft;
  unsigned int ntick = pfft->FFTSize();
  CLHEP::RandFlat flat(*m_pran);
  // Create noise spectrum in frequency.
  unsigned nbin = ntick/2 + 1;
  std::vector<TComplex> noiseFrequency(nbin, 0.);
  double pval = 0.;
  double lofilter = 0.;
  double phase = 0.;
  double rnd[2] = {0.};
  // width of frequencyBin in kHz
  double binWidth = 1.0/(ntick*sampleRate*1.0e-6);
  bool flatAtLowFreq = aLowCutoff < 0.0;
  double frqCutoff = flatAtLowFreq ? -aLowCutoff : 0.0;
  for ( unsigned int i=0; i<ntick/2+1; ++i ) {
    double frq = double(i)*binWidth;
    // For aLowCutoff < 0, we have constant noise below -aLowCutoff.
    // For positive values, we keep the old albeit questionable behavior.
    if ( flatAtLowFreq ) {
      lofilter = 1.0;
      if ( frq < frqCutoff ) frq = frqCutoff;
    } else {
      lofilter = 1.0/(1.0+exp(-(i-aLowCutoff/binWidth)/0.5));
    }
    // exponential noise spectrum 
    pval = aNoiseNorm*exp(-frq/aNoiseWidth);
    // low frequency cutoff     
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
