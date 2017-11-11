// ProtoDUNEChannelNoiseService.cxx

#include "dune/DetSim/Service/ProtoDUNEChannelNoiseService.h"
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
#include "TF1.h"
#include "TMath.h"

using std::cout;
using std::ostream;
using std::endl;
using std::string;
using std::ostringstream;
using rndm::NuRandomService;
using CLHEP::HepJamesRandom;

//**********************************************************************

ProtoDUNEChannelNoiseService::
ProtoDUNEChannelNoiseService(fhicl::ParameterSet const& pset)
: fRandomSeed(0), fLogLevel(1),
  fNoiseHistZ(nullptr), fNoiseHistU(nullptr), fNoiseHistV(nullptr),
  fNoiseChanHist(nullptr),
  m_pran(nullptr) {
  const string myname = "ProtoDUNEChannelNoiseService::ctor: ";
  fLowCutoffZ        = pset.get<double>("LowCutoffZ");
  fLowCutoffU        = pset.get<double>("LowCutoffU");
  fLowCutoffV        = pset.get<double>("LowCutoffV");
  fWhiteNoiseZ       = pset.get<double>("WhiteNoiseZ");
  fWhiteNoiseU       = pset.get<double>("WhiteNoiseU");
  fWhiteNoiseV       = pset.get<double>("WhiteNoiseV");
  fNoiseArrayPoints  = pset.get<unsigned int>("NoiseArrayPoints");
  fOldNoiseIndex     = pset.get<bool>("OldNoiseIndex");
  fWirelengthZ       = pset.get<double>("WireLengthZ");
  fWirelengthU       = pset.get<double>("WireLengthU");
  fWirelengthV       = pset.get<double>("WireLengthV");
  fENOB              = pset.get<double>("EffectiveNBits");
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
  // Assign a unique name for the random number engine ProtoDUNEChannelNoiseServiceVIII
  // III = for each instance of this class.
  string rname = "ProtoDUNEChannelNoiseService";
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
    generateNoise(fWirelengthZ, fENOB, fLowCutoffZ, fNoiseZ[isam], fNoiseHistZ);
    generateNoise(fWirelengthU, fENOB, fLowCutoffU, fNoiseU[isam], fNoiseHistU);
    generateNoise(fWirelengthV, fENOB, fLowCutoffV, fNoiseV[isam], fNoiseHistV);
  }
  if ( fLogLevel > 1 ) print() << endl;
}

//**********************************************************************

ProtoDUNEChannelNoiseService::
ProtoDUNEChannelNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: ProtoDUNEChannelNoiseService(pset) { }

//**********************************************************************

ProtoDUNEChannelNoiseService::~ProtoDUNEChannelNoiseService() {
  const string myname = "ProtoDUNEChannelNoiseService::dtor: ";
  if ( fLogLevel > 0 ) {
    cout << myname << "Deleting random engine with seed " << m_pran->getSeed() << endl;
  }
  delete m_pran;
}

//**********************************************************************

int ProtoDUNEChannelNoiseService::addNoise(Channel chan, AdcSignalVector& sigs) const {
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

ostream& ProtoDUNEChannelNoiseService::print(ostream& out, string prefix) const {
  out << prefix << "ProtoDUNEChannelNoiseService: " << endl;
  out << prefix << "        LowCutoffZ: " << fLowCutoffZ << endl;
  out << prefix << "        LowCutoffU: " << fLowCutoffU << endl;
  out << prefix << "       WireLengthZ: " << fWirelengthZ  << endl;
  out << prefix << "       WireLengthU: " << fWirelengthU  << endl;
  out << prefix << "       WireLengthV: " << fWirelengthV  << endl;
  out << prefix << "       EffectiveNBits: " << fENOB  << endl;
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

void ProtoDUNEChannelNoiseService::
generateNoise(float wirelength, float ENOB, float aLowCutoff,
              AdcSignalVector& noise, TH1* aNoiseHist) const {
  const string myname = "ProtoDUNEChannelNoiseService::generateNoise: ";
  if ( fLogLevel > 1 ) {
    cout << myname << "Generating noise." << endl;
    if ( fLogLevel > 2 ) {
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
//  double lofilter = 0.;
  double phase = 0.;
  double rnd[3] = {0.};
  
  ////////////////////////////// MicroBooNE noise model/////////////////////////////////
  // vars
  double params[1] = {0.};
  double fitpar[9] = {0.};
  double wldparams[2] = {0.};
  
  // calculate FFT magnitude of noise from ENOB
  double baseline_noise = std::sqrt(ntick*1.0/12)*std::pow(2, 12 - fENOB);
  // wire length dependence function 
  TF1* _wld_f = new TF1("_wld_f", "[0] + [1]*x", 0.0, 1000);
  // custom poisson  
  TF1* _poisson = new TF1("_poisson", "[0]**(x) * exp(-[0]) / ROOT::Math::tgamma(x+1.)", 0, 30);
  // gain function in MHz
  TF1* _pfn_f1 = new TF1("_pfn_f1", "([0]*1/(x*[8]/2 + 10) + ([1]*exp(-0.5*(((x*[8]/2)-[2])/[3])**2)*exp(-0.5*pow(x*[8]/(2*[4]),[5])))*[6])+[7]", 0.0, (double)ntick);
  
  // set data-driven parameters
  // poisson mean
  params[0] = 3.30762;
  //wire length dependence parameters
  wldparams[0] = 0.395;
  wldparams[1] = 0.001304;
  _wld_f->SetParameters(wldparams);
  double wldValue = _wld_f->Eval(wirelength);
  fitpar[0] = 9.27790e+02;
  fitpar[1] = 1.20284e+07;
  fitpar[2] = 4.93692e+03;
  fitpar[3] = 1.03438e+03;
  fitpar[4] = 2.33306e+02;
  fitpar[5] = 1.36605e+00;
  fitpar[6] = wldValue;
  fitpar[7] = baseline_noise;
  fitpar[8] = ntick;
  _pfn_f1->SetParameters(fitpar);
  _poisson->SetParameters(params);
  	
  // width of frequencyBin in kHz
  double binWidth = 1.0/(ntick*sampleRate*1.0e-6);
  for ( unsigned int i=0; i<ntick/2+1; ++i ) {
    //MicroBooNE noise model
    double pfnf1val = _pfn_f1->Eval((double)i*binWidth/1000); //convert to MHz
    // define FFT parameters
    double randomizer = _poisson->GetRandom()/params[0];
    pval = pfnf1val * randomizer;
    // low frequency cutoff     
    //lofilter = 1.0/(1.0+exp(-(i-aLowCutoff/binWidth)/0.5));
    // randomize 10%
    flat.fireArray(2, rnd, 0, 1);
    //pval *= lofilter*(0.9 + 0.2*rnd[0]);
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
  for ( unsigned int itck=0; itck<noise.size(); ++itck ) {
    noise[itck] = /*sqrt(ntick)**/tmpnoise[itck];
  }
  for ( unsigned int itck=0; itck<noise.size(); ++itck ) {
    aNoiseHist->Fill(noise[itck]);
  }
}

//**********************************************************************

void ProtoDUNEChannelNoiseService::generateNoise() {
  fNoiseZ.resize(fNoiseArrayPoints);
  fNoiseU.resize(fNoiseArrayPoints);
  fNoiseV.resize(fNoiseArrayPoints);
  for ( unsigned int inch=0; inch<fNoiseArrayPoints; ++inch ) {
    generateNoise(fWirelengthZ, fENOB, fLowCutoffZ, fNoiseZ[inch], fNoiseHistZ);
    generateNoise(fWirelengthU, fENOB, fLowCutoffU, fNoiseZ[inch], fNoiseHistU);
    generateNoise(fWirelengthV, fENOB, fLowCutoffV, fNoiseZ[inch], fNoiseHistV);
  }
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(ProtoDUNEChannelNoiseService, ChannelNoiseService)

//**********************************************************************
