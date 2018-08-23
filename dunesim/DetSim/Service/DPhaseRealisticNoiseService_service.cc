////////////////////////////////////////////////////////////////////////////////
//
// Implement a service including realistic noise using 3x1x1 data
//
// A realistic nose model can be read from .root file, containing the FFT transform
// for each channel. In the fft array is stored average frequency
// on the chosen plane.
//
//
// NB: Hardcoded rotated geometry!
//
////////////////////////////////////////////////////////////////////////////////

#include "dune/DetSim/Service/DPhaseRealisticNoiseService.h"
#include <sstream>
#include <string>
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/LArFFT.h"
#include "larcore/Geometry/Geometry.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "TProfile.h"
#include "TFile.h"
#include "TKey.h"
#include "TF1.h"
#include "TH2F.h"
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

DPhaseRealisticNoiseService::DPhaseRealisticNoiseService(fhicl::ParameterSet const& pset)
: fRandomSeed(0), fLogLevel(1),
  fNoiseHistX(nullptr), fNoiseHistY(nullptr),
  fNoiseChanHist(nullptr),
  m_pran(nullptr) {
  const string myname = "DPhaseRealisticNoiseService::ctor: ";
  fNoiseModel        = pset.get<string>("NoiseModel");
  fRandomizeX        = pset.get<double>("RandomizeX");
  fRandomizeY        = pset.get<double>("RandomizeY");
  fSmooth            = pset.get<double>("Smooth");
  fSetFirst0         = pset.get<bool>("SetFirst0");
  fSetBaseline          = pset.get<bool>("SetBaseline");
  fNoiseArrayPoints  = pset.get<unsigned int>("NoiseArrayPoints");
  fOldNoiseIndex     = pset.get<bool>("OldNoiseIndex");
  bool haveSeed = pset.get_if_present<int>("RandomSeed", fRandomSeed);
  if ( fRandomSeed == 0 ) haveSeed = false;
  pset.get_if_present<int>("LogLevel", fLogLevel);
  fNoiseX.resize(fNoiseArrayPoints);
  fNoiseY.resize(fNoiseArrayPoints);
  int seed = fRandomSeed;
  art::ServiceHandle<art::TFileService> tfs;
  fNoiseHistX = tfs->make<TH1F>("xnoise", ";X Noise [ADC counts];", 1000,   -0.1, 0.1);
  fNoiseHistY = tfs->make<TH1F>("ynoise", ";Y Noise [ADC counts];", 1000,   -0.1, 0.1);
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
    generateNoise(fNoiseModelFrequenciesX,  fNoiseX[isam], fNoiseHistX, fRandomizeX);
    generateNoise(fNoiseModelFrequenciesY,  fNoiseY[isam], fNoiseHistY, fRandomizeY);
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

int DPhaseRealisticNoiseService::addNoise(Channel chan, AdcSignalVector& sigs) const{

  CLHEP::RandFlat flat(*m_pran);

  //define the baseline fluctuations
  //done here because chan is needed to tune the phase
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  unsigned int ntick = detprop->NumberTimeSamples();
  double dt = 1./detprop->SamplingRate();
  unsigned int model_tick = GetModelSize();
  std::map<Channel, double>  fPhaseChannelMap;
  Chan2Phase(fPhaseChannelMap);

  //define a short osclillating baseline function
  double params[3] = {0.};
  TF1 * _sin = new TF1("_sin", "[0]*TMath::Sin( [1]*x + [2] )", 0.,
                                                             (double)ntick*dt);
  params[0] = 1.; //in adc
  params[1] = (2*TMath::Pi())/(2*model_tick*dt);
  params[2] = fPhaseChannelMap[chan];

  _sin->SetParameters(params);

  //add one of the pregenerated random noise vectors to the signal vector
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

  //sigs.resize( detprop->NumberTimeSamples() );

  for ( unsigned int itck=0; itck<ntick; ++itck ) {
    double tnoise = 0.0;
    if ( view==geo::kY ) {
      tnoise = fNoiseY[noisechan][itck];
    } else if ( view==geo::kZ ) {
      tnoise = fNoiseX[noisechan][itck];
    } else {
      tnoise = fNoiseX[noisechan][itck];
    }

    if(fSetBaseline){
      sigs[itck] += _sin->Eval((double)itck*dt);
    }

    sigs[itck] += tnoise;
  }

  return 0;
}

//**********************************************************************

ostream& DPhaseRealisticNoiseService::print(ostream& out, string prefix) const {

  out << prefix << "DPhaseRealisticNoiseService:  " <<  endl;
  out << prefix << " Noise model source file:     " << fNoiseModel   << endl;
  out << prefix << " RandomizationX:              " << fRandomizeX   << endl;
  out << prefix << " RandomizationY:              " << fRandomizeY   << endl;
  out << prefix << " Smooth                       " << fSmooth       << endl;
  out << prefix << " fSetFirst0                   " << fSetFirst0    << endl;
  out << prefix << " fSetBaseline                 " << fSetBaseline  << endl;
  out << prefix << " RandomSeed:                  " << fRandomSeed   << endl;
  out << prefix << " LogLevel:                    " << fLogLevel     << endl;
  out << prefix << " Actual random seed:          " << m_pran->getSeed();

  return out;
}

//**********************************************************************

void DPhaseRealisticNoiseService::importNoiseModel(std::string noiseModel,
std::vector<double> & frequencyArrayX, std::vector<double> & frequencyArrayY) const {

  //import noise model
  art::ServiceHandle<geo::Geometry> geo;

  TFile *fin = TFile::Open(noiseModel.c_str(), "READ");
  if(!fin->IsOpen()){
      cout << "ERROR!: Can't open file: " << noiseModel << endl;
      return;
  }
  else{
      cout << "File: " << noiseModel << " successfully opened!" << endl;
  }
  //get the histogram
  TIter next( fin->GetListOfKeys() );
  TKey *key;
  TProfile *inputHist = new TProfile(); //default initialisazion

  while( (key = (TKey*)next()) ){

    string name( key->GetName() );
    std::string keyName( key->GetClassName() ); //parse char to string
    if( keyName == "TProfile"){
      inputHist = (TProfile*)fin->Get(key->GetName());
    }
    else{
      std::cout << "ERROR! Object: " << keyName << " in file " << noiseModel
                                     << "has not the right format!" << std::endl;
      fin->Close();
      return;
    }

    //loop over frequencies: eliminate the last bin of the fft
    //if geometry is 3x1x1dp, then use both models, otherwise just take 0

    frequencyArrayX.resize( inputHist->GetNbinsX() );
    frequencyArrayY.resize( inputHist->GetNbinsX() );

    for(size_t f=0; f<(size_t)inputHist->GetNbinsX() ; f++){
        if (name.find("_0")!=string::npos){

          //NB: use just the 3m strips view for Far Detector sim
          frequencyArrayY.at(f) = inputHist->GetBinContent(f);
          frequencyArrayX.at(f) = inputHist->GetBinContent(f);
        }
        else if(name.find("_1")!=string::npos){
          //frequencyArrayX.at(f) = inputHist->GetBinContent(f);
          continue;
        }
        else{
          cout << "not valid view " << endl;
          continue;
        }
    } //end f loop
  } //end event loop

  //one may want to set the first bin of the model to 0
  if(fSetFirst0){
    frequencyArrayX[0]=0;
    frequencyArrayY[0]=0;
  }

  fin->Close();
  return;
}

//**********************************************************************

unsigned int DPhaseRealisticNoiseService::GetModelSize() const{
  return fModelsize;
}

void DPhaseRealisticNoiseService::SetModelSize(unsigned int size){
  fModelsize = size;
}

//**********************************************************************

void DPhaseRealisticNoiseService::Chan2Phase(std::map<Channel, double> &PhaseChannelMap) const{
  //create an association between channel and phase.

  CLHEP::RandFlat flat(*m_pran);
  art::ServiceHandle<geo::Geometry> geo;

  double phase = flat.fire();
  double dph = (TMath::Pi()*0.5)/62; //phase small increment. (May be improved)

  for(Channel chan=0; chan< geo->Nchannels() ; chan++){

    if(chan % 64 == 0){
      phase = flat.fire();
    }else if(chan % 32){
      dph = -dph;
      phase += dph;
    }else{
      phase += dph;
    }
    PhaseChannelMap[chan] = phase*2*TMath::Pi();
  }//end for
  return;
}

//**********************************************************************

double DPhaseRealisticNoiseService::GetShift(AdcSignalVector time_vector,
                                                       int window_length) const{
  //get the shift requeired for mirroring the waveform

  int size = time_vector.size();
  double sum=0.;
  for(int i = size - window_length; i< size; i++){ sum+=time_vector[i]; }
  double shift = sum/window_length;

  return shift *=2;
}

//**********************************************************************

void DPhaseRealisticNoiseService::mirrorWaveform(AdcSignalVector& noise,
                                                         int TimeSamples) const{

  //extend by mirroring the waveform length from the  model to match
  // the one of the detector

  int ArraySize = noise.size();
  int n_window = fSmooth; //tick window to smooth the mirroring
  double shift =0.;

  //initialize first shift value from the mean in the time window
  shift = GetShift(noise, n_window);

  noise.resize(2*ArraySize);

  //do first mirror
  for(int s = ArraySize; s<2*ArraySize; s++){
    if(s>TimeSamples){ break; }
    noise[s] = -noise[2*ArraySize -s - 1 ] + shift;
  }

  //if after the first pass the detector time window is not reached, do more...
  if(TimeSamples > 2*ArraySize){

    int n = ceil(TimeSamples/ArraySize);
    int n_pass = ceil(log(n)/log(2))+1;

    for(int pass=1; pass<n_pass; pass++){

      if(noise.size() > (size_t)TimeSamples){ break; } //no need for an extra pass

      shift = GetShift(noise, n_window);

      noise.resize( pow(2,(pass+1))*ArraySize );

      for(int s = pow(2,pass)*ArraySize; s<pow(2,(pass+1))*ArraySize; s++){
        if(s > TimeSamples){ break; }
        noise[s] = -noise[pow(2,(pass+1))*ArraySize -s - 1 ] + shift;
      }
    }//end for pass
  }//end if

  //Correct the waveform inclination
  double y, sx = 0, sy = 0, sxy = 0, sx2 = 0, sy2 = 0;

  for (size_t s = 0; s < (size_t)TimeSamples; ++s)
  {
      y = noise[s]; sx += s; sy += y; sxy += s*y; sx2 += s*s; sy2 += y*y;
  }
  double ssx = sx2 - ((sx * sx) / TimeSamples);
  double c = sxy - ((sx * sy) / TimeSamples);
  double mx = sx / TimeSamples;
  double my = sy / TimeSamples;
  double b = my - ((c / ssx) * mx);
  double a = c / ssx;

  for (size_t s = 0; s < (size_t)TimeSamples; ++s) { noise[s] -= (a*s + b); }

  return;
}//end mirrorWaveform

//**********************************************************************

void DPhaseRealisticNoiseService::generateNoise(std::vector<double> frequencyVector,
                AdcSignalVector& noise, TH1* aNoiseHist, double customRandom){

  const string myname = "DPRealisticNoiseService::generateNoise: ";
  if ( fLogLevel > 1 ) {
    cout << myname << "Generating noise." << endl;
  }

  CLHEP::RandFlat flat(*m_pran);
  CLHEP::RandGauss gaus(*m_pran);

  // rayleigh
  //TF1* _rayleigh = new TF1("_rayleigh", "[0]*( x/([1]*[1]) )*TMath::Exp(-( (x*x)/(2*[1]*[1]) ))", 0, 200);
  //_rayleigh->SetParameter(0, 7.80e+5);

  unsigned int inputFreqSize = frequencyVector.size();
  unsigned int pointFFT = 2*(inputFreqSize-1);

  art::ServiceHandle<util::LArFFT> pfft;
  pfft->ReinitializeFFT( pointFFT ," ", pfft->FFTFitBins() );

  unsigned int model_tick = pfft->FFTSize(); //ticks of the time model
  unsigned nbin = model_tick/2 + 1;
  SetModelSize( model_tick );

  std::vector<TComplex> noiseFrequency(nbin, 0.);
  double pval = 0.;
  double phase = 0.;

  for ( unsigned int i=0; i<model_tick/2+1; ++i ) {

    pval = frequencyVector[i];

    //Randomize amplitude (Gaussian)
    pval += customRandom*gaus.fire();

    //Randomize phase (Constant)
    phase = flat.fire()*2.*TMath::Pi();
    TComplex tc(pval*cos(phase),pval*sin(phase));
    noiseFrequency[i] += tc;
  }

  //Perform the InvFFT of the randomized noise frequency model
  noise.clear();
  noise.resize(model_tick, 0.0);
  std::vector<double> tmpnoise(noise.size());
  pfft->DoInvFFT(noiseFrequency, tmpnoise);
  noiseFrequency.clear();

  //model is not normalized
  for ( unsigned int itck=0; itck<noise.size(); ++itck ) {
    noise[itck] = tmpnoise[itck];
  }

  //here the signal is mirrored until it match the number of time samples for
  //a given detector geometry
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  int ntick = detprop->NumberTimeSamples();

  if(noise.size()<(size_t)ntick){
    mirrorWaveform(noise, ntick);
  }

  for ( unsigned int itck=0; itck<noise.size(); ++itck ) {
    aNoiseHist->Fill(noise[itck]);
  }

  //restore fft with the same number of points as it was before
  pfft->ReinitializeFFT( ntick ," ", pfft->FFTFitBins() );

  return;
}//end generateNoise

//**********************************************************************

void DPhaseRealisticNoiseService::generateNoise() {
  fNoiseX.resize(fNoiseArrayPoints);
  fNoiseY.resize(fNoiseArrayPoints);
  for ( unsigned int inch=0; inch<fNoiseArrayPoints; ++inch ) {
    generateNoise(fNoiseModelFrequenciesX, fNoiseX[inch], fNoiseHistX, fRandomizeX);
    generateNoise(fNoiseModelFrequenciesY, fNoiseY[inch], fNoiseHistY, fRandomizeY);
  }
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(DPhaseRealisticNoiseService, ChannelNoiseService)

//**********************************************************************
