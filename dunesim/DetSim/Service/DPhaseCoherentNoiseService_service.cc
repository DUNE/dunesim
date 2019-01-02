////////////////////////////////////////////////////////////////////////////////
//
// Implement a service including realistic noise using 3x1x1 data
//
// A realistic nose model can be read from .root file, containing the channel vs.
// frequency map for the 3x1x1 detector.
//
// Frequency are injected as sinusoidal functions, with realistic amplitudes and
// phases randomized to respect the correlation patterns inside the detector
//
// NB: kZ = view 1; kY = view 0
//
// mailto:andrea.scarpelli@cern.ch
//
////////////////////////////////////////////////////////////////////////////////

#include "dune/DetSim/Service/DPhaseCoherentNoiseService.h"
#include <sstream>
#include <string>
//#include "art/canvas/Utilities/Exception.h"
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
#include "TProfile2D.h"
#include "TRandom3.h"

using std::cout;
using std::ostream;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::ostringstream;
using rndm::NuRandomService;
using CLHEP::HepJamesRandom;


//**********************************************************************

DPhaseCoherentNoiseService::DPhaseCoherentNoiseService(fhicl::ParameterSet const& pset)
: fRandomSeed(0), fLogLevel(1),
  m_pran(nullptr)
{
  const string myname = "DPhaseCoherentNoiseService::ctor: ";
  fNoiseModel = pset.get<string>("NoiseModel");
  fAmplitudeCut = pset.get<double>("AmplitudeCut");
  fNormalization = pset.get<int>("Normalization");
  fRandomize = pset.get< vector< float > >("Randomize");
  fPhaseShift = pset.get< vector< float > >("PhaseShift");
  fChannelGroup = pset.get< vector<int> >("ChannelGroup");
  fInchoerentNoise = pset.get< vector< float > >("InchoerentNoise");
  fNumberOfPhases = pset.get<int>("NumberOfPhases");
  fLogLevel = pset.get<int>("LogLevel");
  bool haveSeed = pset.get_if_present<int>("RandomSeed", fRandomSeed);
  if ( fRandomSeed == 0 ) haveSeed = false;
  pset.get_if_present<int>("LogLevel", fLogLevel);

  int seed = fRandomSeed;
  art::ServiceHandle<art::TFileService> tfs;

  // Assign a unique name for the random number engine ExponentialChannelNoiseServiceVIII
  // III = for each instance of this class.
  string rname = "DPhaseCoherentNoiseService";
  if ( haveSeed )
  {
    if ( fLogLevel > 0 ) cout << myname << "WARNING: Using hardwired seed." << endl;
    m_pran = new HepJamesRandom(seed);
  }
  else
  {
    if ( fLogLevel > 0 ) cout << myname << "Using NuRandomService." << endl;
    art::ServiceHandle<NuRandomService> seedSvc;
    m_pran = new HepJamesRandom;
    if ( fLogLevel > 0 ) cout << myname << "    Initial seed: " << m_pran->getSeed() << endl;
    seedSvc->registerEngine(NuRandomService::CLHEPengineSeeder(m_pran), rname);
  }
  if ( fLogLevel > 0 ) cout << myname << "  Registered seed: " << m_pran->getSeed() << endl;

  importNoiseModel(fNoiseModel, fChFrequencyMap, fChAmplitudeMap, fAmplitudeCut, fNormalization);

  //sanity checks for the imported maps and find max frequency array length
  int max =0;

  art::ServiceHandle<geo::Geometry> geo;
  for(Channel chan=0; chan< geo->Nchannels() ; chan++){

    vector<float> amplitudeArray;
    vector<float> frequencyArray;

    getFrequencyArray( chan, frequencyArray);
    getAmplitudeArray( chan, amplitudeArray);

    if( frequencyArray.size() != amplitudeArray.size() ){
      if ( fLogLevel > 0 ){ cout << myname << "frequency array and amplitude array have not the same size in chan: " << chan << endl; }
      //throw art::Exception("DPhaseCoherentNoiseService") << "frequency array and amplitude array have not the same size";
    }

    if( max < (int)frequencyArray.size() ){ max = frequencyArray.size(); }
  }

  fMaxFrequencySize = max;

  //pregenerated random phases arrays
  fPhaseMap.resize(fNumberOfPhases);

  for ( int isam=0; isam<fNumberOfPhases; ++isam ) {
    makePhaseMap(fPhaseMap[isam], fMaxFrequencySize, fPhaseShift[0], fPhaseShift[1]);
  }

  if ( fLogLevel > 1 ) print() << endl;
} //m_pran(nullptr)

//**********************************************************************

DPhaseCoherentNoiseService::
DPhaseCoherentNoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: DPhaseCoherentNoiseService(pset) { }

//**********************************************************************

DPhaseCoherentNoiseService::~DPhaseCoherentNoiseService() {
  const string myname = "DPhaseCoherentNoiseService::dtor: ";
  if ( fLogLevel > 0 ) {
    cout << myname << "Deleting random engine with seed " << m_pran->getSeed() << endl;
  }
  delete m_pran;
}

//**********************************************************************

int DPhaseCoherentNoiseService::addNoise(Channel chan, AdcSignalVector& sigs) const {

  const string myname = "DPhaseCoherentNoiseService::addNoise: ";
  if ( fLogLevel > 0 ) {
    cout << myname << " Processing channel: " << chan  << endl;
  }

  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  unsigned int ntick = detprop->NumberTimeSamples();

  art::ServiceHandle<geo::Geometry> geo;
  const geo::View_t view = geo->View(chan);

  //get the map associated to channel and phase array associated to channel
  int num = getNumber( chan );
  Map phaseMap = fPhaseMap.at(num);
  vector<float> phaseArray = phaseMap.at(chan);

  if ( fLogLevel > 0 ) {
    cout << myname << " Map number: " << num  << endl;
  }

  //get amplitude and frequecny
  vector<float> amplitudeArray;
  vector<float> frequencyArray;

  getFrequencyArray( chan, frequencyArray);
  getAmplitudeArray( chan, amplitudeArray);

  //Test if the channel was in the fft model
  if( frequencyArray.size() == 0 && amplitudeArray.size() ==0 ){
    vector<float> nullArray;
    nullArray.resize( ntick );
    sigs = nullArray;

    if ( fLogLevel > 0 ) {
      cout << myname << " No freq and amplitude arrays for " << chan  << endl;
    }

    return 0;
  }


  //resize the signal vector
  sigs.resize( ntick );

  //choose the correct plane for the amplitude randomization
  float randAmp=0;

  if ( view==geo::kY )
  {
    randAmp = fRandomize[0];
  }
  else if (  view==geo::kZ ) {
    randAmp = fRandomize[1];
  }
  else
  {
    if ( fLogLevel > 0 ) {
      cout << myname << "Invalid plane"  << endl;
    }
  }

  //build the noise signal
  getNoiseArray( sigs, amplitudeArray, frequencyArray, phaseArray, randAmp);

  //add incoherent noise if set
  if( fInchoerentNoise.size() > 0 ){

    CLHEP::RandGauss gaus(*m_pran);

    if ( view==geo::kY )
    {
      for ( unsigned int itck=0; itck<ntick; ++itck )
        sigs.at(itck) += gaus.fire( fInchoerentNoise.at(0), fInchoerentNoise.at(1) );
    }
    else if (  view==geo::kZ )
    {
      for ( unsigned int itck=0; itck<ntick; ++itck )
        sigs.at(itck) += gaus.fire( fInchoerentNoise.at(2), fInchoerentNoise.at(3) );
    }
    else {
      if ( fLogLevel > 0 ) {
        cout << myname << "Invalid plane"  << endl;
      }
    }
  }

  if ( fLogLevel > 1 ) {
    cout << myname << " All done "  << endl;
  }

  return 0;

}

//**********************************************************************

ostream& DPhaseCoherentNoiseService::print(ostream& out, string prefix) const {

  out << prefix << "DPhaseCoherentNoiseService:  " <<  endl;
  out << prefix << " Noise model source file:     " << fNoiseModel   << endl;
  out << prefix << " RandomSeed:                  " << fRandomSeed   << endl;
  out << prefix << " LogLevel:                    " << fLogLevel     << endl;
  out << prefix << " Actual random seed:          " << m_pran->getSeed();

  return out;
}

//******************************************************************************

int DPhaseCoherentNoiseService::getNumber( Channel chan ) const{

  //return the correct random number to select the PhaseMap making sure it is
  //correct for the given event

  int num=0;

  if( chan == 0 )
  {
    CLHEP::RandFlat flat(*m_pran);
    num = flat.fire()*fNumberOfPhases;
    fNum=num;
  }
  else
  {
      num = fNum;
  }

  return num;
}

//******************************************************************************

void DPhaseCoherentNoiseService::getAmplitudeArray( Channel chan, vector<float> & array ) const {

  //test if there is the amplitude array for channel chan

  if ( fChAmplitudeMap.find(chan) != fChAmplitudeMap.end() ) {
    array = fChAmplitudeMap.at(chan);
  }

  return;
}

//******************************************************************************

void DPhaseCoherentNoiseService::getFrequencyArray( Channel chan, vector<float> & array ) const {

  //test if there is the amplitude array for channel chan

  if ( fChFrequencyMap.find(chan) != fChFrequencyMap.end() ) {
      array = fChFrequencyMap.at(chan);
  }

  return;
}

//******************************************************************************

void DPhaseCoherentNoiseService::importNoiseModel(std::string noiseModel,
        Map & chFrequencyMap, Map & chAmplitudeMap, double cut, int &normalization ) const  {
  /*
  Import fft from the TProfile2D and store into the frequency and amplitude
  maps if the frequencies are above the cut
  */

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
  TProfile2D *inputHist = new TProfile2D(); //default initialisazion

  while( (key = (TKey*)next()) ){

    string name( key->GetName() );
    string keyName( key->GetClassName() ); //parse char to string
    if( keyName == "TProfile2D"){
      inputHist = (TProfile2D*)fin->Get(key->GetName());
    }
    else{
      std::cout << "ERROR! Object: " << keyName << " in file " << noiseModel
                                   << "has not the right format!" << std::endl;
      fin->Close();
      return;
    }
  }

  //if the model was successfully build the frequency/channel map
  for(int xx=1; xx<inputHist->GetNbinsX()+1; xx++){
    for(int yy=1; yy<inputHist->GetNbinsY()+1; yy++){

      double amplitude = inputHist->GetBinContent( xx, yy );
      double ch = (int)inputHist->GetXaxis()->GetBinLowEdge(xx)+1;  //channel numbered fom 1 to 1280
      double frequency = inputHist->GetYaxis()->GetBinCenter(yy);   //frequencies fom 0 to 1.25 MHz

      if(amplitude >= cut){
        chFrequencyMap[ch].push_back( frequency );
        chAmplitudeMap[ch].push_back( amplitude );
      }
    }
  }

  //get the normalization
  //normalization = inputHist->GetNbinsY();

  return;
}

//******************************************************************************

void DPhaseCoherentNoiseService::makePhaseMap( Map & phaseMap ,int size,
                                               float minShift, float maxShift ){

  //Simpy generates an map with phases from 0 to 2*pi with the same size of size

  CLHEP::RandFlat flat(*m_pran);

  vector<float> phaseVector;   //<<Phase shifts per group of channels
  vector<float> phaseArray;    //<<Phases common to all channels
  phaseVector.resize(size);
  phaseArray.resize(size);

  //make a new phase vector
  for( int f=0; f<size; f++ ){
    double phase = flat.fire()*2.*TMath::Pi();
    phaseArray.at(f) = phase;

  }

  //loop over all channels
  int index=0; //keep track of the groups already checked;
  art::ServiceHandle<geo::Geometry> geo;
  for(Channel chan=0; chan< geo->Nchannels() ; chan++){

    //assign each channel to a group
    if( chan >= (Channel)fChannelGroup.at(index) ){

        index++; //now we can change group

        //make a new phase vector
        for( int f=0; f<size; f++ ){
          double phase = minShift + flat.fire()*( maxShift - minShift );
          phaseVector.at(f) = phase + phaseArray.at(f);

        }
    }//end if channel group

    //add the vector to the map
    phaseMap[chan] = phaseVector;

  }//end for channels

  return;
}

//******************************************************************************

void DPhaseCoherentNoiseService::getNoiseArray(  vector< float > & noiseArray,
   vector< float > ampArray, vector< float > freqArray, vector< float > phaseArray, float randAmp ) const {

  //Sum up all the frequencies and amplitude. Make noise waveform
  const string myname = "DPhaseCoherentNoiseService::getNoiseArray: ";

  //random number generator
  CLHEP::RandGauss gaus(*m_pran);

  //detector properties service
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  //check if frequency vector, and amplitude vector have the same size.
  if( (ampArray.size() != freqArray.size()) || (ampArray.size() < freqArray.size()) ){

    const string myname = "DPhaseCoherentNoiseService::getNoiseArray: ";
    if( fLogLevel > 0 ){
      cout << myname << "ERROR: amplitude array and frequency array have not the same size." << endl;
    }
    return;
  }

  //make the noiseArray: loop over time...
  for( size_t t=0; t<noiseArray.size(); t++ ){

    //...and loop over frequencies
    for( size_t f=0; f<freqArray.size(); f++ ){

      //randomize amplitude
      double amp = ampArray.at(f) + gaus.fire( 0, randAmp  ); //<< Randomization with the expected rms fluctuation calculated from the model

      //make signal for that frequency
      double argument = 2*TMath::Pi()*detprop->SamplingRate()*(1.e-3)*freqArray.at(f)*t + phaseArray.at(f);

      noiseArray.at(t) += ( ((float)1/(float)fNormalization)*amp*sin( argument ) );
    }
  }

  //return the noise array
  return;
}

//******************************************************************************

void DPhaseCoherentNoiseService::getPhases() {

  //pregenerated random phases arrays
  fPhaseMap.resize(fNumberOfPhases);
  for ( int isam=0; isam<fNumberOfPhases; ++isam ) {
    makePhaseMap(fPhaseMap[isam], fMaxFrequencySize, fPhaseShift[0], fPhaseShift[1]);
  }

}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(DPhaseCoherentNoiseService, ChannelNoiseService)
