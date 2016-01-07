// SimWireDUNE_module.cc
//
// David Adams
// December 2015
//
// SimWireDUNE class designed to simulate signal on a wire in the TPC
//
// Developed from  SimWireDUNE35t_module.cc.

#include <vector>
#include <string>
#include <sstream>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "RawData/raw.h"
#include "Geometry/Geometry.h"
#include "Simulation/sim.h"
#include "Simulation/SimChannel.h"
#include "RawData/RawDigit.h"
#include "Utilities/DetectorProperties.h"

#include "dune/DuneInterface/AdcTypes.h"
#include "dune/DuneInterface/AdcSuppressService.h"
#include "dune/DuneInterface/AdcCompressService.h"
#include "dune/DuneInterface/SimChannelExtractService.h"
#include "dune/DuneInterface/ChannelNoiseService.h"
#include "dune/DuneInterface/PedestalAdditionService.h"
#include "dune/DuneInterface/AdcDistortionService.h"

using std::ostringstream;
using std::endl;

//**********************************************************************

// Base class for creation of raw signals on wires. 
class SimWireDUNE : public art::EDProducer {
    
public:
        
  explicit SimWireDUNE(fhicl::ParameterSet const& pset); 
  virtual ~SimWireDUNE();
    
  // read/write access to event
  void produce (art::Event& evt);
  void beginJob();
  void endJob();
  void reconfigure(fhicl::ParameterSet const& p);

private:

  std::string            fDriftEModuleLabel; ///< module making the ionization electrons
  
  // Flags.
  bool fNoiseOn;           ///< noise turned on or off for debugging; default is on
  bool fPedestalOn;        ///< switch for simulation of nonzero pedestals
  bool fDistortOn;         ///< switch for simulation of stuck bits
  bool fSuppressOn;        ///< switch for simulation of zero suppression
  bool fSaveEmptyChannel;  ///< switch for saving channels with all zero entries

  // Services.
  art::ServiceHandle<geo::Geometry> m_pgeo;
  art::ServiceHandle<AdcSuppressService> m_pzs;
  art::ServiceHandle<AdcCompressService> m_pcmp;
  art::ServiceHandle<SimChannelExtractService> m_pscx;
  art::ServiceHandle<ChannelNoiseService> m_pcns;
  art::ServiceHandle<AdcDistortionService> m_pdis;
  art::ServiceHandle<PedestalAdditionService> m_ppad;

}; // class SimWireDUNE

DEFINE_ART_MODULE(SimWireDUNE)

//**********************************************************************

SimWireDUNE::SimWireDUNE(fhicl::ParameterSet const& pset) {
  reconfigure(pset);
  produces< std::vector<raw::RawDigit> >();
}

//**********************************************************************

SimWireDUNE::~SimWireDUNE() { }

//**********************************************************************

void SimWireDUNE::reconfigure(fhicl::ParameterSet const& p) {
  fDriftEModuleLabel = p.get<std::string>("DriftEModuleLabel");
  fNoiseOn           = p.get<bool>("NoiseOn");
  fPedestalOn        = p.get<bool>("PedestalOn");  
  fDistortOn         = p.get<bool>("DistortOn");  
  fSuppressOn        = p.get<bool>("SuppressOn");  
  fSaveEmptyChannel  = p.get< bool >("SaveEmptyChannel");  

  ostringstream out;
  out << "  Compression service:";
  m_pcmp->print(out, "    ");
  if ( fNoiseOn ) {
    out << "  Channel noise service:" << endl;;
    m_pcns->print(out, "    ");
  } else {
    out << "  Channel noise is off.";
  }
  out << endl;
  if ( fPedestalOn ) {
    out << "  Pedestal addition service:" << endl;;
    m_ppad->print(out, "    ");
  } else {
    out << "  Pedestal addition is off.";
  }
  out << endl;
  if ( fSuppressOn ) {
    out << "  ADC suppression service:" << endl;
    m_pzs->print(out, "    ");
  } else {
    out << "ADC suppression is off.";
  }
  out << endl;
  if ( fDistortOn ) {
    out << "  ADC distortion service:" << endl;;
    m_pdis->print(out, "    ");
  } else {
    out << "  ADC distortion bits is off.";
  }
  mf::LogInfo("SimWireDUNE::reconfigure") << out.str();

  return;
}

//**********************************************************************

void SimWireDUNE::beginJob() { }

//**********************************************************************

void SimWireDUNE::endJob() { }

//**********************************************************************

void SimWireDUNE::produce(art::Event& evt) {

  // Make a vector of const sim::SimChannel* that has same number
  // of entries as the number of channels in the detector
  // and set the entries for the channels that have signal on them
  // using the chanHandle
  std::vector<const sim::SimChannel*> chanHandle;
  std::vector<const sim::SimChannel*> simChannels(m_pgeo->Nchannels(), nullptr);
  evt.getView(fDriftEModuleLabel, chanHandle);
  for ( size_t c=0; c<chanHandle.size(); ++c ) {
    simChannels[chanHandle[c]->Channel()] = chanHandle[c];
  }
    
  // make an unique_ptr of sim::SimDigits that allows ownership of the produced
  // digits to be transferred to the art::Event after the put statement below
  std::unique_ptr<std::vector<raw::RawDigit>>  digcol(new std::vector<raw::RawDigit>);
          
  // Fetch the number of ticks to write out for each channel.
  art::ServiceHandle<util::DetectorProperties> detprop;
  unsigned int nTickReadout  = detprop->ReadOutWindowSize();

  // Loop over channels.
  std::map<int,double>::iterator mapIter;      
  for ( unsigned int chan = 0; chan<m_pgeo->Nchannels(); ++chan ) {    

    // Get the SimChannel for this channel
    const sim::SimChannel* psc = simChannels[chan];
    const geo::View_t view = m_pgeo->View(chan);
    if (view != geo::kU && view != geo::kV && view != geo::kZ) {
      mf::LogError("SimWireDUNE") << "ERROR: CHANNEL NUMBER " << chan << " OUTSIDE OF PLANE";
    }

    // Create vector that holds the floating ADC count for each tick.
    std::vector<AdcSignal> fChargeWork;

    // Extract the floating ADC count from the SimChannel for each tick in the channel.
    m_pscx->extract(psc, fChargeWork);

    // Add noise to each tick in the channel.
    if ( fNoiseOn ) {              
      m_pcns->addNoise(chan, fChargeWork);
    }

    // Add pedestal.
    float pedval = 0.0; // Pedestal to be recorded in RawDigits collection
    float pedrms = 0.0; // Pedestal RMS to be recorded in RawDigits collection
    if ( fPedestalOn ) {
      m_ppad->addPedestal(chan, fChargeWork, pedval, pedrms);
    }

    // Convert floating ADC to integral ADC count.
    std::vector<short> adcvec(fChargeWork.size(), 0);        
    const short adcmax = 4095;
    for ( unsigned int itck=0; itck<fChargeWork.size(); ++itck ) {
      AdcSignal adcsig = fChargeWork[itck];
      short adc = 0;
      if ( adcsig > 0 ) adc = (short) (adcsig + 0.5);
      if ( adc > adcmax ) adc = adcmax;
      adcvec[itck] = adc;
    }
    // Resize to the correct number of time samples, dropping extra samples.
    adcvec.resize(nTickReadout);
    
    // Add stuck bits.
    if ( fDistortOn ) {
      m_pdis->modify(chan, adcvec);
    }
    
    // Zero suppress and compress.
    AdcFilterVector keep(adcvec.size(), true);
    if ( fSuppressOn ) {
      m_pzs->filter(adcvec, chan, pedval, keep);
    }
    int nkeep = 0;
    for ( bool kept : keep ) if ( kept ) ++nkeep;
    raw::Compress_t comp = raw::kNone;
    m_pcmp->compress(adcvec, keep, pedval, comp);

    // Create and store raw digit.
    raw::RawDigit rd(chan, nTickReadout, adcvec, comp);
    rd.SetPedestal(pedval, pedrms);
    digcol->push_back(rd);            // add this digit to the collection

  }  // end loop over channels      

  evt.put(std::move(digcol));
 
}
