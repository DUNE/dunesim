// SimWireDUNE_module.cc
//
// David Adams
// December 2015
//
// SimWireDUNE class designed to simulate signal on a wire in the TPC
//
// Developed from the now-obsolete SimWireDUNE35t_module.cc. This implementation
// follows the TSI model where most of the algorithmic code is moved to
// services accessed via service interfaces.
//
// For configuration parameters, see the "Flags" block in the module class
// definition below. Remove the leading "f" to get the parameter name.
//
// There is no flag for compression because this must always be invoked
// to apply zero suppression. Use ReplaceCompressService (prolog cmpreplace)
// to effectively skip the compression while retaining the supression.
//
// The interface names for the accessed services are listed in the "Services"
// block in the header below.
//
// Some useful module and service configurations may be found in
// detsimmodules_dune.fcl, e.g. cmpreplace to skip compression.

#include <vector>
#include <string>
#include <sstream>
#include <iostream>

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/raw.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "dune/ArtSupport/DuneToolManager.h"

#include "dune/DuneInterface/AdcTypes.h"
#include "dune/DuneInterface/AdcSimulator.h"
#include "dune/DuneInterface/AdcSuppressService.h"
#include "dune/DuneInterface/AdcCompressService.h"
#include "dune/DuneInterface/SimChannelExtractService.h"
#include "dune/DuneInterface/ChannelNoiseService.h"
#include "dune/DuneInterface/PedestalAdditionService.h"
#include "dune/DuneInterface/AdcDistortionService.h"

using std::ostringstream;
using std::cout;
using std::endl;
using std::string;

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

  std::string fSimChannelLabel; ///< Data product holding the ionization electrons

  // Flags.
  bool fNoiseOn;           ///< noise turned on or off for debugging; default is on
  bool fPedestalOn;        ///< switch for simulation of nonzero pedestals
  bool fDistortOn;         ///< switch for simulation of stuck bits
  bool fSuppressOn;        ///< switch for simulation of zero suppression
  bool fKeepEmptyChannels; ///< Write out empty channels iff true.

  // Tools.
  std::string fAdcSimulatorName;
  std::unique_ptr<AdcSimulator> m_pads;

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
  string myname = "SimWireDUNE::reconfigure: ";
  string myprefix = myname + "    ";
  fSimChannelLabel   = p.get<std::string>("SimChannelLabel");
  fNoiseOn           = p.get<bool>("NoiseOn");
  fPedestalOn        = p.get<bool>("PedestalOn");
  fDistortOn         = p.get<bool>("DistortOn");
  fSuppressOn        = p.get<bool>("SuppressOn");
  fKeepEmptyChannels = p.get<bool>("KeepEmptyChannels");
  fAdcSimulatorName = p.get<string>("AdcSimulator");
  DuneToolManager* pdtm = DuneToolManager::instance();
  m_pads = pdtm == nullptr ? nullptr : pdtm->getPrivate<AdcSimulator>(fAdcSimulatorName);
  ostringstream out;
  out << myname << "Tools:" << endl;
  out << "  AdcSimulator: " << bool(m_pads) << endl;
  out << myname << "Accessed services:" << endl;
  out << myname << "  SimChannel extraction service:" << endl;
  m_pscx->print(out, myprefix) << endl;
  if ( fNoiseOn ) {
    out << myname << "  Channel noise service:" << endl;;
    m_pcns->print(out, myprefix);
  } else {
    out << myname << "  Channel noise is off.";
  }
  out << endl;
  if ( fPedestalOn ) {
    out << myname << "  Pedestal addition service:" << endl;;
    m_ppad->print(out, myprefix);
  } else {
    out << myname << "  Pedestal addition is off.";
  }
  out << endl;
  if ( fDistortOn ) {
    out << myname << "  ADC distortion service:" << endl;;
    m_pdis->print(out, myprefix);
  } else {
    out << myname << "  ADC distortion is off.";
  }
  out << endl;
  if ( fSuppressOn ) {
    out << myname << "  ADC suppression service:" << endl;
    m_pzs->print(out, myprefix);
  } else {
    out << myname << "  ADC suppression is off.";
  }
  out << endl;
  out << myname << "  Compression service:" << endl;
  out << endl;
  m_pcmp->print(out, myprefix);
  out << myname << "  KeepEmptyChannels:" << fKeepEmptyChannels << endl;
  mf::LogInfo("SimWireDUNE::reconfigure") << out.str();

  return;
}

//**********************************************************************

void SimWireDUNE::beginJob() { }

//**********************************************************************

void SimWireDUNE::endJob() { }

//**********************************************************************

void SimWireDUNE::produce(art::Event& evt) {
  const string myname = "SimWireDUNE::produce: ";

  // Make a vector of const sim::SimChannel* that has same number
  // of entries as the number of channels in the detector
  // and set the entries for the channels that have signal on them
  // using the chanHandle
  std::vector<const sim::SimChannel*> chanHandle;
  std::vector<const sim::SimChannel*> simChannels(m_pgeo->Nchannels(), nullptr);
  evt.getView(fSimChannelLabel, chanHandle);
  for ( size_t c=0; c<chanHandle.size(); ++c ) {
    simChannels[chanHandle[c]->Channel()] = chanHandle[c];
  }

  // make an unique_ptr of sim::SimDigits that allows ownership of the produced
  // digits to be transferred to the art::Event after the put statement below
  std::unique_ptr<std::vector<raw::RawDigit>>  digcol(new std::vector<raw::RawDigit>);

  // Fetch the number of ticks to write out for each channel.
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  unsigned int nTickReadout  = detprop->ReadOutWindowSize();

  // Loop over channels.
  std::map<int,double>::iterator mapIter;
  for ( unsigned int chan = 0; chan<m_pgeo->Nchannels(); ++chan ) {

    // Get the SimChannel for this channel
    const sim::SimChannel* psc = simChannels[chan];

    // Create vector that holds the floating ADC count for each tick.
    std::vector<AdcSignal> fChargeWork;

    // Extract the floating ADC count from the SimChannel for each tick in the channel.
    m_pscx->extract(psc, fChargeWork);

    // Add noise to each tick in the channel.
    if ( fNoiseOn ) {
      m_pcns->addNoise(chan, fChargeWork);
    }

    // Option to display signals before adding pedestals.
    // E.g. logsig = chan==1000.
    bool logsig = false;
    if ( logsig ) {
      cout << myname << "Signals after noise:" << endl;
      for ( unsigned int itck=0; itck<fChargeWork.size(); ++itck ) {
        if(fChargeWork[itck] > 0 )
          cout << myname << " " << itck << ": chg=" << fChargeWork[itck] << endl;
      }
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
      AdcSignal adcin = fChargeWork[itck];
      short adc = 0;
      bool useOldAdc = false;
      // New ADC calculation (with tool).
      if ( m_pads ) {
        if(adcin > 0)
          //cout << " adcin " << adcin << endl;
        adc = m_pads->count(adcin, chan);
      } else {
        //cout << myname << "WARNING: AdcSimulator not found." << endl;
        useOldAdc = true;
      }
      // Old ADC calculation.
      if ( useOldAdc ) {
        short adc1 = 0;
        if ( adcin > 0 ) adc1 = (short) (adcin + 0.5);
        if ( adc1 > adcmax ) adc1 = adcmax;
        bool show = m_pads && adc1 != adc;

        static int ndump = 1000;
        if ( ndump && show ) {
          cout << myname << "  ADC: " << adc1 << " --> " << adc << " (" << adcin << ")" << endl;
          --ndump;
        }
        adc = adc1;
      }
      // Record the ADC value.
      adcvec[itck] = adc;
    }
    // Resize to the correct number of time samples, dropping extra samples.
    adcvec.resize(nTickReadout);

    // Add stuck bits.
    if ( fDistortOn ) {
      m_pdis->modify(chan, adcvec);
    }

    // Zero suppress.
    AdcCountSelection acs(adcvec, chan, pedval);
    if ( fSuppressOn ) {
      m_pzs->filter(acs);
    }
    int nkeep = 0;
    for ( bool kept : acs.filter ) if ( kept ) ++nkeep;

    // If flag is not set and channel is empty, skip it.
    if ( ! fKeepEmptyChannels && nkeep==0 ) continue;

    // Compress.
    raw::Compress_t comp = raw::kNone;
    m_pcmp->compress(adcvec, acs.filter, pedval, comp);

    // Create and store raw digit.
    raw::RawDigit rd(chan, nTickReadout, adcvec, comp);
    rd.SetPedestal(pedval, pedrms);
    digcol->push_back(rd);            // add this digit to the collection

  }  // end loop over channels

  evt.put(std::move(digcol));

}
