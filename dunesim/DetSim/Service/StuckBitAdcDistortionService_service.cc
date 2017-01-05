// StuckBitAdcDistortionService_service.cc

#include "dune/DetSim/Service/StuckBitAdcDistortionService.h"
#include "fhiclcpp/ParameterSet.h"
#include "nutools/RandomUtils/NuRandomService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "TFile.h"
#include "TString.h"
#include "TProfile.h"

using std::string;
using std::cout;
using std::ostream;
using std::endl;
using rndm::NuRandomService;
using CLHEP::HepJamesRandom;

//**********************************************************************

StuckBitAdcDistortionService::
StuckBitAdcDistortionService(const fhicl::ParameterSet& pset, art::ActivityRegistry&)
: m_RandomSeed(0), m_LogLevel(1),
  m_pran(nullptr) {
  const string myname = "StuckBitAdcDistortionService::ctor: ";
  // Fetch parameters.
  fStuckBitsProbabilitiesFname = pset.get<string>("StuckBitsProbabilitiesFname");
  fStuckBitsOverflowProbHistoName = pset.get<string>("StuckBitsOverflowProbHistoName");
  fStuckBitsUnderflowProbHistoName = pset.get<string>("StuckBitsUnderflowProbHistoName");
  bool haveSeed = pset.get_if_present<int>("RandomSeed", m_RandomSeed);
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  if ( m_RandomSeed == 0 ) haveSeed = false;
  // Create random number engine.
  if ( haveSeed ) {
    if ( m_LogLevel > 0 ) cout << myname << "WARNING: Using hardwired seed." << endl;
    m_pran = new HepJamesRandom(m_RandomSeed);
  } else {
    string rname = "StuckBitAdcDistortionService";
    if ( m_LogLevel > 0 ) cout << myname << "Using NuRandomService." << endl;
    art::ServiceHandle<NuRandomService> seedSvc;
    m_pran = new HepJamesRandom;
    if ( m_LogLevel > 0 ) cout << myname << "    Initial seed: " << m_pran->getSeed() << endl;
    seedSvc->registerEngine(NuRandomService::CLHEPengineSeeder(m_pran), rname);
  }
  if ( m_LogLevel > 0 ) cout << myname << "  Registered seed: " << m_pran->getSeed() << endl;
  // Fetch the probabilities.
  mf::LogInfo("SimWireDUNE") << " using ADC stuck code probabilities from .root file " ;
  std::string fname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(fStuckBitsProbabilitiesFname, fname);
  std::unique_ptr<TFile> fin(new TFile(fname.c_str(), "READ"));
  if ( !fin->IsOpen() )
    throw art::Exception(art::errors::NotFound)
          << "Could not find the ADC stuck code probabilities file " << fname;
  TString iOverflowHistoName = Form( "%s", fStuckBitsOverflowProbHistoName.c_str());
  TProfile *overflowtemp = (TProfile*) fin->Get( iOverflowHistoName );
  if ( !overflowtemp )
    throw art::Exception(art::errors::NotFound)
          << "Could not find the ADC code overflow probabilities histogram "
          << fStuckBitsOverflowProbHistoName;
  if ( overflowtemp->GetNbinsX() != 64 )
    throw art::Exception(art::errors::InvalidNumber)
          << "Overflow ADC stuck code probability histograms must have 64 bins.";
  TString iUnderflowHistoName = Form( "%s", fStuckBitsUnderflowProbHistoName.c_str());
  TProfile *underflowtemp = (TProfile*) fin->Get(iUnderflowHistoName);
  if ( !underflowtemp )
    throw art::Exception( art::errors::NotFound )
          << "Could not find the ADC code underflow probabilities histogram "
          << fStuckBitsUnderflowProbHistoName;
  if ( underflowtemp->GetNbinsX() != 64 )
    throw art::Exception(art::errors::InvalidNumber)
          << "Underflow ADC stuck code probability histograms must have 64 bins.";
  for ( unsigned int cellnumber=0; cellnumber < 64; ++cellnumber ) {
    fOverflowProbs[cellnumber] = overflowtemp->GetBinContent(cellnumber+1);
    fUnderflowProbs[cellnumber] = underflowtemp->GetBinContent(cellnumber+1);
  }
  fin->Close();
}
  
//**********************************************************************

StuckBitAdcDistortionService::~StuckBitAdcDistortionService() {
  const string myname = "StuckBitAdcDistortionService:dtor: ";
  if ( m_LogLevel > 0 ) {
    cout << myname << "Deleting random engine with seed " << m_pran->getSeed() << endl;
  }
  delete m_pran;
}

//**********************************************************************

int StuckBitAdcDistortionService::modify(Channel, AdcCountVector& adcvec) const {
  CLHEP::RandFlat stuck_flat(*m_pran);
  for ( size_t itck = 0; itck<adcvec.size(); ++itck ) {
    double rnd = stuck_flat.fire(0,1);
    const unsigned int zeromask = 0xffc0;
    const unsigned int onemask = 0x003f;
    unsigned int sixlsbs = adcvec[itck] & onemask;
    int probability_index = (int)sixlsbs;
    if ( rnd < fUnderflowProbs[probability_index] ) {
      adcvec[itck] = adcvec[itck] | onemask; // 6 LSBs are stuck at 3F
      adcvec[itck] -= 64; // correct 1st MSB value by subtracting 64
    } else if ( rnd > fUnderflowProbs[probability_index] &&
              rnd < fUnderflowProbs[probability_index] + fOverflowProbs[probability_index] ) {
      adcvec[itck] = adcvec[itck] & zeromask; // 6 LSBs are stuck at 0
      adcvec[itck] += 64; // correct 1st MSB value by adding 64
    }
  }
  return 0;
}

//**********************************************************************

ostream& StuckBitAdcDistortionService::print(ostream& out, string prefix) const {
  out << prefix << "StuckBitAdcDistortionService:" << endl;
  out << prefix << "       StuckBitsProbabilitiesFname: " << fStuckBitsProbabilitiesFname << endl;
  out << prefix << "   fStuckBitsOverflowProbHistoName: " << fStuckBitsOverflowProbHistoName << endl;
  out << prefix << "  fStuckBitsUnderflowProbHistoName: " << fStuckBitsUnderflowProbHistoName << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(StuckBitAdcDistortionService, AdcDistortionService)

//**********************************************************************
