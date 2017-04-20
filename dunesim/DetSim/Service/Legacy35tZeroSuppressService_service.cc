// Legacy35tZeroSuppressService.cxx

#include "dune/DetSim/Service/Legacy35tZeroSuppressService.h"
#include <cmath>
#include "fhiclcpp/ParameterSet.h"
#include "dune/DetSim/Utility/AdcCodeHelper.h"

using std::string;
using std::ostream;
using std::endl;

namespace {

string sbool(bool arg) {
  if ( arg ) return "true";
  return "false";
}

}  // end unnamed namespace.

typedef unsigned int Index;

//**********************************************************************

Legacy35tZeroSuppressService::
Legacy35tZeroSuppressService(const fhicl::ParameterSet& pset, art::ActivityRegistry&) {
  m_AdcThreshold       = pset.get<float>("AdcThreshold");
  m_TickRange          = pset.get<unsigned int>("TickRange");
  m_MinTickGap         = pset.get<unsigned int>("MinTickGap");
  m_SuppressStickyBits = pset.get<bool>("SuppressStickyBits");
}
  
//**********************************************************************

Legacy35tZeroSuppressService::
Legacy35tZeroSuppressService(float aAdcThreshold,
                             unsigned int aTickRange,
                             unsigned int aMinTickGap,
                             bool aSuppressStickyBits)
: m_AdcThreshold(aAdcThreshold),
  m_TickRange(aTickRange),
  m_MinTickGap(aMinTickGap),
  m_SuppressStickyBits(aSuppressStickyBits) { }

//**********************************************************************

int Legacy35tZeroSuppressService::
filter(const AdcCountVector& sigs, Channel, AdcPedestal ped, AdcFilterVector& keep) const {
  const unsigned int nsig = sigs.size();
  keep.clear();
  keep.resize(nsig, false);
  if ( nsig == 0 ) return 0;
  AdcCodeHelper ach(64);
  for ( unsigned int isig=0; isig<nsig; ++isig ) {
    AdcCount rawsig = sigs[isig];
    AdcSignal pedsig = ach.subtract(ped, rawsig);
    if ( m_SuppressStickyBits ) {
      if ( ach.hasStickyBits(rawsig) && ach.isSmall(pedsig) ) {
        std::cout << "Legacy35tZeroSuppressService::filter: Suppressing sticky ADC value " << pedsig
             << " (raw=" << rawsig << ", ped=" << ped << ")" << endl;
        pedsig = 0.0;
      }
    }
    if ( fabs(pedsig) > m_AdcThreshold ) {
      Index isig1 = 0;
      if ( m_TickRange < isig ) isig1 = isig - m_TickRange;
      Index isig2 = isig + m_TickRange;
      if ( isig2 > nsig ) isig2 = nsig;
      for ( Index jsig=isig1; jsig<=isig2; ++jsig ) keep[jsig] = true;
    }
  }
  // Make gaps between unsupressed blocks < MinTickGap.
  bool inGap = false;
  unsigned int jsig1 = 0;   // Start of gap.
  unsigned int jsig2 = 0;   // End of gap.
  for ( unsigned int isig=0; isig<nsig; ++isig ) {
    // End of a gap;
    if ( inGap && (isig+1==nsig || keep[isig]) ) {
      inGap = false;
      if ( (jsig2 - jsig1 + 1) < m_MinTickGap ) {
        for ( unsigned int jsig=jsig1; jsig<=jsig2; ++jsig ) {
          keep[jsig] = true;
        }
      }
    // In a gap.
    } else if ( inGap ) {
      jsig2 = isig;
    // Start of gap.
    } else if ( !inGap && !keep[isig] ) {
      inGap = true;
      jsig1 = isig;
      jsig2 = jsig1;
    }
  }
  return 0;
}

//**********************************************************************

ostream& Legacy35tZeroSuppressService::print(ostream& out, string prefix) const {
  out << prefix << "Legacy35tZeroSuppressService:"               << endl;
  out << prefix << "       AdcThreshold: " << m_AdcThreshold       << endl;
  out << prefix << "         MinTickGap: " << m_MinTickGap         << endl;
  out << prefix << "          TickRange: " << m_TickRange          << endl;
  out << prefix << " SuppressStickyBits: " << sbool(m_SuppressStickyBits);
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(Legacy35tZeroSuppressService, AdcSuppressService)

//**********************************************************************
