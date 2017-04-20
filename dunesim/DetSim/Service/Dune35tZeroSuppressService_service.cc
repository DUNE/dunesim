// Dune35tZeroSuppressService.cxx

#include "dune/DetSim/Service/Dune35tZeroSuppressService.h"
#include <iomanip>
#include "dune/DetSim/Utility/SignalTypeConverter.h"
#include "fhiclcpp/ParameterSet.h"

using std::string;
using std::ostream;
using std::cout;
using std::endl;
using std::setw;

typedef Dune35tZeroSuppressService::Index Index;

//**********************************************************************

// Signal state.
//    OUT - outside the signal (may be in the tail of an earlier signal)
//   LIVE - Live region of the signal (above TL or after)
//   DEAD - Dead region of the signal (below TD)
//    END - End of dead region (followed by tail)
enum SigState { OUT, LIVE, DEAD, END };

//**********************************************************************

namespace {

// Calculate sum for nsig bins before isig.
// If skipStuck, the stuck bins do not contribute to sum or count.
// Bins with sig <= thresh, do not contribute to sum.
struct RunningSum {
  RunningSum(const AdcCountVector& sigs, AdcPedestal ped, Index isig, Index ns, AdcCount thresh, bool skipStuck);
  AdcCount sigsum;
  Index count;
  SignalTypeConverter sigcon;
};

RunningSum::
RunningSum(const AdcCountVector& sigs, AdcPedestal ped, Index isig, Index nsig, AdcCount thresh, bool skipStuck) {
  sigsum = 0;
  count = 0;
  Index jsig1 = 0;
  if ( isig > nsig ) jsig1 = isig - nsig;
  Index jsig2 = isig;
  if ( jsig2 > sigs.size() ) jsig2 = sigs.size();
  for ( Index jsig=jsig1; jsig<jsig2; ++jsig ) {
    AdcCount rawsig = sigs[jsig];
    AdcCount sig = sigcon.convert<AdcCount>(rawsig - ped);
    if ( skipStuck ) {
      Index lsb = rawsig & 0x3f;
      if ( lsb == 0 ) continue;
      if ( lsb == 0x3f ) continue;
    }
    ++count;
    if ( abs(sig) > thresh ) {
      sigsum +=  sig;
    }
  }
}

// Convert state to string.
string sstate(SigState state) {
  if ( state == OUT ) return "OUT";
  if ( state == LIVE ) return "LIVE";
  if ( state == DEAD ) return "DEAD";
  if ( state == END ) return "END";
  return "NONE";
}

}  // End unnamed namespace.

//**********************************************************************

Dune35tZeroSuppressService::
Dune35tZeroSuppressService(AdcCount ts, AdcCount tl, AdcCount td,
                       Index ns, Index nl, Index nd, Index nt)
: m_ts(ts), m_tl(tl), m_td(td), m_ns(ns), m_nl(nl), m_nd(nd), m_nt(nt),
  m_LogLevel(1) { }

//**********************************************************************

Dune35tZeroSuppressService::
Dune35tZeroSuppressService(const fhicl::ParameterSet& pset, art::ActivityRegistry&)
: m_LogLevel(1) {
  const string myname = "Dune35tZeroSuppressService::ctor: ";
  m_ts = pset.get<AdcCount>("TS");
  m_tl = pset.get<AdcCount>("TL");
  m_td = pset.get<AdcCount>("TD");
  m_ns = pset.get<Index>("NS");
  m_nl = pset.get<Index>("NL");
  m_nd = pset.get<Index>("ND");
  m_nt = pset.get<Index>("NT");
  pset.get_if_present<int>("LogLevel", m_LogLevel);
  if ( m_LogLevel >= 1 ) print(cout, myname);
}
  
//**********************************************************************

int Dune35tZeroSuppressService::
filter(const AdcCountVector& sigs, Channel chan, AdcPedestal ped, AdcFilterVector& keep) const {
  const string myname = "ZeroSuppress35t::filter: ";
  if ( m_LogLevel >= 2 ) cout << myname << "Filtering channel " << chan << " with tick count " << sigs.size() << endl;
  bool m_skipStuck = false;
  AdcCount m_ts = 0; 
  unsigned int nsig = sigs.size();
  keep.clear();
  keep.resize(nsig, false);
  // Loop over signals.
  SigState state = OUT;
  unsigned int nlow = 0;
  for ( unsigned int isig=0; isig<nsig; ++isig ) {
    AdcCount sig = sigs[isig];
    // Evaluate a running signal sum of the preceding m_nl signals.
    RunningSum rs(sigs, ped, isig, m_ns, m_ts, m_skipStuck);
    AdcCount asigsum = std::abs(rs.sigsum);
    if ( m_LogLevel >= 3 ) cout << myname << setw(6) << isig << setw(6) << sig << setw(5) << sstate(state) << endl;
    // Last tick is outside a signal.
    if ( state == OUT || state == END ) {
      // If this tick is above TL, we are in the live region of a signal.
      // Keep the NL preceding signals.
      AdcCount sumthresh = m_tl*rs.count;
      bool keepit = asigsum > sumthresh;
      if ( m_LogLevel >= 3 ) cout << myname << " RS sum/thresh=" << setw(3) << rs.sigsum << "/" << setw(3) << sumthresh << endl;
      if ( keepit ) {
        if ( m_LogLevel == 2 ) cout << myname << setw(6) << isig << ": RS sum/thresh=" << setw(3) << rs.sigsum
                                    << "/" << setw(3) << sumthresh << endl;
        state = LIVE;
        unsigned int jsig1 = 0;
        if ( isig > m_nl ) jsig1 = isig - m_nl;
        unsigned int jsig2 = isig; 
        for ( unsigned int jsig=jsig1; jsig<jsig2; ++jsig) {
          keep[jsig] = true;
        }
      } else {
        state = OUT;
      }
    } else {
      // Last tick is is in the live region of a signal.
      AdcCount sumthresh = m_td*rs.count;
      if ( m_LogLevel >= 3 ) cout << myname << " RS sum/thresh=" << setw(3) << rs.sigsum << "/" << setw(3) << sumthresh << endl;
        if ( state == LIVE ) {
        // If this tick is below TD, we are in the dead region of a signal.
        if ( asigsum <= sumthresh ) {
          state = DEAD; 
          nlow = 1;
        }
      // Last tick is is in the dead region of a signal.
      } else if ( state == DEAD ) {
        // If signal is above TD, we are back in the live region.
        if ( asigsum > sumthresh ) {
          state = LIVE;
          nlow = 0;
        // If this is the ND'th consecutive signal in the dead region, we
        // have reached the end of the signal.
        // Keep this signal and a tail.
        } else if ( ++nlow >= m_nd ) {
          state = END;
          nlow = 0;
          // Protect the tail.
          unsigned int jsig1 = isig + 1;
          unsigned int jsig2 = jsig1 + m_nt;
          if ( jsig2 > nsig ) jsig2 = nsig;
          for ( unsigned int jsig=jsig1; jsig<jsig2; ++jsig) {
            keep[jsig] = true;
          }
        }
      } else {
        assert(false);
      }
    }
    if ( state != OUT ) keep[isig] = true;
  }
  return 0;
}

//**********************************************************************

ostream& Dune35tZeroSuppressService::print(ostream& out, string prefix) const {
  out << prefix << "    TS = " << m_ts << endl;
  out << prefix << "    TL = " << m_tl << endl;
  out << prefix << "    TD = " << m_td << endl;
  out << prefix << "    NS = " << m_ns << endl;
  out << prefix << "    NL = " << m_nl << endl;
  out << prefix << "    ND = " << m_nd << endl;
  out << prefix << "    NT = " << m_nt << endl;
  out << prefix << "LogLevel: " << m_LogLevel << endl;
  return out;
}

//**********************************************************************

void Dune35tZeroSuppressService::setDebug(int dbg) {
  m_LogLevel = dbg;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(Dune35tZeroSuppressService, AdcSuppressService)

//**********************************************************************
