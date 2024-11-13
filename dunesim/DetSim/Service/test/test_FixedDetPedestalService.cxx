// test_FixedDetPedestalService.cxx

// David Adams
// March 2015
//
// Test FixedDetPedestalService.

#include "../FixedDetPedestalService.h"
#include "dunecore/ArtSupport/ArtServiceHelper.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using std::vector;
using std::setw;
using art::ServiceHandle;
using geo::Geometry;
using lariov::DetPedestalService;
using lariov::DetPedestalProvider;

#undef NDEBUG
#include <cassert>

int test_FixedDetPedestalService() {
  const string myname = "test_FixedDetPedestalService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Create top-level FCL." << endl;
  std::ostringstream oss;
  oss << "#include \"pedestals_dune.fcl\"" << endl;
  oss << "#include \"geometry_dune.fcl\"" << endl;
  oss << "services.DetPedestalService: @local::dune_fixedpeds" << endl;
  oss << "services.Geometry: @local::dune35t_geo" << endl;
  oss << "services.WireReadout: @local::dune35t_wire_readout" << endl;
  ArtServiceHelper::load_services(oss.str());

  cout << myname << line << endl;
  cout << myname << "Fetch geometry service." << endl;
  ServiceHandle<geo::Geometry> hgeo;
  auto const& wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();
  cout << myname << "Detector: " << hgeo->DetectorName() << endl;

  cout << myname << line << endl;
  cout << myname << "Fetch pedestal service." << endl;
  ServiceHandle<DetPedestalService> hdps;
  const DetPedestalProvider* pdpp = hdps->provider();
  assert( pdpp != nullptr );

  cout << myname << line << endl;
  cout << myname << "Checking pedestals for selected channels." << endl;
  vector<int> chans = {0, 1, 200, 300, 400, 500, 600, 700, 800, 900};
  int hchan = 10;
  int hview = 10;
  int hmean = 12;
  cout << setw(hchan) << "Channel"
       << setw(hview) << "View"
       << setw(hmean) << "Mean"
       << setw(hmean) << "RMS"
       << endl;
  for ( int chan : chans ) {
    cout << setw(hchan) << chan
         << setw(hview) << wireReadout.View(chan)
         << setw(hmean) << pdpp->PedMean(chan)
         << setw(hmean) << pdpp->PedRms(chan)
         << endl;
    float meanExp = 1800;
    if ( wireReadout.View(chan) == geo::kZ ) {
      meanExp = 500.0;
    }
    assert( pdpp->PedMean(chan) == meanExp );
  }

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

int main(int argc, char* argv[]) {
  return test_FixedDetPedestalService();
}
