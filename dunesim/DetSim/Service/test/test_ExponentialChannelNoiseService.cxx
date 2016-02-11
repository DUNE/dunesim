// test_ExponentialChannelNoiseService.cxx

// David Adams
// February 2015
//
// Test ExponentialChannelNoiseService.

#include "../ExponentialChannelNoiseService.h"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/TimeService.h"
#include "Utilities/LArFFT.h"
#include "RawData/raw.h"

using std::string;
using std::cout;
using std::endl;
using std::istringstream;
using std::ofstream;

typedef std::vector<short> AdcVector;

#undef NDEBUG
#include <cassert>

int test_ExponentialChannelNoiseService() {
  const string myname = "test_ExponentialChannelNoiseService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  cout << myname << "Fetch art service helper." << endl;
  ArtServiceHelper& ash = ArtServiceHelper::instance();
  ash.print();

  cout << line << endl;
  cout << myname << "Add TFileService." << endl;
  string scfg = "fileName: \"test.root\" service_type: \"TFileService\"";
  assert( ash.addService("TFileService", scfg) == 0 );

  cout << myname << line << endl;
  cout << myname << "Add RandomNumberGenerator service." << endl;
  scfg = "";
  cout << myname << "Configuration: " << scfg << endl;
  assert( ash.addService("RandomNumberGenerator", scfg) == 0 );

  cout << myname << line << endl;
  cout << myname << "Add DetectorProperties service." << endl;
  scfg = "ElectronsToADC: 6.8906513e-3 NumberTimeSamples: 3200 ReadOutWindowSize: 3200 TimeOffsetU: 0 TimeOffsetV: 0 TimeOffsetZ: 0";
  cout << myname << "Configuration: " << scfg << endl;
  assert( ash.addService("DetectorProperties", scfg) == 0 );

  cout << myname << line << endl;
  cout << myname << "Add the TimeService service." << endl;
  scfg = "ClockSpeedExternal: 3.125e1 ClockSpeedOptical: 128 ClockSpeedTPC: 2 ClockSpeedTrigger: 16 DefaultBeamTime: 0 DefaultTrigTime: 0 FramePeriod: 1600 G4RefTime: 0 InheritClockConfig: false TrigModuleName: \"\" TriggerOffsetTPC: 0 service_type: \"TimeService\"";
  cout << myname << "Configuration: " << scfg << endl;
  assert( ash.addService("TimeService", scfg) == 0 );

  cout << myname << line << endl;
  cout << myname << "Add the signal shaping service." << endl;
  bool isFile = false;
  // This doesn't work because defn is in a PROLOG.
  scfg = "services_dune.fcl";
  // This doesn't work because text must be preprocessed.
  scfg = "#include \"services_dune.fcl\"\n";
  scfg += "services: {\n";
  scfg += "user: @local::dune35t_services\n";
  scfg += "}";
  {
    ofstream ofile("test_SignalShapingServiceDUNE35t_services.fcl");
    ofile << scfg << endl;
  }
  // This works if the file contains something like the above text.
  scfg = "test_SignalShapingServiceDUNE35t_services.fcl";
  isFile = true;
  cout << myname << "Configuration: " << scfg << endl;
  assert( ash.addService("SignalShapingServiceDUNE35t", scfg, isFile) == 0 );

  cout << myname << line << endl;
  cout << myname << "Add the FFT service." << endl;
  assert( ash.addService("LArFFT", scfg, isFile) == 0 );

  cout << myname << line << endl;
  cout << myname << "Load services." << endl;
  assert( ash.loadServices() == 1 );
  ash.print();

  cout << line << endl;
  cout << "Create noise service." << endl;
  fhicl::ParameterSet pset;
  pset.put("NoiseFactZ",     0.05);
  pset.put("NoiseFactU",     0.05);
  pset.put("NoiseFactV",     0.05);
  pset.put("NoiseWidthZ", 2000.0);
  pset.put("NoiseWidthU", 2000.0);
  pset.put("NoiseWidthV", 2000.0);
  pset.put("LowCutoffZ",     7.5);
  pset.put("LowCutoffU",     7.5);
  pset.put("LowCutoffV",     7.5);
  pset.put("NoiseArrayPoints", 1000);
  pset.put("OldNoiseIndex", true);
  ExponentialChannelNoiseService noise(pset);

  cout << myname << "Done." << endl;
  return 0;
}

int main(int argc, char* argv[]) {
  int logLevel = 1;
  if ( argc > 1 ) {
    istringstream ssarg(argv[1]);
    ssarg >> logLevel;
  }
  test_ExponentialChannelNoiseService();
  return 0;
}
