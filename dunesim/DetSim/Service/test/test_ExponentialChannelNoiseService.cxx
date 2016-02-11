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
#include <iomanip>
#include "dune/ArtSupport/ArtServiceHelper.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/TimeService.h"
#include "Utilities/LArFFT.h"
#include "Geometry/Geometry.h"
#include "RawData/raw.h"

using std::string;
using std::cout;
using std::endl;
using std::istringstream;
using std::ostringstream;
using std::ofstream;
using std::setw;
using art::ServiceHandle;
using util::LArFFT;

typedef std::vector<short> AdcVector;

#undef NDEBUG
#include <cassert>

void load_services(int fftsize =0) {
  const string myname = "load_services: ";
  const string line = "-----------------------------";

  static bool loaded = false;
  if ( loaded ) return;
  loaded = true;

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

#if 0
  cout << myname << line << endl;
  cout << myname << "Add the signal shaping service." << endl;
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
#endif

  // Work around defect in larsoft: https://cdcvs.fnal.gov/redmine/issues/10618.
  TH1::AddDirectory(kFALSE);

  cout << myname << line << endl;
  ostringstream sscfg;
  sscfg << "FFTOption: \"\" FFTSize: " << fftsize << " FitBins: 5";
  scfg = sscfg.str();
  cout << myname << "Add the FFT service." << endl;
  //assert( ash.addService("LArFFT", scfg, isFile) == 0 );
  assert( ash.addService("LArFFT", scfg) == 0 );

  cout << myname << line << endl;
  string gname = "dune35t4apa_v5";
  cout << myname << "Add Geometry service." << endl;
  scfg = "DisableWiresInG4: true GDML: \"dune35t4apa_v5.gdml\" Name: \"" + gname +
         "\" ROOT: \"" + gname + "\" SortingParameters: { DetectorVersion: \"" + gname +
         "\" } SurfaceY: 0";
  cout << myname << "Configuration: " << scfg << endl;
  assert( ash.addService("Geometry", scfg) == 0 );

  cout << myname << line << endl;
  cout << myname << "Add the DUNE geometry helper service (required to load DUNE geometry)." << endl;
  scfg = "service_provider: \"DUNEGeometryHelper\"";
  cout << myname << "Configuration: " << scfg << endl;
  assert( ash.addService("ExptGeoHelperInterface", scfg) == 0 );

  cout << myname << line << endl;
  cout << myname << "Load services." << endl;
  assert( ash.loadServices() == 1 );
  ash.print();

  return;
}

int test_ExponentialChannelNoiseService(unsigned int ntick) {
  const string myname = "test_ExponentialChannelNoiseService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  const string line = "-----------------------------";

  load_services();

  cout << myname << line << endl;
  cout << myname << "Create noise service." << endl;
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
  ExponentialChannelNoiseService noisvc(pset);

  ServiceHandle<LArFFT> hfftsvc;

  AdcSignal ped = 1000;
  cout << myname << line << endl;
  cout << myname << "Create signal with ped=" << ped << " and " << ntick << " ticks." << endl;
  cout << myname << "FFT size: " << hfftsvc->FFTSize() << endl;
  AdcSignalVector sigs(ntick, ped);
  cout << myname << "Add noise." << endl;
  assert( noisvc.addNoise(0, sigs) == 0 );
  for ( unsigned int isig=0; isig<20; ++isig ) {
  }
  double sum = 0.0;
  double sumsq = 0.0;
  for ( unsigned int isig=0; isig<ntick; ++isig ) {
    if ( isig < 20 || ntick-isig < 20 ) {
      cout << myname << setw(5) << isig << ": " << setw(8) << sigs[isig] << endl;
    }
    if ( isig == 20 && ntick-isig > 20 ) {
      cout << myname << setw(5) << "..." << endl;
    }
    float sig = sigs[isig];
    sum += sig;
    sumsq += sig*sig;
  }
  float mean = sum/ntick;
  float rms = sqrt(sumsq/ntick - mean*mean);
  cout << myname << "  Mean: " << mean << endl;
  cout << myname << "   RMS: " << rms << endl;

  cout << myname << "Done." << endl;
  return 0;
}

int main(int argc, char* argv[]) {
  int fftsize = 0;
  if ( argc > 1 ) {
    istringstream ssarg(argv[1]);
    ssarg >> fftsize;
  }
  load_services(fftsize);
  test_ExponentialChannelNoiseService(100);
  test_ExponentialChannelNoiseService(1000);
  test_ExponentialChannelNoiseService(10000);
  return 0;
}
