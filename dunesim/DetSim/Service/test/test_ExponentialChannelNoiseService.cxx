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
#include "lardata/DetectorInfo/DetectorProperties.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/LArFFT.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/RawData/raw.h"

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

int test_ExponentialChannelNoiseService(unsigned int ntick, unsigned int maxchan =100) {
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
  pset.put("NoiseNormZ",     3.16);
  pset.put("NoiseNormU",     3.16);
  pset.put("NoiseNormV",     3.16);
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
  ServiceHandle<geo::Geometry> hgeosvc;

  AdcSignal ped = 1000;
  cout << myname << line << endl;
  cout << myname << "Create signal with ped=" << ped << " and " << ntick << " ticks." << endl;
  cout << myname << "FFT size: " << hfftsvc->FFTSize() << endl;
  cout << myname << "Detector: " << hgeosvc->DetectorName() << endl;
  AdcSignalVector sigsin(ntick, ped);
  for ( unsigned int isig=0; isig<20; ++isig ) {
  }
  double sum = 0.0;
  double sumsq = 0.0;
  double sumdsq = 0.0;
  double sumcol = 0.0;
  double sumsqcol = 0.0;
  double sumdsqcol = 0.0;
  double sumind = 0.0;
  double sumsqind = 0.0;
  double sumdsqind = 0.0;
  unsigned int nchan = 0;
  unsigned int nchancol = 0;
  unsigned int nchanind = 0;
  unsigned int count = 0;
  unsigned int countcol = 0;
  unsigned int countind = 0;
  // Loop over channels.
  cout << myname << "Looping over channels." << endl;
  for ( unsigned int icha=0; icha<1000000; ++icha ) {
    geo::SigType_t sigtyp = hgeosvc->SignalType(icha);
    if ( sigtyp == geo::kMysteryType ) break;
    if ( icha >= maxchan ) break;
    ++nchan;
    bool firstcol = false;
    bool firstind = false;
    string labtyp = "unknown";
    if      ( sigtyp == geo::kCollection ) {
      firstcol = nchancol == 0;
      labtyp = "collection";
      ++nchancol;
    } else if ( sigtyp == geo::kInduction ) {
      firstind = nchanind == 0;
      labtyp = "induction";
      ++nchanind;
    } else abort();
    AdcSignalVector sigs = sigsin;
    assert( noisvc.addNoise(icha, sigs) == 0 );
    if ( firstcol || firstind ) {
      cout << myname << "First " << labtyp << " channel: " << icha << endl;
    }
    for ( unsigned int isig=0; isig<ntick; ++isig ) {
      if ( firstcol || firstind ) {
        if ( isig < 20 || ntick-isig < 20 ) {
          cout << myname << setw(5) << isig << ": " << setw(8) << sigs[isig] << endl;
        }
        if ( isig == 20 && ntick-isig > 20 ) {
          cout << myname << setw(5) << "..." << endl;
        }
      }
      float sig = sigs[isig];
      double dif = sig - ped;
      sum += sig;
      sumsq += sig*sig;
      sumdsq += dif*dif;
      ++count;
      if ( sigtyp == geo::kCollection ) {
        sumcol += sig;
        sumsqcol += sig*sig;
        sumdsqcol += dif*dif;
        ++countcol;
      } else if ( sigtyp == geo::kInduction ) {
        sumind += sig;
        sumsqind += sig*sig;
        sumdsqind += dif*dif;
        ++countind;
      }
    }
  }
  cout << myname << "           # channels: " << nchan << endl;
  cout << myname << "# collection channels: " << nchancol << endl;
  cout << myname << " # induction channels: " << nchanind << endl;
  double nticktot = double(nchan)*ntick;
  double nticktotcol = double(nchancol)*ntick;
  double nticktotind = double(nchanind)*ntick;
  cout << myname << "            # tick total: " << nticktot << endl;
  cout << myname << "                   Count: " << count << endl;
  cout << myname << " Collection # tick total: " << nticktotcol << endl;
  cout << myname << " Collection        Count: " << countcol << endl;
  cout << myname << "  Induction # tick total: " << nticktotind << endl;
  cout << myname << "  Induction        Count: " << countind << endl;
  assert( nchan > 0 );
  float mean = sum/nticktot;
  float rms = sqrt(sumsq/nticktot - mean*mean);
  float rms2 = sqrt(sumdsq/nticktot);
  float meancol = sumcol/nticktotcol;
  float rmscol = sqrt(sumsqcol/nticktotcol - meancol*meancol);
  float rms2col = sqrt(sumdsqcol/nticktotcol);
  float meanind = sumind/nticktotind;
  float rmsind = sqrt(sumsqind/nticktotind - meanind*meanind);
  float rms2ind = sqrt(sumdsqind/nticktotind);
  cout << myname << "            Mean: " << mean << endl;
  cout << myname << "             RMS: " << rms << endl;
  cout << myname << "            RMS2: " << rms2 << endl;
  cout << myname << " Collection Mean: " << meancol << endl;
  cout << myname << " Collection  RMS: " << rmscol << endl;
  cout << myname << " Collection RMS2: " << rms2col << endl;
  cout << myname << "  Induction Mean: " << meanind << endl;
  cout << myname << "  Induction  RMS: " << rmsind << endl;
  cout << myname << "  Induction RMS2: " << rms2ind << endl;

  cout << myname << "Done." << endl;
  return 0;
}

int main(int argc, char* argv[]) {
  int fftsize = 0;
  int maxchan = 2000;
  int ntick = 1000;
  if ( argc > 1 ) {
    istringstream ssarg(argv[1]);
    ssarg >> fftsize;
  }
  if ( argc > 2 ) {
    istringstream ssarg(argv[2]);
    ssarg >> maxchan;
  }
  load_services(fftsize);
  test_ExponentialChannelNoiseService(ntick, maxchan);
  // We must close TFileService to write out histograms.
  cout << "Close services." << endl;
  ArtServiceHelper::close();
  return 0;
}
