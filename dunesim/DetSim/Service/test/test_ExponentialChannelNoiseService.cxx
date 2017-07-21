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
#include "lardata/Utilities/LArFFT.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/raw.h"

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

// service configuration strings (definitions are at the bottom)
extern const std::string LArPropertiesServiceConfigurationString;
extern const std::string DetectorPropertiesServiceConfigurationString;

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
  string gname = "dune35t4apa_v6";
  cout << myname << "Add Geometry service." << endl;
  scfg = "DisableWiresInG4: true GDML: \"dune35t4apa_v6.gdml\" Name: \"" + gname +
         "\" ROOT: \"" + gname + "\""
         " SortingParameters: { DetectorVersion: \"" + gname + "\""
         " ChannelsPerOpDet: 12" +
         "} SurfaceY: 0";
  cout << myname << "Configuration: " << scfg << endl;
  assert( ash.addService("Geometry", scfg) == 0 );

  cout << myname << line << endl;
  cout << myname << "Add the DUNE geometry helper service (required to load DUNE geometry)." << endl;
  scfg = "service_provider: \"DUNEGeometryHelper\"";
  cout << myname << "Configuration: " << scfg << endl;
  assert( ash.addService("ExptGeoHelperInterface", scfg) == 0 );

  cout << myname << line << endl;
  cout << myname << "Add LArPropertiesService service." << endl;
  scfg = LArPropertiesServiceConfigurationString;
  cout << myname << "Configuration: " << scfg << endl;
  assert( ash.addService("LArPropertiesService", scfg) == 0 );

  cout << myname << line << endl;
  cout << myname << "Add DetectorPropertiesService service." << endl;
  scfg = DetectorPropertiesServiceConfigurationString;
  cout << myname << "Configuration: " << scfg << endl;
  assert( ash.addService("DetectorPropertiesService", scfg) == 0 );

  cout << myname << line << endl;
  cout << myname << "Add the DetectorClocksService service." << endl;
  scfg = "ClockSpeedExternal: 3.125e1 ClockSpeedOptical: 128 ClockSpeedTPC: 2 ClockSpeedTrigger: 16 DefaultBeamTime: 0 DefaultTrigTime: 0 FramePeriod: 1600 G4RefTime: 0 InheritClockConfig: false TrigModuleName: \"\" TriggerOffsetTPC: 0 service_type: \"DetectorClocksService\" service_provider: \"DetectorClocksServiceStandard\" ";
  cout << myname << "Configuration: " << scfg << endl;
  assert( ash.addService("DetectorClocksService", scfg) == 0 );

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
  pset.put("WhiteNoiseZ",    0.0);
  pset.put("WhiteNoiseU",    0.0);
  pset.put("WhiteNoiseV",    0.0);
  pset.put("NoiseArrayPoints", 1000);
  pset.put("OldNoiseIndex", true);
  pset.put("RandomSeed", 54321);
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


//******************************************************************************
//*** Big configuration chunks
//***

// this is a copy of the "standard" LArProperties configuration
extern const std::string LArPropertiesServiceConfigurationString { R"cfg(
 service_provider: "LArPropertiesServiceStandard"

 # For following parameters, see http://pdg.lbl.gov/AtomicNuclearProperties/
 RadiationLength:  19.55   # g/cm^2
 AtomicNumber:     18      # Ar atomic number.
 AtomicMass:       39.948  # Ar atomic mass (g/mol).
 ExcitationEnergy: 188.0   # Ar mean excitation energy (eV).

# realistic Argon 39 decays
# Argon39DecayRate: 0.00141 # decays per cm^3 per second.  Assumes 1.01 Bq/kg and a density of 1.396 g/cc
# switch them off for faster simulation
 Argon39DecayRate: 0.0

 # Optical properties	
 # Fast and slow scintillation emission spectra, from [J Chem Phys vol 91 (1989) 1469]
 FastScintEnergies:    [ 6.0,  6.7,  7.1,  7.4,  7.7, 7.9,  8.1,  8.4,  8.5,  8.6,  8.8,  9.0,  9.1,  9.4,  9.8,  10.4,  10.7]
 SlowScintEnergies:    [ 6.0,  6.7,  7.1,  7.4,  7.7, 7.9,  8.1,  8.4,  8.5,  8.6,  8.8,  9.0,  9.1,  9.4,  9.8,  10.4,  10.7]
 FastScintSpectrum:    [ 0.0,  0.04, 0.12, 0.27, 0.44, 0.62, 0.80, 0.91, 0.92, 0.85, 0.70, 0.50, 0.31, 0.13, 0.04,  0.01, 0.0]
 SlowScintSpectrum:    [ 0.0,  0.04, 0.12, 0.27, 0.44, 0.62, 0.80, 0.91, 0.92, 0.85, 0.70, 0.50, 0.31, 0.13, 0.04,  0.01, 0.0]
 ScintResolutionScale: 1.     # resolution factor used by G4 scintillation
 ScintFastTimeConst:   6.     # fast scintillation time constant (ns)
 ScintSlowTimeConst:   1590.  # slow scintillation time constant (ns)
 ScintBirksConstant:   0.069  # birks constant for LAr (1/MeV cm)
 ScintYield:           24000. # total scintillation yield (ph/Mev)         
 ScintPreScale:        0.03   # Scale factor to reduce number of photons simulated
                              # Later QE should be corrected for this scale
 ScintYieldRatio:      0.3    # fast / slow scint ratio (needs revisitting)
 ScintByParticleType:  false  # whether to use different yields and
                              # quenching per particle in fast op sim
 EnableCerenkovLight: true    # whether to switch on cerenkov light (slow)




 # Scintillation yields and fast/slow ratios per particle type 
 MuonScintYield:          24000
 MuonScintYieldRatio:     0.23
 PionScintYield:          24000
 PionScintYieldRatio:     0.23 
 ElectronScintYield:      20000
 ElectronScintYieldRatio: 0.27
 KaonScintYield:          24000
 KaonScintYieldRatio:     0.23
 ProtonScintYield:        19200
 ProtonScintYieldRatio:   0.29
 AlphaScintYield:         16800
 AlphaScintYieldRatio:    0.56


 # Refractive index as a function of energy (eV) from arXiv:1502.04213v1
 RIndexEnergies: [ 1.900,  2.934,  3.592,  5.566,  6.694,  7.540,  8.574,  9.044,  9.232,  9.420,  9.514,  9.608,  9.702,  9.796,  9.890,  9.984,  10.08,  10.27,  10.45,  10.74,  10.92 ]
 RIndexSpectrum: [ 1.232,  1.236,  1.240,  1.261,  1.282,  1.306,  1.353,  1.387,  1.404,  1.423,  1.434,  1.446,  1.459,  1.473,  1.488,  1.505,  1.524,  1.569,  1.627,  1.751,  1.879 ]

 # absorption length as function of energy
 AbsLengthEnergies: [ 4,     5,     6,     7,     8,     9,     10,    11   ]       
 AbsLengthSpectrum: [ 2000., 2000., 2000., 2000., 2000., 2000., 2000., 2000.] 

 # Rayleigh scattering length (cm) @ 90K as a function of energy (eV) from arXiv:1502.04213
 RayleighEnergies: [   2.80,   3.00,   3.50,   4.00,  5.00,  6.00,  7.00,  8.00,  8.50,  9.00,  9.20,  9.40,  9.50,  9.60,  9.70,  9.80,  9.90,  10.0,  10.2,  10.4,  10.6, 10.8 ]
 RayleighSpectrum: [ 47923., 35981., 18825., 10653., 3972., 1681., 750.9, 334.7, 216.8, 135.0, 109.7, 88.06, 78.32, 69.34, 61.06, 53.46, 46.50, 40.13, 28.91, 19.81, 12.61, 7.20 ]

 # Surface reflectivity data - vector of energy spectrum per
 #   surface type
 ReflectiveSurfaceEnergies:           [ 7, 9, 10 ]             
 ReflectiveSurfaceNames:            [ "STEEL_STAINLESS_Fe7Cr2Ni" ]  
 ReflectiveSurfaceReflectances:     [ [ 0.25, 0.25, 0.25 ] ]        
 ReflectiveSurfaceDiffuseFractions: [ [ 0.5,  0.5,  0.5  ] ]        

 # Information related with the simulation of the Wavelength Shifter (TPB)
 LoadExtraMatProperties: false

 # TPB - WLS
 TpbTimeConstant: 2.5 #wls time constant in s J. Lumin 81(1999) 285

 # WLS - TPB properties original tpb [0.0, 0.0, 0.0, 0.0588,0.235, 0.853, 1.0,1.0,0.9259,0.704,0.0296,0.011, 0.0,0.0, 0.]
 TpbEmmisionEnergies: [0.05,1.0,1.5, 2.25, 2.481, 2.819, 2.952,2.988,3.024, 3.1, 3.14,3.1807, 3.54, 5.5, 50.39]
 TpbEmmisionSpectrum: [0.0, 0.0, 0.0, 0.0588,0.235, 0.853, 1.0,1.0,0.9259,0.704,0.0296,0.011, 0.0,0.0, 0.]
 TpbAbsorptionEnergies: [0.05,1.77,2.0675, 7.42, 7.75, 8.16, 8.73, 9.78,10.69, 50.39]
 TpbAbsorptionSpectrum: [100000.0,100000.0, 100000.0,0.001,0.00000000001,0.00000000001, 0.00000000001, 0.00000000001, 0.00000000001, 0.00000000001]
)cfg"};


// this is a copy of the "standard" DetectorProperties configuration
extern const std::string DetectorPropertiesServiceConfigurationString { R"cfg(
service_provider: "DetectorPropertiesServiceStandard"

# Drift properties 
SternheimerA:     0.1956  # Ar Sternheimer parameter a.
SternheimerK:     3.0000  # Ar Sternheimer parameter k.
SternheimerX0:    0.2000  # Ar Sternheimer parameter x0.
SternheimerX1:    3.0000  # Ar Sternheimer parameter x0.
SternheimerCbar:  5.2146  # Ar Sternheimer parameter Cbar.

Temperature:       87
Electronlifetime:  3.0e3
Efield:           [0.5,0.666,0.8]  #(predicted for microBooNE)
ElectronsToADC:    6.8906513e-3 # 1fC = 43.008 ADC counts for DUNE fd
NumberTimeSamples: 4492         # drift length/drift velocity*sampling rate = (359.4 cm)/(0.16 cm/us)*(2 MHz)
ReadOutWindowSize: 4492         # drift length/drift velocity*sampling rate = (359.4 cm)/(0.16 cm/us)*(2 MHz)
TimeOffsetU:       0.
TimeOffsetV:       0.
TimeOffsetZ:       0.

SimpleBoundaryProcess: true  #enable opticalBoundaryProcessSimple instead of G4 default

)cfg"};
