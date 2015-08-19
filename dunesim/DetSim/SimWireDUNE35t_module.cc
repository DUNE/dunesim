////////////////////////////////////////////////////////////////////////
//
// SimWireDUNE35t class designed to simulate signal on a wire in the TPC
//
//
// jti3@fnal.gov
// - Revised to use sim::RawDigit instead of rawdata::RawDigit, and to
// - save the electron clusters associated with each digit.
//
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <bitset>

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art extensions
#include "artextensions/SeedService/SeedService.hh"

#include "Utilities/LArFFT.h"
#include "RawData/raw.h"
#include "Utilities/LArProperties.h"
#include "dune/Utilities/SignalShapingServiceDUNE35t.h"
#include "Geometry/Geometry.h"

#include "Simulation/sim.h"
#include "Simulation/SimChannel.h"
#include "RawData/RawDigit.h"
#include "Utilities/DetectorProperties.h"

#include "TMath.h"
#include "TComplex.h"
#include "TString.h"
#include "TH2.h"
#include "TH1D.h"
#include "TFile.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

///Detector simulation of raw signals on wires 
namespace detsim {

  // tye used for passing messages for simulating gaps

    typedef enum {
      NONACTIVE, UCOMB, VCOMB, ACTIVE, HORIZGAP, VERTGAP
    } GapType_t;


  // Base class for creation of raw signals on wires. 
  class SimWireDUNE35t : public art::EDProducer {
    
  public:
        
    explicit SimWireDUNE35t(fhicl::ParameterSet const& pset); 
    virtual ~SimWireDUNE35t();
    
    // read/write access to event
    void produce (art::Event& evt);
    void beginJob();
    void endJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    void         GenNoise(std::vector<float>& array);

    std::string            fDriftEModuleLabel;///< module making the ionization electrons
    raw::Compress_t        fCompression;      ///< compression type to use
    unsigned int           fNoiseOn;          ///< noise turned on or off for debugging; default is on
    unsigned int           fNoiseModel;          ///< noise model>
    float                  fNoiseFact;        ///< noise scale factor
    float                  fNoiseWidth;       ///< exponential noise width (kHz)
    float                  fLowCutoff;        ///< low frequency filter cutoff (kHz)
    float                  fNoiseFactZ;        ///< noise scale factor for Z (collection) plane
    float                  fNoiseWidthZ;       ///< exponential noise width (kHz)  for Z (collection) plane
    float                  fLowCutoffZ;        ///< low frequency filter cutoff (kHz) for Z (collection) plane
    float                  fNoiseFactU;        ///< noise scale factor  for U plane
    float                  fNoiseWidthU;       ///< exponential noise width (kHz)   for U plane
    float                  fLowCutoffU;        ///< low frequency filter cutoff (kHz)  for U plane
    float                  fNoiseFactV;        ///< noise scale factor   for V plane
    float                  fNoiseWidthV;       ///< exponential noise width (kHz)   for V plane
    float                  fLowCutoffV;        ///< low frequency filter cutoff (kHz)  for V plane
    unsigned int           fZeroThreshold;    ///< Zero suppression threshold
    int                    fNearestNeighbor;  ///< Maximum distance between hits above threshold before they are separated into different blocks
    unsigned int           fNeighboringChannels; ///< Number of neighboring channels on either side allowed to influence zero suppression
    int                    fNTicks;           ///< number of ticks of the clock
    double                 fSampleRate;       ///< sampling rate in ns
    unsigned int           fNSamplesReadout;  ///< number of ADC readout samples in 1 readout frame
    unsigned int           fNTimeSamples;     ///< number of ADC readout samples in all readout frames (per event)
    unsigned int           fNoiseArrayPoints; ///< number of  points in randomly generated noise array
  
    std::vector<double>    fChargeWork;
    //std::vector< std::vector<float> > fNoise;///< noise on each channel for each time
    std::vector< std::vector<float> > fNoiseZ;///< noise on each channel for each time for Z (collection) plane
    std::vector< std::vector<float> > fNoiseU;///< noise on each channel for each time for U plane
    std::vector< std::vector<float> > fNoiseV;///< noise on each channel for each time for V plane
    
    TH1D*                fNoiseDist;          ///< distribution of noise counts
    

    // variables for simulating the charge deposition in gaps and charge drifting over the comb materials.

    uint32_t               fFirstCollectionChannel;

    // variables for finding the first and last channel numbers on each plane

    std::vector< uint32_t > fFirstChannelsInPlane;
    std::vector< uint32_t > fLastChannelsInPlane;

    //define max ADC value - if one wishes this can
    //be made a fcl parameter but not likely to ever change
    const float adcsaturation = 4095;
    float                  fCollectionPed;    ///< ADC value of baseline for collection plane
    float                  fCollectionPedRMS;    ///< ADC value of baseline RMS for collection plane
    float                  fInductionPed;     ///< ADC value of baseline for induction plane
    float                  fInductionPedRMS;     ///< ADC value of baseline RMS for induction plane
    bool                   fPedestalOn;          // switch for simulation of nonzero pedestals

    // input fcl parameters

    bool                   fSimCombs;          // switch for simulation of the combs
    bool                   fSimStuckBits;      // switch for simulation of stuck bits
    bool                   fSaveEmptyChannel;  // switch for saving channels with all zero entries
    float                  fFractUUCollect;    // fraction of charge that collects on U (non-transparency) when charge drifts over the comb holding U wires
    float                  fFractUVCollect;    // fraction of charge that collects on U (non-transparency) when charge drifts over the comb holding V wires
    float                  fFractVUCollect;    // fraction of charge that collects on V (non-transparency) when charge drifts over the comb holding U wires
    float                  fFractVVCollect;    // fraction of charge that collects on V (non-transparency) when charge drifts over the comb holding V wires
    float                  fFractUUMiss;       // fraction of charge that gets missed on U when charge drifts over the comb holding U
    float                  fFractUVMiss;       // fraction of charge that gets missed on U when charge drifts over the comb holding V
    float                  fFractVUMiss;       // fraction of charge that gets missed on V when charge drifts over the comb holding U
    float                  fFractVVMiss;       // fraction of charge that gets missed on V when charge drifts over the comb holding V
    float                  fFractZUMiss;       // fraction of charge that gets missed on Z (collection)  when charge drifts over the comb holding U
    float                  fFractZVMiss;       // fraction of charge that gets missed on Z (collection)  when charge drifts over the comb holding V
    float                  fFractHorizGapUMiss;     // fraction of charge in the horizontal gap that is missing on U (and not collected)
    float                  fFractVertGapUMiss;     // fraction of charge in the horizontal gaps that is missing on U
    float                  fFractHorizGapVMiss;     // fraction of charge in the horizontal gap that is missing on V
    float                  fFractVertGapVMiss;     // fraction of charge in the horizontal gaps that is missing on V
    float                  fFractHorizGapZMiss;     // fraction of charge in the horizontal gap that is missing on Z (collection)
    float                  fFractVertGapZMiss;     // fraction of charge in the horizontal gaps that is missing on Z (collection
    float                  fFractHorizGapUCollect;     // fraction of charge in the horizontal gap that collects on U
    float                  fFractVertGapUCollect;     // fraction of charge in the horizontal gaps that collects on U
    float                  fFractHorizGapVCollect;     // fraction of charge in the horizontal gap that collects on V
    float                  fFractVertGapVCollect;     // fraction of charge in the horizontal gaps that collects on V

    // boundaries of the combs -- cached here for speed

    double zcomb1,zcomb2,zcomb3,zcomb4,zcomb5,zcomb6;
    double zcomb7,zcomb8,zcomb9,zcomb10,zcomb11,zcomb12;
    double zcomb13,zcomb14,zcomb15,zcomb16,zcomb17,zcomb18;
    double ycomb1,ycomb2,ycomb3,ycomb4,ycomb5,ycomb6;
    double ycomb7,ycomb8,ycomb9,ycomb10,ycomb11,ycomb12;
    double ycomb13,ycomb14,ycomb15,ycomb16,ycomb17,ycomb18;

    GapType_t combtest35t(double x, double y, double z);


    double               fZeroprobs[32];       ///< array of probabilities of 6 LSF bits getting stuck at 000000
    double               fThreefprobs[32];     ///< array of probabilities of 6 LSF bits getting stuck at 111111

    TH1D*                fZeroProbabilities;          ///< histogram for probabilities of 6 LSF bits getting stuck at 000000
    TH1D*                fThreefProbabilities;        ///< histogram for probabilities of 6 LSF bits getting stuck at 111111


  }; // class SimWireDUNE35t

  DEFINE_ART_MODULE(SimWireDUNE35t)

  //-------------------------------------------------
  SimWireDUNE35t::SimWireDUNE35t(fhicl::ParameterSet const& pset)
  {

    this->reconfigure(pset);

    produces< std::vector<raw::RawDigit>   >();

    fCompression = raw::kNone;
    TString compression(pset.get< std::string >("CompressionType"));
    if(compression.Contains("Huffman",TString::kIgnoreCase)) fCompression = raw::kHuffman;    
    if(compression.Contains("ZeroSuppression",TString::kIgnoreCase)) fCompression = raw::kZeroSuppression;

// create a default random engine; obtain the random seed from SeedService,
// unless overridden in configuration with key "Seed"
    art::ServiceHandle<artext::SeedService>()
      ->createEngine(*this, pset, "Seed");

  }

  //-------------------------------------------------
  SimWireDUNE35t::~SimWireDUNE35t()
  {

    fChargeWork.clear();
 
    for(unsigned int i = 0; i < fNoiseZ.size(); ++i) fNoiseZ[i].clear();
    fNoiseZ.clear();
   
    for(unsigned int i = 0; i < fNoiseU.size(); ++i) fNoiseU[i].clear();
    fNoiseU.clear();
   
    for(unsigned int i = 0; i < fNoiseV.size(); ++i) fNoiseV[i].clear();
    fNoiseV.clear();

  }

  //-------------------------------------------------
  void SimWireDUNE35t::reconfigure(fhicl::ParameterSet const& p) 
  {
    fDriftEModuleLabel= p.get< std::string         >("DriftEModuleLabel");


    fNoiseFactZ        = p.get< double              >("NoiseFactZ");
    fNoiseWidthZ       = p.get< double              >("NoiseWidthZ");
    fLowCutoffZ        = p.get< double              >("LowCutoffZ");
    fNoiseFactU        = p.get< double              >("NoiseFactU");
    fNoiseWidthU       = p.get< double              >("NoiseWidthU");
    fLowCutoffU        = p.get< double              >("LowCutoffU");
    fNoiseFactV        = p.get< double              >("NoiseFactV");
    fNoiseWidthV       = p.get< double              >("NoiseWidthV");
    fLowCutoffV        = p.get< double              >("LowCutoffV");
    fZeroThreshold    = p.get< unsigned int         >("ZeroThreshold");
    fNearestNeighbor         = p.get< int           >("NearestNeighbor");
    fNeighboringChannels   = p.get< unsigned int    >("NeighboringChannels");
    fNoiseArrayPoints = p.get< unsigned int         >("NoiseArrayPoints");
    fNoiseOn           = p.get< unsigned int        >("NoiseOn");
    fNoiseModel           = p.get< unsigned int     >("NoiseModel");
    fCollectionPed    = p.get< float                >("CollectionPed");
    fCollectionPedRMS = p.get< float                >("CollectionPedRMS");
    fInductionPed     = p.get< float                >("InductionPed");
    fInductionPedRMS  = p.get< float                >("InductionPedRMS");
    fPedestalOn       = p.get< bool                 >("PedestalOn");  
    art::ServiceHandle<util::DetectorProperties> detprop;
    fSampleRate       = detprop->SamplingRate();
    fNSamplesReadout  = detprop->ReadOutWindowSize();
    fNTimeSamples  = detprop->NumberTimeSamples();
    
    fSimCombs            = p.get< bool >("SimCombs");  
    fSimStuckBits        = p.get< bool >("SimStuckBits");   
    fSaveEmptyChannel    = p.get< bool >("SaveEmptyChannel");  
    fFractUUCollect      = p.get< float >("FractUUCollect");
    fFractUVCollect      = p.get< float >("FractUVCollect");
    fFractVUCollect      = p.get< float >("FractVUCollect");
    fFractVVCollect      = p.get< float >("FractVVCollect");
    fFractUUMiss         = p.get< float >("FractUUMiss");
    fFractUVMiss         = p.get< float >("FractUVMiss");
    fFractVUMiss         = p.get< float >("FractVUMiss");
    fFractVVMiss         = p.get< float >("FractVVMiss");
    fFractZUMiss         = p.get< float >("FractZUMiss");
    fFractZVMiss         = p.get< float >("FractZVMiss");
    fFractHorizGapUMiss  = p.get< float >("FractHorizGapUMiss");
    fFractVertGapUMiss   = p.get< float >("FractVertGapUMiss");
    fFractHorizGapVMiss  = p.get< float >("FractHorizGapVMiss");
    fFractVertGapVMiss   = p.get< float >("FractVertGapVMiss");
    fFractHorizGapZMiss  = p.get< float >("FractHorizGapZMiss");
    fFractVertGapZMiss   = p.get< float >("FractVertGapZMiss");
    fFractHorizGapUCollect  = p.get< float >("FractHorizGapUCollect");
    fFractVertGapUCollect   = p.get< float >("FractVertGapUCollect");
    fFractHorizGapVCollect  = p.get< float >("FractHorizGapVCollect");
    fFractVertGapVCollect   = p.get< float >("FractVertGapVCollect");

    return;
  }

  //-------------------------------------------------
  void SimWireDUNE35t::beginJob() 
  { 

    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;

    fNoiseDist  = tfs->make<TH1D>("Noise", ";Noise  (ADC);", 1000,   -10., 10.);


    art::ServiceHandle<util::LArFFT> fFFT;
    fNTicks = fFFT->FFTSize();

    if ( fNTicks%2 != 0 ) 
      LOG_DEBUG("SimWireDUNE35t") << "Warning: FFTSize not a power of 2. "
				  << "May cause issues in (de)convolution.\n";

    if ( (int)fNSamplesReadout > fNTicks ) 
      mf::LogError("SimWireDUNE35t") << "Cannot have number of readout samples "
				     << "greater than FFTSize!";

    fChargeWork.resize(fNTicks, 0.);
    art::ServiceHandle<geo::Geometry> geo;

    bool foundfirstcollectionchannel = false;
    fFirstChannelsInPlane.push_back(0);
    unsigned int currentPlaneNumber = geo->ChannelToWire(0).at(0).Plane; // ID of current wire plane
    unsigned int currentTPCNumber = geo->ChannelToWire(0).at(0).TPC; // ID of current wire plane

    for (uint32_t ichan=0;ichan<geo->Nchannels();ichan++)
      {

	if(!foundfirstcollectionchannel)
	  {
	    const geo::View_t view = geo->View(ichan);
	    if (view == geo::kZ)
	      {
		foundfirstcollectionchannel = true;
		fFirstCollectionChannel = ichan;
		//break;
	      }
	  }

	const unsigned int thisPlaneNumber = geo->ChannelToWire(ichan).at(0).Plane;
	const unsigned int thisTPCNumber = geo->ChannelToWire(ichan).at(0).TPC;
		
	if(thisPlaneNumber != currentPlaneNumber || (thisPlaneNumber == geo::kZ && thisTPCNumber != currentTPCNumber))
	  {
	    fLastChannelsInPlane.push_back(ichan-1);
	    fFirstChannelsInPlane.push_back(ichan); 
	    currentPlaneNumber = thisPlaneNumber;
	    currentTPCNumber = thisTPCNumber;
	  }

      } 
    if (!foundfirstcollectionchannel)
      {
	throw  cet::exception("SimWireDUNE35t  BeginJob") << " Could not find any collection channels\n";
      }
    
    fLastChannelsInPlane.push_back(geo->Nchannels()-1);

     
    // //Check starting and ending channels for each wire plane
    // for(size_t ip = 0; ip < fFirstChannelsInPlane.size(); ++ip){
    //   std::cout << "First channel in plane " << ip << " is " << fFirstChannelsInPlane.at(ip) << std::endl;
    // }

    // for(size_t ip = 0; ip < fLastChannelsInPlane.size(); ++ip){
    //   std::cout << "Last channel in plane " << ip << " is " << fLastChannelsInPlane.at(ip) << std::endl;
    // }

    //Generate noise if selected to be on
    if(fNoiseOn && fNoiseModel==1){

      //fNoise.resize(geo->Nchannels());
      fNoiseZ.resize(fNoiseArrayPoints);
      fNoiseU.resize(fNoiseArrayPoints);
      fNoiseV.resize(fNoiseArrayPoints);
      
      // GenNoise() will further resize each channel's 
      // fNoise vector to fNoiseArrayPoints long.
      
      for(unsigned int p = 0; p < fNoiseArrayPoints; ++p){
	
	fNoiseFact = fNoiseFactZ;
	fNoiseWidth = fNoiseWidthZ;
	fLowCutoff = fLowCutoffZ;

	GenNoise(fNoiseZ[p]);
	for(int i = 0; i < fNTicks; ++i)
	  fNoiseDist->Fill(fNoiseZ[p][i]);
	
	fNoiseFact = fNoiseFactU;
	fNoiseWidth = fNoiseWidthU;
	fLowCutoff = fLowCutoffU;

	GenNoise(fNoiseU[p]);
	for(int i = 0; i < fNTicks; ++i)	 
	  fNoiseDist->Fill(fNoiseU[p][i]);


	fNoiseFact = fNoiseFactV;
	fNoiseWidth = fNoiseWidthV;
	fLowCutoff = fLowCutoffV;
 
    
	GenNoise(fNoiseV[p]);
	for(int i = 0; i < fNTicks; ++i)
	  fNoiseDist->Fill(fNoiseV[p][i]);
	
      }// end loop over wires
    } 

    // initialize the comb test positions.  This is clumsy here mainly due to the irregular geometry
    // should write something more systematic for the FD.  There is also some duplication as the
    // vertical positions of APA's 0 and 3 are assumed to be the same.  Could think about either adding
    // an exception if they're not, or defining more y positions to hold different APA positions if we want
    // them to be different at a later time.  Simulation may always be perfect though.

    // WireEndPoints takes cryostat, tpc, plane, wire, as ints and returns endpoints
    //geo->WireEndPoints(c,t,p,w,xyzbeg,xyzend);

    // wire endpoints are at the places where the wire hits the comb that supports it.  Bits of
    // wire running over the comb are beyond the endpoints.  So we need to extrapolate.

    double xyzbeg[3],xyzend[3];
    int lastwire = 0;

    // Numbers in comments are from Geometry V3 for debugging purposes.

    // APA 0

    geo->WireEndPoints(0,0,0,0,xyzbeg,xyzend);  // first U wire in TPC 0. 
    zcomb2 = xyzbeg[2];  // 0.0
    ycomb5 = xyzend[1];  // 113.142

    lastwire = geo->Nwires(0,0,0)-1;  // 358 in v3
    geo->WireEndPoints(0,0,0,lastwire,xyzbeg,xyzend);  // last U wire in TPC 0.
    zcomb5 = xyzend[2];  // 50.8929
    ycomb2 = xyzbeg[1];  // -82.9389

    geo->WireEndPoints(0,0,1,0,xyzbeg,xyzend);  // first V wire in TPC 0.  
    zcomb4 = xyzend[2];  //  50.5774
    ycomb4 = xyzbeg[1];  //  113.142

    lastwire = geo->Nwires(1,0,0)-1;  // 344 in v3
    geo->WireEndPoints(0,0,1,lastwire,xyzbeg,xyzend);  // last V wire in TPC 0.  
    zcomb3 = xyzbeg[2];  //  0.3155
    ycomb3 = xyzend[1];  //  -82.6234

    // the collection wires appear to end where they meet their comb.
    //geo->WireEndPoints(0,0,2,0,xyzbeg,xyzend);  // first collection wire in TPC 0
    //ycomb3 = xyzbeg[2];  // -82.308
    //ycomb4 = xyzend[2];  // 113.142

    // need to get zcomb1, zcomb6, ycomb1, and ycomb6 -- extrapolate

    zcomb1 = zcomb2 - (zcomb3 - zcomb2);
    zcomb6 = zcomb5 + (zcomb5 - zcomb4);
    ycomb1 = ycomb2 - (ycomb3 - ycomb2);
    ycomb6 = ycomb5 + (ycomb5 - ycomb4);


    // APA 1

    geo->WireEndPoints(0,2,0,0,xyzbeg,xyzend);  // first U wire in TPC 2. 
    zcomb11 = xyzend[2];  // 102.817
    ycomb8 = xyzbeg[1];  // -85.221

    lastwire = geo->Nwires(0,2,0)-1;  // 194 in v3
    geo->WireEndPoints(0,2,0,lastwire,xyzbeg,xyzend);  // last U wire in TPC 2.
    zcomb8 = xyzbeg[2];  // 51.924
    ycomb11 = xyzend[1];  // -0.831

    geo->WireEndPoints(0,2,1,0,xyzbeg,xyzend);  // first V wire in TPC 2.  
    zcomb9 = xyzbeg[2];  //  52.2395 
    ycomb9 = xyzend[1];  //  -85.222

    lastwire = geo->Nwires(1,2,0)-1;  // 188 in v3
    geo->WireEndPoints(0,2,1,lastwire,xyzbeg,xyzend);  // last V wire in TPC 2.  
    zcomb10 = xyzend[2];  //  102.501
    ycomb10 = xyzbeg[1];  //  -1.14655

    //geo->WireEndPoints(0,2,2,0,xyzbeg,xyzend);  // first collection wire in TPC 2
    //ycombx = xyzbeg[2];  // -85.222   edges of the combs
    //ycombx = xyzend[2];  // -1.46205

    // need to get zcomb7, zcomb12, ycomb7, and ycomb12 -- extrapolate

    zcomb7 = zcomb8 - (zcomb9 - zcomb8);
    zcomb12 = zcomb11 + (zcomb11 - zcomb10);
    ycomb7 = ycomb8 - (ycomb9 - ycomb8);
    ycomb12 = ycomb11 + (ycomb11 - ycomb10);

    // APA 2

    geo->WireEndPoints(0,4,0,0,xyzbeg,xyzend);  // first U wire in TPC 4.
    zcomb8 = xyzbeg[2]; // 51.924 -- same as above
    ycomb17 = xyzend[1];  // 113.142 -- same as above 

    lastwire = geo->Nwires(0,4,0)-1;  // 235 in v3
    geo->WireEndPoints(0,4,0,lastwire,xyzbeg,xyzend);  // last U wire in TPC 4.
    zcomb11 = xyzend[2];  // 102.817 -- same as above 
    ycomb14 = xyzbeg[1];  // 0.83105 

    geo->WireEndPoints(0,4,1,0,xyzbeg,xyzend);  // first V wire in TPC 4.  
    zcomb10 = xyzend[2];  //   102.501 -- same as above
    ycomb16 = xyzbeg[1];  //  113.142 -- everything ends here in y

    lastwire = geo->Nwires(1,4,0)-1;  // 227 in v3
    geo->WireEndPoints(0,4,1,lastwire,xyzbeg,xyzend);  // last V wire in TPC 4.  
    zcomb9 = xyzbeg[2];  //  52.2395  -- same as above
    ycomb15 = xyzend[1];  //  1.14655

    //geo->WireEndPoints(0,4,2,0,xyzbeg,xyzend);  // first collection wire in TPC 1
    //ycombx = xyzbeg[2];  // 52.2234   edges of the combs -- not what we want
    //ycombx = xyzend[2];  // 113.142   for this

    // need to get zcomb7, zcomb12, ycomb13, and ycomb18 -- extrapolate
    // the z's are just recalculations of the numbers above

    zcomb7 = zcomb8 - (zcomb9 - zcomb8);
    zcomb12 = zcomb11 + (zcomb11 - zcomb10);
    ycomb13 = ycomb14 - (ycomb15 - ycomb14);
    ycomb18 = ycomb17 + (ycomb17 - ycomb16);

    // APA 3 -- a lot like APA 0

    geo->WireEndPoints(0,6,0,0,xyzbeg,xyzend);  // first U wire in TPC 6.
    zcomb14 = xyzbeg[2];  // 103.84
    ycomb5 = xyzend[1];  //  113.142 -- same as above

    lastwire = geo->Nwires(0,6,0)-1;  // 358 in v3
    geo->WireEndPoints(0,6,0,lastwire,xyzbeg,xyzend);  // last U wire in TPC 6.
    zcomb17 = xyzend[2];  // 154.741
    ycomb2 = xyzbeg[1];  // -82.9389 -- same as above

    geo->WireEndPoints(0,6,1,0,xyzbeg,xyzend);  // first V wire in TPC 6.  
    zcomb16 = xyzend[2];  //  154.425
    ycomb4 = xyzbeg[1];  //  113.142 -- same as above

    lastwire = geo->Nwires(1,6,0)-1;  // 344 in v3
    geo->WireEndPoints(0,6,1,lastwire,xyzbeg,xyzend);  // last V wire in TPC 6.  
    zcomb15 = xyzbeg[2];  //  104.164
    ycomb3 = xyzend[1];  //  -82.6234 -- same as above

    // the collection wires appear to end where they meet their comb.
    //geo->WireEndPoints(0,6,2,0,xyzbeg,xyzend);  // first collection wire in TPC 0
    //ycomb3 = xyzbeg[2];  // -82.308
    //ycomb4 = xyzend[2];  // 113.142

    // need to get zcomb13, zcomb18, ycomb1, and ycomb6 -- extrapolate
    // the ycomb1 and ycomb6 are just copies.

    zcomb13 = zcomb14 - (zcomb15 - zcomb14);
    zcomb18 = zcomb17 + (zcomb17 - zcomb16);
    ycomb1 = ycomb2 - (ycomb3 - ycomb2);
    ycomb6 = ycomb5 + (ycomb5 - ycomb4);

    if(fSimStuckBits){
  
      fZeroprobs[0] = 0.38;
      fZeroprobs[1] = 0.38;
      fZeroprobs[2] = 0.245;
      fZeroprobs[3] = 0.155;
      fZeroprobs[4] = 0.12;
      fZeroprobs[5] = 0.082;
      fZeroprobs[6] = 0.07;
      fZeroprobs[7] = 0.075;
      fZeroprobs[8] = 0.078;
      fZeroprobs[9] = 0.08;
      fZeroprobs[10] = 0.082;
      fZeroprobs[11] = 0.105;
      fZeroprobs[12] = 0.11;
      fZeroprobs[13] = 0.112;
      fZeroprobs[14] = 0.11;
      fZeroprobs[15] = 0.075;
      fZeroprobs[16] = 0.065;
      fZeroprobs[17] = 0.125;
      fZeroprobs[18] = 0.145;
      fZeroprobs[19] = 0.145;
      fZeroprobs[20] = 0.1;
      fZeroprobs[21] = 0.145;
      fZeroprobs[22] = 0.175;
      fZeroprobs[23] = 0.21;
      fZeroprobs[24] = 0.245;
      fZeroprobs[25] = 0.25;
      fZeroprobs[26] = 0.24;
      fZeroprobs[27] = 0.265;
      fZeroprobs[28] = 0.305;
      fZeroprobs[29] = 0.26;
      fZeroprobs[30] = 0.275;
      fZeroprobs[31] = 0.27;
      
      fThreefprobs[0] = 0.27;
      fThreefprobs[1] = 0.2175;
      fThreefprobs[2] = 0.23;
      fThreefprobs[3] = 0.18;
      fThreefprobs[4] = 0.12;
      fThreefprobs[5] = 0.178;
      fThreefprobs[6] = 0.125;
      fThreefprobs[7] = 0.11;
      fThreefprobs[8] = 0.10;
      fThreefprobs[9] = 0.097;
      fThreefprobs[10] = 0.05;
      fThreefprobs[11] = 0.055;
      fThreefprobs[12] = 0.05;
      fThreefprobs[13] = 0.054;
      fThreefprobs[14] = 0.02;
      fThreefprobs[15] = 0.02;
      fThreefprobs[16] = 0.02;
      fThreefprobs[17] = 0.018;
      fThreefprobs[18] = 0.025;
      fThreefprobs[19] = 0.017;
      fThreefprobs[20] = 0.018;
      fThreefprobs[21] = 0.017;
      fThreefprobs[22] = 0.025;
      fThreefprobs[23] = 0.0252;
      fThreefprobs[24] = 0.02;
      fThreefprobs[25] = 0.035;
      fThreefprobs[26] = 0.05;
      fThreefprobs[27] = 0.065;
      fThreefprobs[28] = 0.07;
      fThreefprobs[29] = 0.12;
      fThreefprobs[30] = 0.165;
      fThreefprobs[31] = 0.38;

      fZeroProbabilities  = tfs->make<TH1D>("ZeroProbabilities", "Probability of 6 LSBs sticking to 0x00; Input voltage (ADC modulo 64);", 32,   0., 64.);
      fThreefProbabilities  = tfs->make<TH1D>("ThreefProbabilities", "Probability of 6 LSBs sticking to 0x3f; Input voltage (ADC modulo 64);", 32,   0., 64.);

      for(int i=0;i<32;++i){
	fZeroProbabilities->SetBinContent(i+1,fZeroprobs[i]);
	fThreefProbabilities->SetBinContent(i+1,fThreefprobs[i]);
      }

    }
    return;

  }

  //-------------------------------------------------
  void SimWireDUNE35t::endJob() 
  {
  }

  //-------------------------------------------------
  void SimWireDUNE35t::produce(art::Event& evt)
  {
    // get the geometry to be able to figure out signal types and chan -> plane mappings
    art::ServiceHandle<geo::Geometry> geo;
    unsigned int signalSize = fNTicks;

    // vectors for working
    std::vector<short>    adcvec(signalSize, 0);	
    std::vector<const sim::SimChannel*> chanHandle;
    evt.getView(fDriftEModuleLabel,chanHandle);

    //Get fIndShape and fColShape from SignalShapingService, on the fly
    art::ServiceHandle<util::SignalShapingServiceDUNE35t> sss;

    // make a vector of const sim::SimChannel* that has same number
    // of entries as the number of channels in the detector
    // and set the entries for the channels that have signal on them
    // using the chanHandle
    std::vector<const sim::SimChannel*> channels(geo->Nchannels());
    for(size_t c = 0; c < chanHandle.size(); ++c){
      channels[chanHandle[c]->Channel()] = chanHandle[c];
    }
    
    // make an unique_ptr of sim::SimDigits that allows ownership of the produced
    // digits to be transferred to the art::Event after the put statement below
    std::unique_ptr< std::vector<raw::RawDigit>   >  digcol(new std::vector<raw::RawDigit>);
	  
    unsigned int chan = 0; 
    fChargeWork.clear();
    fChargeWork.resize(fNTicks, 0.);
	  
 
    std::vector<double> fChargeWorkCollInd;

    art::ServiceHandle<util::LArFFT> fFFT;

    // Add all channels  
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat flat(engine);

    std::map<int,double>::iterator mapIter;      


    // make ring buffer to hold neighboring channels in order to enable nearest neighbor-influenced zero suppression

    boost::circular_buffer<std::vector<short>> adcvec_neighbors(fNeighboringChannels*2+1); 

    // make vector of adc vectors to hold first channels in induction plane in order to wrap around to the start for nearest neighbor-influenced zero suppression

    std::vector<std::vector<short>> adcvec_inductionplanestart; 

    unsigned int plane_number = 0;

    for(chan = 0; chan < geo->Nchannels(); chan++) {    
      
 
      fChargeWork.clear();    
      //      fChargeWork.resize(fNTicks, 0.);    
      fChargeWork.resize(fNTimeSamples, 0.);    
      if (fSimCombs)
	{
	  fChargeWorkCollInd.clear();
	  fChargeWorkCollInd.resize(fNTimeSamples, 0.);
	}

      // get the sim::SimChannel for this channel
      const sim::SimChannel* sc = channels[chan];
      const geo::View_t view = geo->View(chan);


      if( sc ){      
	// loop over the tdcs and grab the number of electrons for each
	for(size_t t = 0; t < fChargeWork.size(); ++t) 
	  if (fSimCombs)
	    {
	      const std::vector<sim::IDE> ides = sc->TrackIDsAndEnergies(t,t);
	      for (auto const &ide : ides)
		{
		  GapType_t gaptype = combtest35t(ide.x,ide.y,ide.z);
		  switch (gaptype)
		    {
		    case ACTIVE:
		      {
			fChargeWork[t] += ide.numElectrons;
			break;
		      }
		    case UCOMB:
		      {
			switch (view)
			  {
			  case geo::kU:
			    {
			      fChargeWork[t] += ide.numElectrons * (1.0-fFractUUCollect-fFractUUMiss);
			      fChargeWorkCollInd[t] += ide.numElectrons * fFractUUCollect;
			      break;
			    }
			  case geo::kV:
			    {
			      fChargeWork[t] += ide.numElectrons * (1.0-fFractVUCollect-fFractUUCollect-fFractVUMiss);
			      fChargeWorkCollInd[t] += ide.numElectrons * fFractVUCollect;
			      break;
			    }
			  case geo::kZ:
			    {
			      fChargeWork[t] += ide.numElectrons * (1.0-fFractVUCollect-fFractUUCollect-fFractZUMiss);
			      break;
			    }
			  default:
			    {
			      throw cet::exception("SimWireDUNE35t") << "ILLEGAL VIEW Type: " << view <<"\n";
			    }
			  }
			break;
		      }
		    case VCOMB:
		      {
			switch (view)
			  {
			  case geo::kU:
			    {
			      fChargeWork[t] += ide.numElectrons * (1.0-fFractUVCollect-fFractUVMiss);
			      fChargeWorkCollInd[t] += ide.numElectrons * fFractUVCollect;
			      break;
			    }
			  case geo::kV:
			    {
			      fChargeWork[t] += ide.numElectrons * (1.0-fFractUVCollect-fFractVVCollect-fFractVVMiss);
			      fChargeWorkCollInd[t] += ide.numElectrons * fFractVVCollect;
			      break;
			    }
			  case geo::kZ:
			    {
			      fChargeWork[t] += ide.numElectrons * (1.0-fFractVVCollect-fFractUVCollect-fFractZVMiss);
			      break;
			    }
			  default:
			    {
			      throw cet::exception("SimWireDUNE35t") << "ILLEGAL VIEW Type: " << view <<"\n";
			    }
			  }
			break;
		      }
		    case HORIZGAP:
		      {
			switch (view)
			  {
			  case geo::kU:
			    {
			      fChargeWork[t] += ide.numElectrons * (1.0-fFractHorizGapUMiss-fFractHorizGapUCollect);
			      fChargeWorkCollInd[t] += ide.numElectrons * fFractHorizGapUCollect;
			      break;
			    }
			  case geo::kV:
			    {
			      fChargeWork[t] += ide.numElectrons * (1.0-fFractHorizGapVMiss-fFractHorizGapUCollect-fFractHorizGapVCollect);
			      fChargeWorkCollInd[t] += ide.numElectrons * fFractHorizGapVCollect;
			      break;
			    }
			  case geo::kZ:
			    {
			      fChargeWork[t] += ide.numElectrons * (1.0-fFractHorizGapZMiss-fFractHorizGapUCollect-fFractHorizGapVCollect);
			      break;
			    }
			  default:
			    {
			      throw cet::exception("SimWireDUNE35t") << "ILLEGAL VIEW Type: " << view <<"\n";
			    }
			  }
			break;
		      }
		    case VERTGAP:
		      {
			switch (view)
			  {
			  case geo::kU:
			    {
			      fChargeWork[t] += ide.numElectrons * (1.0-fFractVertGapUMiss-fFractVertGapUCollect);
			      fChargeWorkCollInd[t] += ide.numElectrons * fFractVertGapUCollect;
			      break;
			    }
			  case geo::kV:
			    {
			      fChargeWork[t] += ide.numElectrons * (1.0-fFractVertGapVMiss-fFractVertGapUCollect-fFractVertGapVCollect);
			      fChargeWorkCollInd[t] += ide.numElectrons * fFractVertGapVCollect;
			      break;
			    }
			  case geo::kZ:
			    {
			      fChargeWork[t] += ide.numElectrons * (1.0-fFractVertGapZMiss-fFractVertGapUCollect-fFractVertGapVCollect);
			      break;
			    }
			  default:
			    {
			      throw cet::exception("SimWireDUNE35t") << "ILLEGAL VIEW Type: " << view <<"\n";
			    }
			  }
			break;
		      }
		    case NONACTIVE: 
		      { 
			break;
		      }
		    }		  
		}
	      // the line all this replaced.
	      // fChargeWork[t] = sc->Charge(t); 
	    }
	  else
	    {
	      fChargeWork[t] = sc->Charge(t);
	      //if (chan == 180 ) std::cout << "Xin1: " << t << " " << fChargeWork[t] << std::endl;
	    }      

        // Convolve charge with appropriate response function 

	fChargeWork.resize(fNTicks,0);
	sss->Convolute(chan,fChargeWork);
	// if (chan == 180 ) {
	//   for(size_t t = 0; t < fChargeWork.size(); ++t) {
	//     std::cout << "Xin2: " << t << " " << fChargeWork[t] << std::endl;
	//   }
	// }

	fChargeWorkCollInd.resize(fNTicks,0);
        sss->Convolute(fFirstCollectionChannel,fChargeWorkCollInd); 

      }

      float ped_mean = fCollectionPed;
      float ped_rms = fCollectionPedRMS;
      geo::SigType_t sigtype = geo->SignalType(chan);
      if (sigtype == geo::kInduction){
        ped_mean = fInductionPed;
	ped_rms = fInductionPedRMS;
      }
      else if (sigtype == geo::kCollection){
        ped_mean = fCollectionPed;
	ped_rms = fCollectionPedRMS;
      }

      // noise was already generated for each wire in the event
      // raw digit vec is already in channel order
      // pick a new "noise channel" for every channel  - this makes sure    
      // the noise has the right coherent characteristics to be on one channel

      int noisechan = nearbyint(flat.fire()*(1.*(fNoiseArrayPoints-1)+0.1));
      // optimize for speed -- access vectors as arrays 

      double *fChargeWork_a=0;
      double *fChargeWorkCollInd_a=0;
      short *adcvec_a=0;
      float *noise_a_U=0;
      float *noise_a_V=0;
      float *noise_a_Z=0;

      if (signalSize>0)	{
	fChargeWork_a = fChargeWork.data();
	fChargeWorkCollInd_a = fChargeWorkCollInd.data();
	adcvec_a = adcvec.data();
	if (fNoiseOn && fNoiseModel==1) {
          noise_a_U=(fNoiseU[noisechan]).data();
	  noise_a_V=(fNoiseV[noisechan]).data();
	  noise_a_Z=(fNoiseZ[noisechan]).data();
	}
      }

      float tmpfv=0;  // this is here so we do our own rounding from floats to short ints (saves CPU time)
      float tnoise=0;

      if (view != geo::kU && view != geo::kV && view != geo::kZ) {
	mf::LogError("SimWireDUNE35t") << "ERROR: CHANNEL NUMBER " << chan << " OUTSIDE OF PLANE";
      }

      if(fNoiseOn && fNoiseModel==1) {	      
	for(unsigned int i = 0; i < signalSize; ++i){
	  if(view==geo::kU)       { tnoise = noise_a_U[i]; }
	  else if (view==geo::kV) { tnoise = noise_a_V[i]; }
	  else                    { tnoise = noise_a_Z[i]; }
          tmpfv = tnoise + fChargeWork_a[i] ;
	  if (fSimCombs)  tmpfv += fChargeWorkCollInd_a[i];
	  //allow for ADC saturation
	  if ( tmpfv > adcsaturation - ped_mean)
	    tmpfv = adcsaturation- ped_mean;
	  //don't allow for "negative" saturation
	  if ( tmpfv < 0 - ped_mean)
	    tmpfv = 0- ped_mean;

	  adcvec_a[i] = (tmpfv >=0) ? (short) (tmpfv+0.5) : (short) (tmpfv-0.5); 
	}
      }else if (fNoiseOn && fNoiseModel==2){

	float fASICGain      = sss->GetASICGain(chan);  
	
	double fShapingTime   = sss->GetShapingTime(chan);
	std::map< double, int > fShapingTimeOrder;
	fShapingTimeOrder = { {0.5, 0}, {1.0, 1}, {2.0, 2}, {3.0, 3} };
	DoubleVec              fNoiseFactVec;

	//

	auto tempNoiseVec = sss->GetNoiseFactVec();

	if ( fShapingTimeOrder.find( fShapingTime ) != fShapingTimeOrder.end() ){
	  size_t i = 0;
	  fNoiseFactVec.resize(2);
	  for (auto& item : tempNoiseVec) {
	    fNoiseFactVec[i]   = item.at( fShapingTimeOrder.find( fShapingTime )->second );
	    fNoiseFactVec[i] *= fASICGain/4.7;
	    ++i;
	  }
	}
	else {//Throw exception...
	  throw cet::exception("SimWireDUNE35t")
	    << "\033[93m"
	    << "Shaping Time received from signalservices_dune.fcl is not one of allowed values"
	    << std::endl
	    << "Allowed values: 0.5, 1.0, 2.0, 3.0 usec"
	    << "\033[00m"
	    << std::endl;
	}
	//std::cout << "Xin " << fASICGain << " " << fShapingTime << " " << fNoiseFactVec[0] << " " << fNoiseFactVec[1] << std::endl;

	art::ServiceHandle<art::RandomNumberGenerator> rng;
	CLHEP::HepRandomEngine &engine = rng->getEngine();
	CLHEP::RandGaussQ rGauss_Ind(engine, 0.0, fNoiseFactVec[0]);
	CLHEP::RandGaussQ rGauss_Col(engine, 0.0, fNoiseFactVec[1]);


	for(unsigned int i = 0; i < signalSize; ++i){
	  if(view==geo::kU)       { tnoise = rGauss_Ind.fire(); }
	  else if (view==geo::kV) { tnoise = rGauss_Ind.fire(); }
	  else                    { tnoise = rGauss_Col.fire(); }
          tmpfv = tnoise + fChargeWork_a[i] ;
	  if (fSimCombs)  tmpfv += fChargeWorkCollInd_a[i];
	  //allow for ADC saturation
	  if ( tmpfv > adcsaturation - ped_mean)
	    tmpfv = adcsaturation- ped_mean;
	  //don't allow for "negative" saturation
	  if ( tmpfv < 0 - ped_mean)
	    tmpfv = 0- ped_mean;
	  adcvec_a[i] = (tmpfv >=0) ? (short) (tmpfv+0.5) : (short) (tmpfv-0.5); 
	}
      }else {   // no noise, so just round the values to nearest short ints and store them
	for(unsigned int i = 0; i < signalSize; ++i){
	  tmpfv = fChargeWork_a[i];
	  if (fSimCombs) tmpfv += fChargeWorkCollInd_a[i] ;
	  //allow for ADC saturation
	  if ( tmpfv > adcsaturation - ped_mean)
	    tmpfv = adcsaturation- ped_mean;
	  //don't allow for "negative" saturation
	  if ( tmpfv < 0 - ped_mean)
	    tmpfv = 0- ped_mean;
	  adcvec_a[i] = (tmpfv >=0) ? (short) (tmpfv+0.5) : (short) (tmpfv-0.5); 
	}
      }

      // resize the adcvec to be the correct number of time samples, 
      // just drop the extra samples


      adcvec.resize(fNSamplesReadout);

      // add pedestal values
      if(fPedestalOn)
	{
	  art::ServiceHandle<art::RandomNumberGenerator> rng;
	  CLHEP::HepRandomEngine &engine = rng->getEngine();
	  CLHEP::RandGaussQ rGauss_Ped(engine, 0.0, ped_rms);
	  for(unsigned int i = 0; i < signalSize; ++i){
	    
	    float ped_variation = rGauss_Ped.fire();
	    tmpfv = adcvec_a[i] + ped_mean + ped_variation;
	  
	    adcvec_a[i] = (short) tmpfv; 
	  
	  }
	}

 


      if(fSimStuckBits)//
      	{



	  for(size_t i = 0; i < adcvec.size(); ++i){

	    art::ServiceHandle<art::RandomNumberGenerator> rng;
	    CLHEP::HepRandomEngine &engine = rng->getEngine();
	    CLHEP::RandFlat flat(engine);
	    
	    
	    double rnd[2] = {0.};
	    flat.fireArray(2,rnd,0,1);

	    unsigned int zeromask = 0xffc0;
	    unsigned int onemask = 0x003f;

	    unsigned int sixlsbs = adcvec_a[i] & onemask;

	    int probability_index = (int)sixlsbs/2.0;

	    if(adcvec_a[i] != 0){
	      if(rnd[0] < rnd[1]){
		if(rnd[0] < fZeroprobs[probability_index]){  // 6 LSBs are stuck at 0
		  
		  adcvec_a[i] = adcvec_a[i] & zeromask;
		  
		}
	      }
	      
	      else{ 
		if(rnd[1] < fThreefprobs[probability_index]){  // 6 LSBs are stuck at 3F
		  
		  adcvec_a[i] = adcvec_a[i] | onemask;
		  
		}	      
	      }
	    }
	  }
 
      	}

      if(fNeighboringChannels==0){ // case where neighboring channels are disregarded in zero suppression


	// compress the adc vector using the desired compression scheme,
	// if raw::kNone is selected nothing happens to adcvec
	// This shrinks adcvec, if fCompression is not kNone.
	
	
	raw::Compress(adcvec, fCompression, fZeroThreshold, fNearestNeighbor); 


	raw::RawDigit rd(chan, fNSamplesReadout, adcvec, fCompression);
	rd.SetPedestal(ped_mean,ped_rms);

	adcvec.resize(signalSize);        // Then, resize adcvec back to full length.  Do not initialize to zero (slow)
	if(fSaveEmptyChannel || adcvec[1]>0)
	  digcol->push_back(rd);            // add this digit to the collection
	
	
      }
      else{ //case where contribution of neighboring channels is included in zero suppression


	if(sigtype == geo::kCollection)
	  {
	      
	
	    // push the adc vector to the ring buffer to enable zero suppression with neighboring channel contributions
      
	    adcvec_neighbors.push_back(adcvec);


	    if(!(adcvec_neighbors.full()))
	      {
		if(adcvec_neighbors.size()>fNeighboringChannels){ // apply zero suppression to entries at start of collection plane, once ring buffer is full enough with neighbors
		  
		  adcvec = adcvec_neighbors.at(adcvec_neighbors.size()-fNeighboringChannels);
		  raw::Compress(adcvec_neighbors, adcvec, fCompression, fZeroThreshold, fNearestNeighbor); // apply zero suppression to entries at start of collection plane, once ring buffer is full enough with neighbors
		  
		  raw::RawDigit rd(chan-fNeighboringChannels, fNSamplesReadout, adcvec, fCompression);
		  rd.SetPedestal(ped_mean,ped_rms);
	
		  if(fSaveEmptyChannel || adcvec[1]>0)
		    digcol->push_back(rd);            // add this digit to the collection

		}
	      }	    
	    else{ // ring buffer is full
	      
	      // compress the adc vector using the desired compression scheme,
	      // if raw::kNone is selected nothing happens to adcvec
	      // This shrinks adcvec, if fCompression is not kNone.
	      
	      adcvec = adcvec_neighbors.at(fNeighboringChannels);
	      raw::Compress(adcvec_neighbors, adcvec, fCompression, fZeroThreshold, fNearestNeighbor); // apply zero suppression to entry in middle of ring buffer
	      
	      raw::RawDigit rd(chan-fNeighboringChannels, fNSamplesReadout, adcvec, fCompression);
	      rd.SetPedestal(ped_mean,ped_rms);

	      if(fSaveEmptyChannel || adcvec[1]>0)
		digcol->push_back(rd);            // add this digit to the collection

	      if(chan == fLastChannelsInPlane.at(plane_number)) // End of collection plane is reached, so apply zero suppression to last entries as well
		{
		  

		  unsigned int channel_number = chan - fNeighboringChannels;
		  //std::cout << "We have reached the end of a collection plane!" << std::endl;
		  for(size_t lastentries = 0; lastentries < fNeighboringChannels; ++lastentries)
		    { 
		      ++channel_number;
		   
		      adcvec_neighbors.pop_front(); // remove first entry from ring buffer, to exclude it from nearest neighboring wire checks in zero
		  	
		      adcvec = adcvec_neighbors.at(fNeighboringChannels);
		      raw::Compress(adcvec_neighbors, adcvec, fCompression, fZeroThreshold, fNearestNeighbor); // apply zero suppression to entry in middle of ring buffer
		      
		      raw::RawDigit rd2(channel_number, fNSamplesReadout, adcvec, fCompression);
		      rd2.SetPedestal(ped_mean,ped_rms);

		      if(fSaveEmptyChannel || adcvec[1]>0)
			digcol->push_back(rd2);            // add this digit to the collection

		    }
		  
		  adcvec_neighbors.clear(); // clear ring buffer for next wire plane
		}
	    } // ring buffer is full

	  } // collection plane
	else if(sigtype == geo::kInduction)
	  {
	    // push the adc vector to the ring buffer to enable zero suppression with neighboring channel contributions
      

	    adcvec_neighbors.push_back(adcvec);
	    if(!(adcvec_neighbors.full()))
	      {
		adcvec_inductionplanestart.push_back(adcvec);


	      }
	    else //ring buffer is full
	      {

		// compress the adc vector using the desired compression scheme,
		// if raw::kNone is selected nothing happens to adcvec
		// This shrinks adcvec, if fCompression is not kNone.
		
		adcvec = adcvec_neighbors.at(fNeighboringChannels);
		raw::Compress(adcvec_neighbors, adcvec, fCompression, fZeroThreshold, fNearestNeighbor); // apply zero suppression to entry in middle of ring buffer
		
		raw::RawDigit rd(chan-fNeighboringChannels, fNSamplesReadout, adcvec, fCompression);
		rd.SetPedestal(ped_mean,ped_rms);

		if(fSaveEmptyChannel || adcvec[1]>0)
		  digcol->push_back(rd);            // add this digit to the collection
		
	      }
	    if(chan==fLastChannelsInPlane.at(plane_number)){ // reached the last channel of the induction plane

	      unsigned int channel_number = chan-fNeighboringChannels;

	      for(size_t lastentries = 0; lastentries < 2*fNeighboringChannels; ++lastentries)
		{ 
		  
		  ++channel_number;

		  if(channel_number > fLastChannelsInPlane.at(plane_number))
		    channel_number = fFirstChannelsInPlane.at(plane_number); // set channel number back to start of induction plane after looping around

		  //std::cout << "Channel number of looping-around sections of induction plane = " << channel_number << std::endl;
		  
		  adcvec_neighbors.push_back(adcvec_inductionplanestart.at(lastentries)); // push channel from start of plane onto ring buffer
		  	
		  adcvec = adcvec_neighbors.at(fNeighboringChannels);	

		  raw::Compress(adcvec_neighbors, adcvec, fCompression, fZeroThreshold, fNearestNeighbor); // apply zero suppression to entry in middle of ring buffer
		  
		  raw::RawDigit rd2(channel_number, fNSamplesReadout, adcvec, fCompression);
		  rd2.SetPedestal(ped_mean,ped_rms);

		  if(fSaveEmptyChannel || adcvec[1]>0)
		    digcol->push_back(rd2);            // add this digit to the collection

		}
	      
	      adcvec_neighbors.clear(); // clear ring buffer for next wire plane
	      adcvec_inductionplanestart.clear(); // clear vector of starting adc vectors for next induction plane
	    }
	    
	  } // induction plane

      } // zero suppression with nearest neighboring wire influence complete
    
      adcvec.resize(signalSize);        // Then, resize adcvec back to full length.  Do not initialize to zero (slow)
 
      if(chan==fLastChannelsInPlane.at(plane_number))
	++plane_number;

    }// end loop over channels      

    evt.put(std::move(digcol));

    return;
  }

  //-------------------------------------------------
  void SimWireDUNE35t::GenNoise(std::vector<float>& noise)
  {
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine();
    CLHEP::RandFlat flat(engine);

    noise.clear();
    noise.resize(fNTicks,0.0);
    // noise in frequency space
    std::vector<TComplex> noiseFrequency(fNTicks/2+1, 0.);

    double pval = 0.; 
    double lofilter = 0.;
    double phase = 0.;
    double rnd[2] = {0.};

    // width of frequencyBin in kHz
    double binWidth = 1.0/(fNTicks*fSampleRate*1.0e-6);
    for(int i=0; i< fNTicks/2+1; ++i){
      // exponential noise spectrum 
      pval = fNoiseFact*exp(-(double)i*binWidth/fNoiseWidth);
      // low frequency cutoff     
      lofilter = 1.0/(1.0+exp(-(i-fLowCutoff/binWidth)/0.5));
      // randomize 10%
      flat.fireArray(2,rnd,0,1);
      pval *= lofilter*(0.9+0.2*rnd[0]);
      // random pahse angle
      phase = rnd[1]*2.*TMath::Pi();

      TComplex tc(pval*cos(phase),pval*sin(phase));
      noiseFrequency[i] += tc;
    }
  
   
    // inverse FFT MCSignal
    art::ServiceHandle<util::LArFFT> fFFT;
    fFFT->DoInvFFT(noiseFrequency, noise);

    noiseFrequency.clear();

    // multiply each noise value by fNTicks as the InvFFT 
    // divides each bin by fNTicks assuming that a forward FFT
    // has already been done.
    for(unsigned int i = 0; i < noise.size(); ++i) noise[i] *= 1.*fNTicks;

    return;
  }



  //-------------------------------------------------


  // see the ASCII cartoon of APA's at the bottom of this file for a picture of what all the boundaries are

  //-------------------------------------------------
  GapType_t SimWireDUNE35t::combtest35t(double x, double y, double z)
  {

    if (z<zcomb1) return VERTGAP;  // off to the side of the first APA -- kind of like being in a vertical gap
    if (z<zcomb2) return UCOMB;  // over U comb
    if (z<zcomb3) return VCOMB;  // over V comb
    if (z<zcomb4) 
      {
	if (y<ycomb1) return HORIZGAP; // below the bottom
	if (y<ycomb2) return UCOMB; // over U comb
	if (y<ycomb3) return VCOMB; // over V comb
	if (y<ycomb4) return ACTIVE; // active volume
	if (y<ycomb5) return VCOMB; // over V comb
	if (y<ycomb6) return UCOMB; // over U comb
	return HORIZGAP; // outside top edge

      }
    if (z<zcomb5) return VCOMB;  // over V comb
    if (z<zcomb6) return UCOMB;  // over U comb

    if (z<zcomb7) return VERTGAP; // in gap
    if (z<zcomb8) return UCOMB; // over U comb
    if (z<zcomb9) return VCOMB; // over V comb
    if (z<zcomb10) 
      {
	if (y<ycomb7) return HORIZGAP; // off the bottom
	if (y<ycomb8) return UCOMB; // over U comb
	if (y<ycomb9) return VCOMB; // over V comb
	if (y<ycomb10) return ACTIVE; // active
	if (y<ycomb11) return VCOMB; // over V comb
	if (y<ycomb12) return UCOMB; // over U comb
	if (y<ycomb13) return HORIZGAP; // over gap
	if (y<ycomb14) return UCOMB; // over U comb
	if (y<ycomb15) return VCOMB; // over V comb
	if (y<ycomb16) return ACTIVE; // active volume
	if (y<ycomb17) return VCOMB; // over V comb
	if (y<ycomb18) return UCOMB; // over U comb
	return HORIZGAP;  // above the top edge
      }
    if (z<zcomb11) return VCOMB;  // over V comb
    if (z<zcomb12) return UCOMB;  // over U comb

    if (z<zcomb13) return VERTGAP;  // outside first APA
    if (z<zcomb14) return UCOMB;  // over U comb
    if (z<zcomb15) return VCOMB;  // over V comb
    if (z<zcomb16) 
      {
	if (y<ycomb1) return HORIZGAP; // below the bottom
	if (y<ycomb2) return UCOMB; // over U comb
	if (y<ycomb3) return VCOMB; // over V comb
	if (y<ycomb4) return ACTIVE; // active volume
	if (y<ycomb5) return VCOMB; // over V comb
	if (y<ycomb6) return UCOMB; // over U comb
	return HORIZGAP; // outside top edge
      }
    if (z<zcomb17) return VCOMB;  // over V comb
    if (z<zcomb18) return UCOMB;  // over U comb
    return VERTGAP; // off the end in Z.

  }

}


/* -------------------------------------------------
   APA Cartoons for the combtest35t method

   z->

   ^
   |
   y


   zcomb1                                       zcomb6
    zcomb2                                    zcomb5
     zcomb3                                  zcomb4
   ______________________________________________  ycomb6
   |____________________________________________|  ycomb5
   ||__________________________________________||  ycomb4
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                 APA0                   |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   ||__________________________________________||  ycomb3
   |____________________________________________|  ycomb2
   ______________________________________________  ycomb1


   z->

   ^
   |
   y


   zcomb7                                       zcomb12
    zcomb8                                    zcomb11
     zcomb9                                  zcomb10
   ______________________________________________  ycomb18
   |____________________________________________|  ycomb17
   ||__________________________________________||  ycomb16
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||               APA2                     |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   ||__________________________________________||  ycomb15
   |____________________________________________|  ycomb14
   ______________________________________________  ycomb13

   ______________________________________________  ycomb12
   |____________________________________________|  ycomb11
   ||__________________________________________||  ycomb10
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||              APA1                      |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   ||__________________________________________||  ycomb9
   |____________________________________________|  ycomb8
   ______________________________________________  ycomb7


   APA 0 Cartoon:

   z->

   ^
   |
   y


   zcomb13                                      zcomb18
    zcomb14                                   zcomb17
     zcomb15                                 zcomb16
   ______________________________________________  ycomb6
   |____________________________________________|  ycomb5
   ||__________________________________________||  ycomb4
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||         APA3                           |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   |||                                        |||
   ||__________________________________________||  ycomb3
   |____________________________________________|  ycomb2
   ______________________________________________  ycomb1



*/
  
