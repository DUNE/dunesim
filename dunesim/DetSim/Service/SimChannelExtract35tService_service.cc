// SimChannelExtract35tService.cxx

#include "SimChannelExtract35tService.h"
#include <string>
#include "Geometry/Geometry.h"
#include "Simulation/SimChannel.h"

using std::string;

//**********************************************************************

SimChannelExtract35tService::
SimChannelExtract35tService(fhicl::ParameterSet const& pset, art::ActivityRegistry&)
: m_ntick(0),
  fFirstCollectionChannel(9999999) {
  fFractUUCollect         = pset.get< float >("FractUUCollect");
  fFractUVCollect         = pset.get< float >("FractUVCollect");
  fFractVUCollect         = pset.get< float >("FractVUCollect");
  fFractVVCollect         = pset.get< float >("FractVVCollect");
  fFractUUMiss            = pset.get< float >("FractUUMiss");
  fFractUVMiss            = pset.get< float >("FractUVMiss");
  fFractVUMiss            = pset.get< float >("FractVUMiss");
  fFractVVMiss            = pset.get< float >("FractVVMiss");
  fFractZUMiss            = pset.get< float >("FractZUMiss");
  fFractZVMiss            = pset.get< float >("FractZVMiss");
  fFractHorizGapUMiss     = pset.get< float >("FractHorizGapUMiss");
  fFractVertGapUMiss      = pset.get< float >("FractVertGapUMiss");
  fFractHorizGapVMiss     = pset.get< float >("FractHorizGapVMiss");
  fFractVertGapVMiss      = pset.get< float >("FractVertGapVMiss");
  fFractHorizGapZMiss     = pset.get< float >("FractHorizGapZMiss");
  fFractVertGapZMiss      = pset.get< float >("FractVertGapZMiss");
  fFractHorizGapUCollect  = pset.get< float >("FractHorizGapUCollect");
  fFractVertGapUCollect   = pset.get< float >("FractVertGapUCollect");
  fFractHorizGapVCollect  = pset.get< float >("FractHorizGapVCollect");
  fFractVertGapVCollect   = pset.get< float >("FractVertGapVCollect");
  init();
}

//**********************************************************************

int SimChannelExtract35tService::
extract(const sim::SimChannel* psc, AdcSignalVector& fChargeWork) const {
  fChargeWork.clear();
  fChargeWork.resize(m_ntick, 0.0);
  if ( psc == nullptr ) return 0;
  AdcSignalVector fChargeWorkCollInd(m_ntick, 0.0);
  string fname = "SimChannelExtract35tService::extract";
  art::ServiceHandle<geo::Geometry> geo;
  unsigned int chan = psc->Channel();
  const geo::View_t view = geo->View(chan);
  fChargeWorkCollInd.clear();
  fChargeWorkCollInd.resize(fChargeWork.size(), 0.);
  for ( size_t t=0; t<fChargeWork.size(); ++t ) {
    const std::vector<sim::IDE> ides = psc->TrackIDsAndEnergies(t,t);
    for ( auto const &ide : ides ) {
      GapType_t gaptype = combtest35t(ide.x,ide.y,ide.z);
      switch ( gaptype ) {
        case ACTIVE:
          fChargeWork[t] += ide.numElectrons;
          break;
        case UCOMB:
          switch (view) {
            case geo::kU:
              fChargeWork[t] += ide.numElectrons * (1.0-fFractUUCollect-fFractUUMiss);
              fChargeWorkCollInd[t] += ide.numElectrons * fFractUUCollect;
              break;
            case geo::kV:
              fChargeWork[t] += ide.numElectrons * (1.0-fFractVUCollect-fFractUUCollect-fFractVUMiss);
              fChargeWorkCollInd[t] += ide.numElectrons * fFractVUCollect;
              break;
            case geo::kZ:
              fChargeWork[t] += ide.numElectrons * (1.0-fFractVUCollect-fFractUUCollect-fFractZUMiss);
              break;
            default:
              throw cet::exception(fname) << "ILLEGAL VIEW Type: " << view <<"\n";
          }
          break;
        case VCOMB:
          switch (view) {
            case geo::kU:
              fChargeWork[t] += ide.numElectrons * (1.0-fFractUVCollect-fFractUVMiss);
              fChargeWorkCollInd[t] += ide.numElectrons * fFractUVCollect;
              break;
            case geo::kV:
              fChargeWork[t] += ide.numElectrons * (1.0-fFractUVCollect-fFractVVCollect-fFractVVMiss);
              fChargeWorkCollInd[t] += ide.numElectrons * fFractVVCollect;
              break;
            case geo::kZ:
              fChargeWork[t] += ide.numElectrons * (1.0-fFractVVCollect-fFractUVCollect-fFractZVMiss);
              break;
            default:
              throw cet::exception(fname) << "ILLEGAL VIEW Type: " << view <<"\n";
          }
          break;
        case HORIZGAP:
          switch (view) {
            case geo::kU:
              fChargeWork[t] += ide.numElectrons * (1.0-fFractHorizGapUMiss-fFractHorizGapUCollect);
              fChargeWorkCollInd[t] += ide.numElectrons * fFractHorizGapUCollect;
              break;
            case geo::kV:
              fChargeWork[t] += ide.numElectrons * (1.0-fFractHorizGapVMiss-fFractHorizGapUCollect-fFractHorizGapVCollect);
              fChargeWorkCollInd[t] += ide.numElectrons * fFractHorizGapVCollect;
              break;
            case geo::kZ:
              fChargeWork[t] += ide.numElectrons * (1.0-fFractHorizGapZMiss-fFractHorizGapUCollect-fFractHorizGapVCollect);
              break;
            default:
              throw cet::exception(fname) << "ILLEGAL VIEW Type: " << view <<"\n";
          }
          break;
        case VERTGAP:
          switch (view) {
            case geo::kU:
              fChargeWork[t] += ide.numElectrons * (1.0-fFractVertGapUMiss-fFractVertGapUCollect);
              fChargeWorkCollInd[t] += ide.numElectrons * fFractVertGapUCollect;
              break;
            case geo::kV:
              fChargeWork[t] += ide.numElectrons * (1.0-fFractVertGapVMiss-fFractVertGapUCollect-fFractVertGapVCollect);
              fChargeWorkCollInd[t] += ide.numElectrons * fFractVertGapVCollect;
              break;
            case geo::kZ:
              fChargeWork[t] += ide.numElectrons * (1.0-fFractVertGapZMiss-fFractVertGapUCollect-fFractVertGapVCollect);
              break;
            default:
              throw cet::exception(fname) << "ILLEGAL VIEW Type: " << view <<"\n";
          }
          break;
        case NONACTIVE:
          break;
      }  // end switch gaptype
    }  // end loop over ides
  }  // end loop over ticks
  // Convolute and combine charges.
  m_psss->Convolute(chan, fChargeWork);
  m_psss->Convolute(fFirstCollectionChannel, fChargeWorkCollInd);
  for ( unsigned int itck=0; itck<m_ntick; ++itck ) {
    fChargeWork[itck] += fChargeWorkCollInd[itck];
  }
  return 0;
}

//**********************************************************************

void SimChannelExtract35tService::init() {

  if ( m_init ) return;

  // Fetch the number of ticks.
  m_ntick = m_pfft->FFTSize();
  if ( m_ntick%2 != 0 )
    throw cet::exception("SimChannelExtractAllService")
          << "FFT size is not a power of 2.";
  art::ServiceHandle<geo::Geometry> geo;

  // Find the first collection channel.
  for ( uint32_t ichan=0; ichan<geo->Nchannels(); ++ichan ) {
    if ( geo->View(ichan) == geo::kZ ) {
      fFirstCollectionChannel = ichan;
      break;
    }
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

  // need to get zcomb13, zcomb18, ycomb1, and ycomb6 -- extrapolate
  // the ycomb1 and ycomb6 are just copies.

  zcomb13 = zcomb14 - (zcomb15 - zcomb14);
  zcomb18 = zcomb17 + (zcomb17 - zcomb16);
  ycomb1 = ycomb2 - (ycomb3 - ycomb2);
  ycomb6 = ycomb5 + (ycomb5 - ycomb4);

  m_init = true;

}

//**********************************************************************

// see the ASCII cartoon of APA's at the bottom of this file for a picture of what all the boundaries are

SimChannelExtract35tService::GapType_t 
SimChannelExtract35tService::combtest35t(double x, double y, double z) const {
  if (z<zcomb1) return VERTGAP;  // off to the side of the first APA -- kind of like being in a vertical gap
  if (z<zcomb2) return UCOMB;  // over U comb
  if (z<zcomb3) return VCOMB;  // over V comb
  if (z<zcomb4) {
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
  if (z<zcomb10) {
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
  if (z<zcomb16) {
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

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(SimChannelExtract35tService, SimChannelExtractService)

//**********************************************************************
