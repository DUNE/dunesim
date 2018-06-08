#include <iostream>
#include "dune/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEBeamInstrument.h"

namespace sim{

  //--------------------------constructors-----------------------------------------------
  ProtoDUNEBeamInstrument::ProtoDUNEBeamInstrument(){
    fInstrumentName = "default";
    fX = 0.0;
    fY = 0.0;
    fZ = 0.0;
    fT = 0.0;
    fPx = 0.0;
    fPy = 0.0;
    fPz = 0.0;
    fPDGid = 0;
    fEventID = 0;
    fTrackID = 0;
  }

  //-------------------------------default destructor------------------------------------------
  ProtoDUNEBeamInstrument::~ProtoDUNEBeamInstrument(){

  }

  ProtoDUNEBeamInstrument::ProtoDUNEBeamInstrument(std::string name,
      Float_t x,
      Float_t y,
      Float_t z,
      Float_t t,
      Float_t Px,
      Float_t Py,
      Float_t Pz,
      Int_t PDGid,
      Int_t EventID,
      Int_t TrackID
      ){
    fInstrumentName = name;
    fX = x; fY = y; fZ = z; fT = t;
    fPx = Px; fPy = Py; fPz = Pz;
    fPDGid = PDGid;
    fEventID = EventID;
    fTrackID = TrackID;
  }

  /// Vector-based constructor
  ProtoDUNEBeamInstrument::ProtoDUNEBeamInstrument(std::string name,
      std::vector<Float_t> position,
      Float_t t,
      std::vector<Float_t> momentum,
      Int_t PDGid,
      Int_t EventID,
      Int_t TrackID
      ){
    fInstrumentName = name;
    fX = position[0]; fY = position[1]; fZ = position[2]; fT = t;
    fPx = momentum[0]; fPy = momentum[1]; fPz = momentum[2];
    fPDGid = PDGid;
    fEventID = EventID;
    fTrackID = TrackID;
  }

  // TVector3-based constructor
  ProtoDUNEBeamInstrument::ProtoDUNEBeamInstrument(std::string name,
      TVector3 position,
      Float_t t,
      TVector3 momentum,
      Int_t PDGid,
      Int_t EventID,
      Int_t TrackID
      ){
    fInstrumentName = name;
    fX = position.X(); fY = position.Y(); fZ = position.Z(); fT = t;
    fPx = momentum.X(); fPy = momentum.Y(); fPz = momentum.Z();
    fPDGid = PDGid;
    fEventID = EventID;
    fTrackID = TrackID;

  }

  ProtoDUNEBeamInstrument::ProtoDUNEBeamInstrument(const ProtoDUNEBeamInstrument& rhs){
    fX = rhs.fX; fY = rhs.fY; fZ = rhs.fZ; fT = rhs.fT;
    fPx = rhs.fPx; fPy = rhs.fPy; fPz = rhs.fPz;
    fPDGid = rhs.fPDGid;
    fEventID = rhs.fEventID;
    fTrackID = rhs.fTrackID;
    fInstrumentName = rhs.fInstrumentName;
  }

}

