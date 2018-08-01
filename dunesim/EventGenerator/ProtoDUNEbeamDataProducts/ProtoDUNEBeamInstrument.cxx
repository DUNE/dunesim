//// Created and Modified by Pablo and Leigh H. Howard, 
////Smear important variables of beam monitors, mimic Cherenkov monitors response
////and store some histograms through art utilities
//// pablo.fer@cern.ch
//// July 2018
/////////////////////////////////////////////////////////////

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
    fSmearedVar1 = 0.0;
    fSmearedVar2 = 0.0;
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
    if (name == "TOF1" || name == "TRIG2"){
      Float_t delta_t = 5.;
      Float_t resolt = 1.0;
      srand (static_cast <unsigned> (time(0)));
      if (name == "TOF1") srand (static_cast <unsigned> (time(0)*9));
      Float_t t_test = t -delta_t + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(delta_t)));
      Float_t p_test = static_cast <float> (rand()) /( static_cast <float> (RAND_MAX));
      Float_t p_gauss = ProtoDUNEBeamInstrument::UnitGauss(t,t_test,resolt);
      while (p_test > p_gauss){
        p_test = static_cast <float> (rand()) /( static_cast <float> (RAND_MAX));
        t_test = t -delta_t + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(delta_t)));
        p_gauss = ProtoDUNEBeamInstrument::UnitGauss(t,t_test,resolt);
      }
      fSmearedVar1 = t_test;
  }
    if (name == "BPROFEXT" || name == "BPROF4"){
      Float_t delta_x = 10.;
      Float_t resolx = 2.0/pow(12,0.5);
      srand (static_cast <unsigned> (time(0)*99));
      Float_t x_test = x -delta_x + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(delta_x)));
      Float_t p_test = static_cast <float> (rand()) /( static_cast <float> (RAND_MAX));
      Float_t p_gauss = ProtoDUNEBeamInstrument::UnitGauss(x,x_test,resolx);
      while (p_test > p_gauss){
        p_test = static_cast <float> (rand()) /( static_cast <float> (RAND_MAX)); 
        x_test = x -delta_x + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(delta_x)));
        p_gauss = ProtoDUNEBeamInstrument::UnitGauss(x,x_test,resolx);
      }
      fSmearedVar1 = x_test;

      Float_t delta_y = 10.;
      Float_t resoly = 2.0/pow(12,0.5);
      srand (static_cast <unsigned> (time(0)*7));
      Float_t y_test = y -delta_y + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(delta_y)));
      p_test = static_cast <float> (rand()) /( static_cast <float> (RAND_MAX));
      p_gauss = ProtoDUNEBeamInstrument::UnitGauss(y,y_test,resoly);
      while (p_test > p_gauss){
        p_test = static_cast <float> (rand()) /( static_cast <float> (RAND_MAX));
        y_test = y -delta_y + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(delta_y)));
        p_gauss = ProtoDUNEBeamInstrument::UnitGauss(y,y_test,resoly);
      }
      fSmearedVar2 = y_test;
}
    if (name == "CHERENKOV1"){
      Float_t Ptot = pow(pow(Px,2)+pow(Py,2)+pow(Pz,2),0.5)/1000.;
      if (Ptot <= 2.0){
        if (abs(PDGid) == 2212 || abs(PDGid) == 310) fSmearedVar1 = 0;
        if (abs(PDGid) == 11) fSmearedVar1 = 1;
        if (abs(PDGid) == 211) fSmearedVar1 = 1;
}
      if (Ptot <= 3.0 && Ptot >2.0){
        if (abs(PDGid) == 2212 || abs(PDGid) == 310) fSmearedVar1 = 0;
        if (abs(PDGid) == 11) fSmearedVar1 = 1;
        if (abs(PDGid) == 211) fSmearedVar1 = 1;
}
      if (Ptot <= 5.0 && Ptot > 3.0){
        if (abs(PDGid) == 2212 || abs(PDGid) == 310) fSmearedVar1 = 0;
        if (abs(PDGid) == 11) fSmearedVar1 = 1;
        if (abs(PDGid) == 211) fSmearedVar1 = 1;
}
      if (Ptot > 5.0){
        if (abs(PDGid) == 2212 || abs(PDGid) == 310) fSmearedVar1 = 0;
        if (abs(PDGid) == 11) fSmearedVar1 = 0;
        if (abs(PDGid) == 211) fSmearedVar1 = 1;
}
}
    if (name == "CHERENKOV2"){
      Float_t Ptot = pow(pow(Px,2)+pow(Py,2)+pow(Pz,2),0.5)/1000.;
      if (Ptot <= 2.0){
        if (abs(PDGid) == 2212 || abs(PDGid) == 310) fSmearedVar1 = 0;
        if (abs(PDGid) == 11) fSmearedVar1 = 0;
        if (abs(PDGid) == 211) fSmearedVar1 = 1;
}
      if (Ptot <= 3.0 && Ptot >2.0){
        if (abs(PDGid) == 2212 || abs(PDGid) == 310) fSmearedVar1 = 0;
        if (abs(PDGid) == 11) fSmearedVar1 = 0;
        if (abs(PDGid) == 211) fSmearedVar1 = 1;
}
      if (Ptot <= 5.0 && Ptot > 3.0){
        if (abs(PDGid) == 2212 || abs(PDGid) == 310) fSmearedVar1 = 0;
        if (abs(PDGid) == 11) fSmearedVar1 = 0;
        if (abs(PDGid) == 211) fSmearedVar1 = 1;
}
      if (Ptot > 5.0){
        if (abs(PDGid) == 2212 || abs(PDGid) == 310) fSmearedVar1 = 0;
        if (abs(PDGid) == 11) fSmearedVar1 = 0;
        if (abs(PDGid) == 211) fSmearedVar1 = 1;
}
}
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
    fSmearedVar1 = rhs.fSmearedVar1;
    fSmearedVar2 = rhs.fSmearedVar2;
  }

}

