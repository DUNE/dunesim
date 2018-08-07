// ProtoDUNEBeamInstrument.h
//
// Leigh Whitehead
// leigh.howard.whitehead@cern.ch
// June 2018
//
// This class contains all of the information from a given
// simulated beam line component

#ifndef PROTODUNEBEAMINSTRUMENT_H
#define PROTODUNEBEAMINSTRUMENT_H

#include <string>
#include <vector>
#include "TVector3.h"

namespace sim {
  class ProtoDUNEBeamInstrument  {

    public:
      ProtoDUNEBeamInstrument(); //constructor
      ~ProtoDUNEBeamInstrument(); //destructor

      ProtoDUNEBeamInstrument(std::string name,
        Float_t x,
        Float_t y,
        Float_t z,
        Float_t t,
        Float_t Px,
        Float_t Py,
        Float_t Pz,
        Int_t PDGid,
        Int_t EventID,
        Int_t TrackID,
        Float_t Resolution
      );

      /// Vector-based constructor
      ProtoDUNEBeamInstrument(std::string name, 
        std::vector<Float_t> position,
        Float_t t,
        std::vector<Float_t> momentum,
        Int_t PDGid,
        Int_t EventID,
        Int_t TrackID,
        Float_t Resolution
      );      

      // TVector3-based constructor
      ProtoDUNEBeamInstrument(std::string name, 
        TVector3 position,
        Float_t t,
        TVector3 momentum,
        Int_t PDGid,
        Int_t EventID,
        Int_t TrackID,
        Float_t Resolution
      );      

      // Copy constructor
      ProtoDUNEBeamInstrument(const ProtoDUNEBeamInstrument& rhs);

      std::string GetInstrumentName() const {return fInstrumentName;};
      Float_t GetX() const {return fX;};
      Float_t GetY() const {return fY;};
      Float_t GetZ() const {return fZ;};
      Float_t GetT() const {return fT;};
      Float_t GetPx() const {return fPx;};
      Float_t GetPy() const {return fPy;};
      Float_t GetPz() const {return fPz;};
      Int_t GetPDGid() const {return fPDGid;};
      Int_t GetEventID() const {return fEventID;}; 
      Int_t GetTrackID() const {return fTrackID;};
      Float_t GetSmearedVar1() const {return fSmearedVar1;};
      Float_t GetSmearedVar2() const {return fSmearedVar2;};
      Float_t GetResolution() const {return fResolution;};

      void SetInstrumentName(std::string name) {fInstrumentName = name;};
      void SetX(Float_t val) {fX = val;};
      void SetY(Float_t val) {fY = val;};
      void SetZ(Float_t val) {fZ = val;};
      void SetT(Float_t val) {fT = val;};
      void SetPx(Float_t val) {fPx = val;};
      void SetPy(Float_t val) {fPy = val;};
      void SetPz(Float_t val) {fPz = val;};
      void SetPDGid(Int_t val) {fPDGid = val;};
      void SetEventID(Int_t val) {fEventID = val;};
      void SetTrackID(Int_t val) {fTrackID = val;};
      void SetSmearedVar1(Float_t val) {fSmearedVar1 = val;};
      void SetSmearedVar2(Float_t val) {fSmearedVar2 = val;};
      void SetResolution(Float_t val) {fResolution = val;};


      Float_t UnitGauss(Float_t mean, Float_t value, Float_t sigma){
      Float_t p = exp(-0.5*pow((mean-value),2)/(pow(sigma,2)));
      return p;}



    private:

      std::string fInstrumentName;

      Float_t fX;
      Float_t fY;
      Float_t fZ;
      Float_t fT;
      Float_t fPx;
      Float_t fPy;
      Float_t fPz;
      Int_t fPDGid;
      Int_t fEventID;
      Int_t fTrackID;
      Float_t fSmearedVar1;
      Float_t fSmearedVar2;
      Float_t fResolution;
  };
}

////////////////////////////////////////////////////////////////////////
#endif // PROTODUNEBEAMINSTRUMENT_H
