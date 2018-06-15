/////////////////////////////////////////////////////////
///////// ProtoDUNEbeamsim.cxx///////////////////////////
///////// Caroline Zhang carolineligezhang@gmail.com/////
///////// August 2017 ///////////////////////////////////
/////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include "dune/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEbeamsim.h"
#include <string.h>
#include <time.h>
#include <cmath>
#include <iomanip>

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace sim{

  //--------------------------constructors-----------------------------------------------
  ProtoDUNEbeamsim::ProtoDUNEbeamsim(){
/*     BPROFEXT_x=0;
     BPROFEXT_y=0;
     BPROFEXT_z=0;
     BPROFEXT_Px=0;
     BPROFEXT_Py=0;
     BPROFEXT_Pz=0;
     BPROFEXT_PDGid=0;
     BPROFEXT_EventID=0;
     BPROFEXT_TrackID=0;
     BPROF4_x=0;
     BPROF4_y=0;
     BPROF4_z=0;
     BPROF4_Px=0;
     BPROF4_Py=0;
     BPROF4_Pz=0;
     BPROF4_PDGid=0;
     BPROF4_EventID=0;
     BPROF4_TrackID=0;
     TOF1_t=0;
     TOF1_x=0;
     TOF1_y=0;
     TOF1_z=0;
     TOF1_Px=0;
     TOF1_Py=0;
     TOF1_Pz=0;
     TOF1_PDGid=0;
     TRIG2_x=0;
     TRIG2_y=0;
     TRIG2_z=0;
     TRIG2_Px=0;
     TRIG2_Py=0;
     TRIG2_Pz=0;
     TRIG2_PDGid=0;
     TRIG2_EventID=0;
     TRIG2_TrackID=0;
     NP04FieldCage_x=0;
     NP04FieldCage_y=0;
     NP04FieldCage_z=0;
     NP04FieldCage_Px=0;
     NP04FieldCage_Py=0;
     NP04FieldCage_Pz=0;
     NP04FieldCage_EventID=0;
     NP04FieldCage_TrackID=0; */
  }
  
  //-------------------------------default destructor------------------------------------------
  ProtoDUNEbeamsim::~ProtoDUNEbeamsim(){ }


  //------------------------alternate constructor-------------------------------------------------
/*  ProtoDUNEbeamsim::ProtoDUNEbeamsim(
    Float_t         vBPROFEXT_x,
    Float_t         vBPROFEXT_y,
    Float_t         vBPROFEXT_z,
    Float_t         vBPROFEXT_Px,
    Float_t         vBPROFEXT_Py,
    Float_t         vBPROFEXT_Pz,
    Float_t         vBPROFEXT_PDGid,
    Float_t         vBPROFEXT_EventID,
    Float_t         vBPROFEXT_TrackID,
    Float_t         vBPROF4_x,
    Float_t         vBPROF4_y,
    Float_t         vBPROF4_z,
    Float_t         vBPROF4_Px,
    Float_t         vBPROF4_Py,
    Float_t         vBPROF4_Pz,
    Float_t         vBPROF4_PDGid,
    Float_t         vBPROF4_EventID,
    Float_t         vBPROF4_TrackID,
    Float_t         vTOF1_t,
    Float_t         vTOF1_x,
    Float_t         vTOF1_y,
    Float_t         vTOF1_z,
    Float_t         vTOF1_Px,
    Float_t         vTOF1_Py,
    Float_t         vTOF1_Pz,
    Float_t         vTOF1_PDGid,
    Float_t         vTRIG2_x,
    Float_t         vTRIG2_y,
    Float_t         vTRIG2_z,
    Float_t         vTRIG2_Px,
    Float_t         vTRIG2_Py,
    Float_t         vTRIG2_Pz,
    Float_t         vTRIG2_EventID,
    Float_t         vTRIG2_TrackID,
    Float_t         vNP04FieldCage_x,
    Float_t         vNP04FieldCage_y,
    Float_t         vNP04FieldCage_z,
    Float_t         vNP04FieldCage_Px,
    Float_t         vNP04FieldCage_Py,
    Float_t         vNP04FieldCage_Pz,
    Float_t         vNP04FieldCage_EventID,
    Float_t         vNP04FieldCage_TrackID){
     BPROFEXT_x=vBPROFEXT_x;
     BPROFEXT_y=vBPROFEXT_y;
     BPROFEXT_z=vBPROFEXT_z;
     BPROFEXT_Px=vBPROFEXT_Px;
     BPROFEXT_Py=vBPROFEXT_Py;
     BPROFEXT_Pz=vBPROFEXT_Pz;
     BPROFEXT_PDGid=vBPROFEXT_PDGid;
     BPROFEXT_EventID=vBPROFEXT_EventID;
     BPROFEXT_TrackID=vBPROFEXT_TrackID;
     BPROF4_x=vBPROF4_x;
     BPROF4_y=vBPROF4_y;
     BPROF4_z=vBPROF4_z;
     BPROF4_Px=vBPROF4_Px;
     BPROF4_Py=vBPROF4_Py;
     BPROF4_Pz=vBPROF4_Pz;
     BPROF4_PDGid=vBPROF4_PDGid;
     BPROF4_EventID=vBPROF4_EventID;
     BPROF4_TrackID=vBPROF4_TrackID;
     TOF1_t=vTOF1_t;
     TOF1_x=vTOF1_x;
     TOF1_y=vTOF1_y;
     TOF1_z=vTOF1_z;
     TOF1_Px=vTOF1_Px;
     TOF1_Py=vTOF1_Py;
     TOF1_Pz=vTOF1_Pz;
     TOF1_PDGid=vTOF1_PDGid;
     TRIG2_x=vTRIG2_x;
     TRIG2_y=vTRIG2_y;
     TRIG2_z=vTRIG2_z;
     TRIG2_Px=vTRIG2_Px;
     TRIG2_Py=vTRIG2_Py;
     TRIG2_Pz=vTRIG2_Pz;
     TRIG2_EventID=vTRIG2_EventID;
     TRIG2_TrackID=vTRIG2_TrackID;
     NP04FieldCage_x=vNP04FieldCage_x;
     NP04FieldCage_y=vNP04FieldCage_y;
     NP04FieldCage_z=vNP04FieldCage_z;
     NP04FieldCage_Px=vNP04FieldCage_Px;
     NP04FieldCage_Py=vNP04FieldCage_Py;
     NP04FieldCage_Pz=vNP04FieldCage_Pz;
     NP04FieldCage_EventID=vNP04FieldCage_EventID;
     NP04FieldCage_TrackID=vNP04FieldCage_TrackID;
  }
*/
  // Leigh - I have encapsulated each beam instrument into a ProtoDUNEBeamInstrument object.
  // The below provides the interface for these.
  void ProtoDUNEbeamsim::AddInstrument(ProtoDUNEBeamInstrument newInst){

    bool alreadyExists = false;
    for(auto const &inst : fAllInstruments){
      if(newInst.GetInstrumentName() == inst.GetInstrumentName()){
        alreadyExists = true;
        break;
      }
    }

    if(!alreadyExists){
      fAllInstruments.push_back(newInst);
    }
    else{
      mf::LogError("ProtoDUNEbeamsim") << "Beam Instrument " << newInst.GetInstrumentName() << " already exists." << std::endl;
    }

  }

  ProtoDUNEBeamInstrument ProtoDUNEbeamsim::GetInstrument(std::string name) const{
    ProtoDUNEBeamInstrument temp;
    bool found = false;
    for(auto const inst : fAllInstruments){
      if(name == inst.GetInstrumentName()){
        temp = inst;
        found = true;
        break;
      }
    }

    if(!found){
      mf::LogWarning("ProtoDUNEbeamsim") << "Beam Instrument " << name << " not found, returning empty object." << std::endl;
    }

    return temp;
  }

}// namespace






