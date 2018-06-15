// ProtoDUNEbeamsim.h
//
// Caroline Zhang
// carolineligezhang@gmail.com
// August 2017
//

#ifndef PROTODUNEBEAMSIM_H
#define PROTODUNEBEAMSIM_H

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TDirectory.h"
#include "TStyle.h"
#include "TMath.h"

#include "dune/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEBeamInstrument.h"

namespace sim {
  class ProtoDUNEbeamsim  {
    
  public:
    ProtoDUNEbeamsim(); //constructor
    //each particle should have a event ID, PDG ID, Momentum/Position @ BP4, TRIG2, LAG_ENTRY

    ~ProtoDUNEbeamsim(); //destructor
    
  private:

    //values I need to contain
/*    Float_t         BPROFEXT_x;
    Float_t         BPROFEXT_y;
    Float_t         BPROFEXT_z;
    Float_t         BPROFEXT_Px;
    Float_t         BPROFEXT_Py;
    Float_t         BPROFEXT_Pz;
    Float_t         BPROFEXT_PDGid;
    Float_t         BPROFEXT_EventID;
    Float_t         BPROFEXT_TrackID;
    Float_t         BPROF4_x;
    Float_t         BPROF4_y;
    Float_t         BPROF4_z;
    Float_t         BPROF4_Px;
    Float_t         BPROF4_Py;
    Float_t         BPROF4_Pz;
    Float_t         BPROF4_PDGid;
    Float_t         BPROF4_EventID;
    Float_t         BPROF4_TrackID;
    Float_t         TOF1_t;
    Float_t         TOF1_x;
    Float_t         TOF1_y;
    Float_t         TOF1_z;
    Float_t         TOF1_Px;
    Float_t         TOF1_Py;
    Float_t         TOF1_Pz;
    Float_t         TOF1_PDGid;
    Float_t         TRIG2_x;
    Float_t         TRIG2_y;
    Float_t         TRIG2_z;
    Float_t         TRIG2_Px;
    Float_t         TRIG2_Py;
    Float_t         TRIG2_Pz;
    Float_t         TRIG2_PDGid;
    Float_t         TRIG2_EventID;
    Float_t         TRIG2_TrackID;
    Float_t         NP04FieldCage_x;
    Float_t         NP04FieldCage_y;
    Float_t         NP04FieldCage_z;
    Float_t         NP04FieldCage_Px;
    Float_t         NP04FieldCage_Py;
    Float_t         NP04FieldCage_Pz;
    Float_t         NP04FieldCage_PDGid;
    Float_t         NP04FieldCage_EventID;
    Float_t         NP04FieldCage_TrackID; */
   
    // Leigh: I think it would make much more sense to store a ProtoDUNEBeamInstrument object
    // for each of the instruments
    std::vector<ProtoDUNEBeamInstrument> fAllInstruments;

  public:

    /*ProtoDUNEbeamsim(  Float_t         BPROFEXT_x,
                         Float_t         BPROFEXT_y,
                         Float_t         BPROFEXT_z,
                         Float_t         BPROFEXT_Px,
                         Float_t         BPROFEXT_Py,
                         Float_t         BPROFEXT_Pz,
                         Float_t         BPROFEXT_PDGid,
                         Float_t         BPROFEXT_EventID,
                         Float_t         BPROFEXT_TrackID,
			 Float_t         BPROF4_x,
			 Float_t         BPROF4_y,
			 Float_t         BPROF4_z,
			 Float_t         BPROF4_Px,
			 Float_t         BPROF4_Py,
			 Float_t         BPROF4_Pz,
			 Float_t         BPROF4_PDGid,
			 Float_t         BPROF4_EventID,
			 Float_t         BPROF4_TrackID,
			 Float_t         TOF1_t,
			 Float_t         TOF1_x,
                         Float_t         TOF1_y,
                         Float_t         TOF1_z,
                         Float_t         TOF1_Px,
                         Float_t         TOF1_Py,
                         Float_t         TOF1_Pz,
                         Float_t         TOF1_PDGid,
			 Float_t         TRIG2_x,
			 Float_t         TRIG2_y,
			 Float_t         TRIG2_z,
			 Float_t         TRIG2_Px,
			 Float_t         TRIG2_Py,
			 Float_t         TRIG2_Pz,
			 Float_t         TRIG2_EventID,
			 Float_t         TRIG2_TrackID,
			 Float_t         NP04FieldCage_x,
                         Float_t         NP04FieldCage_y,
                         Float_t         NP04FieldCage_z,
                         Float_t         NP04FieldCage_Px,
                         Float_t         NP04FieldCage_Py,
                         Float_t         NP04FieldCage_Pz,
                         Float_t         NP04FieldCage_EventID,
                         Float_t         NP04FieldCage_TrackID);*/

    void AddInstrument(ProtoDUNEBeamInstrument newInst);
    ProtoDUNEBeamInstrument GetInstrument(std::string name) const;
    unsigned short NInstruments() {return fAllInstruments.size();};
    
  };
  
  
}

////////////////////////////////////////////////////////////////////////
#endif // PROTODUNESIM_H
