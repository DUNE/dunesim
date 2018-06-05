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
    Float_t         BPROFEXT_x;
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
    Float_t         NP04FieldCage_TrackID; 
   
    // Leigh: I thibk it would make much more sense to store a ProtoDUNEBeamInstrument object
    // for each of the instruments
    std::vector<ProtoDUNEBeamInstrument> fAllInstruments;

  public:

    ProtoDUNEbeamsim(    Float_t         BPROFEXT_x,
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
                         Float_t         NP04FieldCage_TrackID);

//inspectors 
		Float_t         get_BPROFEXT_x() const;
                Float_t         get_BPROFEXT_y() const;
                Float_t         get_BPROFEXT_z() const;
                Float_t         get_BPROFEXT_Px() const;
                Float_t         get_BPROFEXT_Py() const;
                Float_t         get_BPROFEXT_Pz() const;
                Float_t         get_BPROFEXT_PDGid() const;
                Float_t         get_BPROFEXT_EventID() const;
                Float_t         get_BPROFEXT_TrackID() const;
		Float_t         get_BPROF4_x() const;
	    	Float_t         get_BPROF4_y() const;
	    	Float_t         get_BPROF4_z() const;
	    	Float_t         get_BPROF4_Px() const;
	    	Float_t         get_BPROF4_Py() const;
	    	Float_t         get_BPROF4_Pz() const;
	   	Float_t         get_BPROF4_PDGid() const;
	    	Float_t         get_BPROF4_EventID() const;
	    	Float_t         get_BPROF4_TrackID() const;
                Float_t         get_TOF1_t() const;
                Float_t         get_TOF1_x() const;
                Float_t         get_TOF1_y() const;
                Float_t         get_TOF1_z() const;
                Float_t         get_TOF1_Px() const;
                Float_t         get_TOF1_Py() const;
                Float_t         get_TOF1_Pz() const;
                Float_t         get_TOF1_PDGid() const;
	    	Float_t         get_TRIG2_x() const;
	    	Float_t         get_TRIG2_y() const;
	    	Float_t         get_TRIG2_z() const;
	    	Float_t         get_TRIG2_Px() const;
	    	Float_t         get_TRIG2_Py() const;
	    	Float_t         get_TRIG2_Pz() const;
	    	Float_t         get_TRIG2_EventID() const;
	    	Float_t         get_TRIG2_TrackID() const;
                Float_t         get_NP04FieldCage_x() const;
                Float_t         get_NP04FieldCage_y() const;
                Float_t         get_NP04FieldCage_z() const;
                Float_t         get_NP04FieldCage_Px() const;
                Float_t         get_NP04FieldCage_Py() const;
                Float_t         get_NP04FieldCage_Pz() const;
                Float_t         get_NP04FieldCage_EventID() const;
                Float_t         get_NP04FieldCage_TrackID() const;
 
  //mutators 

		 void SetBPROFEXT_x(Float_t val);
                void SetBPROFEXT_y(Float_t val);
                void SetBPROFEXT_z(Float_t val);
                void SetBPROFEXT_Px(Float_t val);
                void SetBPROFEXT_Py(Float_t val);
                void SetBPROFEXT_Pz(Float_t val);
                void SetBPROFEXT_PDGid(Float_t val);
                void SetBPROFEXT_EventID(Float_t val);
                void SetBPROFEXT_TrackID(Float_t val);
		void SetBPROF4_x(Float_t val);
	    	void SetBPROF4_y(Float_t val);
	    	void SetBPROF4_z(Float_t val);
	    	void SetBPROF4_Px(Float_t val);
	    	void SetBPROF4_Py(Float_t val);
	    	void SetBPROF4_Pz(Float_t val);
	   	void SetBPROF4_PDGid(Float_t val);
	    	void SetBPROF4_EventID(Float_t val);
	    	void SetBPROF4_TrackID(Float_t val);
                void SetTOF1_t(Float_t val);
                void SetTOF1_x(Float_t val);
                void SetTOF1_y(Float_t val);
                void SetTOF1_z(Float_t val);
                void SetTOF1_Px(Float_t val);
                void SetTOF1_Py(Float_t val);
                void SetTOF1_Pz(Float_t val);
                void SetTOF1_PDGid(Float_t val);
	    	void SetTRIG2_x(Float_t val);
	    	void SetTRIG2_y(Float_t val);
	    	void SetTRIG2_z(Float_t val);
	    	void SetTRIG2_Px(Float_t val);
	    	void SetTRIG2_Py(Float_t val);
	    	void SetTRIG2_Pz(Float_t val);
	    	void SetTRIG2_PDGid(Float_t val);
	    	void SetTRIG2_EventID(Float_t val);
	    	void SetTRIG2_TrackID(Float_t val);
                void SetNP04FieldCage_x(Float_t val);
                void SetNP04FieldCage_y(Float_t val);
                void SetNP04FieldCage_z(Float_t val);
                void SetNP04FieldCage_Px(Float_t val);
                void SetNP04FieldCage_Py(Float_t val);
                void SetNP04FieldCage_Pz(Float_t val);
                void SetNP04FieldCage_PDGid(Float_t val);
                void SetNP04FieldCage_EventID(Float_t val);
                void SetNP04FieldCage_TrackID(Float_t val);
  

    void AddInstrument(ProtoDUNEBeamInstrument newInst);
    ProtoDUNEBeamInstrument GetInstrument(std::string name) const;
    
  };
  
  
}

////////////////////////////////////////////////////////////////////////
#endif // PROTODUNESIM_H
