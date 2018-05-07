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

namespace sim {
  class ProtoDUNEbeamsim  {
    
  public:
    ProtoDUNEbeamsim(); //constructor
    //each particle should have a event ID, PDG ID, Momentum/Position @ BP4, TRIG2, LAG_ENTRY

    ~ProtoDUNEbeamsim(); //destructor
    
  private:

    //values I need to contain
    Float_t         BPROF4_x;
    Float_t         BPROF4_y;
    Float_t         BPROF4_z;
    Float_t         BPROF4_Px;
    Float_t         BPROF4_Py;
    Float_t         BPROF4_Pz;
    Float_t         BPROF4_PDGid;
    Float_t         BPROF4_EventID;
    Float_t         BPROF4_TrackID;
    Float_t         TRIG2_x;
    Float_t         TRIG2_y;
    Float_t         TRIG2_z;
    Float_t         TRIG2_Px;
    Float_t         TRIG2_Py;
    Float_t         TRIG2_Pz;
    Float_t         TRIG2_PDGid;
    Float_t         TRIG2_EventID;
    Float_t         TRIG2_TrackID;
    Float_t         Lag_ENTRY_x;
    Float_t         Lag_ENTRY_y;
    Float_t         Lag_ENTRY_z;
    Float_t         Lag_ENTRY_Px;
    Float_t         Lag_ENTRY_Py;
    Float_t         Lag_ENTRY_Pz;
    Float_t         Lag_ENTRY_PDGid;
    Float_t         Lag_ENTRY_EventID;
    Float_t         Lag_ENTRY_TrackID;
    
    
  public:

    ProtoDUNEbeamsim(    Float_t         BPROF4_x,
			 Float_t         BPROF4_y,
			 Float_t         BPROF4_z,
			 Float_t         BPROF4_Px,
			 Float_t         BPROF4_Py,
			 Float_t         BPROF4_Pz,
			 Float_t         BPROF4_PDGid,
			 Float_t         BPROF4_EventID,
			 Float_t         BPROF4_TrackID,
			 Float_t         TRIG2_x,
			 Float_t         TRIG2_y,
			 Float_t         TRIG2_z,
			 Float_t         TRIG2_Px,
			 Float_t         TRIG2_Py,
			 Float_t         TRIG2_Pz,
			 Float_t         TRIG2_EventID,
			 Float_t         TRIG2_TrackID,
			 Float_t         Lag_ENTRY_x,
			 Float_t         Lag_ENTRY_y,
			 Float_t         Lag_ENTRY_z,
			 Float_t         Lag_ENTRY_Px,
			 Float_t         Lag_ENTRY_Py,
			 Float_t         Lag_ENTRY_Pz,
			 Float_t         Lag_ENTRY_EventID,
		         Float_t         Lag_ENTRY_TrackID);

//inspectors 
		Float_t         get_BPROF4_x() const;
	    	Float_t         get_BPROF4_y() const;
	    	Float_t         get_BPROF4_z() const;
	    	Float_t         get_BPROF4_Px() const;
	    	Float_t         get_BPROF4_Py() const;
	    	Float_t         get_BPROF4_Pz() const;
	   	Float_t         get_BPROF4_PDGid() const;
	    	Float_t         get_BPROF4_EventID() const;
	    	Float_t         get_BPROF4_TrackID() const;
	    	Float_t         get_TRIG2_x() const;
	    	Float_t         get_TRIG2_y() const;
	    	Float_t         get_TRIG2_z() const;
	    	Float_t         get_TRIG2_Px() const;
	    	Float_t         get_TRIG2_Py() const;
	    	Float_t         get_TRIG2_Pz() const;
	    	Float_t         get_TRIG2_EventID() const;
	    	Float_t         get_TRIG2_TrackID() const;
	    	Float_t         get_Lag_ENTRY_x() const;
	    	Float_t         get_Lag_ENTRY_y() const;
	    	Float_t         get_Lag_ENTRY_z() const;
	    	Float_t         get_Lag_ENTRY_Px() const;
	    	Float_t         get_Lag_ENTRY_Py() const;
	    	Float_t         get_Lag_ENTRY_Pz() const;
	    	Float_t         get_Lag_ENTRY_EventID() const;
	    	Float_t         get_Lag_ENTRY_TrackID() const;
    
  //mutators 

		void SetBPROF4_x(Float_t val);
	    	void SetBPROF4_y(Float_t val);
	    	void SetBPROF4_z(Float_t val);
	    	void SetBPROF4_Px(Float_t val);
	    	void SetBPROF4_Py(Float_t val);
	    	void SetBPROF4_Pz(Float_t val);
	   	void SetBPROF4_PDGid(Float_t val);
	    	void SetBPROF4_EventID(Float_t val);
	    	void SetBPROF4_TrackID(Float_t val);
	    	void SetTRIG2_x(Float_t val);
	    	void SetTRIG2_y(Float_t val);
	    	void SetTRIG2_z(Float_t val);
	    	void SetTRIG2_Px(Float_t val);
	    	void SetTRIG2_Py(Float_t val);
	    	void SetTRIG2_Pz(Float_t val);
	    	void SetTRIG2_PDGid(Float_t val);
	    	void SetTRIG2_EventID(Float_t val);
	    	void SetTRIG2_TrackID(Float_t val);
	    	void SetLag_ENTRY_x(Float_t val);
	    	void SetLag_ENTRY_y(Float_t val);
	    	void SetLag_ENTRY_z(Float_t val);
	    	void SetLag_ENTRY_Px(Float_t val);
	    	void SetLag_ENTRY_Py(Float_t val);
	    	void SetLag_ENTRY_Pz(Float_t val);
	    	void SetLag_ENTRY_PDGid(Float_t val);
	    	void SetLag_ENTRY_EventID(Float_t val);
	    	void SetLag_ENTRY_TrackID(Float_t val);
  
    
  };
  
  
}

////////////////////////////////////////////////////////////////////////
#endif // PROTODUNESIM_H
