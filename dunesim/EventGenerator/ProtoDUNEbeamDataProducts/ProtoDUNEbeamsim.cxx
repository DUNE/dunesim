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
namespace sim{

  //--------------------------constructors-----------------------------------------------
  ProtoDUNEbeamsim::ProtoDUNEbeamsim(){
     BPROF4_x=0;
     BPROF4_y=0;
     BPROF4_z=0;
     BPROF4_Px=0;
     BPROF4_Py=0;
     BPROF4_Pz=0;
     BPROF4_PDGid=0;
     BPROF4_EventID=0;
     BPROF4_TrackID=0;
     TRIG2_x=0;
     TRIG2_y=0;
     TRIG2_z=0;
     TRIG2_Px=0;
     TRIG2_Py=0;
     TRIG2_Pz=0;
     TRIG2_PDGid=0;
     TRIG2_EventID=0;
     TRIG2_TrackID=0;
     Lag_ENTRY_x=0;
     Lag_ENTRY_y=0;
     Lag_ENTRY_z=0;
     Lag_ENTRY_Px=0;
     Lag_ENTRY_Py=0;
     Lag_ENTRY_Pz=0;
     Lag_ENTRY_EventID=0;
     Lag_ENTRY_TrackID=0;
   
  }
  
  //-------------------------------default destructor------------------------------------------
  ProtoDUNEbeamsim::~ProtoDUNEbeamsim(){ }


  //------------------------alternate constructor-------------------------------------------------
  ProtoDUNEbeamsim::ProtoDUNEbeamsim(
    Float_t         vBPROF4_x,
    Float_t         vBPROF4_y,
    Float_t         vBPROF4_z,
    Float_t         vBPROF4_Px,
    Float_t         vBPROF4_Py,
    Float_t         vBPROF4_Pz,
    Float_t         vBPROF4_PDGid,
    Float_t         vBPROF4_EventID,
    Float_t         vBPROF4_TrackID,
    Float_t         vTRIG2_x,
    Float_t         vTRIG2_y,
    Float_t         vTRIG2_z,
    Float_t         vTRIG2_Px,
    Float_t         vTRIG2_Py,
    Float_t         vTRIG2_Pz,
    Float_t         vTRIG2_EventID,
    Float_t         vTRIG2_TrackID,
    Float_t         vLag_ENTRY_x,
    Float_t         vLag_ENTRY_y,
    Float_t         vLag_ENTRY_z,
    Float_t         vLag_ENTRY_Px,
    Float_t         vLag_ENTRY_Py,
    Float_t         vLag_ENTRY_Pz,
    Float_t         vLag_ENTRY_EventID,
    Float_t         vLag_ENTRY_TrackID){
     BPROF4_x=vBPROF4_x;
     BPROF4_y=vBPROF4_y;
     BPROF4_z=vBPROF4_z;
     BPROF4_Px=vBPROF4_Px;
     BPROF4_Py=vBPROF4_Py;
     BPROF4_Pz=vBPROF4_Pz;
     BPROF4_PDGid=vBPROF4_PDGid;
     BPROF4_EventID=vBPROF4_EventID;
     BPROF4_TrackID=vBPROF4_TrackID;
     TRIG2_x=vTRIG2_x;
     TRIG2_y=vTRIG2_y;
     TRIG2_z=vTRIG2_z;
     TRIG2_Px=vTRIG2_Px;
     TRIG2_Py=vTRIG2_Py;
     TRIG2_Pz=vTRIG2_Pz;
     TRIG2_EventID=vTRIG2_EventID;
     TRIG2_TrackID=vTRIG2_TrackID;
     Lag_ENTRY_x=vLag_ENTRY_x;
     Lag_ENTRY_y=vLag_ENTRY_y;
     Lag_ENTRY_z=vLag_ENTRY_z;
     Lag_ENTRY_Px=vLag_ENTRY_Px;
     Lag_ENTRY_Py=vLag_ENTRY_Py;
     Lag_ENTRY_Pz=vLag_ENTRY_Pz;
     Lag_ENTRY_EventID=vLag_ENTRY_EventID;
     Lag_ENTRY_TrackID=vLag_ENTRY_TrackID;


  }
  
  
  
  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_BPROF4_x() const
  {
    return BPROF4_x;
  }
  
  //-------------------------------------------------------------------------
 
  Float_t ProtoDUNEbeamsim::get_BPROF4_y() const
  {
    return BPROF4_y;
  }

  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_BPROF4_z() const
  {
    return BPROF4_z;
  }
  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_BPROF4_Px() const
  {
    return BPROF4_Px;
  }
    //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_BPROF4_Py() const
  {
    return BPROF4_Py;
  }
  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_BPROF4_Pz() const
  {
    return BPROF4_Pz;
  }

  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_BPROF4_PDGid() const
  {
    return BPROF4_PDGid;
  }

  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_BPROF4_EventID() const
  {
    return BPROF4_EventID;
  }
  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_BPROF4_TrackID() const
  {
    return BPROF4_TrackID;
  }
  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_TRIG2_x() const
  {
    return TRIG2_x;
  }
  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_TRIG2_y() const
  {
    return TRIG2_y;
  }

  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_TRIG2_z() const
  {
    return TRIG2_z;
  }

  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_TRIG2_Px() const
  {
    return TRIG2_Px;
  }

  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_TRIG2_Py() const
  {
    return TRIG2_Py;
  }

  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_TRIG2_Pz() const
  {
    return TRIG2_Pz;
  }

  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_TRIG2_EventID() const
  {
    return TRIG2_EventID;
  }
  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_TRIG2_TrackID() const
  {
    return TRIG2_TrackID;
  }
  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_Lag_ENTRY_x() const
  {
    return Lag_ENTRY_x;
  }
  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_Lag_ENTRY_y() const
  {
    return Lag_ENTRY_y;
  }
  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_Lag_ENTRY_z() const
  {
    return Lag_ENTRY_z;
  }
  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_Lag_ENTRY_Px() const
  {
    return Lag_ENTRY_Px;
  }
  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_Lag_ENTRY_Py() const
  {
    return Lag_ENTRY_Py;
  }
  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_Lag_ENTRY_Pz() const
  {
    return Lag_ENTRY_Pz;
  }

  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_Lag_ENTRY_EventID() const
  {
    return Lag_ENTRY_EventID;
  }
  //-------------------------------------------------------------------------
  Float_t ProtoDUNEbeamsim::get_Lag_ENTRY_TrackID() const
  {
    return Lag_ENTRY_TrackID;
  }
 //-----------MUTATORS-------------------//

  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetBPROF4_x(Float_t val)
  { 
    BPROF4_x= val;
  }

  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetBPROF4_y(Float_t val)
  { 
    BPROF4_y= val;
  }

  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetBPROF4_z(Float_t val)
  { 
    BPROF4_z= val;
  }


  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetBPROF4_Px(Float_t val)
  { 
    BPROF4_Px= val;
  }
  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetBPROF4_Py(Float_t val)
  { 
    BPROF4_Py= val;
  }
  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetBPROF4_Pz(Float_t val)
  { 
    BPROF4_Pz= val;
  }

  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetBPROF4_PDGid(Float_t val)
  { 
    BPROF4_PDGid= val;
  }

  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetBPROF4_EventID(Float_t val)
  { 
    BPROF4_EventID= val;
  }
  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetBPROF4_TrackID(Float_t val)
  { 
    BPROF4_TrackID= val;
  }

  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetTRIG2_x(Float_t val)
  { 
    TRIG2_x= val;
  }

  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetTRIG2_y(Float_t val)
  { 
    TRIG2_y= val;
  }

  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetTRIG2_z(Float_t val)
  { 
    TRIG2_z= val;
  }


  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetTRIG2_Px(Float_t val)
  { 
    TRIG2_Px= val;
  }
  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetTRIG2_Py(Float_t val)
  { 
    TRIG2_Py= val;
  }
  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetTRIG2_Pz(Float_t val)
  { 
    TRIG2_Pz= val;
  }

  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetTRIG2_PDGid(Float_t val)
  { 
    TRIG2_PDGid= val;
  }

  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetTRIG2_EventID(Float_t val)
  { 
    TRIG2_EventID= val;
  }
  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetTRIG2_TrackID(Float_t val)
  { 
    TRIG2_TrackID= val;
  }

// SECTION FOR Lag_ENTRY

  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetLag_ENTRY_x(Float_t val)
  { 
    Lag_ENTRY_x= val;
  }

  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetLag_ENTRY_y(Float_t val)
  { 
    Lag_ENTRY_y= val;
  }

  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetLag_ENTRY_z(Float_t val)
  { 
    Lag_ENTRY_z= val;
  }


  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetLag_ENTRY_Px(Float_t val)
  { 
    Lag_ENTRY_Px= val;
  }
  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetLag_ENTRY_Py(Float_t val)
  { 
    Lag_ENTRY_Py= val;
  }
  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetLag_ENTRY_Pz(Float_t val)
  { 
    Lag_ENTRY_Pz= val;
  }

  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetLag_ENTRY_PDGid(Float_t val)
  { 
    Lag_ENTRY_PDGid= val;
  }

  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetLag_ENTRY_EventID(Float_t val)
  { 
    Lag_ENTRY_EventID= val;
  }
  //-------------------------------------------------------------------------
  void  ProtoDUNEbeamsim::SetLag_ENTRY_TrackID(Float_t val)
  { 
    Lag_ENTRY_TrackID= val;
  }

}// namespace






