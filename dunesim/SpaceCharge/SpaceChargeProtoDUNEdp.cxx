////////////////////////////////////////////////////////////////////////
// \file SpaceChargeProtoDUNE-DP.cxx
//
// \brief implementation of class for storing/accessing space charge distortions for ProtoDUNE-DP
//
// based on\author mrmooney@colostate.edu
// \author jdawson@in2p3.fr
//  \author zambelli@lapp.in2p3.fr
////////////////////////////////////////////////////////////////////////
// C++ language includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include "stdio.h"
// LArSoft includes
#include "dunesim/SpaceCharge/SpaceChargeProtoDUNEdp.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#
// Framework includes
#include "cetlib_except/exception.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
// ROOT includes
#include "TFile.h"
#include "TH3.h"
#include "TTree.h"
#include "TLeaf.h"
//-----------------------------------------------
bool isWithinHist(geo::Point_t const& point, TH3F* hist);
bool isWithinHistOuter(geo::Point_t const& point, TH3F* hist);

double ExtrapolateAtEdge(TH3F *hist, double x, double y, double z);
bool IsBorder(TAxis *ax, double v);
double GetSymmetric(double x, double xc);
//-----------------------------------------------
spacecharge::SpaceChargeProtoDUNEdp::SpaceChargeProtoDUNEdp(
  fhicl::ParameterSet const& pset
)
{
  driftcoordinate=0;
  //Configure(pset);
}
//------------------------------------------------
bool spacecharge::SpaceChargeProtoDUNEdp::Configure(fhicl::ParameterSet const& pset,
                                                    detinfo::DetectorPropertiesData const& detProp)
{

  fEnableSimSpatialSCE = pset.get<bool>("EnableSimSpatialSCE");
  fEnableSimEfieldSCE = pset.get<bool>("EnableSimEfieldSCE");
  fEnableCalSpatialSCE = pset.get<bool>("EnableCalSpatialSCE");
  fEnableCalEfieldSCE = pset.get<bool>("EnableCalEfieldSCE");

  fEfield = detProp.Efield();
  std::cout<<"Efield : "<<fEfield<<std::endl;

  art::ServiceHandle<geo::Geometry> geom;
  //  auto const* geom = lar::providerFrom<geo::GeometryCore>();
  driftcoordinate = geom->TPC().DetectDriftDirection();
  if( driftcoordinate==1 || driftcoordinate==2 )
   {
        std::cout<<" drift coordinate: "<<driftcoordinate<<std::endl;
   }else{ throw cet::exception("CRTGen") << "unknown drift coordinate "<< driftcoordinate << " \n"; }


     /*    +1: positive x
       +2: positive y
       +3: positive z
       -1: negative x
       -2: negative y
       -3: negative z
       0: other (or algorithm failed)
     */


  if((fEnableSimSpatialSCE == true) || (fEnableSimEfieldSCE == true))
  {

    //x,y, z (drift is x)

    fInputFilename = pset.get<std::string>("InputFilename");


    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fInputFilename,fname);
    std::cout<<"spacecharge dualphase protodune: Find file: "<<fInputFilename<<" "<<fname<<std::endl;

    std::unique_ptr<TFile> infile(new TFile(fname.c_str(), "READ"));
    if(!infile->IsOpen()) throw cet::exception("SpaceChargeProtoDUNEdp") << "Could not find the space charge effect file '" << fname << "'!\n";
    // units are cm and V/cm
    //   if (fRepresentationType == "Voxelized_TH3") {
    //set up for driftcoordinate=1, driftcoordinate=2 is y-axis drift
    // drift coordinate uses -> hDeltaLength
    TH3F* hDx_sim_orig;
    TH3F* hDy_sim_orig;
    TH3F* hEndPoint_sim_orig; //was X
    if (driftcoordinate==1){
      hDx_sim_orig= (TH3F*)infile->Get("hDeltaLength");
      hDy_sim_orig = (TH3F*)infile->Get("hEndPointDY");
      hEndPoint_sim_orig = (TH3F*)infile->Get("hEndPointX");
      Anodebin = hEndPoint_sim_orig->GetXaxis()->GetXmax();

    }else if (driftcoordinate==2){
      hDx_sim_orig= (TH3F*)infile->Get("hEndPointDX");
      hDy_sim_orig = (TH3F*)infile->Get("hDeltaLength");
      hEndPoint_sim_orig = (TH3F*)infile->Get("hEndPointY");
      Anodebin = hEndPoint_sim_orig->GetYaxis()->GetXmax();//?
    }else{
      throw cet::exception("SpaceChargeProtoDUNEdp") << "Driftcoordinate "<<driftcoordinate<<" unknown\n";
    }

    TH3F* hDz_sim_orig = (TH3F*)infile->Get("hEndPointDZ");
    TH3F* hEx_sim_orig = (TH3F*)infile->Get("hEx");
    TH3F* hEy_sim_orig = (TH3F*)infile->Get("hEy");
    TH3F* hEz_sim_orig = (TH3F*)infile->Get("hEz");



    TH3F* hEndPointDrift_sim = (TH3F*) hEndPoint_sim_orig->Clone("hEndPointDrift");

    TH3F* hDx_sim = (TH3F*)hDx_sim_orig->Clone("hDx_sim");
    TH3F* hDy_sim = (TH3F*)hDy_sim_orig->Clone("hDy_sim");
    TH3F* hDz_sim = (TH3F*)hDz_sim_orig->Clone("hDz_sim");
    TH3F* hEx_sim = (TH3F*)hEx_sim_orig->Clone("hEx_sim");
    TH3F* hEy_sim = (TH3F*)hEy_sim_orig->Clone("hEy_sim");
    TH3F* hEz_sim = (TH3F*)hEz_sim_orig->Clone("hEz_sim");

    hEndPointDrift_sim->SetDirectory(0);

    hDx_sim->SetDirectory(0);
    hDy_sim->SetDirectory(0);
    hDz_sim->SetDirectory(0);
    hEx_sim->SetDirectory(0);
    hEy_sim->SetDirectory(0);
    hEz_sim->SetDirectory(0);



    SCEhistograms = {hDx_sim, hDy_sim, hDz_sim, hEx_sim, hEy_sim, hEz_sim, hEndPointDrift_sim};

      //   } //other representations not included yet..

  }
  fEnableCalEfieldSCE = false;
  fEnableCalSpatialSCE = false;

   if((fEnableCalSpatialSCE == true) || (fEnableCalEfieldSCE == true))
  {


    fCalInputFilename = pset.get<std::string>("CalibrationInputFilename");

    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fCalInputFilename,fname);
    std::cout<<"spacecharge dualphase protodune: Find file: "<<fCalInputFilename<<" "<<fname<<std::endl;

    std::unique_ptr<TFile> calinfile(new TFile(fname.c_str(), "READ"));
    if(!calinfile->IsOpen()) throw cet::exception("SpaceChargeProtoDUNEdp") << "Could not find the space charge effect calibration file '" << fname << "'!\n";

    TH3F* hDx_cal_orig = (TH3F*)calinfile->Get("hoffsetX");
    TH3F* hDy_cal_orig = (TH3F*)calinfile->Get("hoffsetY");
    TH3F* hDz_cal_orig = (TH3F*)calinfile->Get("hoffsetZ");

    TH3F* hEx_cal_orig = (TH3F*)calinfile->Get("hEx");
    TH3F* hEy_cal_orig = (TH3F*)calinfile->Get("hEy");
    TH3F* hEz_cal_orig = (TH3F*)calinfile->Get("hEz");



    TH3F* hDx_cal = (TH3F*)hDx_cal_orig->Clone("hDx_cal");
    TH3F* hDy_cal = (TH3F*)hDy_cal_orig->Clone("hDy_cal");
    TH3F* hDz_cal = (TH3F*)hDz_cal_orig->Clone("hDz_cal");
    TH3F* hEx_cal = (TH3F*)hEx_cal_orig->Clone("hEx_cal");
    TH3F* hEy_cal = (TH3F*)hEy_cal_orig->Clone("hEy_cal");
    TH3F* hEz_cal = (TH3F*)hEz_cal_orig->Clone("hEz_cal");

    hDx_cal->SetDirectory(0);
    hDy_cal->SetDirectory(0);
    hDz_cal->SetDirectory(0);
    hEx_cal->SetDirectory(0);
    hEy_cal->SetDirectory(0);
    hEz_cal->SetDirectory(0);



    CalSCEhistograms = {hDx_cal, hDy_cal, hDz_cal, hEx_cal, hEy_cal, hEz_cal};





   }


  return true;
}
//------------------------------------------------
bool spacecharge::SpaceChargeProtoDUNEdp::Update(uint64_t ts)
{
  if (ts == 0) return false;
  return true;
}
//----------------------------------------------------------------------------
/// Return boolean indicating whether or not to turn simulation of SCE on for
/// spatial distortions
bool spacecharge::SpaceChargeProtoDUNEdp::EnableSimSpatialSCE() const
{
  return fEnableSimSpatialSCE;
}
//----------------------------------------------------------------------------
/// Return boolean indicating whether or not to turn simulation of SCE on for
/// E field distortions
bool spacecharge::SpaceChargeProtoDUNEdp::EnableSimEfieldSCE() const
{
  return fEnableSimEfieldSCE;
}
//----------------------------------------------------------------------------
/// Return boolean indicating whether or not to apply SCE corrections
//bool spacecharge::SpaceChargeProtoDUNE::EnableCorrSCE() const
//{
//  return fEnableCorrSCE;
//}

/// Return boolean indicating whether or not to apply SCE corrections
bool spacecharge::SpaceChargeProtoDUNEdp::EnableCalSpatialSCE() const
{
  return fEnableCalSpatialSCE;
}

/// Return boolean indicating whether or not to apply SCE corrections
bool spacecharge::SpaceChargeProtoDUNEdp::EnableCalEfieldSCE() const
{
  return fEnableCalEfieldSCE;
}
//----------------------------------------------------------------------------
/// Primary working method of service that provides position offsets to be
/// used in ionization electron drift
geo::Vector_t spacecharge::SpaceChargeProtoDUNEdp::GetPosOffsets(geo::Point_t const& tmp_point) const
{
  //x [-300,300], y [-300,300], z[0,600]
  std::vector<double> thePosOffsets;
  geo::Point_t point = {tmp_point.X(), tmp_point.Y(), tmp_point.Z()};


  //  if (fRepresentationType=="Voxelized_TH3"){

        thePosOffsets = GetOffsetsVoxel(point, SCEhistograms.at(0), SCEhistograms.at(1), SCEhistograms.at(2));
        //drift offset = bended trajectory  - straight, LArVoxelReadout requires (straight-bended), so flip


        //}  else thePosOffsets.resize(3,0.0);

        //CHECK EndPointX (does this point get to the anode?), if not make offsets so large point goes outside volume...
        //	std::cout<<"pos: "<< tmp_point.X()<<" "<< tmp_point.Y()<<std::endl;

        if(isWithinHist(point, SCEhistograms.at(6)))
          {
            // std::cout<<"within..?"<<SCEhistograms.at(6)->Interpolate(point.X(), point.Y(), point.Z())<<" "<<SCEhistograms.at(6)->GetXaxis()->GetXmax()<<std::endl;
            if (SCEhistograms.at(6)->Interpolate(point.X(), point.Y(), point.Z()) >= Anodebin )
              {

                //	              std::cout<<"hits anode: "<< tmp_point.X()<<" "<< tmp_point.Y()<<" "<< tmp_point.Z()<<" "<<thePosOffsets[0]<<" "<< thePosOffsets[1]<<" "<< thePosOffsets[2]<<std::endl;
                if(driftcoordinate==1){
                  return {-thePosOffsets[0], thePosOffsets[1], thePosOffsets[2] };
                }else if (driftcoordinate==2)
                   { //is Voxel ready for y-shift?
                     return {thePosOffsets[0], -thePosOffsets[1], thePosOffsets[2] };
                   }else{
                  throw cet::exception("SpaceChargeProtoDUNEdp") << "Driftcoordinate "<<driftcoordinate<<" unknown\n";
                }
            }
          }

          return {9999999,9999999,9999999 };

          }

//----------------------------------------------------------------------------
/// Primary working method of service that provides position offsets to be
/// used in calibration of space charge
geo::Vector_t spacecharge::SpaceChargeProtoDUNEdp::GetCalPosOffsets(geo::Point_t const& tmp_point, int const& TPCid) const
{

  std::vector<double> thePosOffsets;
   geo::Point_t point = tmp_point;


  //  if (fRepresentationType == "Voxelized_TH3"){

      thePosOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(0), CalSCEhistograms.at(1), CalSCEhistograms.at(2));

      // } else thePosOffsets.resize(3,0.0);

  return { thePosOffsets[0], thePosOffsets[1], thePosOffsets[2] };
}
//Laura Zambelli's code from 3x1x1 sim
bool IsBorder(TAxis *ax, double v) {

       /*
              In the external half of the first/last bin or outside the range
               */
    if( ax->GetNbins()==1) return true;
    if( v<= ax->GetBinCenter(1) ) return true;
    //a->GetNbins() or GetNbins()-1..?
    if( v>= ax->GetBinCenter(ax->GetNbins()) ) return true;

    return false;
  }
double GetSymmetric(double x, double xc){
     return xc + (xc-x);
   }

double ExtrapolateAtEdge(TH3F *hist, double x, double y, double z){
  /*
     Want to get the value at x, close to the wall
    |----+----+--+--+---->
         x    c  a  b
     c is the center of the voxel (interpolation limit)
     b is the symetric of x wrt to c
     a is between c and b
     ->make linear interpolation from B, A to X

    */

     int binx = hist->GetXaxis()->FindBin(x);
     int biny = hist->GetYaxis()->FindBin(y);
     int binz = hist->GetZaxis()->FindBin(z);


     //get coordinate of the center of the voxel
     double xcenter = IsBorder(hist->GetXaxis(),x) ? hist->GetXaxis()->GetBinCenter(binx) : x;
     double ycenter = IsBorder(hist->GetYaxis(),y) ? hist->GetYaxis()->GetBinCenter(biny) : y;
     double zcenter = IsBorder(hist->GetZaxis(),z) ? hist->GetZaxis()->GetBinCenter(binz) : z;

     int Nx = hist->GetXaxis()->GetNbins();
     int Ny = hist->GetYaxis()->GetNbins();
     int Nz = hist->GetZaxis()->GetNbins();

     /*     if (Nz==1){ // when the GAr map is used (only one Z bin at the moment)
        return InterpolateAtEdge(x,y,z);// LZ
        }*/

      //special case when the starting point is at the center of the voxel
     if(x==xcenter && IsBorder(hist->GetXaxis(),x)){
         int binx2 = (binx==Nx) ? binx-1 : 2;
         xcenter =  0.5*(hist->GetXaxis()->GetBinCenter(binx2)+x);
        }
     if(y==ycenter&& IsBorder(hist->GetYaxis(),y)){
       int biny2 = (biny==Ny) ? biny-1 : 2;
       ycenter =  0.5*(hist->GetYaxis()->GetBinCenter(biny2)+y);
     }
     if(z==zcenter&& IsBorder(hist->GetZaxis(), z)){
            int binz2 = (binz==Nz) ? binz-1 : 2;
         zcenter =  0.5*(hist->GetZaxis()->GetBinCenter(binz2)+z);
       }


      //Get the symmetric of edgy point wrt to center
     double xb = IsBorder(hist->GetXaxis(),x) ? GetSymmetric(x,xcenter) : x;
     double yb = IsBorder(hist->GetYaxis(),y) ? GetSymmetric(y,ycenter) : y;
     double zb = IsBorder(hist->GetZaxis(),z) ? GetSymmetric(z,zcenter) : z;

       //Get middle point of [symmetric, center]
     double xa = IsBorder(hist->GetXaxis(),x) ?  0.5*(xcenter+xb) : x;
     double ya = IsBorder(hist->GetYaxis(),y) ?  0.5*(ycenter+yb) : y;
     double za = IsBorder(hist->GetZaxis(),z) ?  0.5*(zcenter+zb) : z;

     double valb = hist->Interpolate(xb, yb, zb); //symmetric point value
     double vala =  hist->Interpolate(xa, ya, za);
     double val = 0.;

     double Nedges = 0.;
     if(IsBorder(hist->GetXaxis(),x)) Nedges = Nedges+1.;
     if(IsBorder(hist->GetYaxis(),y)) Nedges = Nedges+1.;
     if(IsBorder(hist->GetZaxis(),z)) Nedges = Nedges+1.;



     double delta = vala - valb;//

     double dx = IsBorder(hist->GetXaxis(),x) ? xa-xb : 1;
     double dy = IsBorder(hist->GetYaxis(),y) ? ya-yb : 1;
     double dz = IsBorder(hist->GetZaxis(),z) ? za-zb : 1;

     val = x*delta/dx + (xa*valb - xb*vala)/dx
         +  y*delta/dy + (ya*valb - yb*vala)/dy
         +  z*delta/dz + (za*valb - zb*vala)/dz;


     val = val/Nedges;
     return val;

  }
//----------------------------------------------------------------------------
bool isWithinHistOuter(geo::Point_t const& point, TH3F* hist)
{
  if(point.X()> hist->GetXaxis()->GetBinLowEdge(0) &&
     point.X()< hist->GetXaxis()->GetBinLowEdge(hist->GetXaxis()->GetNbins()+1) &&
     point.Y()> hist->GetYaxis()->GetBinLowEdge(0) &&
     point.Y()< hist->GetYaxis()->GetBinLowEdge(hist->GetYaxis()->GetNbins()+1) &&
     point.Z()> hist->GetZaxis()->GetBinLowEdge(0) &&
     point.Z()< hist->GetZaxis()->GetBinLowEdge(hist->GetZaxis()->GetNbins()+1)  ){

    return true;}
  return false;

}

bool isWithinHist(geo::Point_t const& point, TH3F* hist)
{
  if(point.X()> hist->GetXaxis()->GetBinCenter(1)  &&
     point.X()< hist->GetXaxis()->GetBinCenter(hist->GetXaxis()->GetNbins()) &&
     point.Y()> hist->GetYaxis()->GetBinCenter(1) &&
     point.Y()< hist->GetYaxis()->GetBinCenter(hist->GetYaxis()->GetNbins()) &&
     point.Z()> hist->GetZaxis()->GetBinCenter(1) &&
     point.Z()< hist->GetZaxis()->GetBinCenter(hist->GetZaxis()->GetNbins()) ){

    return true;}
  return false;

}
/// Provides position offsets using voxelized interpolation
std::vector<double> spacecharge::SpaceChargeProtoDUNEdp::GetOffsetsVoxel
  (geo::Point_t const& point, TH3F* hX, TH3F* hY, TH3F* hZ) const
{

  if (isWithinHist(point, hX) && isWithinHist(point, hY) && isWithinHist(point, hZ)) {

    //need to check, X Y and Z on all three histograms..

    return {
        hX->Interpolate(point.X(),point.Y(),point.Z()),
        hY->Interpolate(point.X(),point.Y(),point.Z()),
        hZ->Interpolate(point.X(),point.Y(),point.Z())
    };
  }
  if (isWithinHistOuter(point, hX) && isWithinHistOuter(point, hY) && isWithinHistOuter(point, hZ)) {
  //use Laura's extrapolation method
    //    std::cout<<"extrap: ";

  return {ExtrapolateAtEdge(hX, point.X(),point.Y(),point.Z()),
          ExtrapolateAtEdge(hY,point.X(),point.Y(),point.Z()),
         ExtrapolateAtEdge(hZ,point.X(),point.Y(),point.Z())};
  }
  //  std::cout<<"outside: ";
      return {0.0,0.0,0.0}; //no offset and no efield!



}


//----------------------------------------------------------------------------

/// Primary working method of service that provides E field offsets to be
/// used in charge/light yield calculation (e.g.)
geo::Vector_t spacecharge::SpaceChargeProtoDUNEdp::GetEfieldOffsets(geo::Point_t const& tmp_point) const
{

  std::vector<double> theEfieldOffsets;
  geo::Point_t point = {tmp_point.X(), tmp_point.Y(), tmp_point.Z()};



       theEfieldOffsets = GetOffsetsVoxel(point, SCEhistograms.at(3), SCEhistograms.at(4), SCEhistograms.at(5));
       //THESE are the absoluted values of the Efield
       //ISCalculation expects difference from expected field
       /*
 geo::Point_t midPoint
      { ( step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->Ge\
tPosition() ) * 0.5/CLHEP::cm };
    auto EfieldDelta = fEfield * SCE->GetEfieldOffsets(midPoint);
    geo::Vector_t EfieldVec
      = { fEfield + EfieldDelta.X(), EfieldDelta.Y(), EfieldDelta.Z() };
      std::cout<<"ISCalculation : "<<EfieldVec.R()<<std::endl;    return EfieldVe\
c.R();
        */
       //     fEfield = detprop->Efield(); kV/cm..
       // Efield in larsoft is in kV/cm
       //std::cout<<"SpaceCharge efield: "<<pow( pow(theEfieldOffsets[0],2) + pow(theEfieldOffsets[1],2) +pow(theEfieldOffsets[2],2) ,0.5)/1000.0<<std::endl;


       //These maps are in: V/cm
       theEfieldOffsets[0]= theEfieldOffsets[0]/(-1000.0);
       theEfieldOffsets[1]= theEfieldOffsets[1]/(-1000.0);
       theEfieldOffsets[2]= theEfieldOffsets[2]/(-1000.0);

       if(driftcoordinate==1){
       theEfieldOffsets[0]=theEfieldOffsets[0]-fEfield;

       }else if(driftcoordinate==2){
         theEfieldOffsets[1]=theEfieldOffsets[1]-fEfield;
       }else{
         std::cout<<" Problem with drift field coordiate system for Efield calculation "<<std::endl;
       }
       theEfieldOffsets[0]= theEfieldOffsets[0]/(fEfield);
       theEfieldOffsets[1]= theEfieldOffsets[1]/(fEfield);
       theEfieldOffsets[2]= theEfieldOffsets[2]/(fEfield);



       //       std::cout<<theEfieldOffsets[0]<<" "<<theEfieldOffsets[1]<<" "<<
       //theEfieldOffsets[2]<<std::endl;


     //  }  else theEfieldOffsets.resize(3,0.0);

   return { theEfieldOffsets[0], theEfieldOffsets[1], theEfieldOffsets[2] };
   // }


}
//----------------------------------------------------------------------------
/// Primary working method of service that provides E field offsets to be
/// used in charge/light yield calculation (e.g.) for calibration
geo::Vector_t spacecharge::SpaceChargeProtoDUNEdp::GetCalEfieldOffsets(geo::Point_t const& tmp_point, int const& TPCid) const
{ //not implemented yet..
  std::vector<double> theEfieldOffsets;
  geo::Point_t point = tmp_point;


  // Efieldoffset or ratio of Efield?

  //if (fRepresentationType == "Voxelized_TH3"){

    theEfieldOffsets = GetOffsetsVoxel(point, CalSCEhistograms.at(3), CalSCEhistograms.at(4), CalSCEhistograms.at(5));
      //   }else
      //theEfieldOffsets.resize(3,0.0);

  return { -theEfieldOffsets[0], -theEfieldOffsets[1], -theEfieldOffsets[2] };
}

//----------------------------------------------------------------------------
