////////////////////////////////////////////////////////////////////////
// \file SpaceCharge3x1x1dphase.cxx
//
// \brief implementation of class for storing/accessing space charge distortions for the 3x1x1 detector.
// \Adapted from SpaceChargeProtoDUNE.cxx
//
// \author kevin.fusshoeller@cern.ch
// 
////////////////////////////////////////////////////////////////////////

// C++ language includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "math.h"
#include "stdio.h"

// LArSoft includes
#include "dune/SpaceCharge/SpaceCharge3x1x1dphase.h"

// Framework includes
#include "cetlib_except/exception.h"

//-----------------------------------------------
spacecharge::SpaceCharge3x1x1dphase::SpaceCharge3x1x1dphase(
  fhicl::ParameterSet const& pset
)
{
  Configure(pset);
}

//------------------------------------------------
bool spacecharge::SpaceCharge3x1x1dphase::Configure(fhicl::ParameterSet const& pset)
{  
  fEnableSimSpatialSCE = pset.get<bool>("EnableSimSpatialSCE");
  fEnableSimEfieldSCE = pset.get<bool>("EnableSimEfieldSCE");
  fEnableCorrSCE = pset.get<bool>("EnableCorrSCE");

  if((fEnableSimSpatialSCE == true) || (fEnableSimEfieldSCE == true))
  {
    fRepresentationType = pset.get<std::string>("RepresentationType");
    fInputFilename = pset.get<std::string>("InputFilename");

    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(fInputFilename,fname);

    std::unique_ptr<TFile> infile(new TFile(fname.c_str(), "READ"));
    if(!infile->IsOpen()) throw cet::exception("SpaceCharge3x1x1dphase") << "Could not find the space charge effect file '" << fname << "'!\n";

    if(fRepresentationType == "Parametric")
    {      
      for(int i = 0; i < 5; i++)
      {
        g1_x[i] = (TGraph*)infile->Get(Form("deltaX/g1_%d",i));
        g2_x[i] = (TGraph*)infile->Get(Form("deltaX/g2_%d",i));
        g3_x[i] = (TGraph*)infile->Get(Form("deltaX/g3_%d",i));   
        g4_x[i] = (TGraph*)infile->Get(Form("deltaX/g4_%d",i));
        g5_x[i] = (TGraph*)infile->Get(Form("deltaX/g5_%d",i));

        g1_y[i] = (TGraph*)infile->Get(Form("deltaY/g1_%d",i));
        g2_y[i] = (TGraph*)infile->Get(Form("deltaY/g2_%d",i));
        g3_y[i] = (TGraph*)infile->Get(Form("deltaY/g3_%d",i));   
        g4_y[i] = (TGraph*)infile->Get(Form("deltaY/g4_%d",i));
        g5_y[i] = (TGraph*)infile->Get(Form("deltaY/g5_%d",i));
        g6_y[i] = (TGraph*)infile->Get(Form("deltaY/g6_%d",i));

        g1_z[i] = (TGraph*)infile->Get(Form("deltaZ/g1_%d",i));
        g2_z[i] = (TGraph*)infile->Get(Form("deltaZ/g2_%d",i));
        g3_z[i] = (TGraph*)infile->Get(Form("deltaZ/g3_%d",i));   
        g4_z[i] = (TGraph*)infile->Get(Form("deltaZ/g4_%d",i));

	g1_Ex[i] = (TGraph*)infile->Get(Form("deltaExOverE/g1_%d",i));
	g2_Ex[i] = (TGraph*)infile->Get(Form("deltaExOverE/g2_%d",i));
	g3_Ex[i] = (TGraph*)infile->Get(Form("deltaExOverE/g3_%d",i));
	g4_Ex[i] = (TGraph*)infile->Get(Form("deltaExOverE/g4_%d",i));
	g5_Ex[i] = (TGraph*)infile->Get(Form("deltaExOverE/g5_%d",i));

	g1_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g1_%d",i));
	g2_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g2_%d",i));
	g3_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g3_%d",i));
	g4_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g4_%d",i));
	g5_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g5_%d",i));
	g6_Ey[i] = (TGraph*)infile->Get(Form("deltaEyOverE/g6_%d",i));

	g1_Ez[i] = (TGraph*)infile->Get(Form("deltaEzOverE/g1_%d",i));
	g2_Ez[i] = (TGraph*)infile->Get(Form("deltaEzOverE/g2_%d",i));
	g3_Ez[i] = (TGraph*)infile->Get(Form("deltaEzOverE/g3_%d",i));
	g4_Ez[i] = (TGraph*)infile->Get(Form("deltaEzOverE/g4_%d",i));
      }

      g1_x[5] = (TGraph*)infile->Get("deltaX/g1_5");
      g2_x[5] = (TGraph*)infile->Get("deltaX/g2_5");
      g3_x[5] = (TGraph*)infile->Get("deltaX/g3_5");   
      g4_x[5] = (TGraph*)infile->Get("deltaX/g4_5");
      g5_x[5] = (TGraph*)infile->Get("deltaX/g5_5");

      g1_y[5] = (TGraph*)infile->Get("deltaY/g1_5");
      g2_y[5] = (TGraph*)infile->Get("deltaY/g2_5");
      g3_y[5] = (TGraph*)infile->Get("deltaY/g3_5");   
      g4_y[5] = (TGraph*)infile->Get("deltaY/g4_5");
      g5_y[5] = (TGraph*)infile->Get("deltaY/g5_5");
      g6_y[5] = (TGraph*)infile->Get("deltaY/g6_5");
      
      g1_x[6] = (TGraph*)infile->Get("deltaX/g1_6");
      g2_x[6] = (TGraph*)infile->Get("deltaX/g2_6");
      g3_x[6] = (TGraph*)infile->Get("deltaX/g3_6");
      g4_x[6] = (TGraph*)infile->Get("deltaX/g4_6");
      g5_x[6] = (TGraph*)infile->Get("deltaX/g5_6");

      g1_Ex[5] = (TGraph*)infile->Get("deltaExOverE/g1_5");
      g2_Ex[5] = (TGraph*)infile->Get("deltaExOverE/g2_5");
      g3_Ex[5] = (TGraph*)infile->Get("deltaExOverE/g3_5");
      g4_Ex[5] = (TGraph*)infile->Get("deltaExOverE/g4_5");
      g5_Ex[5] = (TGraph*)infile->Get("deltaExOverE/g5_5");

      g1_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g1_5");
      g2_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g2_5");
      g3_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g3_5");
      g4_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g4_5");
      g5_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g5_5");
      g6_Ey[5] = (TGraph*)infile->Get("deltaEyOverE/g6_5");

      g1_Ex[6] = (TGraph*)infile->Get("deltaExOverE/g1_6");
      g2_Ex[6] = (TGraph*)infile->Get("deltaExOverE/g2_6");
      g3_Ex[6] = (TGraph*)infile->Get("deltaExOverE/g3_6");
      g4_Ex[6] = (TGraph*)infile->Get("deltaExOverE/g4_6");
      g5_Ex[6] = (TGraph*)infile->Get("deltaExOverE/g5_6");
    }

    infile->Close();
  }

  if(fEnableCorrSCE == true)
  {
    // Grab other parameters from pset  
  }

  return true;
}

//------------------------------------------------
bool spacecharge::SpaceCharge3x1x1dphase::Update(uint64_t ts) 
{
  if (ts == 0) return false;

  return true;
}

//----------------------------------------------------------------------------
/// Return boolean indicating whether or not to turn simulation of SCE on for
/// spatial distortions
bool spacecharge::SpaceCharge3x1x1dphase::EnableSimSpatialSCE() const
{
  return fEnableSimSpatialSCE;
}

//----------------------------------------------------------------------------
/// Return boolean indicating whether or not to turn simulation of SCE on for
/// E field distortions
bool spacecharge::SpaceCharge3x1x1dphase::EnableSimEfieldSCE() const
{
  return fEnableSimEfieldSCE;
}

//----------------------------------------------------------------------------
/// Return boolean indicating whether or not to apply SCE corrections
bool spacecharge::SpaceCharge3x1x1dphase::EnableCorrSCE() const
{
  return fEnableCorrSCE;
}

//----------------------------------------------------------------------------
/// Primary working method of service that provides position offsets to be
/// used in ionization electron drift
geo::Vector_t spacecharge::SpaceCharge3x1x1dphase::GetPosOffsets(geo::Point_t const& point) const
{
  std::vector<double> thePosOffsets;

  if(IsInsideBoundaries(point.X(), point.Y(), point.Z()) == false)
  {
    thePosOffsets.resize(3,0.0);
  }
  else
  {
    if(fRepresentationType == "Parametric")
      thePosOffsets = GetPosOffsetsParametric(point.X(), point.Y(), point.Z());
    else
      thePosOffsets.resize(3,0.0);
  }

  return { thePosOffsets[0], thePosOffsets[1], thePosOffsets[2] };
}

//----------------------------------------------------------------------------
/// Provides position offsets using a parametric representation
std::vector<double> spacecharge::SpaceCharge3x1x1dphase::GetPosOffsetsParametric(double xVal, double yVal, double zVal) const
{
  std::vector<double> thePosOffsetsParametric;

  double xValNew = TransformX(xVal);
  double yValNew = TransformY(yVal);
  double zValNew = TransformZ(zVal);

  thePosOffsetsParametric.push_back(GetOnePosOffsetParametric(xValNew,yValNew,zValNew,"X"));
  thePosOffsetsParametric.push_back(GetOnePosOffsetParametric(xValNew,yValNew,zValNew,"Y"));
  thePosOffsetsParametric.push_back(GetOnePosOffsetParametric(xValNew,yValNew,zValNew,"Z"));

  return thePosOffsetsParametric;
}

//----------------------------------------------------------------------------
/// Provides one position offset using a parametric representation, for a given
/// axis
double spacecharge::SpaceCharge3x1x1dphase::GetOnePosOffsetParametric(double xValNew, double yValNew, double zValNew, std::string axis) const
{      
  double parA[6][7];
  double parB[6];
  
  for(int j = 0; j < 6; j++)
  {
    for(int i = 0; i < 7; i++)
      parA[j][i] = 0.0;
  
    parB[j] = 0.0;
  }
  
  if(axis == "X")
  {
    for(int j = 0; j < 7; j++)
    {
      parA[0][j] = g1_x[j]->Eval(zValNew);
      parA[1][j] = g2_x[j]->Eval(zValNew);
      parA[2][j] = g3_x[j]->Eval(zValNew);
      parA[3][j] = g4_x[j]->Eval(zValNew);
      parA[4][j] = g5_x[j]->Eval(zValNew);
    }
  
    f1_x->SetParameters(parA[0]);
    f2_x->SetParameters(parA[1]);
    f3_x->SetParameters(parA[2]);
    f4_x->SetParameters(parA[3]);
    f5_x->SetParameters(parA[4]);
  }
  else if(axis == "Y")
  {
    for(int j = 0; j < 6; j++)
    {
      parA[0][j] = g1_y[j]->Eval(zValNew);
      parA[1][j] = g2_y[j]->Eval(zValNew);
      parA[2][j] = g3_y[j]->Eval(zValNew);
      parA[3][j] = g4_y[j]->Eval(zValNew);
      parA[4][j] = g5_y[j]->Eval(zValNew);
      parA[5][j] = g6_y[j]->Eval(zValNew);
    }
  
    f1_y->SetParameters(parA[0]);
    f2_y->SetParameters(parA[1]);
    f3_y->SetParameters(parA[2]);
    f4_y->SetParameters(parA[3]);
    f5_y->SetParameters(parA[4]);
    f6_y->SetParameters(parA[5]);
  }
  else if(axis == "Z")
  {
    for(int j = 0; j < 5; j++)
    {
      parA[0][j] = g1_z[j]->Eval(zValNew);
      parA[1][j] = g2_z[j]->Eval(zValNew);
      parA[2][j] = g3_z[j]->Eval(zValNew);
      parA[3][j] = g4_z[j]->Eval(zValNew);
    }
  
    f1_z->SetParameters(parA[0]);
    f2_z->SetParameters(parA[1]);
    f3_z->SetParameters(parA[2]);
    f4_z->SetParameters(parA[3]);
  }
  
  double aValNew;
  double bValNew;
  
  if(axis == "Y")
  {
    aValNew = xValNew;
    bValNew = yValNew;
  }
  else
  {
    aValNew = yValNew;
    bValNew = xValNew;
  }
  
  double offsetValNew = 0.0;
  if(axis == "X")
  {
    parB[0] = f1_x->Eval(aValNew);
    parB[1] = f2_x->Eval(aValNew);
    parB[2] = f3_x->Eval(aValNew);
    parB[3] = f4_x->Eval(aValNew);
    parB[4] = f5_x->Eval(aValNew);
  
    fFinal_x->SetParameters(parB);
    offsetValNew = fFinal_x->Eval(bValNew);
  }
  else if(axis == "Y")
  {
    parB[0] = f1_y->Eval(aValNew);
    parB[1] = f2_y->Eval(aValNew);
    parB[2] = f3_y->Eval(aValNew);
    parB[3] = f4_y->Eval(aValNew);
    parB[4] = f5_y->Eval(aValNew);
    parB[5] = f6_y->Eval(aValNew);
  
    fFinal_y->SetParameters(parB);
    offsetValNew = fFinal_y->Eval(bValNew);
  }
  else if(axis == "Z")
  {
    parB[0] = f1_z->Eval(aValNew);
    parB[1] = f2_z->Eval(aValNew);
    parB[2] = f3_z->Eval(aValNew);
    parB[3] = f4_z->Eval(aValNew);
  
    fFinal_z->SetParameters(parB);
    offsetValNew = fFinal_z->Eval(bValNew);
  }
   
  return offsetValNew;
}

//----------------------------------------------------------------------------
/// Primary working method of service that provides E field offsets to be
/// used in charge/light yield calculation (e.g.)
geo::Vector_t spacecharge::SpaceCharge3x1x1dphase::GetEfieldOffsets(geo::Point_t const& point) const
{
  std::vector<double> theEfieldOffsets;

  if(IsInsideBoundaries(point.X(), point.Y(), point.Z()) == false)
  {
    theEfieldOffsets.resize(3,0.0);
  }
  else
  {
    if(fRepresentationType == "Parametric")
      theEfieldOffsets = GetEfieldOffsetsParametric(point.X(), point.Y(), point.Z());
    else
      theEfieldOffsets.resize(3,0.0);
  }

  return { -theEfieldOffsets[0], -theEfieldOffsets[1], -theEfieldOffsets[2] };
}

//----------------------------------------------------------------------------
/// Provides E field offsets using a parametric representation
std::vector<double> spacecharge::SpaceCharge3x1x1dphase::GetEfieldOffsetsParametric(double xVal, double yVal, double zVal) const
{
  std::vector<double> theEfieldOffsetsParametric;

  double xValNew = TransformX(xVal);
  double yValNew = TransformY(yVal);
  double zValNew = TransformZ(zVal);

  theEfieldOffsetsParametric.push_back(GetOneEfieldOffsetParametric(xValNew,yValNew,zValNew,"X"));
  theEfieldOffsetsParametric.push_back(GetOneEfieldOffsetParametric(xValNew,yValNew,zValNew,"Y"));
  theEfieldOffsetsParametric.push_back(GetOneEfieldOffsetParametric(xValNew,yValNew,zValNew,"Z"));

  return theEfieldOffsetsParametric;
}

//----------------------------------------------------------------------------
/// Provides one E field offset using a parametric representation, for a given
/// axis
double spacecharge::SpaceCharge3x1x1dphase::GetOneEfieldOffsetParametric(double xValNew, double yValNew, double zValNew, std::string axis) const
{      
  double parA[6][7];
  double parB[6];
  
  for(int j = 0; j < 6; j++)
  {
    for(int i = 0; i < 7; i++)
      parA[j][i] = 0.0;
  
    parB[j] = 0.0;
  }
  
  if(axis == "X")
  {
    for(int j = 0; j < 7; j++)
    {
      parA[0][j] = g1_Ex[j]->Eval(zValNew);
      parA[1][j] = g2_Ex[j]->Eval(zValNew);
      parA[2][j] = g3_Ex[j]->Eval(zValNew);
      parA[3][j] = g4_Ex[j]->Eval(zValNew);
      parA[4][j] = g5_Ex[j]->Eval(zValNew);
    }
  
    f1_Ex->SetParameters(parA[0]);
    f2_Ex->SetParameters(parA[1]);
    f3_Ex->SetParameters(parA[2]);
    f4_Ex->SetParameters(parA[3]);
    f5_Ex->SetParameters(parA[4]);
  }
  else if(axis == "Y")
  {
    for(int j = 0; j < 6; j++)
    {
      parA[0][j] = g1_Ey[j]->Eval(zValNew);
      parA[1][j] = g2_Ey[j]->Eval(zValNew);
      parA[2][j] = g3_Ey[j]->Eval(zValNew);
      parA[3][j] = g4_Ey[j]->Eval(zValNew);
      parA[4][j] = g5_Ey[j]->Eval(zValNew);
      parA[5][j] = g6_Ey[j]->Eval(zValNew);
    }
  
    f1_Ey->SetParameters(parA[0]);
    f2_Ey->SetParameters(parA[1]);
    f3_Ey->SetParameters(parA[2]);
    f4_Ey->SetParameters(parA[3]);
    f5_Ey->SetParameters(parA[4]);
    f6_Ey->SetParameters(parA[5]);
  }
  else if(axis == "Z")
  {
    for(int j = 0; j < 5; j++)
    {
      parA[0][j] = g1_Ez[j]->Eval(zValNew);
      parA[1][j] = g2_Ez[j]->Eval(zValNew);
      parA[2][j] = g3_Ez[j]->Eval(zValNew);
      parA[3][j] = g4_Ez[j]->Eval(zValNew);
    }
  
    f1_Ez->SetParameters(parA[0]);
    f2_Ez->SetParameters(parA[1]);
    f3_Ez->SetParameters(parA[2]);
    f4_Ez->SetParameters(parA[3]);
  }
  
  double aValNew;
  double bValNew;
  
  if(axis == "Y")
  {
    aValNew = xValNew;
    bValNew = yValNew;
  }
  else
  {
    aValNew = yValNew;
    bValNew = xValNew;
  }
  
  double offsetValNew = 0.0;
  if(axis == "X")
  {
    parB[0] = f1_Ex->Eval(aValNew);
    parB[1] = f2_Ex->Eval(aValNew);
    parB[2] = f3_Ex->Eval(aValNew);
    parB[3] = f4_Ex->Eval(aValNew);
    parB[4] = f5_Ex->Eval(aValNew);
  
    fFinal_Ex->SetParameters(parB);
    offsetValNew = fFinal_Ex->Eval(bValNew);
  }
  else if(axis == "Y")
  {
    parB[0] = f1_Ey->Eval(aValNew);
    parB[1] = f2_Ey->Eval(aValNew);
    parB[2] = f3_Ey->Eval(aValNew);
    parB[3] = f4_Ey->Eval(aValNew);
    parB[4] = f5_Ey->Eval(aValNew);
    parB[5] = f6_Ey->Eval(aValNew);
  
    fFinal_Ey->SetParameters(parB);
    offsetValNew = fFinal_Ey->Eval(bValNew);
  }
  else if(axis == "Z")
  {
    parB[0] = f1_Ez->Eval(aValNew);
    parB[1] = f2_Ez->Eval(aValNew);
    parB[2] = f3_Ez->Eval(aValNew);
    parB[3] = f4_Ez->Eval(aValNew);
  
    fFinal_Ez->SetParameters(parB);
    offsetValNew = fFinal_Ez->Eval(bValNew);
  }
  
  return offsetValNew;
}

//----------------------------------------------------------------------------
/// Transform X to SCE X coordinate:  [0.0,3.6] --> [0.0,3.6]
double spacecharge::SpaceCharge3x1x1dphase::TransformX(double xVal) const
{
  double xValNew = xVal;

  return xValNew;
}

//----------------------------------------------------------------------------
/// Transform Y to SCE Y coordinate:  [0.0,6.08] --> [0.0,6.0]
double spacecharge::SpaceCharge3x1x1dphase::TransformY(double yVal) const
{
  double yValNew =  yVal;

  return yValNew;
}

//----------------------------------------------------------------------------
/// Transform Z to SCE Z coordinate:  [0.0,6.97] --> [0.0,7.2]
double spacecharge::SpaceCharge3x1x1dphase::TransformZ(double zVal) const
{
  double zValNew = zVal;

  return zValNew;
}

//----------------------------------------------------------------------------
/// Check to see if point is inside boundaries of map (allow to go slightly out of range)
bool spacecharge::SpaceCharge3x1x1dphase::IsInsideBoundaries(double xVal, double yVal, double zVal) const
{
  bool isInside = true;

  if((xVal < -360.0) || (xVal > 360.0) || (yVal < -5.0) || (yVal > 615.0) || (zVal < -5.0) || (zVal > 705.0))
  {
    isInside = false;
  }

  return isInside;
}
