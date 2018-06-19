/////////////////////////////////////////////////////////
///////////// ProtoDUNEbeammatch.h //////////////////////
///////////// Pablo F. pablo.fer@cern.ch ////////////////
///////////// June 2018 /////////////////////////////////
/////////////////////////////////////////////////////////

#ifndef PROTODUNEBEAMMATCH_H
#define PROTODUNEBEAMMATCH_H

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

//#include "dune/EventGenerator/ProtoDUNEbeamTPCmatching/ProtoDUNEBeamDirection.h"
#include "dune/EventGenerator/ProtoDUNEbeamTPCmatching/ProtoDUNEBeamToF.h"

namespace match {
 class ProtoDUNEbeammatch  {

  public:

    ProtoDUNEbeammatch(); //constructor

    ~ProtoDUNEbeammatch(); //destructor

  private:

  public:

//    ProtoDUNEBeamDirection GetDir() const;  
    ProtoDUNEBeamToF GetToF() const;



};
}//namespace
#endif // PROTODUNEMATCH_H

