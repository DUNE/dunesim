/////////////////////////////////////////////////////////
/////////// ProtoDUNEbeammatch.cxx //////////////////////
/////////// Pablo F. pablo.fer@cern.ch //////////////////
/////////// June 2018 ///////////////////////////////////
/////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include "dune/EventGenerator/ProtoDUNEbeamTPCMatching/ProtoDUNEbeammatch.h"
#include "dune/EventGenerator/ProtoDUNEbeamTPCMatching/ProtoDUNEbeamToF.h"
#include <string.h>
#include <time.h>
#include <cmath>
#include <iomanip>

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace match{

// constructor
  ProtoDUNEbeammatch::ProtoDUNEbeammatch(){}

// destructor
  ProtoDUNEbeammatch::~ProtoDUNEbeammatch(){}

// Get monitors for direction (BPROFEXT and BPROF4) and ToF (TOF1 and TRIG2)
//  ProtoDUNEbeamDirection ProtoDUNEbeammatch::GetDir() const{} 
  ProtoDUNEbeamToF ProtoDUNEbeammatch::GetToF() const{} 




}//namespace
