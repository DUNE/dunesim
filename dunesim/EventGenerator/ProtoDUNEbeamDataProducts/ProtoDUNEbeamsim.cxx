/////////////////////////////////////////////////////////
///////// ProtoDUNEbeamsim.cxx///////////////////////////
///////// Caroline Zhang carolineligezhang@gmail.com/////
///////// August 2017 ///////////////////////////////////
/////////////////////////////////////////////////////////
// Modified by Pablo and Leigh H. Howard, 
// pablo.fer@cern.ch
// July 2018
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
  ProtoDUNEbeamsim::ProtoDUNEbeamsim(){}
  
  //-------------------------------default destructor------------------------------------------
  ProtoDUNEbeamsim::~ProtoDUNEbeamsim(){ }

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






