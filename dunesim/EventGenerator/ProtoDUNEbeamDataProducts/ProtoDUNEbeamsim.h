// ProtoDUNEbeamsim.h
//
// Caroline Zhang
// carolineligezhang@gmail.com
// August 2017
//
// Modified by Pablo and Leigh H. Howard, 
// // pablo.fer@cern.ch
// // July 2018
// /////////////////////////////////////////////////////////

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

    // Leigh: I think it would make much more sense to store a ProtoDUNEBeamInstrument object
    // for each of the instruments
    std::vector<ProtoDUNEBeamInstrument> fAllInstruments;

  public:

    void AddInstrument(ProtoDUNEBeamInstrument newInst);
    ProtoDUNEBeamInstrument GetInstrument(std::string name) const;
    unsigned short NInstruments() const {return fAllInstruments.size();};
    
  };
}
////////////////////////////////////////////////////////////////////////
#endif // PROTODUNESIM_H
