////////////////////////////////////////////////////////////////////////
// Class:       SimCounter35t
// Module Type: producer
// File:        SimCounter35t_module.cc
//
// Generated at Wed Mar 18 05:42:26 2015 by Matthew Thiesse using artmod
// from cetpkgsupport v1_08_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "artextensions/SeedService/SeedService.hh"

#include "RawData/raw.h"
#include "RawData/ExternalTrigger.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/TimeService.h"
#include "Simulation/AuxDetSimChannel.h"

#include "CLHEP/Random/RandFlat.h"

#include <memory>
#include <iostream>
#include <sstream>

#include "TTree.h"

namespace detsim {
  class SimCounter35t;
}

class detsim::SimCounter35t : public art::EDProducer {
public:
  explicit SimCounter35t(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SimCounter35t(SimCounter35t const &) = delete;
  SimCounter35t(SimCounter35t &&) = delete;
  SimCounter35t & operator = (SimCounter35t const &) = delete;
  SimCounter35t & operator = (SimCounter35t &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // chanTick represents a single time "tick" for a single counter.
  // multiple hits during the same time window on the same counter are possible.
  // using this struct is an alternative to initializing an array of 50000 values for 
  // each AuxDetSimChannel where most of the entries are zero. This way, information is
  // stored only if there is nonzero energy deposition.
  struct chanTick {
    unsigned int tick;
    double eDep;
    unsigned int auxDetID;
    int numHits;
  
    chanTick(int t, int adid, double ed) 
      : tick(t), eDep(ed), auxDetID(adid), numHits(1) {
    }

    void Add(double ed) { 
      eDep += ed;
      ++numHits;
    }
  };

  // container to hold info for ticks with nonzero energy deposition
  std::vector<chanTick> tickv;

  // fhicl parameters
  std::string fLArG4ModuleLabel;
  bool fMakeTree;
  double fBSUTriggerThreshold; // MeV
  double fTSUTriggerThreshold; // MeV
  double fTriggerEfficiency;
  double fClockSpeedCounter;
  unsigned int fCombinedTimeDelay; // ns

  // Tree containing aux det hit info. copied from dune/LArG4/CheckAuxDet_module.cc
  TTree *fTree;
  int run;
  int subrun;
  int event;
  int nauxdets;
  uint32_t auxdetid;
  int ntrkids;
  double entryx;
  double entryy;
  double entryz;
  double entryt;
  double exitx;
  double exity;
  double exitz;
  double exitt;
  double exitpx;
  double exitpy;
  double exitpz;
  int trackid;
  double energy;
};

/////////////////////////////////////////////////////////////////////////////////

detsim::SimCounter35t::SimCounter35t(fhicl::ParameterSet const & p)
  :
  fLArG4ModuleLabel(p.get<std::string>("LArGeantModuleLabel", "largeant")),
  fMakeTree(p.get<bool>("MakeTree",false)),
  fBSUTriggerThreshold(p.get<double>("BSUTriggerThreshold",0.5)),// MeV
  fTSUTriggerThreshold(p.get<double>("TSUTriggerThreshold",0.25)),// MeV
  fTriggerEfficiency(p.get<double>("TriggerEfficiency",1.)),
  fClockSpeedCounter(p.get<double>("ClockSpeedCounter",31.25)), // MHz
  fCombinedTimeDelay(p.get<double>("CombinedTimeDelay",160)) // ns
{
  art::ServiceHandle<artext::SeedService>()->createEngine(*this, "HepJamesRandom", "rand");
 
  produces< std::vector< raw::ExternalTrigger > >();
}

/////////////////////////////////////////////////////////////////////////////////

void detsim::SimCounter35t::produce(art::Event & e)
{
  int skippedHitsIneff = 0;
  int skippedHitsOutRange = 0;

  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine &engine = rng->getEngine("rand");
  CLHEP::RandFlat flat(engine,0,1);

  art::ServiceHandle<util::TimeService> ts;
  art::ServiceHandle<util::DetectorProperties> detprop;

  // make unique_ptr that allows ownership of the produced triggers to be transferred to the art::Event after the put statement below
  std::unique_ptr<std::vector<raw::ExternalTrigger>> trigcol(new std::vector<raw::ExternalTrigger>);

  run = e.run();
  subrun = e.subRun();
  event = e.id().event();

  // get all auzDetSimChannel information containing AuxDetIDEs
  std::vector<const sim::AuxDetSimChannel*> fAuxDetSimChannels;
  e.getView(fLArG4ModuleLabel, fAuxDetSimChannels);

  nauxdets = fAuxDetSimChannels.size();

  // loop over all AuxDetSimChannels
  for (size_t i=0; i<fAuxDetSimChannels.size(); ++i) {
    const sim::AuxDetSimChannel* c = fAuxDetSimChannels[i];
    auxdetid = c->AuxDetID();

    // get AuxDetIDEs associated with sim channel
    const std::vector<sim::AuxDetIDE>& setOfIDEs = c->AuxDetIDEs();
    ntrkids = setOfIDEs.size();

    // loop over all AuxDetIDEs
    for (size_t j=0; j<setOfIDEs.size(); ++j) {

      // fill tree
      if (fMakeTree) {
        entryx = setOfIDEs[j].entryX;
        entryy = setOfIDEs[j].entryY;
        entryz = setOfIDEs[j].entryZ;
        entryt = setOfIDEs[j].entryT;
        exitx = setOfIDEs[j].exitX;
        exity = setOfIDEs[j].exitY;
        exitz = setOfIDEs[j].exitZ;
        exitt = setOfIDEs[j].exitT;
        exitpx = setOfIDEs[j].exitMomentumX;
        exitpy = setOfIDEs[j].exitMomentumY;
        exitpz = setOfIDEs[j].exitMomentumZ;
        energy = setOfIDEs[j].energyDeposited;
        trackid = setOfIDEs[j].trackID;
        fTree->Fill();
      }

      double randEff = flat.fire();
      if (randEff > fTriggerEfficiency) {
	++skippedHitsIneff;
	continue;
      }

      // calculate the time length of one window
      double triggerOffsetTPC = ts->TriggerOffsetTPC()*1.e3; // ns
      double readoutWindowSizeTPC = detprop->ReadOutWindowSize(); // tpc ticks
      double clockSpeedTPC = ts->TPCClock().Frequency()/1.e6; // MHz
      double windowLength = readoutWindowSizeTPC/clockSpeedTPC; // us

      // get information from AuxDetIDE
      double time = setOfIDEs[j].entryT+fCombinedTimeDelay-triggerOffsetTPC; // ns
      if (time<0 || time>windowLength*1000) {
	++skippedHitsOutRange;
	continue;
      }
      uint32_t tickIDE = time*fClockSpeedCounter/1000;
      double edepIDE = setOfIDEs[j].energyDeposited*1000;//MeV
      
      // loop over tickv to add eDep
      std::vector<chanTick>::iterator it;
      for (it = tickv.begin(); it != tickv.end(); ++it) {
	chanTick* ct = &*it;
	if (ct->tick == tickIDE && ct->auxDetID == auxdetid) {
	  ct->Add(edepIDE);
	  break;
	}
      } 

      // if the chanTick doesn't exist already, make one
      if (it==tickv.end()) {
	// initialise chanTick with IDE values
	tickv.push_back(chanTick(tickIDE,auxdetid,edepIDE));
      }
    }
  }    

  // fill collection of raw::ExternalTriggers, if eDep is above threshold
  // Build triggers from individual hits
  uint iRM=0; uint iCL=0; uint iNU=0; uint iNL=0; uint iSU=0; uint iSL=0; uint iWU=0; uint iEL=0;
  std::vector<chanTick>::iterator it;
  for (it = tickv.begin(); it != tickv.end(); ++it) {
    chanTick ct = *it;
    if ( (ct.auxDetID >= 44 && ct.auxDetID <= 91 && ct.eDep > fBSUTriggerThreshold) || 
	 (ct.auxDetID >= 0 && ct.auxDetID <=43 && ct.eDep > fTSUTriggerThreshold) ||
	 (ct.auxDetID >= 92 && ct.eDep > 1.e-6) )
      trigcol->push_back(raw::ExternalTrigger(ct.auxDetID,ct.tick));
    if (ct.auxDetID<6) iSL=ct.tick; 
    if (ct.auxDetID>5 &&ct.auxDetID<16) iEL=ct.tick; 
    if (ct.auxDetID>15 &&ct.auxDetID<22) iNL=ct.tick; 
    if (ct.auxDetID>21 &&ct.auxDetID<28) iNU=ct.tick; 
    if (ct.auxDetID>27 &&ct.auxDetID<38) iWU=ct.tick; 
    if (ct.auxDetID>43 &&ct.auxDetID<57) iCL=ct.tick; 
    if (ct.auxDetID>66 &&ct.auxDetID<83) iRM=ct.tick; 
  }

  if (iRM>0 && iCL>0 && (fabs(iRM-iCL)<5)) trigcol->push_back(raw::ExternalTrigger(110,0.5*(iRM+iCL)));
  if (iEL>0 && iWU>0 && (fabs(iEL-iWU)<5)) trigcol->push_back(raw::ExternalTrigger(111,0.5*(iEL+iWU)));
  if (iNU>0 && iSL>0 && (fabs(iNU-iSL)<5)) trigcol->push_back(raw::ExternalTrigger(112,0.5*(iNU+iSL)));
  if (iSU>0 && iNL>0 && (fabs(iSU-iNL)<5)) trigcol->push_back(raw::ExternalTrigger(113,0.5*(iSU+iNL)));


  // put ExternalTrigger collection into event record
  e.put(std::move(trigcol));


  // output hit information for event
  std::ostringstream out;
  for (std::vector<chanTick>::iterator it = tickv.begin(); it != tickv.end(); ++it) {
    chanTick ct = *it;
    if ( (ct.auxDetID >= 44 && ct.auxDetID <= 91 && ct.eDep > fBSUTriggerThreshold) ||
         (ct.auxDetID >= 0 && ct.auxDetID <=43 && ct.eDep > fTSUTriggerThreshold) ||
	 (ct.auxDetID >= 92 && ct.eDep > 1.e-6) )
      out << "AuxDet " << ct.auxDetID << " had " << ct.numHits << " hits at readout tick " << ct.tick << ". Total eDep = " << ct.eDep << " MeV.\n";
  }
  if (skippedHitsIneff) out << skippedHitsIneff << " hits were skipped due to counter inefficiency.\n";
  if (skippedHitsOutRange) out << skippedHitsOutRange << " hits were skipped for being out of TPC window range.\n";
  mf::LogInfo("SimCounter35t") << out.str();

  tickv.clear();

}

/////////////////////////////////////////////////////////////////////////////////

void detsim::SimCounter35t::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  if (fMakeTree) {
    fTree = tfs->make<TTree>("SimCounter35t","SimCounter35t");
    fTree->Branch("run",&run,"run/I");
    fTree->Branch("subrun",&subrun,"subrun/I");
    fTree->Branch("event",&event,"event/I");
    fTree->Branch("nauxdets",&nauxdets,"nauxdets/I");
    fTree->Branch("auxdetid",&auxdetid,"auxdetid/I");
    fTree->Branch("ntrkids",&ntrkids,"ntrkids/I");
    fTree->Branch("entryx",&entryx,"entryx/D");
    fTree->Branch("entryy",&entryy,"entryy/D");
    fTree->Branch("entryz",&entryz,"entryz/D");
    fTree->Branch("entryt",&entryt,"entryt/D");
    fTree->Branch("exitx",&exitx,"exitx/D");
    fTree->Branch("exity",&exity,"exity/D");
    fTree->Branch("exitz",&exitz,"exitz/D");
    fTree->Branch("exitt",&exitt,"exitt/D");
    fTree->Branch("exitpx",&exitpx,"exitpx/D");
    fTree->Branch("exitpy",&exitpy,"exitpy/D");
    fTree->Branch("exitpz",&exitpz,"exitpz/D");
    fTree->Branch("trackid",&trackid,"trackid/D");
    fTree->Branch("energy",&energy,"energy/D");
  }
}

DEFINE_ART_MODULE(detsim::SimCounter35t)
