
////////////////////////////////////////////////////////////////////////
// Class:       ProtoDUNEBeam
// Module Type: producer
// File:        ProtoDUNEBeam_module.cc
//
// Generated at Thu Nov 17 11:20:31 2016 by Leigh Howard Whitehead,42 3-039,+41227672470, using artmod
// from cetpkgsupport v1_11_00.
////////////////////////////////////////////////////////////////////////
// Modified by Caroline Zhang with inputs from Karl Wharburton  
// August 2017 for beam simulation storage
// Email: carolineligezhang@gmail.com
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSetRegistry.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include <memory>
#include <string>
#include <map>
#include <utility>
#include <vector>
// Added for ProtoDUNE beam simulation storage
//#include "lardataobj/Simulation/ProtoDUNEbeamsim.h"
#include "dune/EventGenerator/ProtoDUNEbeamDataProducts/ProtoDUNEbeamsim.h"
#include "lardata/Utilities/AssociationUtil.h"
// art extensions
#include "nutools/RandomUtils/NuRandomService.h"
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include "TSystem.h"
#include "CLHEP/Random/RandFlat.h"
#include "ifdh.h"
#include <sys/stat.h>

namespace evgen{
    class ProtoDUNEBeam;
    
    class ProtoDUNEBeam : public art::EDProducer {
    public:
        explicit ProtoDUNEBeam(fhicl::ParameterSet const & p);
        // The destructor generated by the compiler is fine for classes
        // without bare pointers or other resource use.
        ~ProtoDUNEBeam();
        
        // Plugins should not be copied or assigned.
        ProtoDUNEBeam(ProtoDUNEBeam const &) = delete;
        ProtoDUNEBeam(ProtoDUNEBeam &&) = delete;
        ProtoDUNEBeam & operator = (ProtoDUNEBeam const &) = delete;
        ProtoDUNEBeam & operator = (ProtoDUNEBeam &&) = delete;
        
        // Required functions.
        void produce(art::Event & e) override;
        void beginJob() override;
        void beginRun(art::Run& run) override;
        void endJob() override;
        
    private:
        
        // We need to make a map of good particle event numbers and all
        // matching entries in the overlay events in the main particle list.
        std::map<int,std::vector<std::pair<int,std::vector<int> > > > fEventParticleMap;
        
        // A second map storing the trigger time of the good particle.
        std::map<int,float> fGoodParticleTriggerTime;
        
        // Track ID of the good particle.
        std::map<int,int> fGoodParticleTrackID;
        
        // A list of good events and an index for it.
        unsigned int fCurrentGoodEvent;
        std::vector<int> fGoodEventList;
        
        // Calculate how many overlay events we need.
        void CalculateNOverlays();
        
        // Check if a given beam event is close enough to a good particle event to be useful.
        int IsOverlayEvent(int event, int nOverlay);
        std::vector<int> GetAllOverlays(int event, int nOverlay);
        
        // Fill the above maps and vector.
        void FillParticleMaps();
        
        // Generate a true event based on a single entry from the input tree.
        void GenerateTrueEvent(simb::MCTruth &mcTruth, std::vector<sim::ProtoDUNEbeamsim> &beamsimcol);
        
        // Handle root files from beam instrumentation group
        void OpenInputFile();
        
        // Generate a TLorentzVector for position making sure we get the
        // coordinates as we need them.
        TLorentzVector ConvertCoordinates(float x, float y, float z, float t);
        
        // Make the momentum vector, rotating as required.
        TLorentzVector MakeMomentumVector(float px, float py, float pz, int pdg);
        
        std::string fFileName;
        std::string fGoodParticleTreeName;
        std::string fAllParticlesTreeName;
        
        // The current event number. Ideally this could be an unsigned int,
        // but we will need to compare it to some ints later on.
        int fEventNumber;
        
        // Let the user define the event to start at
        int fStartEvent;
        
        TFile* fInputFile;
        // Input file provides a TTree that we need to read.
        TTree* fGoodParticleTree;
        TTree* fAllParticlesTree;
        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        //////////////////// Good Particle Tree Variables added by Caroline for beam simulation storage
        ////////////////////////////////////////////////////////////////////////////////////////////////
        //For the LAG_ENTRY part
        Float_t fGoodLag_ENTRY_x;
        Float_t fGoodLag_ENTRY_y;
        Float_t fGoodLag_ENTRY_z;
        Float_t fGoodLagPDGID;
        Float_t fGoodLagPx;
        Float_t fGoodLagPy;
        Float_t fGoodLagPz;
        
        // add more lines for the trig 2 part
        Float_t fGoodTRIG2_x;
        Float_t fGoodTRIG2_y;
        Float_t fGoodTRIG2_z;
        Float_t fGoodTRIG2_PDGid;
        Float_t fGoodTRIG2_Px;
        Float_t fGoodTRIG2_Py;
        Float_t fGoodTRIG2_Pz;
        Float_t fGoodTRIG2_EventID;
        Float_t fGoodTRIG2_TrackID;
        
        //add more lines for the BPROF4 part
        Float_t fGoodBPROF4_x;
        Float_t fGoodBPROF4_y;
        Float_t fGoodBPROF4_z;
        Float_t fGoodBPROF4_PDGid;
        Float_t fGoodBPROF4_Px;
        Float_t fGoodBPROF4_Py;
        Float_t fGoodBPROF4_Pz;
        Float_t fGoodBPROF4_EventID;
        Float_t fGoodBPROF4_TrackID;
        //////////////////////////////////
        
        
        // Members we need to extract from the tree
        float fX, fY, fZ;
        float fPx, fPy, fPz;
        float fPDG; // Input tree has all floats
        float fBeamEvent;
        float fTrackID;
        // We need two times: the trigger time, and the time at the entry point
        // to the TPC where we generate the event.
        float fEntryT, fTriggerT;
        
        // Define the coordinate transform from the beam frame to the detector frame
        float fBeamX;
        float fBeamY;
        float fBeamZ;
        float fRotateXZ;
        float fRotateYZ;
        
        // Parameters from the .fcl file to deal with overlaying events
        float fIntensity; // Number of interactions on the secondary target per SPS spill
        float fReadoutWindow; // Readout window (needs to match the values used in the simulation) in milliseconds
        float fBeamSpillLength; // The SPS spill length in seconds
        
        // Number of beam interactions to overlay.
        int fOverlays;
        
        // In the absense of the full shielding in the beam simulation, add a momentum cut away from the beam pipe.
        float fBackgroundMomentumCut; // Momentum cut in GeV/c
        float fBackgroundRadiusCut; // Radius away from beam pipe to consider in cm
        
        ifdh_ns::ifdh* fIFDH;
    };
}


//---------------------------------------------------------------------------------
//----------------------------------------constructors-----------------------------
evgen::ProtoDUNEBeam::ProtoDUNEBeam(fhicl::ParameterSet const & pset)
{
    
    // Call appropriate produces<>() functions here.
    produces< std::vector<simb::MCTruth> >();
    produces<std::vector<sim::ProtoDUNEbeamsim>>();
    produces< sumdata::RunData, art::InRun >();
    //produces< art::Assns<sim::ProtoDUNEbeamsim, simb::MCTruth>>();
    // File reading variable initialisations
    fFileName = pset.get< std::string>("FileName");
    fGoodParticleTreeName = pset.get< std::string>("GoodParticleTreeName");
    fAllParticlesTreeName = pset.get< std::string>("AllParticlesTreeName");
    
    // Intensity variables
    fIntensity = pset.get<float>("Intensity");
    fReadoutWindow = pset.get<float>("ReadoutWindow");
    fBeamSpillLength = pset.get<float>("BeamSpillLength");
    
    // Background cut variables
    fBackgroundMomentumCut = pset.get<float>("BkgMomentumCut");
    fBackgroundRadiusCut = pset.get<float>("BkgRadiusCut");
    
    // See if the user wants to start at an event other than zero.
    fStartEvent = pset.get<int>("StartEvent");
    
    // Or maybe there was --nskip specified in the command line or skipEvents in FHiCL?
    for (auto const & p : fhicl::ParameterSetRegistry::get())
    {
        if (p.second.has_key("source.skipEvents"))
        {
            fStartEvent += p.second.get<int>("source.skipEvents");
            break; // take the first occurence
        } // no "else", if parameter not found, then just don't change anything
    }
    // ...and if there is -e option or firstEvent in FHiCL, this add up to the no. of events to skip.
    for (auto const & p : fhicl::ParameterSetRegistry::get())
    {
        if (p.second.has_key("source.firstEvent"))
        {
            int fe = p.second.get<int>("source.firstEvent") - 1; // events base index is 1
            if (fe > 0) fStartEvent += fe;
            break; // take the first occurence
        } // no "else", if parameter not found, then just don't change anything
    }
    mf::LogInfo("ProtoDUNEBeam") << "Skip " << fStartEvent << " first events from the input file.";
    
    fEventNumber = 0;
    
    // Coordinate transform
    fBeamX = pset.get<float>("BeamX");
    fBeamY = pset.get<float>("BeamY");
    fBeamZ = pset.get<float>("BeamZ");
    fRotateXZ = pset.get<float>("RotateXZ");
    fRotateYZ = pset.get<float>("RotateYZ");
    
    // Initialise the input file and tree to be null.
    fInputFile = 0x0;
    fGoodParticleTree = 0x0;
    fAllParticlesTree = 0x0;
    fIFDH = 0;
    
    fCurrentGoodEvent = 0;
    
    // Make sure we use ifdh to open the beam input file.
    OpenInputFile();
    
    // Create the random number generator
    std::string const instanceName = "protoDUNEBeam";
    auto& Seeds = *(art::ServiceHandle<rndm::NuRandomService>());
    
    // declare an engine; NuRandomService associates an (unknown) engine, in
    // the current module and an instance name, with a seed (returned)
    auto const seed = Seeds.declareEngine(instanceName);
    
    // now create the engine (for example, use art); seed will be set
    createEngine(seed, "HepJamesRandom", instanceName);
    
    // finally, complete the registration; seed will be set again
    //	art::ServiceHandle<art::RandomNumberGenerator> RNG;
    //	Seeds.defineEngine(RNG->getEngine(instanceName));
}



//-----------------------------default destructor----------------------------------
//------------------------------------------------------------------------------------
evgen::ProtoDUNEBeam::~ProtoDUNEBeam()
{
    fIFDH->cleanup();
    
}

//-------------------------------------------------------------------------------------
void evgen::ProtoDUNEBeam::beginJob(){
    
    fInputFile = new TFile(fFileName.c_str(),"READ");
    // Check we have the file
    if(fInputFile == 0x0){
        throw cet::exception("ProtoDUNEBeam") << "Input file " << fFileName << " cannot be read.\n";
    }
    
    fGoodParticleTree = (TTree*)fInputFile->Get(fGoodParticleTreeName.c_str());
    // Check we have the tree
    if(fGoodParticleTree == 0x0){
        throw cet::exception("ProtoDUNEBeam") << "Input tree " << fGoodParticleTreeName << " cannot be read.\n";
    }
    
    fAllParticlesTree = (TTree*)fInputFile->Get(fAllParticlesTreeName.c_str());
    // Check we have the tree
    if(fAllParticlesTree == 0x0){
        throw cet::exception("ProtoDUNEBeam") << "Input tree " << fAllParticlesTreeName << " cannot be read.\n";
    }
    
    // Since this is technically an ntuple, all objects are floats
    // Position four-vector components
    fAllParticlesTree->SetBranchAddress("x",&fX);
    fAllParticlesTree->SetBranchAddress("y",&fY);
    fAllParticlesTree->SetBranchAddress("z",&fZ);
    fAllParticlesTree->SetBranchAddress("t",&fEntryT);
    // Momentum components
    fAllParticlesTree->SetBranchAddress("Px",&fPx);
    fAllParticlesTree->SetBranchAddress("Py",&fPy);
    fAllParticlesTree->SetBranchAddress("Pz",&fPz);
    // PDG code
    fAllParticlesTree->SetBranchAddress("PDGid",&fPDG);
    // Event and track number
    fAllParticlesTree->SetBranchAddress("EventID",&fBeamEvent);
    fAllParticlesTree->SetBranchAddress("TrackID",&fTrackID);
    
    // We only need the trigger time and event number from the good particle tree.
    // The good particle tree variable should match the names of the other trees
    std::string namePrefix = fAllParticlesTreeName.substr(fAllParticlesTreeName.find_last_of("\\/")+1,std::string::npos);

    fGoodParticleTree->SetBranchAddress((namePrefix+"_EventID").c_str(),&fBeamEvent);
    fGoodParticleTree->SetBranchAddress((namePrefix+"_TrackID").c_str(),&fTrackID);
 
    ////////************added by Caroline for beam simulation storage for good particles ***************//////////////
    // add lines for the LAG_ENTRY part
    fGoodParticleTree->SetBranchAddress((namePrefix+"_x").c_str(),&fGoodLag_ENTRY_x);
    fGoodParticleTree->SetBranchAddress((namePrefix+"_y").c_str(),&fGoodLag_ENTRY_y);
    fGoodParticleTree->SetBranchAddress((namePrefix+"_z").c_str(),&fGoodLag_ENTRY_z);
    fGoodParticleTree->SetBranchAddress((namePrefix+"_PDGid").c_str(),&fGoodLagPDGID);
    fGoodParticleTree->SetBranchAddress((namePrefix+"_Px").c_str(),&fGoodLagPx);
    fGoodParticleTree->SetBranchAddress((namePrefix+"_Py").c_str(),&fGoodLagPy);
    fGoodParticleTree->SetBranchAddress((namePrefix+"_Pz").c_str(),&fGoodLagPz);
    
    // add more lines for the trig 2 part
    fGoodParticleTree->SetBranchAddress("TRIG2_t",&fTriggerT);
    fGoodParticleTree->SetBranchAddress("TRIG2_x",&fGoodTRIG2_x);
    fGoodParticleTree->SetBranchAddress("TRIG2_y",&fGoodTRIG2_y);
    fGoodParticleTree->SetBranchAddress("TRIG2_z",&fGoodTRIG2_z);
    fGoodParticleTree->SetBranchAddress("TRIG2_PDGid",&fGoodTRIG2_PDGid);
    fGoodParticleTree->SetBranchAddress("TRIG2_Px",&fGoodTRIG2_Px);
    fGoodParticleTree->SetBranchAddress("TRIG2_Py",&fGoodTRIG2_Py);
    fGoodParticleTree->SetBranchAddress("TRIG2_Pz",&fGoodTRIG2_Pz);
    fGoodParticleTree->SetBranchAddress("TRIG2_EventID",&fGoodTRIG2_EventID);
    fGoodParticleTree->SetBranchAddress("TRIG2_TrackID",&fGoodTRIG2_TrackID);
    
    //add more lines for the BPROF4 part
    fGoodParticleTree->SetBranchAddress("BPROF4_x",&fGoodBPROF4_x);
    fGoodParticleTree->SetBranchAddress("BPROF4_y",&fGoodBPROF4_y);
    fGoodParticleTree->SetBranchAddress("BPROF4_z",&fGoodBPROF4_z);
    fGoodParticleTree->SetBranchAddress("BPROF4_PDGid",&fGoodBPROF4_PDGid);
    fGoodParticleTree->SetBranchAddress("BPROF4_Px",&fGoodBPROF4_Px);
    fGoodParticleTree->SetBranchAddress("BPROF4_Py",&fGoodBPROF4_Py);
    fGoodParticleTree->SetBranchAddress("BPROF4_Pz",&fGoodBPROF4_Pz);
    fGoodParticleTree->SetBranchAddress("BPROF4_EventID",&fGoodBPROF4_EventID);
    fGoodParticleTree->SetBranchAddress("BPROF4_TrackID",&fGoodBPROF4_TrackID);
    
    //************************************end of caroline's beam particle tree******************************/////////////
        
    // Calculate the number of events to overlay
    CalculateNOverlays();
    
    // Now we need to fill the particle map
    FillParticleMaps();
}

//----------------------------------------------------------------------------------------

void evgen::ProtoDUNEBeam::beginRun(art::Run& run)
{
    // Grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;
    std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));
    run.put(std::move(runcol));
}

//--------------------------------------------------------------------------------------------

void evgen::ProtoDUNEBeam::endJob(){
    fInputFile->Close();
}


//--------------------------------------------------------------------------------------------
void evgen::ProtoDUNEBeam::produce(art::Event & e)
{
    
    // Define the truth collection for this event.
    auto truthcol = std::make_unique< std::vector<simb::MCTruth> >();
    
    //------------ Added by Caroline for beam simulation storage----------------------------
    std::unique_ptr<std::vector<sim::ProtoDUNEbeamsim>> beamsimcol (new std::vector<sim::ProtoDUNEbeamsim>);
    //std::unique_ptr<art::Assns<sim::ProtoDUNEbeamsim, simb::MCTruth> > beamsimassn (new art::Assns<sim::ProtoDUNEbeamsim, simb::MCTruth>);
    simb::MCTruth truth;
    
    // Fill the MCTruth object
    GenerateTrueEvent(truth, (*beamsimcol) );
    
    //call the event generation fuction to obtain the values
    std::cout<<"the size of *beamsimcol: "<<(*beamsimcol).size()<<std::endl;
    
    ///get your vector of mc particles ...loop over  get the track id inside the loop, loop through data product...track ID -ASSOSCIATED
    
    // Add the MCTruth to the vector
    truthcol->push_back(truth);
    
    //Make the assn                                                                                                  
    // util::CreateAssn(*this, e, *beamsimcol, truth, *beamsimassn);
    // Finally, add the MCTruth to the event
    e.put(std::move(truthcol));
    e.put(std::move(beamsimcol));
    //puts the vector object on to each event
    // We have made our event, increment the event number.
    ++fEventNumber;
}
//--------------------------------------------------------------------------------------

// Fill the particle maps using the input files. This links the events of interest
// to the entry number in fAllParticlesTree.
void evgen::ProtoDUNEBeam::FillParticleMaps(){
    
    // First off, loop over the good particles tree.
    int goodEventCounter = 0;
    for(int i = 0; i < fGoodParticleTree->GetEntries(); ++i){
        // If we want to skip some events, make sure we don't bother reading them in.
        if(fStartEvent > goodEventCounter){
            ++goodEventCounter;
            continue;
        }
        else{
            ++goodEventCounter;
        }
        
        fGoodParticleTree->GetEntry(i);
        int event = (int)fBeamEvent;
        
        // Initialise the event - particle map. This will be filled
        // in the next loop.
        if(fEventParticleMap.find(event) == fEventParticleMap.end()){
            std::vector<int> tempVec; // This will be the vector of track ids for each event near to the good event.
            std::pair<int,std::vector<int> > tempPair = std::make_pair(event,tempVec);
            std::vector<std::pair<int,std::vector<int> > > tempMainVec;
            tempMainVec.push_back(tempPair);
            fEventParticleMap.insert(std::make_pair(event,tempMainVec));
            fGoodEventList.push_back(event);
        }
        
        
        fGoodParticleTriggerTime.insert(std::make_pair(event,fTriggerT));
        
        fGoodParticleTrackID.insert(std::make_pair(event,fTrackID));
    }
    
    // Print a message in case a user starts thinking something has broken.
    mf::LogInfo("ProtoDUNEBeam") << "About to loop over the beam simulation tree, this could take some time.";
    
    // Now we need to loop over the main particle tree
    for(int i = 0; i < fAllParticlesTree->GetEntries(); ++i){
        fAllParticlesTree->GetEntry(i);
        
        if (i%100000==0) std::cout << "Looking at entry " << i << std::endl;
        
        int event = int(fBeamEvent);
       
        // Look at which good events this should be overlaid with
        std::vector<int> goodEventList = GetAllOverlays(event,fOverlays);

        for(auto const goodEvent : goodEventList){
            // Ignore this good event if for some reason it avoided the map
            if(fEventParticleMap.find(goodEvent) == fEventParticleMap.end()) continue;

            // Store the index of this event so that we can quickly access
            // it later when building events
            std::vector<std::pair<int, std::vector<int> > > tracksForEvents = fEventParticleMap[goodEvent];
            bool foundEvent = false;
            unsigned int element = 0;
            for(unsigned int v = 0; v < tracksForEvents.size(); ++v){
                if(tracksForEvents[v].first == event){
                    foundEvent = true;
                    element = v;
                    break;
                }
            }
            if(foundEvent){
                fEventParticleMap[goodEvent][element].second.push_back(i);
            }
            else{
                std::vector<int> newVec;
                newVec.push_back(i);
                std::pair<int,std::vector<int> > newEvent = std::make_pair(event,newVec);
                fEventParticleMap[goodEvent].push_back(newEvent);
            }
        } // End loop over matching events to overlay (this re-uses beam interactions...)
    } // End loop over the main tree.
    
    mf::LogInfo("ProtoDUNEBeam") << "Found " << fGoodEventList.size() << " good events containing " << goodEventCounter << " good particles.";
    mf::LogInfo("ProtoDUNEBeam") << "All maps built, beginning event generation.";
    
}



//-------------------------------------------------------------------------------------------------
//modified by Caroline for beam sim storage
void evgen::ProtoDUNEBeam::GenerateTrueEvent(simb::MCTruth &mcTruth, std::vector<sim::ProtoDUNEbeamsim> &beamsimcol){
    //  std::unique_ptr<std::vector<sim::ProtoDUNEbeamsim>> beamsimcol
    // Check we haven't exceeded the length of the input tree
    if(fEventNumber >= (int)fGoodEventList.size()){
        throw cet::exception("ProtoDUNEBeam") << "Requested entry " << fEventNumber
        << " but tree only has entries 0 to "
        << fGoodEventList.size() - 1 << std::endl;
    } //end of if statement
    
    // Get the list of entries for the current event
    int beamEvent = fGoodEventList[fCurrentGoodEvent];
    
    // Get the random number generator service and make some CLHEP generators
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine("protoDUNEBeam");
    CLHEP::RandFlat flatRnd(engine);
    
    // A single particle seems the most accurate description.
    mcTruth.SetOrigin(simb::kSingleParticle);
    
    // Find the entries that we are interested in.
    //	std::cout << "Finding all particles associated with good particle event " << beamEvent << std::endl;
    for(auto const e : fEventParticleMap[beamEvent]){
        // Is this the event we would have triggered on?
        bool trigEvent = (e.first == beamEvent);
        float baseTime;
        if(trigEvent){
            // Set the base time for the triggered event equal to the negative of the good particle time.
            // This will be corrected later on to set the time to zero, but keep time offsets within the event.
            baseTime = -1.0 * fGoodParticleTriggerTime[beamEvent];
        }
        else{
            // Get a random time from -fReadoutWindow to +fReadoutWindow in ns (fReadoutWindow value is in ms).
            baseTime = (flatRnd.fire() - 0.5)*2.0*(fReadoutWindow*1000.*1000.);
        }
        for(auto const t : e.second){
            // Get the entry from the tree for this event and track.
            fAllParticlesTree->GetEntry(t);
            
            // Convert the pdgCode to an int
            int intPDG = (int)fPDG;
            // We need to ignore nuclei for now...
            if(intPDG > 100000) continue;
            
            // Check to see if this should be a primary beam particle (good particle) or beam background
            std::string process="primaryBackground";
            // If this track is a "good particle", use the usual "primary" tag
            if(trigEvent && (fGoodParticleTrackID[beamEvent] == (int)fTrackID)){
                process="primary";
            }
            
            // Get the position four vector, converting mm to cm
            TLorentzVector pos = ConvertCoordinates(fX/10.,fY/10.,fZ/10.,baseTime + fEntryT);
            // Get momentum four vector, remembering to convert MeV to GeV
            TLorentzVector mom = MakeMomentumVector(fPx/1000.,fPy/1000.,fPz/1000.,intPDG);
            
            // Apply an energy cut to things outside the beam pipe to represent the shielding.
            float r = sqrt(fX*fX + fY*fY);
            r = r / 10.;
            if(r > fBackgroundRadiusCut && mom.Vect().Mag() < fBackgroundMomentumCut){
                continue;
            }
            
            //if(process=="primary") std::cout << " - Got the good particle (" << intPDG << ", " << e.first << ", " << (int)fTrackID << ") with momentum = " << mom.Vect().Mag() << std::endl;
            
            // Track ID needs to be negative for primaries
            int trackID = -1*(mcTruth.NParticles() + 1); //g4trkid in larsoft
            
            // Create the particle and add the starting position and momentum
            simb::MCParticle newParticle(trackID,intPDG,process);
            newParticle.AddTrajectoryPoint(pos,mom);
            
            // Add the MCParticle to the MCTruth for the event.
            mcTruth.Add(newParticle);

            //////------caroline's new code-------
	    //Make the assn 
	    //util::CreateAssn(*this, e, *beamsimcol, newParticle, *beamsimassn);
            if(trigEvent && (fGoodParticleTrackID[beamEvent] == (int)fTrackID)){
                // process="primary";
                
                int EarlierTrackID = fTrackID;
                for (int i =0; i<fGoodParticleTree->GetEntries();++i){
                    fGoodParticleTree->GetEntry(i);
                    if ((int)fTrackID == EarlierTrackID){
                        
                        
                        sim::ProtoDUNEbeamsim temp (fGoodBPROF4_x,fGoodBPROF4_y,fGoodBPROF4_z,fGoodBPROF4_Px,fGoodBPROF4_Py,fGoodBPROF4_Pz,fGoodBPROF4_PDGid,fGoodBPROF4_EventID,fGoodBPROF4_TrackID,fGoodTRIG2_x,fGoodTRIG2_y,fGoodTRIG2_z,fGoodTRIG2_Px,fGoodTRIG2_Py,fGoodTRIG2_Pz,fGoodTRIG2_EventID,fGoodTRIG2_TrackID,fGoodLag_ENTRY_x,fGoodLag_ENTRY_y,fGoodLag_ENTRY_z,fGoodLagPx,fGoodLagPy,fGoodLagPz,fBeamEvent,fTrackID);
                        beamsimcol.push_back(temp);
                        // std::cout<<" test value beam profile monitor: TTREE   "<<fGoodBPROF4_x<<std::endl;
                        std::cout<<"From TTree TRIG2_TRACKID: "<<fGoodTRIG2_TrackID<<std::endl;
                        
                        // std::cout<<"the testing for beam profile monitor information:  "<<fGoodBPROF4_z<<std::endl;
                        
                        
                        std::cout<< "From the data product:  TRIG2TRACKID:   "<<temp.get_TRIG2_TrackID()<<std::endl;
                        //check the last index of the vector
                        sim::ProtoDUNEbeamsim lastelement = beamsimcol.back();
                        
                        std::cout<<"From the vector TRIG2_TRACKID: "<<lastelement.get_TRIG2_TrackID()<<std::endl;
                        
                        
                        //Make the assn                                                                                                                                  
			//util::CreateAssn(*this, e, *beamsimcol, newParticle, *beamsimassn)
                        
                    }
                }
            }
            
            //------------caroline added finish---
            
        } // End loop over interesting tracks for each event
    } // End loop over the vector of interesting events
    
    mf::LogInfo("ProtoDUNEBeam") << "Created event with " << mcTruth.NParticles() << " particles.";
    std::cout << "Created event with " << mcTruth.NParticles() << " particles.";
    
    // Move on the good event iterator
    ++fCurrentGoodEvent;
}

//---------------------------------------------------------------------------------------

// Function written in similar way as "openDBs()" in CORSIKAGen_module.cc
void evgen::ProtoDUNEBeam::OpenInputFile()
{
    // Setup ifdh object
    if (!fIFDH)
    {
        fIFDH = new ifdh_ns::ifdh;
    }
    
    const char* ifdh_debug_env = std::getenv("IFDH_DEBUG_LEVEL");
    if ( ifdh_debug_env )
    {
        mf::LogInfo("ProtoDUNEBeam") << "IFDH_DEBUG_LEVEL: " << ifdh_debug_env<<"\n";
        fIFDH->set_debug(ifdh_debug_env);
    }
    
    std::string path(gSystem->DirName(fFileName.c_str()));
    std::string pattern(gSystem->BaseName(fFileName.c_str()));
    
    auto flist = fIFDH->findMatchingFiles(path,pattern);
    if (flist.empty())
    {
        struct stat buffer;
        if (stat(fFileName.c_str(), &buffer) != 0)
        {
            throw cet::exception("ProtoDUNEBeam") << "No files returned for path:pattern: "<<path<<":"<<pattern<<std::endl;
        }
        else
        {
            mf::LogInfo("ProtoDUNEBeam") << "For "<< fFileName <<"\n";
        }
    }
    else
    {
        std::pair<std::string, long> f = flist.front();
        
        mf::LogInfo("ProtoDUNEBeam") << "For "<< fFileName <<"\n";
        
        // Do the fetching, store local filepaths in locallist
        
        mf::LogInfo("ProtoDUNEBeam")
        << "Fetching: " << f.first << " " << f.second <<"\n";
        std::string fetchedfile(fIFDH->fetchInput(f.first));
        LOG_DEBUG("ProtoDUNEBeam") << " Fetched; local path: " << fetchedfile;
        
        fFileName = fetchedfile;
    }
}


//----------------------------------------------------------------------------------

TLorentzVector evgen::ProtoDUNEBeam::ConvertCoordinates(float x, float y, float z, float t){
    
    float finalX = x + fBeamX;
    float finalY = y + fBeamY;
    float finalZ = (z - z) + fBeamZ; // Just use the z position
    
    TLorentzVector newPos(finalX,finalY,finalZ,t);
    return newPos;
}

//--------------------------------------------------------------------------------

TLorentzVector evgen::ProtoDUNEBeam::MakeMomentumVector(float px, float py, float pz, int pdg){
    
    float rotationXZ = fRotateXZ;
    float rotationYZ = fRotateYZ;
    
    // Make the momentum vector and rotate it
    TVector3 momVec(px,py,pz);
    momVec.RotateY(rotationXZ * TMath::Pi() / 180.);
    momVec.RotateX(rotationYZ * TMath::Pi() / 180.);
    
    // Find the particle mass so we can form the energy
    const TDatabasePDG* databasePDG = TDatabasePDG::Instance();
    const TParticlePDG* definition = databasePDG->GetParticle(pdg);
    float mass = definition->Mass();
    
    float energy = sqrt(mass*mass + momVec.Mag2());
    
    TLorentzVector newMom(momVec,energy);
    return newMom;
}


//-----------------------------------------------------------------------------

void evgen::ProtoDUNEBeam::CalculateNOverlays(){
    
    // The number of events to overlay is as follows:
    // N = Intensity * 2.0 * ReadoutWindow / BeamSpillLength
    fOverlays = fIntensity * (2.0 * fReadoutWindow / 1000.) / fBeamSpillLength;
    
}


//-------------------------------------------------------------------------------
int evgen::ProtoDUNEBeam::IsOverlayEvent(int event, int nOverlay){
    
    // Check if this event lies within nOverlay/2 of each
    for(auto const e : fGoodEventList){
        if(fabs(event - e) < nOverlay/2){
            return e;
        }
    }
    return -1;
}


//---------------------------------------------------------------------------------
std::vector<int> evgen::ProtoDUNEBeam::GetAllOverlays(int event, int nOverlay){
    
    std::vector<int> nMatches;
    for(auto const e : fGoodEventList){
        if(fabs(event - e) < nOverlay/2){
            nMatches.push_back(e);
        }
    }
    return nMatches;
    
}
//----------------------------------------------------------------------------------
DEFINE_ART_MODULE(evgen::ProtoDUNEBeam)

