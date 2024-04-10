#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "dunesim/EventGenerator/Utils/H4BeamFileService.h"

// constructor

dune::H4BeamFileService::H4BeamFileService(
    fhicl::ParameterSet const& p,
    art::ActivityRegistry& areg)
  : fBeamUtils(p) {

  fAllBeamEvents = std::make_unique<std::map<int, evgen::BeamEvent>>();
  fFinalTriggeredEventIDs = std::make_unique<std::vector<int>>();
  fTOF1TreeName      = p.get<std::string>("TOF1TreeName");
  fBPROF1TreeName    = p.get<std::string>("BPROF1TreeName");
  fBPROF2TreeName    = p.get<std::string>("BPROF2TreeName");
  fBPROF3TreeName    = p.get<std::string>("BPROF3TreeName");
  fTRIG1TreeName     = p.get<std::string>("TRIG1TreeName");
  fBPROFEXTTreeName  = p.get<std::string>("BPROFEXTTreeName", "");
  fBPROF4TreeName    = p.get<std::string>("BPROF4TreeName");
  fTRIG2TreeName      = p.get<std::string>("TRIG2TreeName");
  fNP04frontTreeName = p.get<std::string>("NP04frontTreeName");
  fIsNP02 = p.get<bool>("IsNP02", false);

}

// sets the pointer and assumes ownership of it

void dune::H4BeamFileService::OpenFile(const std::string & filename) {

  std::cout << "Opening " << filename << std::endl;
  TFile * inputFile = TFile::Open(filename.c_str());
  std::cout << inputFile << std::endl;
  // Check we have the file
  if(inputFile == 0x0){
      throw cet::exception("H4BeamFileService::OpenFile") <<
          "Input file " << filename << " cannot be read.\n";
  }


  int first_pos = filename.find_last_of("_")+1;
  int last_pos = filename.find_last_of(".");
  std::cout << "Getting run " <<
               last_pos << " " <<
               first_pos << " " <<
               filename.substr(
                 first_pos, last_pos - first_pos) << std::endl;

  fRun = std::stoul(filename.substr(
    first_pos,
    last_pos - first_pos));

  TTree *frontFaceTree = (TTree*)inputFile->Get(fNP04frontTreeName.c_str());
  std::cout << "Front face tree: " << frontFaceTree << std::endl;
  // Fill all potential events from the NP04front tree
  fBeamUtils.FillParticleMaps(frontFaceTree, (*fAllBeamEvents.get()));

  //// Now search for trigger events
  (*fFinalTriggeredEventIDs.get()) = fBeamUtils.FindTriggeredEvents(
      inputFile, fTRIG1TreeName, fTRIG2TreeName, (*fAllBeamEvents.get()));
  std::cout << "Proto trigger list has " << fFinalTriggeredEventIDs->size() <<
               " events" << std::endl;

  // For triggered events, we now need to attach the other instrument information
  std::vector<std::string> otherInstrumentTreeNames;
  otherInstrumentTreeNames.push_back(fTOF1TreeName.c_str());
  otherInstrumentTreeNames.push_back(fBPROF1TreeName.c_str());
  otherInstrumentTreeNames.push_back(fBPROF2TreeName.c_str());
  otherInstrumentTreeNames.push_back(fBPROF3TreeName.c_str());
  if (!fIsNP02)
    otherInstrumentTreeNames.push_back(fBPROFEXTTreeName.c_str());
  otherInstrumentTreeNames.push_back(fBPROF4TreeName.c_str());

  for(const std::string &treeName : otherInstrumentTreeNames){ 
    TTree *instrumentTree = (TTree*)inputFile->Get(treeName.c_str());
    fBeamUtils.FillInstrumentInformation((*fFinalTriggeredEventIDs.get()),
                                         instrumentTree,
                                         (*fAllBeamEvents.get()));
    std::cout << " - Finished adding information from " << treeName <<
                 std::endl;
  }
  std::cout << "Final trigger list has " <<
               fFinalTriggeredEventIDs->size() << " events" << std::endl;



  inputFile->Close();
}

void dune::H4BeamFileService::Reset() {
  if (fAllBeamEvents) 
    fAllBeamEvents->clear();
  if (fFinalTriggeredEventIDs)
    fFinalTriggeredEventIDs->clear();
}


DEFINE_ART_SERVICE(dune::H4BeamFileService)
