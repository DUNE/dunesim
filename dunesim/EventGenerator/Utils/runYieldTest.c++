#include "ProtoDUNETriggeredBeamUtils.h"
#include "BeamMiscUtils.h"
#include "TFile.h"
#include "TTree.h"
#include <thread>

int main(int argc, char ** argv){

  std::string fcl_file;
  std::string output_filename;
  std::string mc_file;
  int nworkers = 1;
  // Options to run
  for (int iArg = 1; iArg < argc; iArg++) {
    if (!strcasecmp(argv[iArg],"-c")) {
     fcl_file = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-o")) {
     output_filename = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-n")) {
     nworkers = std::atoi(argv[++iArg]);
    }
    if (!strcasecmp(argv[iArg],"-h")) {
      std::cout << "Usage: runEventCounter -c fclfile.fcl -o output_file [-n nthreads]" << std::endl;
      return 1;
    }
  }

  //Fcl pars
  auto pset = beammisc::GetPars(fcl_file);
  auto generator = beammisc::GetGenerator(pset);
  auto file_list = pset.get<std::vector<std::string>>("Files");

  beammisc::CountConfig config(generator);

  //Make beam utils instance
  evgen::ProtoDUNETriggeredBeamUtils beam_utils(generator);

  ROOT::EnableThreadSafety();//Needed for root

  return 0;
}

void CountFile(std::string filename, beammisc::CountConfig & config, std::ofstream & output_file, evgen::ProtoDUNETriggeredBeamUtils & beam_utils) {
  std::cout << filename << std::endl;
  if (filename.find("/pnfs") != 0) {
    throw cet::exception("ProtoDUNETriggeredBeam") << "Filename " <<
        filename << " does not start with /pnfs as required for streaming" << std::endl;
  }
  filename.replace(0, 5, "root://fndca1.fnal.gov:1094//pnfs/fnal.gov/usr");

  auto * fIn = TFile::Open(filename.c_str());
  auto * tree = static_cast<TTree*>(fIn->Get("NTuples/GoodParticle"));
  if (tree == nullptr) {
    std::cout << "Error. File " << filename << " tree is malformed. Skipping" <<
                 std::endl;
    return;
  }
  std::cout << tree->GetEntries() << std::endl;

  // Fill all potential events from the NP04front tree
  TTree * frontFaceTree =
      static_cast<TTree*>(fIn->Get(config.fNP04frontTreeName.c_str()));
  std::map<int, evgen::BeamEvent> allBeamEvents;
  beam_utils.FillParticleMaps(frontFaceTree, allBeamEvents);

  std::vector<int> triggeredEventIDs = beam_utils.FindTriggeredEvents(
      fIn, config.fTRIG1TreeName, config.fTRIG2TreeName, allBeamEvents);
  std::cout << "Proto trigger list has " << triggeredEventIDs.size() << " events" << std::endl; 
  // For triggered events, we now need to attach the other instrument information
  std::vector<std::string> otherInstrumentTreeNames;
  otherInstrumentTreeNames.push_back(config.fTOF1TreeName.c_str());
  otherInstrumentTreeNames.push_back(config.fBPROF1TreeName.c_str());
  otherInstrumentTreeNames.push_back(config.fBPROF2TreeName.c_str());
  otherInstrumentTreeNames.push_back(config.fBPROF3TreeName.c_str());
  if (!beam_utils.GetIsNP02())
    otherInstrumentTreeNames.push_back(config.fBPROFEXTTreeName.c_str());
  otherInstrumentTreeNames.push_back(config.fBPROF4TreeName.c_str());

  for(const std::string treeName : otherInstrumentTreeNames){ 
    TTree *instrumentTree = (TTree*)fIn->Get(treeName.c_str());
    beam_utils.FillInstrumentInformation(triggeredEventIDs, instrumentTree, allBeamEvents);
    std::cout << " - Finished adding information from " << treeName << std::endl;
  }
  std::cout << "Final trigger list has " << triggeredEventIDs.size() << " events" << std::endl; 
  WriteOut(output_file, filename, triggeredEventIDs.size());

  fIn->Close();
}
