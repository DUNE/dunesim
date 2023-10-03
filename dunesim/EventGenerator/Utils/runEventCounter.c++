#include "ProtoDUNETriggeredBeamUtils.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include <thread>
//#include <ranges>

struct CountConfig {
  CountConfig(fhicl::ParameterSet & generator)
    : fNP04frontTreeName(generator.get<std::string>("NP04frontTreeName")),
      fTOF1TreeName(generator.get<std::string>("TOF1TreeName")),
      fBPROF1TreeName(generator.get<std::string>("BPROF1TreeName")),
      fBPROF2TreeName(generator.get<std::string>("BPROF2TreeName")),
      fBPROF3TreeName(generator.get<std::string>("BPROF3TreeName")),
      fTRIG1TreeName(generator.get<std::string>("TRIG1TreeName")),
      fBPROFEXTTreeName(generator.get<std::string>("BPROFEXTTreeName", "")),
      fBPROF4TreeName(generator.get<std::string>("BPROF4TreeName")),
      fTRIG2TreeName(generator.get<std::string>("TRIG2TreeName")) {};

  std::string fNP04frontTreeName;
  std::string fTOF1TreeName;
  std::string fBPROF1TreeName;
  std::string fBPROF2TreeName;
  std::string fBPROF3TreeName;
  std::string fTRIG1TreeName;
  std::string fBPROFEXTTreeName;
  std::string fBPROF4TreeName;
  std::string fTRIG2TreeName;
};


fhicl::ParameterSet GetGenerator(fhicl::ParameterSet & pset);
fhicl::ParameterSet GetPars(std::string fcl_file);
void WriteOut(std::ofstream & outputfile, std::string filename, int count);
void CountFile(std::string filename, CountConfig & config, std::ofstream & output_file, evgen::ProtoDUNETriggeredBeamUtils & beam_utils);
void Loop(std::vector<std::string> filelist, CountConfig config,
          std::ofstream & output_file, evgen::ProtoDUNETriggeredBeamUtils beam_utils,
          size_t worker_id, std::vector<size_t> n_files);


std::mutex event_count_mutex;

int main(int argc, char ** argv){

  std::string fcl_file;
  std::string output_filename;
  std::string mc_file;
  std::string data_file;
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
  auto pset = GetPars(fcl_file);
  auto generator = GetGenerator(pset);
  auto file_list = pset.get<std::vector<std::string>>("Files");

  CountConfig config(generator);
  CountConfig config_cp(config);
  std::cout << config_cp.fNP04frontTreeName << std::endl;

  //Make beam utils instance
  evgen::ProtoDUNETriggeredBeamUtils beam_utils(generator);

  std::ofstream output_file(output_filename);

  ROOT::EnableThreadSafety();//Needed for root

  //To split up the file list
  int nfiles = file_list.size() / nworkers;
  int remainder = file_list.size() % nworkers;
  std::vector<size_t> file_list_split(nworkers, nfiles);
  while (remainder > nworkers) {
    nfiles = nfiles/nworkers;
    remainder = nfiles % nworkers;

    for (int i = 0; i < nworkers; ++i) {
      file_list_split[i] += nfiles;
    }
  }
  for (int i = 0; i < remainder; ++i) {
    file_list_split[i] += 1;
  }

  //Make threads
  std::vector<std::thread> workers;
  for (size_t workerid = 0; workerid < static_cast<size_t>(nworkers);
       ++workerid) {
    std::thread worker(
        &Loop, file_list, std::ref(config), std::ref(output_file),
        std::ref(beam_utils), workerid, std::ref(file_list_split)
    );

    workers.emplace_back(std::move(worker));
  }
  for (auto &&worker : workers) { worker.join();}

  output_file.close();

  return 0;
}

void WriteOut(std::ofstream & outputfile, std::string filename, int count) {
  std::cout << "Writing out " << filename << " " << std::to_string(count) << std::endl;
  std::lock_guard<std::mutex> guard(event_count_mutex);
  outputfile << filename << " " << std::to_string(count) << "\n";
}

void Loop(std::vector<std::string> file_list, CountConfig config,
          std::ofstream & output_file, evgen::ProtoDUNETriggeredBeamUtils beam_utils,
          size_t worker_id, std::vector<size_t> n_files) {

  size_t start_file = 0;
  for (size_t i = 0; i < worker_id; ++i) {
    start_file += n_files[i];
  }

  size_t end_file = start_file + n_files[worker_id];
  std::cout << "Loop " << worker_id << " " << start_file << " " << end_file << std::endl;

  //for (std::string & filename : file_list) {
  for (size_t i = start_file; i != end_file; ++i) {
    std::string & filename = file_list[i];
    CountFile(filename, config, output_file, beam_utils);
  } 
}

fhicl::ParameterSet GetPars(std::string fcl_file) {
  std::cout << "Fcl file: " << fcl_file << std::endl;
  fhicl::ParameterSet pset;

  // Configuration file lookup policy.
  char const* fhicl_env = getenv("FHICL_FILE_PATH");
  std::string search_path;

  if (fhicl_env == nullptr) {
    std::cerr << "Expected environment variable FHICL_FILE_PATH is missing or empty: using \".\"\n";
    search_path = ".";
  }
  else {
    search_path = std::string{fhicl_env};
  }

  cet::filepath_first_absolute_or_lookup_with_dot lookupPolicy{search_path};

  //fhicl::make_ParameterSet(fcl_file, lookupPolicy, pset);
  pset = fhicl::ParameterSet::make(fcl_file, lookupPolicy);
  return pset;

}

fhicl::ParameterSet GetGenerator(fhicl::ParameterSet & pset) {
  return pset.get<fhicl::ParameterSet>("physics.producers.generator");
}

void CountFile(std::string filename, CountConfig & config, std::ofstream & output_file, evgen::ProtoDUNETriggeredBeamUtils & beam_utils) {
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

  for(const std::string &treeName : otherInstrumentTreeNames){ 
    TTree *instrumentTree = (TTree*)fIn->Get(treeName.c_str());
    beam_utils.FillInstrumentInformation(triggeredEventIDs, instrumentTree, allBeamEvents);
    std::cout << " - Finished adding information from " << treeName << std::endl;
  }
  std::cout << "Final trigger list has " << triggeredEventIDs.size() << " events" << std::endl; 
  WriteOut(output_file, filename, triggeredEventIDs.size());

  fIn->Close();
}
