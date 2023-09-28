#include "ProtoDUNETriggeredBeamUtils.h"
#include "BeamMiscUtils.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include <thread>

struct tree_buffer {
  int trigger_pdg;
  int trigger_id;
  int interacted;
  double trigger_front_x, trigger_front_y, trigger_front_z;
  double trigger_front_px, trigger_front_py, trigger_front_pz, trigger_front_p;

  std::vector<int> overlays;

  std::vector<int> background_pdgs, background_ids;
  std::vector<double> background_front_x, background_front_y, background_front_z,
                      background_front_px, background_front_py,
                      background_front_pz, background_front_p;
  void Clear() {
    background_pdgs.clear();
    background_ids.clear();
    background_front_x.clear();
    background_front_y.clear();
    background_front_z.clear();
    background_front_px.clear();
    background_front_py.clear();
    background_front_pz.clear();
    background_front_p.clear();
    overlays.clear();
  };

  void AddTrigger(const evgen::BeamParticle & part) {
    trigger_pdg = part.fPDG;
    trigger_id = part.fTrackID;
    trigger_front_x = part.fPosX;
    trigger_front_y = part.fPosY;
    trigger_front_z = part.fPosZ;
    trigger_front_px = part.fMomX;
    trigger_front_py = part.fMomY;
    trigger_front_pz = part.fMomZ;
    trigger_front_p = sqrt(
      part.fMomX*part.fMomX + part.fMomY*part.fMomY + part.fMomZ*part.fMomZ);
  };

  void AddBG(const evgen::BeamParticle & part) {
    background_pdgs.push_back(part.fPDG);
    background_ids.push_back(part.fTrackID);
    background_front_x.push_back(part.fPosX);
    background_front_y.push_back(part.fPosY);
    background_front_z.push_back(part.fPosZ);
    background_front_px.push_back(part.fMomX);
    background_front_py.push_back(part.fMomY);
    background_front_pz.push_back(part.fMomZ);
    background_front_p.push_back(
      sqrt(part.fMomX*part.fMomX + part.fMomY*part.fMomY +
           part.fMomZ*part.fMomZ));
  };
};


void WriteOut(std::ofstream & outputfile, std::string filename, int count);
void CountFile(std::string filename, beammisc::CountConfig & config, std::ofstream & output_file, TTree & tree, tree_buffer & vals, evgen::ProtoDUNETriggeredBeamUtils & beam_utils);
void Loop(std::vector<std::string> filelist, beammisc::CountConfig config,
          std::ofstream & output_file, TTree & tree, evgen::ProtoDUNETriggeredBeamUtils beam_utils,
          size_t worker_id, std::vector<size_t> n_files);


std::mutex event_count_mutex;

int main(int argc, char ** argv){

  std::string fcl_file;
  std::string output_filename;
  std::string output_root_filename;
  std::string mc_file = "";
  int nworkers = 1;
  // Options to run
  for (int iArg = 1; iArg < argc; iArg++) {
    if (!strcasecmp(argv[iArg],"-c")) {
     fcl_file = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-i")) {
     mc_file = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-o")) {
     output_filename = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-n")) {
     nworkers = std::atoi(argv[++iArg]);
    }
    if (!strcasecmp(argv[iArg], "-t")) {
      output_root_filename = argv[++iArg];
    }
    if (!strcasecmp(argv[iArg],"-h")) {
      std::cout << "Usage: runEventCounter -c fclfile.fcl -o output_file -t output_root_file [-n nthreads]" << std::endl;
      return 1;
    }
  }

  //Fcl pars
  auto pset = beammisc::GetPars(fcl_file);
  auto generator = beammisc::GetGenerator(pset);
  std::vector<std::string> file_list;
  if (mc_file != "") {
    file_list.push_back(mc_file);
  }
  else {
    file_list = pset.get<std::vector<std::string>>("Files");
  }

  //bool skip_gamma = pset.get<bool>("SkipGamma", false);

  beammisc::CountConfig config(generator);
  //CountConfig config_cp(config);
  //std::cout << config_cp.fNP04frontTreeName << std::endl;

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

  //Make Output Root File and Tree
  //TFIile root_output_file(output_root_filename.c_str(), "recreate");
  std::vector<TTree *> trees;
  std::vector<TFile *> root_files;

  for (size_t workerid = 0; workerid < static_cast<size_t>(nworkers);
       ++workerid) {
    TFile * root_output_file
        = new TFile(
            TString::Format("%s_%zu.root", output_root_filename.c_str(), workerid),
            "recreate");
    TTree * tree = new TTree("tree", "");
    trees.push_back(tree);
    root_files.push_back(root_output_file);
  }
  

  //Make threads
  std::vector<std::thread> workers;
  for (size_t workerid = 0; workerid < static_cast<size_t>(nworkers);
       ++workerid) {
    std::thread worker(
        &Loop, file_list, std::ref(config), std::ref(output_file),
        std::ref(*trees[workerid]),
        std::ref(beam_utils), workerid, std::ref(file_list_split)
    );

    workers.emplace_back(std::move(worker));
  }
  for (auto &&worker : workers) { worker.join();}

  for (size_t i = 0; i < static_cast<size_t>(nworkers); ++i) {
    root_files[i]->cd();
    trees[i]->Write();
    root_files[i]->Close();
  }

  output_file.close();

  return 0;
}

void WriteOut(std::ofstream & outputfile, std::string filename, int count) {
  std::cout << "Writing out " << filename << " " << std::to_string(count) << std::endl;
  std::lock_guard<std::mutex> guard(event_count_mutex);
  outputfile << filename << " " << std::to_string(count) << "\n";
}

void Loop(std::vector<std::string> file_list, beammisc::CountConfig config,
          std::ofstream & output_file, TTree & tree,
          evgen::ProtoDUNETriggeredBeamUtils beam_utils, size_t worker_id, std::vector<size_t> n_files) {


  size_t start_file = 0;
  for (size_t i = 0; i < worker_id; ++i) {
    start_file += n_files[i];
  }

  size_t end_file = start_file + n_files[worker_id];
  std::cout << "Loop " << worker_id << " " << start_file << " " << end_file << std::endl;


  //int trigger_pdg;
  //double trigger_front_x, trigger_front_y, trigger_front_z;
  //double trigger_front_px, trigger_front_py, trigger_front_pz;
  tree_buffer vals;
  tree.Branch("trigger_pdg", &vals.trigger_pdg, "trigger_pdg/I");
  tree.Branch("trigger_id", &vals.trigger_id, "trigger_id/I");
  tree.Branch("interacted", &vals.interacted, "interacted/I");
  tree.Branch("trigger_front_x", &vals.trigger_front_x, "trigger_front_x/D");
  tree.Branch("trigger_front_y", &vals.trigger_front_y, "trigger_front_y/D");
  tree.Branch("trigger_front_z", &vals.trigger_front_z, "trigger_front_z/D");
  tree.Branch("trigger_front_px", &vals.trigger_front_px, "trigger_front_px/D");
  tree.Branch("trigger_front_py", &vals.trigger_front_py, "trigger_front_py/D");
  tree.Branch("trigger_front_pz", &vals.trigger_front_pz, "trigger_front_pz/D");
  tree.Branch("trigger_front_p", &vals.trigger_front_p, "trigger_front_p/D");
  tree.Branch("overlays", &vals.overlays);

  tree.Branch("background_pdgs", &vals.background_pdgs);
  tree.Branch("background_ids", &vals.background_ids);
  tree.Branch("background_front_x", &vals.background_front_x);
  tree.Branch("background_front_y", &vals.background_front_y);
  tree.Branch("background_front_z", &vals.background_front_z);
  tree.Branch("background_front_px", &vals.background_front_px);
  tree.Branch("background_front_py", &vals.background_front_py);
  tree.Branch("background_front_pz", &vals.background_front_pz);
  tree.Branch("background_front_p", &vals.background_front_p);

  for (size_t i = start_file; i != end_file; ++i) {
    std::string & filename = file_list[i];
    CountFile(filename, config, output_file, tree, vals, beam_utils);
  }

  //tree.Write();
  //root_output_file.Close();
}

void CountFile(std::string filename, beammisc::CountConfig & config,
               std::ofstream & output_file, TTree & output_tree,
               tree_buffer & vals,
               evgen::ProtoDUNETriggeredBeamUtils & beam_utils) {
  std::cout << filename << std::endl;
  if (filename.find("/pnfs") == 0) {
    //throw cet::exception("ProtoDUNETriggeredBeam") << "Filename " <<
    //    filename << " does not start with /pnfs as required for streaming" << std::endl;
    filename.replace(0, 5, "root://fndca1.fnal.gov:1094//pnfs/fnal.gov/usr");
  }
  //filename.replace(0, 5, "root://fndca1.fnal.gov:1094//pnfs/fnal.gov/usr");

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


  //Here -- go through the triggered events
  //Each one will be overlaid with some others
  //The only 'good' particle is the primary triggered one
  //The rest of the front particles + the overlay particles are all
  //background
  for (const int & trigger_id : triggeredEventIDs) {
    std::cout << trigger_id << std::endl;

    auto overlay_event = beam_utils.GenerateOverlaidEvent(
      trigger_id, allBeamEvents, config.fOverlays);

    std::cout << "N overlays: " << overlay_event.fOverlayEventIDs.size() <<
                 std::endl;

    const auto & event = allBeamEvents.at(trigger_id);
    std::cout << "Triggered: " << event.fTriggerID << std::endl;

    evgen::BeamParticle trigParticle;
    if(!event.fHasInteracted){
      trigParticle = event.fParticlesFront.at(event.fTriggerID);
    }
    else{
      const std::string trig2Name = config.fTRIG2TreeName.substr(config.fTRIG2TreeName.find("/")+1);
      trigParticle = event.fTriggeredParticleInfo.at(trig2Name);
    }
    vals.Clear();
    vals.AddTrigger(trigParticle);
    vals.interacted = event.fHasInteracted;
    vals.overlays.insert(
        vals.overlays.end(),
        overlay_event.fOverlayEventIDs.begin(),
        overlay_event.fOverlayEventIDs.end());
    //vals.trigger_pdg = trigParticle.fPDG;
    //vals.trigger_id = event.fTriggerID; 

    for (auto & eventid : overlay_event.fOverlayEventIDs) {
      const auto & event = allBeamEvents.at(eventid);
      auto & front_particles = event.fParticlesFront;
      std::cout << front_particles.size() << std::endl;
      for (const auto & part : front_particles) {
        //if (part.second.fPDG == 22 /*&& skip_gamma*/) continue;
        std::cout << "\t";
        part.second.Print();
        //std::cout << part.second << std::endl;  
        if (part.second.fTrackID == event.fTriggerID) continue;
        vals.AddBG(part.second);
        std::cout << vals.background_pdgs.size() << std::endl;
      }
    }
    output_tree.Fill();
  }

  fIn->Close();
}

