#include <TTree.h>
#include <TFile.h>
#include <iostream>
#include <map>
#include "dunecore/DuneObj/SimpleBeamParticle.h"
bool verbose = false;

// class DetectorHolder {

// public:
//     float x{0};
//     float y{0};
//     float z{0};
//     float Px{0};
//     float Py{0};
//     float Pz{0};
//     float t{0};
//     float PDGid{0};
//     float EventID{0};
//     float TrackID{0};
//     float ParentID{0};
//     int current_event{0};
//     bool done{false};
//     std::map<std::string, float *> branches{
//         {"x", &x},
//         {"y", &y},
//         {"z", &z},
//         {"Px", &Px},
//         {"Py", &Py},
//         {"Pz", &Pz},
//         {"t", &t},
//         {"PDGid", &PDGid},
//         {"EventID", &EventID},
//         {"TrackID", &TrackID},
//         {"ParentID", &ParentID}
//     };
//     // TTree * fTree = 0x0;
//     void SetTree(TTree * tree) {
//         // fTree = tree;
//         for (auto & [key, ptr] : branches) {
//             tree->SetBranchAddress(key.c_str(), ptr);
//         }
//     };

//     friend std::ostream& operator<<(std::ostream& os, const DetectorHolder& obj) {
//         os << "TrackID: " << obj.TrackID << "\n" <<
//               "PDGid: " << obj.PDGid << "\n" <<
//               "EventID: " << obj.EventID << "\n" <<
//               "ParentID: " << obj.ParentID << "\n" <<
//               "x: " << obj.x << "\n" <<
//               "y: " << obj.y << "\n" <<
//               "z: " << obj.z << "\n" <<
//               "t: " << obj.t << "\n" <<
//               "Px: " << obj.Px << "\n" <<
//               "Py: " << obj.Py << "\n" <<
//               "Pz: " << obj.Pz << "\n";
//         return os;
//     }
// };

// struct SimpleBeamParticle {


//     std::vector<float> x, y, z, Px, Py, Pz, t;
//     std::vector<std::string> detectors;
//     int PDGid{0}, EventID{0}, TrackID{0}, ParentID{0};

//     SimpleBeamParticle(DetectorHolder * holder, std::string name) {
//         TrackID = static_cast<int>(holder->TrackID);
//         PDGid = static_cast<int>(holder->PDGid);
//         EventID = static_cast<int>(holder->EventID);
//         ParentID = static_cast<int>(holder->ParentID);

//         detectors.push_back(name);
//         x.push_back(holder->x);
//         y.push_back(holder->y);
//         z.push_back(holder->z);
//         t.push_back(holder->t);
//         Px.push_back(holder->Px);
//         Py.push_back(holder->Py);
//         Pz.push_back(holder->Pz);
//     }

//     void merge(DetectorHolder * holder, std::string name, int overall_event) {
//         if (PDGid != 22 &&
//             ((TrackID != static_cast<int>(holder->TrackID)) ||
//             (PDGid != static_cast<int>(holder->PDGid)) ||
//             (EventID != static_cast<int>(holder->EventID)) || 
//             (ParentID != static_cast<int>(holder->ParentID)))
//             /*std::find(detectors.begin(), detectors.end(), name) != detectors.end()*/) {
//                 std::cerr << "(" << overall_event << ") This part: " << *this << "\nInput: (" << name << ") " << (*holder) << std::endl; 
//                 throw std::runtime_error("ATTEMPTING TO MERGE INVALID PARTS");
//             }
//         detectors.push_back(name);
//         x.push_back(holder->x);
//         y.push_back(holder->y);
//         z.push_back(holder->z);
//         t.push_back(holder->t);
//         Px.push_back(holder->Px);
//         Py.push_back(holder->Py);
//         Pz.push_back(holder->Pz);
//     }

//     friend std::ostream& operator<<(std::ostream& os, const SimpleBeamParticle& obj) {
//         os << "\tTrackID: " << obj.TrackID << "\n" <<
//               "\tPDGid: " << obj.PDGid << "\n" <<
//               "\tEventID: " << obj.EventID << "\n" <<
//               "\tParentID: " << obj.ParentID << "\n";
//               for (size_t i = 0; i < obj.detectors.size(); ++i) {
//                 os << "\t\t" << obj.detectors[i] << " " << obj.x[i] << " " << obj.y[i] << " " << obj.z[i] << " " << obj.t[i] << " " << obj.Px[i] << " " << obj.Py[i] << " " << obj.Pz[i] << "\n";
//               }
//         return os;
//     }

//     SimpleBeamParticle(){}
// };

bool ProcessDetector(TTree * tree, std::map<int, SimpleBeamParticle> & particle_map, DetectorHolder * holder, int overall_event, std::string name) {

    if (holder->done) {
        // std::cout << "Past " << name << std::endl;
        return true;
    }

    int i = holder->current_event;
    for (; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (verbose)
            std::cout << "Got event " << i << " in " << name << std::endl;
        
        if (holder->EventID > overall_event) {
            if (verbose) std::cout << "\t" << holder->EventID << " Past current event " << overall_event << std::endl;
            holder->current_event = i;
            break;
        }
        
        if (verbose) std::cout << *holder << std::endl;
        if (particle_map.find(holder->TrackID) == particle_map.end()) {
            particle_map.insert({holder->TrackID, SimpleBeamParticle(holder, name)});
        }
        else {
            particle_map[holder->TrackID].merge(holder, name, overall_event);
        }
    }
    if (i >= tree->GetEntries()) {
        holder->done = true;
    }


    if (holder->done) {
        // std::cout << "Past " << name << std::endl;
        return true;
    }

    return false;

}

int main(int argc, char ** argv) {
// int reformat_g4bl(std::string input_file, std::string output_file, int n=-1) {
    std::string input_file{""}, output_file{""};
    int start = std::numeric_limits<int>::max();
    int max_events{std::numeric_limits<int>::max()};
    bool verbose = false;
    for (int iarg = 1; iarg < argc; ++iarg) {
        if (strcmp(argv[iarg], "-i") == 0) {
            input_file = argv[++iarg];
            std::cout << "INPUT FILE " << input_file << std::endl;
        }
        else if (strcmp(argv[iarg], "-o") == 0) {
            output_file = argv[++iarg];
            std::cout << "OUTPUT FILE " << output_file << std::endl;
        }
        else if (strcmp(argv[iarg], "-n") == 0) {
            int these_max_events = std::atoi(argv[++iarg]);
            if (these_max_events > 0) max_events = these_max_events;
            std::cout << "Max events: " << max_events << std::endl;
        }
        else if (strcmp(argv[iarg], "-v") == 0) {
            std::cout << "RUNNING VERBOSELY" << std::endl;
            verbose = true;
        }
        else if (strcmp(argv[iarg], "--start") == 0) {
            start = std::atoi(argv[++iarg]);
            std::cout << "START AT " << start << std::endl;
        }
    }

    //Define the set of 'detectors' (aka devices) through which particles may travel
    //std::string first_detector = "BeforeTarget";
    std::vector<std::string> next_detector_list = {
        "VirtualDetector/BeforeTarget",
        "VirtualDetector/AfterTarget",
        "VirtualDetector/TOF1",
        "VirtualDetector/COLL1",
        "VirtualDetector/BPROF1",
        "VirtualDetector/BPROF2",
        "VirtualDetector/BPROF3",
        "VirtualDetector/TRIG1",
        "VirtualDetector/BPROFEXT",
        "VirtualDetector/BPROF4",
        "VirtualDetector/TRIG2",
        "Detector/NP04front"
    };

    std::cout << "Will process\n";
    for (const auto & s : next_detector_list) std::cout << "\t" << s << std::endl;

    //Open the input file and grab the first detector
    TFile * f = TFile::Open(input_file.c_str());
    //TTree * first_tree = (TTree*)f->Get(first_detector.c_str());

    //Make a representation of it and set the branches
    //DetectorHolder before_target;
    //before_target.SetTree(first_tree);
    //first_tree->GetEntry(0);

    //Also make reps for each other tree
    std::map<std::string, TTree*> next_trees;
    std::map<std::string, DetectorHolder*> next_holders;
    
    int overall_event = std::numeric_limits<int>::max();
    int max_tree_count = 0;
    for (const auto & s : next_detector_list) {
        next_trees.insert({s, (TTree*)f->Get(s.c_str())});
        DetectorHolder * holder = new DetectorHolder();
        holder->SetTree(next_trees[s]);
        next_holders.insert({s, holder});
        next_holders[s]->SetTree(next_trees[s]);
        next_trees[s]->GetEntry(0);

        //Get the earliest event
        std::cout << "First event in " << s << " " << next_holders[s]->EventID << std::endl;
        if (next_holders[s]->EventID < overall_event) {
            overall_event = next_holders[s]->EventID;
        }
        int entries = next_trees[s]->GetEntries();
        max_tree_count = (max_tree_count < entries) ? entries : max_tree_count;
    }
    // std::cout << "Max count " << max_tree_count << std::endl;
    max_events = (max_events > max_tree_count) ? max_tree_count : max_events;
    std::cout << "Max events: " << max_events << std::endl;

    if (start != std::numeric_limits<int>::max()) overall_event = start;

    int bt_count = -1;
    // bool first_ibt = true;

    std::map<int, SimpleBeamParticle> particle_map;

    // // Loop over the entries in before target
    // for (int ibt = 0; ibt < first_tree->GetEntries(); ++ibt) {
    //     first_tree->GetEntry(ibt);
    //     int this_event = before_target.EventID;

    //     //We're in a new event in the overall simulation -- go through the current 
    //     if (this_event != overall_event) {
    //         if (!first_ibt) {
    //             if (verbose) {

    //                 std::cout << "!!!!!!!!!!!!!!!!!!!!New event" << std::endl;
                    
    //                 std::cout << "Parts:" << std::endl;
    //                 for (const auto & [id, part] : particle_map) {
    //                     // std::cout << "\t" << id << " " << part.TrackID << " " << part.PDGid << std::endl;
    //                     std::cout << part << std::endl;
    //                 }
    //             }

    //             for (auto & name : next_detector_list) {
    //                 if (verbose)
    //                     std::cout << "Processing " << name << "------------------------------" << std::endl;
    //                 ProcessDetector(next_trees[name], particle_map, next_holders[name], overall_event, name);
    //             }

    //             if (verbose) {
    //                 std::cout << "Parts:" << std::endl;
    //                 for (const auto & [id, part] : particle_map) {
    //                     // std::cout << "\t" << id << " " << part.TrackID << " " << part.PDGid << std::endl;
    //                     std::cout << part << std::endl;
    //                 }
    //             }
    //         }
    //         first_ibt = false;

    //         particle_map.clear();
    //         bt_count++;
    //         if (bt_count >= max_events) {
    //             std::cout << "Hit max events" << std::endl;
    //             break;
    //         }
    //     }
    //     overall_event = this_event;

    //     //Check if track id in particle_map keys
    //     if (particle_map.find(before_target.TrackID) != particle_map.end()) {
    //         throw std::runtime_error("Duplicate track ID in overall event");
    //     }

    //     SimpleBeamParticle this_particle(&before_target, first_detector);       
    //     particle_map.insert({before_target.TrackID, this_particle});
    // }

    // Loop over the entries in before target
    
    TFile fout(output_file.c_str(), "recreate");
    TTree out_tree("tree", "");
    std::vector<SimpleBeamParticle> output_vector;
    std::vector<int> output_TrackID, output_EventID, output_PDGid, output_ParentID;
    std::vector<std::vector<std::string>> output_detectors;
    std::vector<std::vector<float>> output_xs, output_ys, output_zs, output_ts,
                                         output_Pxs, output_Pys, output_Pzs;
    // out_tree.Branch("TrackIDs", &output_TrackID);
    // out_tree.Branch("EventIDs", &output_EventID);
    // out_tree.Branch("PDGids", &output_PDGid);
    // out_tree.Branch("ParentIDs", &output_ParentID);
    // out_tree.Branch("detectors", &output_detectors);
    // out_tree.Branch("xs", &output_xs);
    // out_tree.Branch("ys", &output_ys);
    // out_tree.Branch("zs", &output_zs);
    // out_tree.Branch("Pxs", &output_Pxs);
    // out_tree.Branch("Pys", &output_Pys);
    // out_tree.Branch("Pzs", &output_Pzs);
    // out_tree.Branch("ts", &output_ts);
    // out_tree.Branch("particles", &output_vector);
    out_tree.Branch("particle_map", &particle_map);
    while (bt_count < (max_events-1)) {
                if (verbose) {

                    std::cout << "New event " << bt_count << " " << overall_event << "!!!!!!!!!!!!" << std::endl;
                    
                    // std::cout << "Parts:" << std::endl;
                    // for (const auto & [id, part] : particle_map) {
                    //     // std::cout << "\t" << id << " " << part.TrackID << " " << part.PDGid << std::endl;
                    //     std::cout << part << std::endl;
                    // }
                }

                bool all_done = true;
                for (auto & name : next_detector_list) {
                    if (verbose)
                        std::cout << "Processing " << name << "------------------------------" << std::endl;
                    all_done &= ProcessDetector(next_trees[name], particle_map, next_holders[name], overall_event, name);
                }

                if (verbose) {
                    std::cout << "Parts:" << std::endl;
                    for (const auto & [id, part] : particle_map) {
                        // std::cout << "\t" << id << " " << part.TrackID << " " << part.PDGid << std::endl;
                        std::cout << part << std::endl;
                    }
                }


            //Fill output

            if (particle_map.size() > 0) {
                ++bt_count;
                std::cout << bt_count << " " << particle_map.size() << std::endl;
            }
            else {
                std::cout << overall_event << " Empty event. Moving on and not counting for limit%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << std::endl;  
            }

            output_vector.clear();
            for (auto & [id, part] : particle_map) {
                output_vector.push_back(part);
                if (verbose) std::cout << "Writing part " << part << std::endl;
            }
            out_tree.Fill();
            particle_map.clear();
            ++overall_event;
            if (all_done) {
                std::cout << "Hit limit on all trees" << std::endl;
                break;
            }
    }
    
    out_tree.Write();
    fout.Close();
    f->Close();

    return 0;
}
