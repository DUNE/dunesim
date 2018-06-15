////////////////////////////////////////////////////////////////////////
//
// Event generator for the 311 (Can be also used for ProtoDUNE) with realistic
// muon trigger
//
//
////////////////////////////////////////////////////////////////////////
//#ifndef EVGEN_311Gen_H
//#define EVGEN_311Gen_H

// ROOT includes
#include "TRandom3.h"
#include "TDatabasePDG.h"
#include "TString.h"
#include "TSystem.h" //need BaseName and DirName
#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"
#include "TVector3.h"

// Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// art extensions
#include "nutools/RandomUtils/NuRandomService.h"

// larsoft includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/EventGeneratorBase/evgenbase.h"
#include "larcorealg/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcoreobj/SummaryData/RunData.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "ifdh.h"  //to handle flux files

// c/c++ includes
#include <sqlite3.h>
#include <memory>
#include <stdio.h>
#include <dirent.h>

namespace evgendp{

  class Gen311 : public art::EDProducer {

    public:
      explicit Gen311(fhicl::ParameterSet const & p);
      virtual ~Gen311();

      Gen311(Gen311 const &) = delete;
      Gen311(Gen311 &&) = delete;
      Gen311 & operator = (Gen311 const &) = delete;
      Gen311 & operator = (Gen311 &&) = delete;

      void produce(art::Event & e) override;

      void reconfigure(fhicl::ParameterSet const& p);
      void beginJob() override;
      void beginRun(art::Run& run) override;

      int fShowerInputs=0; ///< Number of shower inputs to process from
      std::vector<double> fNShowersPerEvent; ///< Number of showers to put in each event of duration fSampleTime; one per showerinput
      std::vector<int> fMaxShowers; //< Max number of showers to query, one per showerinput
      double fShowerBounds[6]={0.,0.,0.,0.,0.,0.}; ///< Boundaries of area over which showers are to be distributed
      double fToffset_corsika=0.; ///< Timing offset to account for propagation time through atmosphere, populated from db
      ifdh_ns::ifdh* fIFDH=0; ///< (optional) flux file handling

      //fcl parameters
      double fProjectToHeight=0.; ///< Height to which particles will be projected [cm]
      std::vector< std::string > fShowerInputFiles; ///< Set of CORSIKA shower data files to use
      std::vector< double > fShowerFluxConstants; ///< Set of flux constants to be associated with each shower data file
      double fSampleTime=0.; ///< Duration of sample [s]
      double fToffset=0.; ///< Time offset of sample, defaults to zero (no offset) [s]
      std::vector<double> fBuffBox; ///< Buffer box extensions to cryostat in each direction (6 of them: x_lo,x_hi,y_lo,y_hi,z_lo,z_hi) [cm]
      double fShowerAreaExtension=0.; ///< Extend distribution of corsika particles in x,z by this much (e.g. 1000 will extend 10 m in -x, +x, -z, and +z) [cm]
      sqlite3* fdb[5]; ///< Pointers to sqlite3 database object, max of 5
      double fRandomYZShift=0.; ///< Each shower will be shifted by a random amount in xz so that showers won't repeatedly sample the same space [cm]
      std::vector<double> fActiveVolumeCut; ///< Active volume cut
      std::vector< double > fLeadingParticlesList; /// < List of pdg codes for particles that can be potentially used as triggers
      bool fUseIFDH; ///<< use ifdh protocol

      int fRun;
      int fEvent;

      double fTheta;
      double fPhi;
      double fMomX;
      double fMomY;
      double fMomZ;
      double fMomT;
      double fStartX;
      double fStartY;
      double fStartZ;
      TTree *fTree;

      art::ServiceHandle<geo::Geometry> geom;

    private:

      void openDBs();
      void populateNShowers();
      void populateTOffset();

      bool InTPC(const simb::MCParticle particle);
      double wrapvar( const double var, const double low, const double high);
      double wrapvarBoxNo( const double var, const double low, const double high, int& boxno);
      void GetSample(simb::MCTruth&);

    };

    class Trigger{

      public:
        Trigger();
        ~Trigger();

        //setter
        void AddMuon( simb::MCParticle particle ){
          fMuonList.push_back(particle); return;
        };

        //getters
        int GetTriggerId(){
          return fTriggerID;
        };

        double GetTriggerOffsetX(){
          return fTriggerOffsetX;
        };

        double GetTriggerOffsetY(){
          return fTriggerOffsetY;
        };

        double GetTriggerOffsetZ(){
          return fTriggerOffsetZ;
        };

        double GetTriggerOffsetT(){
          return fTriggerOffsetT;
        };

        TLorentzVector GetTriggerPos(){
          return fTriggerPos;
        };

        double GetTriggerPosX(){
          return fTriggerPosX;
        };

        double GetTriggerPosY(){
          return fTriggerPosY;
        };

        double GetTriggerPosZ(){
          return fTriggerPosZ;
        };

        double GetTriggerPosT(){
          return fTriggerPosT;
        };

        TLorentzVector GetTriggerMom(){
          return fTriggerMu.Momentum();
        };

        double GetTriggerMomX(){
          return fTriggerMu.Momentum().X();
        };

        double GetTriggerMomY(){
          return fTriggerMu.Momentum().Y();
        };

        double GetTriggerMomZ(){
          return fTriggerMu.Momentum().Z();
        };

        double GetTriggerMomT(){
          return fTriggerMu.Momentum().T();
        };

        //Find the trigger particle
        void MakeTrigger();

        //geo utilities
        double GetPhi( const double py, const double pz );
        void GetMatrix( double theta, double phi, double (*p_R)[3][3] );
        void DoRotation( double urv[], double dir[], double theta, double phi );
        bool Intersect(const double x0[], const double dx[], const double bounds[]);
        void ProjectToBoxEdge(	const double 	xyz[],
                                const double 	indxyz[],
                                const double 	xlo,
                                const double 	xhi,
                                const double 	ylo,
                                const double 	yhi,
                                const double 	zlo,
                                const double 	zhi,
                                double xyzout[]
                              );
        void GetTPCSize( double tpc[]);
        void GetCryoSize( double cryo[]);
        void SetCryoBuffer( std::vector<double> buffer ){
          for(int i=0; i<6; i++){ fCryoBuffer.assign( buffer.begin(), buffer.end() ); }
        }
        void SetTPCBuffer( std::vector<double> buffer ){
          for(int i=0; i<6; i++){ fTPCBuffer.assign( buffer.begin(), buffer.end() ); }
        }


      private:

        art::ServiceHandle<geo::Geometry> geom;

        std::vector<simb::MCParticle>  fMuonList;
        simb::MCParticle fTriggerMu;
        TLorentzVector fTriggerPos;
        double fTriggerID = -999;
        double fTriggerPosX;
        double fTriggerPosY;
        double fTriggerPosZ;
        double fTriggerPosT;
        double fTriggerOffsetX;
        double fTriggerOffsetY;
        double fTriggerOffsetZ;
        double fTriggerOffsetT;
        std::vector<double> fTPCBuffer;
        std::vector<double> fCryoBuffer;
    };

}//end namespace


  evgendp::Gen311::Gen311(fhicl::ParameterSet const & p)
    : fProjectToHeight(p.get< double >("ProjectToHeight",0.)),
      fShowerInputFiles(p.get< std::vector< std::string > >("ShowerInputFiles")),
      fShowerFluxConstants(p.get< std::vector< double > >("ShowerFluxConstants")),
      fSampleTime(p.get< double >("SampleTime",0.)),
      fToffset(p.get< double >("TimeOffset",0.)),
      fBuffBox(p.get< std::vector< double > >("BufferBox",{0.0, 0.0, 0.0, 0.0, 0.0, 0.0})),
      fShowerAreaExtension(p.get< double >("ShowerAreaExtension",0.)),
      fRandomYZShift(p.get< double >("RandomYZShift",0.)),
      fActiveVolumeCut(p.get< std::vector< double > >("ActiveVolumeCut")),
      fLeadingParticlesList(p.get< std::vector< double > >("LeadingParticlesList")),
      fUseIFDH(p.get< bool >("UseIFDH"))
   {

    if(fShowerInputFiles.size() != fShowerFluxConstants.size() || fShowerInputFiles.size()==0 || fShowerFluxConstants.size()==0)
      throw cet::exception("Gen311") << "ShowerInputFiles and ShowerFluxConstants have different or invalid sizes!"<<"\n";
    fShowerInputs=fShowerInputFiles.size();

    if(fSampleTime==0.) throw cet::exception("Gen311") << "SampleTime not set!";

    if(fProjectToHeight==0.) mf::LogInfo("Gen311")<<"Using 0. for fProjectToHeight!"
    ;
    // create a default random engine; obtain the random seed from NuRandomService,
    // unless overridden in configuration with key "Seed"
    art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "gen", p, { "Seed", "SeedGenerator" });
    art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, "HepJamesRandom", "pois", p, "SeedPoisson");

    this->reconfigure(p);

    this->openDBs();
    this->populateNShowers();
    this->populateTOffset();

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();
   }

  evgendp::Gen311::~Gen311(){
    for(int i=0; i<fShowerInputs; i++){
      sqlite3_close(fdb[i]);
    }
    //cleanup temp files
    fIFDH->cleanup();
  }

  void evgendp::Gen311::reconfigure(fhicl::ParameterSet const& p){

    return;
  }

  void evgendp::Gen311::beginJob(){


    art::ServiceHandle<art::TFileService> tfs;

    fTree = tfs->make<TTree>("entries", "entries tree");
    fTree->Branch("Run", &fRun, "Run/I");
    fTree->Branch("Event", &fEvent, "Event/I");
    fTree->Branch("Theta", &fTheta, "Theta/D");
    fTree->Branch("Phi", &fPhi, "Phi/D");
    fTree->Branch("MomX", &fMomX, "MomX/D");
    fTree->Branch("MomY", &fMomY, "MomY/D");
    fTree->Branch("MomZ", &fMomZ, "MomZ/D");
    fTree->Branch("MomT", &fMomT, "MomT/D");
    fTree->Branch("StartX", &fStartX, "StartX/D");
    fTree->Branch("StartY", &fStartY, "StartY/D");
    fTree->Branch("StartZ", &fStartZ, "StartZ/D");

  }

  void evgendp::Gen311::beginRun(art::Run& run){
   // grab the geometry object to see what geometry we are using
   art::ServiceHandle<geo::Geometry> geo;
   std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));
   run.put(std::move(runcol));
   return;
  }


  void evgendp::Gen311::openDBs(){
    //choose files based on fShowerInputFiles, copy them with ifdh, open them
    //sqlite3_stmt *statement;
    //get rng engine
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine("gen");
    CLHEP::RandFlat flat(engine);

    if(fUseIFDH){
      //setup ifdh object
      if ( ! fIFDH ) fIFDH = new ifdh_ns::ifdh;
      const char* ifdh_debug_env = std::getenv("IFDH_DEBUG_LEVEL");
      if ( ifdh_debug_env ) {
        mf::LogInfo("CORSIKAGendp") << "IFDH_DEBUG_LEVEL: " << ifdh_debug_env<<"\n";
        fIFDH->set_debug(ifdh_debug_env);
      }
    }

    //setup ifdh object
    if ( ! fIFDH ) fIFDH = new ifdh_ns::ifdh;
    const char* ifdh_debug_env = std::getenv("IFDH_DEBUG_LEVEL");
    if ( ifdh_debug_env ) {
      mf::LogInfo("CORSIKAGendp") << "IFDH_DEBUG_LEVEL: " << ifdh_debug_env<<"\n";
      fIFDH->set_debug(ifdh_debug_env);
    }

    std::vector<std::pair<std::string,long>> selectedflist;
    for(int i=0; i<fShowerInputs; i++){
      if(fShowerInputFiles[i].find("*")==std::string::npos){
        //if there are no wildcards, don't call findMatchingFiles
        selectedflist.push_back(std::make_pair(fShowerInputFiles[i],0));
        mf::LogInfo("Gen311") << "Selected"<<selectedflist.back().first<<"\n";
      }else{
        //read all the file mathching the pattern
        std::vector<std::pair<std::string,long>> flist;
        std::string path(gSystem->DirName(fShowerInputFiles[i].c_str()));
        std::string pattern(gSystem->BaseName(fShowerInputFiles[i].c_str()));

        if(fUseIFDH){
          //acess the files using ifdh
          flist = fIFDH->findMatchingFiles(path,pattern);
        }else{
          //access the files using drent
          auto wildcardPosition = pattern.find("*");
          pattern = pattern.substr( 0, wildcardPosition );

           DIR *dir;
           struct dirent *ent;
           int index=0;
           if ((dir = opendir( path.c_str() )) != NULL) {
               while ((ent = readdir (dir)) != NULL) {
                 index++;
                 std::pair<std::string,long> name;
                 std::string parsename(ent->d_name);

                 if( parsename.substr(0, wildcardPosition) == pattern ){
                   name.first= path+"/"+parsename;
                   name.second = index;
                   flist.push_back( name );
                }
              }
              closedir (dir);
           }else {
               throw cet::exception("Gen311") << "Can't open directory with pattern: "<<path<<":"<<pattern<<std::endl;
           }
        }

        unsigned int selIndex=-1;
        if(flist.size()==1){ //0th element is the search path:pattern
          selIndex=0;
        }else if(flist.size()>1){
          selIndex= (unsigned int) (flat()*(flist.size()-1)+0.5); //rnd with rounding, dont allow picking the 0th element
        }else{
          throw cet::exception("Gen311") << "No files returned for path:pattern: "<<path<<":"<<pattern<<std::endl;
        }
        selectedflist.push_back(flist[selIndex]);
        mf::LogInfo("Gen311") << "For "<<fShowerInputFiles[i]<<":"<<pattern
        <<"\nFound "<< flist.size() << " candidate files"
        <<"\nChoosing file number "<< selIndex << "\n"
        <<"\nSelected "<<selectedflist.back().first<<"\n";
     }
    }

    //open the files in fShowerInputFilesLocalPaths with sqlite3
    std::vector<std::string> locallist;
    for(unsigned int i=0; i<selectedflist.size(); i++){
      mf::LogInfo("Gen311")
        << "Fetching: "<<selectedflist[i].first<<" "<<selectedflist[i].second<<"\n";
      std::string fetchedfile(selectedflist[i].first);
      LOG_DEBUG("Gen311") << "Fetched; local path: "<<fetchedfile;
      locallist.push_back(fetchedfile);
    }

    for(unsigned int i=0; i<locallist.size(); i++){
      //prepare and execute statement to attach db file
      int res=sqlite3_open(locallist[i].c_str(),&fdb[i]);
      if (res!= SQLITE_OK)
        throw cet::exception("Gen311") << "Error opening db: (" <<locallist[i].c_str()<<") ("<<res<<"): " << sqlite3_errmsg(fdb[i]) << "; memory used:<<"<<sqlite3_memory_used()<<"/"<<sqlite3_memory_highwater(0)<<"\n";
      else
        mf::LogInfo("Gen311")<<"Attached db "<< locallist[i]<<"\n";
    }
  }//end openDBs

/*
void evgendp::Gen311::openDBs(){
  //choose files based on fShowerInputFiles, copy them with ifdh, open them
  //sqlite3_stmt *statement;
  //get rng engine
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine &engine = rng->getEngine("gen");
  CLHEP::RandFlat flat(engine);

  //setup ifdh object
  if ( ! fIFDH ) fIFDH = new ifdh_ns::ifdh;
  const char* ifdh_debug_env = std::getenv("IFDH_DEBUG_LEVEL");
  if ( ifdh_debug_env ) {
    mf::LogInfo("CORSIKAGendp") << "IFDH_DEBUG_LEVEL: " << ifdh_debug_env<<"\n";
    fIFDH->set_debug(ifdh_debug_env);
  }

  //get ifdh path for each file in fShowerInputFiles, put into selectedflist
  //if 1 file returned, use that file
  //if >1 file returned, randomly select one file
  //if 0 returned, make exeption for missing files
  std::vector<std::pair<std::string,long>> selectedflist;
  for(int i=0; i<fShowerInputs; i++){
    if(fShowerInputFiles[i].find("*")==std::string::npos){
      //if there are no wildcards, don't call findMatchingFiles
      selectedflist.push_back(std::make_pair(fShowerInputFiles[i],0));
      mf::LogInfo("CorsikaGendp") << "Selected"<<selectedflist.back().first<<"\n";
  }else{
      //use findMatchingFiles
      std::vector<std::pair<std::string,long>> flist;
  std::string path(gSystem->DirName(fShowerInputFiles[i].c_str()));
  std::string pattern(gSystem->BaseName(fShowerInputFiles[i].c_str()));
  //pattern="/"+pattern;


  std::cout << path << " " << pattern << std::endl;

  flist = fIFDH->findMatchingFiles(path,pattern);
  unsigned int selIndex=-1;
  if(flist.size()==1){ //0th element is the search path:pattern
    selIndex=0;
  }else if(flist.size()>1){
    selIndex= (unsigned int) (flat()*(flist.size()-1)+0.5); //rnd with rounding, dont allow picking the 0th element
  }else{
    throw cet::exception("CORSIKAGendp") << "No files returned for path:pattern: "<<path<<":"<<pattern<<std::endl;
  }

  selectedflist.push_back(flist[selIndex]);
  mf::LogInfo("CorsikaGendp") << "For "<<fShowerInputFiles[i]<<":"<<pattern
      <<"\nFound "<< flist.size() << " candidate files"
<<"\nChoosing file number "<< selIndex << "\n"
      <<"\nSelected "<<selectedflist.back().first<<"\n";
   }

  }

  //do the fetching, store local filepaths in locallist
  std::vector<std::string> locallist;
  for(unsigned int i=0; i<selectedflist.size(); i++){
    mf::LogInfo("CorsikaGendp")
      << "Fetching: "<<selectedflist[i].first<<" "<<selectedflist[i].second<<"\n";
    std::string fetchedfile(fIFDH->fetchInput(selectedflist[i].first));
    LOG_DEBUG("CorsikaGendp") << "Fetched; local path: "<<fetchedfile;
    locallist.push_back(fetchedfile);
  }

  //open the files in fShowerInputFilesLocalPaths with sqlite3
  for(unsigned int i=0; i<locallist.size(); i++){
    //prepare and execute statement to attach db file
    int res=sqlite3_open(locallist[i].c_str(),&fdb[i]);
    if (res!= SQLITE_OK)
      throw cet::exception("CORSIKAGendp") << "Error opening db: (" <<locallist[i].c_str()<<") ("<<res<<"): " << sqlite3_errmsg(fdb[i]) << "; memory used:<<"<<sqlite3_memory_used()<<"/"<<sqlite3_memory_highwater(0)<<"\n";
    else
      mf::LogInfo("CORSIKAGendp")<<"Attached db "<< locallist[i]<<"\n";
  }
}
*/

  void evgendp::Gen311::populateTOffset(){
    //populate TOffset_corsika by finding minimum ParticleTime from db file

    sqlite3_stmt *statement;
    const std::string kStatement("select min(t) from particles");
    double t=0.;

    for(int i=0; i<fShowerInputs; i++){
        //build and do query to get run min(t) from each db
        if ( sqlite3_prepare(fdb[i], kStatement.c_str(), -1, &statement, 0 ) == SQLITE_OK ){
          int res=0;
          res = sqlite3_step(statement);
          if ( res == SQLITE_ROW ){
            t=sqlite3_column_double(statement,0);
            mf::LogInfo("Gen311")<<"For showers input "<< i<<" found particles.min(t)="<<t<<"\n";
            if (i==0 || t<fToffset_corsika) fToffset_corsika=t;
          }else{
            throw cet::exception("Gen311") << "Unexpected sqlite3_step return value: (" <<res<<"); "<<"ERROR:"<<sqlite3_errmsg(fdb[i])<<"\n";
          }
        }else{
          throw cet::exception("Gen311") << "Error preparing statement: (" <<kStatement<<"); "<<"ERROR:"<<sqlite3_errmsg(fdb[i])<<"\n";
        }
    }

    mf::LogInfo("Gen311")<<"Found corsika timeoffset [ns]: "<< fToffset_corsika<<"\n";
  }

  void evgendp::Gen311::populateNShowers(){
    //populate vector of the number of showers per event based on:
      //AREA the showers are being distributed over
      //TIME of the event (fSampleTime)
      //flux constants that determine the overall normalizations (fShowerFluxConstants)
      //Energy range over which the sample was generated (ERANGE_*)
      //power spectrum over which the sample was generated (ESLOPE)


    //compute shower area based on the maximal x,z dimensions of cryostat boundaries + fShowerAreaExtension
    //art::ServiceHandle<geo::Geometry> geom;
    for(unsigned int c = 0; c < geom->Ncryostats(); ++c){
      double bounds[6] = {0.};
      geom->CryostatBoundaries(bounds, c);
      for (unsigned int bnd = 0; bnd<6; bnd++){
        mf::LogVerbatim("Gen311")<<"Cryo Boundary: "<<bnd<<"="<<bounds[bnd]<<" ( + Buffer="<<fBuffBox[bnd]<<")\n";
        if(fabs(bounds[bnd])>fabs(fShowerBounds[bnd])){
          fShowerBounds[bnd]=bounds[bnd];
        }
      }
    }
    //add on fShowerAreaExtension without being clever
    fShowerBounds[2] = fShowerBounds[2] - fShowerAreaExtension;
    fShowerBounds[3] = fShowerBounds[3] + fShowerAreaExtension;
    fShowerBounds[4] = fShowerBounds[4] - fShowerAreaExtension;
    fShowerBounds[5] = fShowerBounds[5] + fShowerAreaExtension;

    double showersArea=(fShowerBounds[3]/100-fShowerBounds[2]/100)*(fShowerBounds[5]/100-fShowerBounds[4]/100);

    mf::LogInfo("Gen311")
      <<  "Area extended by : "<<fShowerAreaExtension
      <<"\nShowers to be distributed betweeen: y="<<fShowerBounds[2]<<","<<fShowerBounds[3]
                             <<" & z="<<fShowerBounds[4]<<","<<fShowerBounds[5]
      <<"\nShowers to be distributed with random YZ shift: "<<fRandomYZShift<<" cm"
      <<"\nShowers to be distributed over area: "<<showersArea<<" m^2"
      <<"\nShowers to be distributed over time: "<<fSampleTime<<" s"
      <<"\nShowers to be distributed with time offset: "<<fToffset<<" s"
      <<"\nShowers to be distributed at x: "<<fShowerBounds[1]<<" cm"
      ;

    //db variables
    sqlite3_stmt *statement;
    const std::string kStatement("select erange_high,erange_low,eslope,nshow from input");
    double upperLimitOfEnergyRange=0.,lowerLimitOfEnergyRange=0.,energySlope=0.,oneMinusGamma=0.,EiToOneMinusGamma=0.,EfToOneMinusGamma=0.;

    for(int i=0; i<fShowerInputs; i++){
        //build and do query to get run info from databases
      //  double thisrnd=flat();//need a new random number for each query
        if ( sqlite3_prepare(fdb[i], kStatement.c_str(), -1, &statement, 0 ) == SQLITE_OK ){
          int res=0;
          res = sqlite3_step(statement);
          if ( res == SQLITE_ROW ){
            upperLimitOfEnergyRange=sqlite3_column_double(statement,0);
            lowerLimitOfEnergyRange=sqlite3_column_double(statement,1);
            energySlope = sqlite3_column_double(statement,2);
            fMaxShowers.push_back(sqlite3_column_int(statement,3));
            oneMinusGamma = 1 + energySlope;
            EiToOneMinusGamma = pow(lowerLimitOfEnergyRange, oneMinusGamma);
            EfToOneMinusGamma = pow(upperLimitOfEnergyRange, oneMinusGamma);
            mf::LogVerbatim("Gen311")<<"For showers input "<< i<<" found e_hi="<<upperLimitOfEnergyRange<<", e_lo="<<lowerLimitOfEnergyRange<<", slope="<<energySlope<<", k="<<fShowerFluxConstants[i]<<"\n";
          }else{
            throw cet::exception("Gen311") << "Unexpected sqlite3_step return value: (" <<res<<"); "<<"ERROR:"<<sqlite3_errmsg(fdb[i])<<"\n";
          }
        }else{
          throw cet::exception("Gen311") << "Error preparing statement: (" <<kStatement<<"); "<<"ERROR:"<<sqlite3_errmsg(fdb[i])<<"\n";
        }

      //this is computed, how?
      double NShowers=( M_PI * showersArea * fShowerFluxConstants[i] * (EfToOneMinusGamma - EiToOneMinusGamma) / oneMinusGamma )*fSampleTime;
      fNShowersPerEvent.push_back(NShowers);
      mf::LogVerbatim("Gen311")<<"For showers input "<< i
                               <<" the number of showers per event is "<<(int)NShowers<<"\n";
    }
  }

  double evgendp::Gen311::wrapvar( const double var, const double low, const double high){
    //wrap variable so that it's always between low and high
    return (var - (high - low) * floor(var/(high-low))) + low;
  }

  double evgendp::Gen311::wrapvarBoxNo( const double var, const double low, const double high, int& boxno){
    //wrap variable so that it's always between low and high
    boxno=int(floor(var/(high-low)));
    return (var - (high - low) * floor(var/(high-low))) + low;
  }

  void evgendp::Gen311::GetSample(simb::MCTruth& mctruth){
    //for each input, randomly pull fNShowersPerEvent[i] showers from the Particles table
    //and randomly place them in time (between -fSampleTime/2 and fSampleTime/2)
    //wrap their positions based on the size of the area under consideration
    //based on http://nusoft.fnal.gov/larsoft/doxsvn/html/CRYHelper_8cxx_source.html (Sample)

    //query from sqlite db with select * from particles where shower in (select id from showers ORDER BY substr(id*0.51123124141,length(id)+2) limit 100000) ORDER BY substr(shower*0.51123124141,length(shower)+2);
    //where 0.51123124141 is a random seed to allow randomly selecting rows and should be randomly generated for each query
    //the inner order by is to select randomly from the possible shower id's
    //the outer order by is to make sure the shower numbers are ordered randomly (without this, the showers always come out ordered by shower number
    //and 100000 is the number of showers to be selected at random and needs to be less than the number of showers in the showers table

    //dummy holder for particles object
    std::map< int, std::vector<simb::MCParticle> >  ParticleMap;
    std::map< int, int >  ShowerTrkIDMap;

    //define the trigger object
    Trigger trg;
    trg.SetCryoBuffer( fBuffBox );
    trg.SetTPCBuffer( fActiveVolumeCut );

    //TDatabasePDG is for looking up particle masses
    static TDatabasePDG* pdgt = TDatabasePDG::Instance();

    //db variables
    sqlite3_stmt *statement;
    const TString kStatement("select shower,pdg,px,py,pz,x,z,t,e from particles where shower in (select id from showers ORDER BY substr(id*%f,length(id)+2) limit %d) ORDER BY substr(shower*%f,length(shower)+2)");

    //get rng engine
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine("gen");
    CLHEP::RandFlat flat(engine);

    CLHEP::HepRandomEngine &engine_pois = rng->getEngine("pois");
    CLHEP::RandPoissonQ randpois(engine_pois);

    // get geometry and figure where to project particles to, based on CRYHelper
    //art::ServiceHandle<geo::Geometry> geom;
    double x1, x2;
    double y1, y2;
    double z1, z2;
    geom->WorldBox(&x1, &x2, &y1, &y2, &z1, &z2);

    // make the world box slightly smaller so that the projection to
    // the edge avoids possible rounding errors later on with Geant4
    double fBoxDelta=1.e-5;

    geom->WorldBox(&y1, &y2, &x1, &x2, &z1, &z2);
    y1 += fBoxDelta;
    y2 -= fBoxDelta;
    x1 += fBoxDelta;
    x2 = fProjectToHeight;
    z1 += fBoxDelta;
    z2 -= fBoxDelta;

    //populate mctruth
    int ntotalCtr=0; //count number of particles added to mctruth
    //int lastShower=0; //keep track of last shower id so that t can be randomized on every new shower
    int nShowerCntr=0; //keep track of how many showers are left to be added to mctruth
    int nShowerQry=0; //number of showers to query from db
    int shower,pdg;
    double px,py,pz,etot, y=0,z,t; //tParticleTime,,showerTime=0.,showerTimex=0.,showerTimez=0.,showerXOffset=0.,showerZOffset=0.,t;
    for(int i=0; i<fShowerInputs; i++){
      nShowerCntr=randpois.fire(fNShowersPerEvent[i]);
      mf::LogInfo("Gen311") << " Shower input " << i << " with mean " << fNShowersPerEvent[i] << " generating " << nShowerCntr;

      while(nShowerCntr>0){
        //how many showers should we query?
        if(nShowerCntr>fMaxShowers[i]){
          nShowerQry=fMaxShowers[i]; //take the group size
        }else{
          nShowerQry=nShowerCntr; //take the rest that are needed
        }
        //build and do query to get nshowers
        double thisrnd=flat(); //need a new random number for each query
        TString kthisStatement=TString::Format(kStatement.Data(),thisrnd,nShowerQry,thisrnd);
        LOG_DEBUG("Gen311")<<"Executing: "<<kthisStatement;
        if ( sqlite3_prepare(fdb[i], kthisStatement.Data(), -1, &statement, 0 ) == SQLITE_OK ){
          int res=0;
          //loop over database rows, pushing particles into mctruth object
          while(1){
            res = sqlite3_step(statement);
            if ( res == SQLITE_ROW ){
              shower=sqlite3_column_int(statement,0);

              pdg=sqlite3_column_int(statement,1);
              //get mass for this particle
              double m = 0.; // in GeV
              TParticlePDG* pdgp = pdgt->GetParticle(pdg);
              if (pdgp) m = pdgp->Mass();

              //get momentum components
              //NR rotation from CORSIKA ref system to LArSoft ref system:
              //LArsoft: px, py, pz;
              //uBoone px', py', pz'; << from the sqlite database instances
              // px = py'
              // py = pz'
              // pz = px'
              px = sqlite3_column_double(statement,3);
              py = sqlite3_column_double(statement,4);
              pz = sqlite3_column_double(statement,2);
              etot=sqlite3_column_double(statement,8);
              y = sqlite3_column_double(statement,6);
              z = sqlite3_column_double(statement,5);
              t = sqlite3_column_double(statement,7); //time offset, includes propagation time from top of atmosphere

              TLorentzVector pos(fShowerBounds[1], y, z,t);// time needs to be in ns to match GENIE, etc
              TLorentzVector mom(px,py,pz,etot);

              simb::MCParticle p(ntotalCtr,pdg,"primary",-200,m,1);
              p.AddTrajectoryPoint(pos,mom);

              //if it is muon, fill up the list of possible triggers
              std::vector<double>::iterator pdgIt;
              pdgIt = find( fLeadingParticlesList.begin(), fLeadingParticlesList.end(), p.PdgCode() );
              if( pdgIt != fLeadingParticlesList.end() ){ trg.AddMuon( p ); }

              ParticleMap[ shower ].push_back(p);
              ShowerTrkIDMap[ ntotalCtr ] = shower; //<< use the unique trackID to trace the source shower back

              ntotalCtr++;
            }else if ( res == SQLITE_DONE ){
              break;
            }else{
              throw cet::exception("Gen311") << "Unexpected sqlite3_step return value: (" <<res<<"); "<<"ERROR:"<<sqlite3_errmsg(fdb[i])<<"\n";
            }
          }
        }else{
          throw cet::exception("Gen311") << "Error preparing statement: (" <<kthisStatement<<"); "<<"ERROR:"<<sqlite3_errmsg(fdb[i])<<"\n";
        }
        nShowerCntr=nShowerCntr-nShowerQry;
      }
    } //end loop over showers

    trg.MakeTrigger(); //<<--Select a muon to use as trigger

    double showerTime=0, showerTimey=0, showerTimez=0, showerYOffset=0, showerZOffset=0;
    int boxnoY=0, boxnoZ=0;

    for( auto ParticleMapIt : ParticleMap ){

      if( ParticleMapIt.first == ShowerTrkIDMap[ trg.GetTriggerId() ] ){

        showerTime =1e9*trg.GetTriggerOffsetT(); //converting from s to ns
        showerTimey=1e9*trg.GetTriggerOffsetT();
        showerTimez=1e9*trg.GetTriggerOffsetT();

        //and a random offset in both z and x controlled by the fRandomYZShift parameter
        showerYOffset=1e9*trg.GetTriggerOffsetY();
        showerZOffset=1e9*trg.GetTriggerOffsetZ();

      }else{

        showerTime =1e9*(flat()*fSampleTime); //converting from s to ns
        showerTimey=1e9*(flat()*fSampleTime); //converting from s to ns
        showerTimez=1e9*(flat()*fSampleTime); //converting from s to ns
        //and a random offset in both z and x controlled by the fRandomYZShift parameter
        showerYOffset=flat()*fRandomYZShift - (fRandomYZShift/2);
        showerZOffset=flat()*fRandomYZShift - (fRandomYZShift/2);
      }

      for( auto particleIt : ParticleMapIt.second ){

        simb::MCParticle particle(particleIt.TrackId(),particleIt.PdgCode(),"primary",-200,particleIt.Mass(),1);

        double x0[3]={0.};
        double xyzo[3]={0.};
        double dx[3]={0.};

        if( particleIt.TrackId() == trg.GetTriggerId() ){

          x0[0] = trg.GetTriggerPos().X();
          x0[1] = trg.GetTriggerPos().Y();
          x0[2] = trg.GetTriggerPos().Z();
          dx[0] = trg.GetTriggerMom().X();
          dx[1] = trg.GetTriggerMom().Y();
          dx[2] = trg.GetTriggerMom().Z();

          t = trg.GetTriggerPos().T();

        }else{

          y = wrapvarBoxNo(particleIt.Position().Y()+showerYOffset,fShowerBounds[2],fShowerBounds[3],boxnoY);
          z = wrapvarBoxNo(particleIt.Position().Z()+showerZOffset,fShowerBounds[4],fShowerBounds[5],boxnoZ);

          t=particleIt.Position().T()+showerTime+(1e9*fToffset)-fToffset_corsika + showerTimey*boxnoY + showerTimez*boxnoZ;
          t=wrapvar(t,(1e9*fToffset),1e9*(fToffset+fSampleTime));

          x0[0]= fShowerBounds[1];
          x0[1]= y;
          x0[2]= z;
          dx[0]= particleIt.Momentum().X();
          dx[1]= particleIt.Momentum().Y();
          dx[2]= particleIt.Momentum().Z();

        }

        trg.ProjectToBoxEdge(x0, dx, x1, x2, y1, y2, z1, z2, xyzo);
        TLorentzVector pos(xyzo[0],xyzo[1],xyzo[2], t);
        TLorentzVector mom(dx[0],dx[1],dx[2], particleIt.Momentum().T() );
        particle.AddTrajectoryPoint( pos, mom );

        mctruth.Add( particle );
      } //end loop over particles
    }//loop over showers
  } //end GetSample


//Functions to study the trigger ***********************************************************************************

evgendp::Trigger::Trigger(){
  fMuonList.clear();
}

evgendp::Trigger::~Trigger(){}

double evgendp::Trigger::GetPhi( const double py, const double pz ){

    if( pz !=0 ){
      if( pz > 0 ){
        return atan( py/pz );
      }
      else if( atan( py/pz ) < 0 ){
        return atan( py/pz )+3.1415926;
      }
      else{
        return atan( py/pz )-3.1415926;
      }
    }
    else if( py > 0 ){
      return 3.1415926/2.;
    }
    else{
      return -3.1415926/2.;
    }
}

void evgendp::Trigger::GetMatrix( double theta, double phi, double (*p_R)[3][3] ){

  (*p_R)[0][0]= cos(phi)*cos(theta);
  (*p_R)[0][1]= cos(phi)*sin(theta);
  (*p_R)[0][2]= -sin(phi);
  (*p_R)[1][0]= -sin(theta);
  (*p_R)[1][1]= cos(theta);
  (*p_R)[1][2]= 0;
  (*p_R)[2][0]= sin(phi)*cos(theta);
  (*p_R)[2][1]= sin(phi)*sin(theta);
  (*p_R)[2][2]= cos(phi);

  return;
}

void evgendp::Trigger::DoRotation( double urv[], double dir[], double theta, double phi ){

  //build the rotation matrix;
  double R[3][3];
  double (*p_R)[3][3]=&R;
  GetMatrix(theta, phi, p_R);

  for(int i=0; i<3; i++){
    for(int j=0; j<3; j++){
      dir[i] += (*p_R)[i][j]*urv[j];
      //cout << i << " " << j << " " << urv[j] << " " << (*p_R)[i][j] << " " << dir[i] << endl;
    }
  }

  return;
}

bool evgendp::Trigger::Intersect(const double x0[], const double dx[], const double bounds[]){
  //check if the particle intercept the tpc

  bool intersects_tpc = false;
      for (int bnd=0; bnd!=6; ++bnd) {
        if (bnd<2) {
          double p2[3] = {bounds[bnd],  x0[1] + (dx[1]/dx[0])*(bounds[bnd] - x0[0]), x0[2] + (dx[2]/dx[0])*(bounds[bnd] - x0[0])};
          if ( p2[1] >= bounds[2] && p2[1] <= bounds[3] &&
               p2[2] >= bounds[4] && p2[2] <= bounds[5] ) {
            intersects_tpc = true;
            break;
          }
        }
        else if (bnd>=2 && bnd<4) {
          double p2[3] = {x0[0] + (dx[0]/dx[1])*(bounds[bnd] - x0[1]), bounds[bnd], x0[2] + (dx[2]/dx[1])*(bounds[bnd] - x0[1])};
          if ( p2[0] >= bounds[0] && p2[0] <= bounds[1] &&
               p2[2] >= bounds[4] && p2[2] <= bounds[5] ) {
            intersects_tpc = true;
      break;
          }
        }
        else if (bnd>=4) {
          double p2[3] = {x0[0] + (dx[0]/dx[2])*(bounds[bnd] - x0[2]), x0[1] + (dx[1]/dx[2])*(bounds[bnd] - x0[2]), bounds[bnd]};
          if ( p2[0] >= bounds[0] && p2[0] <= bounds[1] &&
               p2[1] >= bounds[2] && p2[1] <= bounds[3] ) {
            intersects_tpc = true;
      break;
          }
        }
      }

      return intersects_tpc;
}

void evgendp::Trigger::GetTPCSize( double tpc[]){

  //art::ServiceHandle<geo::Geometry> geom;

  double minx=0, miny=0, minz=0;
  double maxx=0, maxy=0, maxz=0;
  for (size_t c = 0; c < geom->Ncryostats(); c++){
         const geo::CryostatGeo& cryostat = geom->Cryostat(c);
         for (size_t t = 0; t < cryostat.NTPC(); t++){
             const geo::TPCGeo& tpcg = cryostat.TPC(t);
             if (tpcg.MinX() < minx) minx = tpcg.MinX();
             if (tpcg.MaxX() > maxx) maxx = tpcg.MaxX();
             if (tpcg.MinY() < miny) miny = tpcg.MinY();
             if (tpcg.MaxY() > maxy) maxy = tpcg.MaxY();
             if (tpcg.MinZ() < minz) minz = tpcg.MinZ();
             if (tpcg.MaxZ() > maxz) maxz = tpcg.MaxZ();
         }
  }

  tpc[0] = minx;
  tpc[1] = maxx;
  tpc[2] = miny;
  tpc[3] = maxy;
  tpc[4] = minz;
  tpc[5] = maxz;

  return;
}

void evgendp::Trigger::GetCryoSize( double cryo[]){

  //;

  double dummy[6] = {0};
  for(unsigned int c = 0; c < geom->Ncryostats(); ++c){
    geom->CryostatBoundaries(dummy, c);
  }

  for(int i=0;i<6;i++){cryo[i] = dummy[i];}

  return;
}

void evgendp::Trigger::ProjectToBoxEdge(	const double 	xyz[],
                                      const double 	indxyz[],
                                      const double 	xlo,
                                      const double 	xhi,
                                      const double 	ylo,
                                      const double 	yhi,
                                      const double 	zlo,
                                      const double 	zhi,
                                      double xyzout[]	 ){


  //we want to project backwards, so take mirror of momentum
  const double dxyz[3]={-indxyz[0],-indxyz[1],-indxyz[2]};

  // Compute the distances to the x/y/z walls
  double dx = 99.E99;
  double dy = 99.E99;
  double dz = 99.E99;
  if      (dxyz[0] > 0.0) { dx = (xhi-xyz[0])/dxyz[0]; }
  else if (dxyz[0] < 0.0) { dx = (xlo-xyz[0])/dxyz[0]; }
  if      (dxyz[1] > 0.0) { dy = (yhi-xyz[1])/dxyz[1]; }
  else if (dxyz[1] < 0.0) { dy = (ylo-xyz[1])/dxyz[1]; }
  if      (dxyz[2] > 0.0) { dz = (zhi-xyz[2])/dxyz[2]; }
  else if (dxyz[2] < 0.0) { dz = (zlo-xyz[2])/dxyz[2]; }


  // Choose the shortest distance
  double d = 0.0;
  if      (dx < dy && dx < dz) d = dx;
  else if (dy < dz && dy < dx) d = dy;
  else if (dz < dx && dz < dy) d = dz;

  // Make the step
  for (int i = 0; i < 3; ++i) {
    xyzout[i] = xyz[i] + dxyz[i]*d;
  }

}

void evgendp::Trigger::MakeTrigger(){

    //get the coordinates of the tpc and an active volume arount it
    double tpc[6] ={0.}; this->GetTPCSize(tpc);
    for(int i=0; i<6; i++){ tpc[i] += fTPCBuffer[i]; }

    //get the coordinates of the cryo + a buffer around it
    double cryo[6] ={0.}; this->GetCryoSize(cryo);
    for(int i=0; i<6; i++){ cryo[i] += fCryoBuffer[i]; }

    //get the random engine:
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine &engine = rng->getEngine("gen");
    CLHEP::RandFlat flat(engine);

    //choose a random muon
    if( fMuonList.size() == 0){
      mf::LogInfo("Gen311") << "Trigger muon not found! Only Background";
      return;
    }

    fTriggerMu = fMuonList.at( (int)flat()*fMuonList.size() );

    double px=0,py=0,pz=0;
    double p=0, theta=0, phi=0;

    px = fTriggerMu.Momentum().X();
    py = fTriggerMu.Momentum().Y();
    pz = fTriggerMu.Momentum().Z();


    p = sqrt(px*px+py*py+pz*pz);
    theta = acos(px/p);
    phi = this->GetPhi( py, pz );

    double xyz[3] ={0.};
    xyz[0]=(tpc[5]-tpc[4])/2-tpc[5];
    xyz[1]=(tpc[1]-tpc[0])/2-tpc[1];
    xyz[2]=(tpc[3]-tpc[2])/2-tpc[2];

    //define direction
    double cos_x=cos(theta);
    double cos_y=sin(phi)*sin(theta);
    double cos_z=cos(phi)*sin(theta);
    double dxyz[3] = { cos_z, cos_x, cos_y };

    //get cr impact size along u and v
    double length[2]={0.}; //0 along u, 1 along v
    double corner[3]={0.};
    double rdm_start[3] ={0.};
    double xyzo[3]={0.};

    //loop over the tpc boundaries and project them on u and v
    for(int i=0; i<2; i++){
      for(int j=2; j<4; j++){
        for(int k=4; k<6; k++){

            double dot=0;
            double proj[3]={0.}; //projection on the plane
            double u=0; double v=0;

            //coordinates of the tpc corner
            corner[0] = tpc[k];
            corner[1] = tpc[i];
            corner[2] = tpc[j];

            //now get the dot product between the corner vector and the plane versor
            for( int p=0; p<3; p++ ){  dot+=(corner[p]-xyz[p])*dxyz[p]; }
            for( int p=0; p<3; p++ ){  proj[p] = corner[p]-(dot*dxyz[p]); }

            proj[0] -= (tpc[5]-tpc[4])/2; //NB: not really elegant...makes things work though

            //build the rotation matrix;
            double R[3][3];
            double (*p_R)[3][3]=&R;
            this->GetMatrix(theta, phi, p_R);

            //project proj along u and pick the longest
            for( int p=0; p<3; p++ ){ u += (*p_R)[p][0]*proj[p]; }
            if(length[0] < u){ length[0] = u; }

            //project proj along v and pick the longest projection
            for( int p=0; p<3; p++ ){ v += (*p_R)[p][2]*proj[p]; }
            if(length[1] < v){ length[1] = v; }

        }//end for k
      }//end for j
    } //end for i

    bool is_inside=false;
    int iteration = 20; //NB: It is harcoded

    while(!is_inside && iteration > 0){

        //We extract a random start over the plane u, v
        double u = flat()*2*length[0]-length[0];
        double v = flat()*2*length[1]-length[1];;

        //now rotate along the r direction
        double urv[3] ={ u,  0,  v };
        double dir[3]={ 0., 0., 0. };
        this->DoRotation( urv, dir, theta, phi );

        //parse the correct direction of our reference system
        rdm_start[0] = (double)dir[1];
        rdm_start[1] = (double)dir[2];
        rdm_start[2] = (double)dir[0]+(tpc[5]-tpc[4])/2;

        double dx[3]={px,py,pz};
        this->ProjectToBoxEdge(rdm_start, dx, cryo[0], cryo[1], cryo[2], cryo[3], cryo[4], cryo[5], xyzo);

        if( this->Intersect(xyzo, dx, tpc) ){
          is_inside=true;
          //mf::LogInfo("Gen311")<< "======> " << iteration;
        }else{
          //mf::LogInfo("Gen311")<< "======> " << iteration;
          iteration--;
        }
    }//end while


    if(iteration==0){
      mf::LogInfo("Gen311") << "Trigger muon not found in 20 iterations. Only Background";
    }

    /*

    double rdm_start[3] ={0.};
    rdm_start[0] = 0;
    rdm_start[1] = flat()*(tpc[3]-tpc[2])-(tpc[3]-tpc[2])/2;
    rdm_start[2] = flat()*(tpc[5]-tpc[4])-(tpc[5]-tpc[4])/2;

    */

    fTriggerID = fTriggerMu.TrackId();
    fTriggerPosX = rdm_start[0];
    fTriggerPosY = rdm_start[1];
    fTriggerPosZ = rdm_start[2];
    fTriggerPosT = 0;
    fTriggerPos.SetXYZT(rdm_start[0],rdm_start[1],rdm_start[2],0);
    fTriggerOffsetX = rdm_start[0] - fTriggerMu.Position().X();
    fTriggerOffsetY = rdm_start[1] - fTriggerMu.Position().Y();
    fTriggerOffsetZ = rdm_start[2] - fTriggerMu.Position().Z();
    fTriggerOffsetT = 0 - fTriggerMu.Position().T();

    return;
  }//end make trigger

//******************************************************************************************************************

  void evgendp::Gen311::produce(art::Event & e){
    //Fetch the particle generated by GetParticle

    fRun = e.id().run();
    fEvent = e.id().event();

    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    art::ServiceHandle<geo::Geometry> geom;

    simb::MCTruth truth;
    truth.SetOrigin(simb::kCosmicRay);

    simb::MCTruth pretruth;

    GetSample(pretruth);
    mf::LogInfo("CORSIKAGendp")<<"GetSample number of particles returned: "<<pretruth.NParticles()<<"\n";
    // loop over particles in the truth object

      //add a buffer box around the cryostat bounds to increase the acceptance and account for scattering
      //By default, the buffer box has zero size

    for(int i = 0; i < pretruth.NParticles(); ++i){
      simb::MCParticle particle = pretruth.GetParticle(i);

      TLorentzVector v4 = particle.Position();
      const TLorentzVector& p4 = particle.Momentum();
      double x0[3] = {v4.X(),  v4.Y(),  v4.Z() };
      double dx[3] = {p4.Px(), p4.Py(), p4.Pz()};
      // now check if the particle goes through any cryostat in the detector
      // if so, add it to the truth object.
      for(unsigned int c = 0; c < geom->Ncryostats(); ++c){
        double bounds[6] = {0.};
        geom->CryostatBoundaries(bounds, c);

        //add a buffer box around the cryostat bounds to increase the acceptance and account for scattering
        //By default, the buffer box has zero size
        for (unsigned int cb=0; cb<6; cb++){
           bounds[cb] = bounds[cb]+fBuffBox[cb];
         }

        //calculate the intersection point with each cryostat surface
        bool intersects_cryo = false;
        for (int bnd=0; bnd!=6; ++bnd) {
          if (bnd<2) {
            double p2[3] = {bounds[bnd],  x0[1] + (dx[1]/dx[0])*(bounds[bnd] - x0[0]), x0[2] + (dx[2]/dx[0])*(bounds[bnd] - x0[0])};
            if ( p2[1] >= bounds[2] && p2[1] <= bounds[3] &&
                 p2[2] >= bounds[4] && p2[2] <= bounds[5] ) {
              intersects_cryo = true;
              break;
            }
          }
          else if (bnd>=2 && bnd<4) {
            double p2[3] = {x0[0] + (dx[0]/dx[1])*(bounds[bnd] - x0[1]), bounds[bnd], x0[2] + (dx[2]/dx[1])*(bounds[bnd] - x0[1])};
            if ( p2[0] >= bounds[0] && p2[0] <= bounds[1] &&
                 p2[2] >= bounds[4] && p2[2] <= bounds[5] ) {
              intersects_cryo = true;
        break;
            }
          }
          else if (bnd>=4) {
            double p2[3] = {x0[0] + (dx[0]/dx[2])*(bounds[bnd] - x0[2]), x0[1] + (dx[1]/dx[2])*(bounds[bnd] - x0[2]), bounds[bnd]};
            if ( p2[0] >= bounds[0] && p2[0] <= bounds[1] &&
                 p2[1] >= bounds[2] && p2[1] <= bounds[3] ) {
              intersects_cryo = true;
        break;
            }
          }
        }

        if (intersects_cryo){
          truth.Add(particle);
          break; //leave loop over cryostats to avoid adding particle multiple times
        }// end if particle goes into a cryostat
      }// end loop over cryostats in the detector

    }// loop on particles

    mf::LogInfo("Gen311")<<"Number of particles from getsample crossing cryostat + bounding box: "<<truth.NParticles()<<"\n";

    truthcol->push_back(truth);
    e.put(std::move(truthcol));

    return;
  }


DEFINE_ART_MODULE(evgendp::Gen311)

//#endif // EVGEN_311Gen_H
