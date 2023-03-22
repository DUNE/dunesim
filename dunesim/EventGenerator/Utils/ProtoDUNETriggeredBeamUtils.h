#include "TTree.h"
#include "TVector3.h"
#include "TFile.h"
#include <iostream>
#include "fhiclcpp/ParameterSetRegistry.h"
#include "fhiclcpp/ParameterSet.h"
// Fill the particle maps using the input files. This links the events of interest
// to the entry number in fAllParticlesTree.
namespace evgen {
// Simple struct to store the information for each particle at the front face
 struct BeamParticle {
   BeamParticle(){
     fTrackID=-999;
     fPDG=-999;
     fParentID=-999;
     fPosX = -999;
     fPosY = -999;
     fPosZ = -999;
     fPosT = -999;
     fMomX = -999;
     fMomY = -999;
     fMomZ = -999; 
   };
   BeamParticle(int trackid, int pdg, int parentid, float posX, float posY, float posZ, float posT,
                float momX, float momY, float momZ){
     fTrackID = trackid;
     fPDG     = pdg;
     fParentID= parentid;
     fPosX    = posX;
     fPosY    = posY;
     fPosZ    = posZ;  
     fPosT    = posT;
     fMomX    = momX;
     fMomY    = momY;
     fMomZ    = momZ;
   };
   void Print(){
     std::cout << "Particle " << fPDG << ": (" << fPosX << "," << fPosY << "," << fPosZ << "," << fPosT << ") "
                              << ": (" << fMomX << "," << fMomY << "," << fMomZ << ") " << std::endl;
   };
   int fTrackID, fPDG, fParentID;
   float fPosX, fPosY, fPosZ, fPosT;
   float fMomX, fMomY, fMomZ;
 };

 // Struct to contain the particles reaching the cryostat wall for each event
 // in the beam simulation files
 struct BeamEvent {
   BeamEvent(){
     fEventID = -999;
     fTriggerID = -999;
     fHasInteracted = false;
   }
   BeamEvent(int eventid){
     fEventID = eventid;
     fTriggerID = -999;
     fHasInteracted = false;
   };

   void AddParticle(BeamParticle particle){
     fParticlesFront.insert(std::make_pair(particle.fTrackID,particle));
   };

   int fEventID;

   // Map of particles to the track ID
   std::map<int,BeamParticle> fParticlesFront;

   int fTriggerID;

   // We need information for each point in the beamline
   std::map<std::string,BeamParticle> fTriggeredParticleInfo;

   // Some events can interact before between TRIG2 and NP04front
   bool fHasInteracted;
   std::vector<int> fSecondaryTrackIDs;
 };

 // Convenience struct to encapsulate all particles that would
 // deposit energy within one readout window of the TPC
 struct OverlaidTriggerEvent {

   OverlaidTriggerEvent(int trigID){
     fTriggerEventID = trigID;
   };

   void AddOverlay(int overlayID){
     fOverlayEventIDs.push_back(overlayID);
   };

   std::vector<int> fOverlayEventIDs;
   int fTriggerEventID;

 };

 class ProtoDUNETriggeredBeamUtils {
  public:
   ProtoDUNETriggeredBeamUtils(fhicl::ParameterSet const & pset);

   // Fill the above maps and vector.
   void FillParticleMaps(TTree * frontFaceTree,
                         std::map<int, BeamEvent> & allBeamEvents);

   // Add the triggered particle information for a given instrument
   void FillInstrumentInformation(
       std::vector<int> &eventIDs, TTree *instrumentTree,
       std::map<int, BeamEvent> & allBeamEvents);

   // Convert to the detector coordinate frame
   void ConvertCoordinates(float & x, float & y, float & z) {
     // Convert to cm and shift to the detector coordinate frame
     x = (x/10.) + fBeamX;
     y = (y/10.) + fBeamY;
     z = fBeamZ; // Just use the z position    
   };

   // Convert the momentum to GeV and rotate as required.
   void ConvertMomentum(float & momX, float & momY, float & momZ);

   // Methods for making beam instrument tracks
   TVector3 ConvertProfCoordinates(double x, double y, double z,
       double zOffset);

   // Find trigger events
   std::vector<int> FindTriggeredEvents(
       TFile * inputFile, std::string trig1TreeName, std::string trig2TreeName,
       std::map<int, BeamEvent> & allBeamEvents);

  void BeamMonitorBasisVectors();
  void RotateMonitorVector(TVector3 &vec);

  bool GetIsNP02() {return fIsNP02;};

  private:
   bool fReduceNP04frontArea;
   float fBeamX, fBeamY, fBeamZ;
   bool fIsNP02;
   bool fNP02XDrift;
   double fNP02Rotation;
   double fBeamPhiShift, fBeamThetaShift;
   float fTRIG2Pos;
   float fRotateMonitorXZ;
   float fRotateMonitorYZ;
   float fNP04frontPos;
   TVector3 fBMBasisX; 
   TVector3 fBMBasisY; 
   TVector3 fBMBasisZ; 
 };
}

