#include "ProtoDUNETriggeredBeamUtils.h"
#include "TVector3.h"

evgen::ProtoDUNETriggeredBeamUtils::ProtoDUNETriggeredBeamUtils(fhicl::ParameterSet const & pset)
 : fReduceNP04frontArea(pset.get<bool>("ReduceNP04frontArea")),
   fBeamX(pset.get<double>("BeamX")),
   fBeamY(pset.get<double>("BeamY")),
   fBeamZ(pset.get<double>("BeamZ")),
   fIsNP02(pset.get<bool>("IsNP02", false)),
   fNP02Rotation(pset.get<double>("NP02Rotation", 0.)),
   fBeamPhiShift(pset.get<double>("BeamPhiShift", 0.)),
   fBeamThetaShift(pset.get<double>("BeamThetaShift", 0.)) {
}

void evgen::ProtoDUNETriggeredBeamUtils::FillParticleMaps(
    TTree * frontFaceTree, std::map<int, BeamEvent> & allBeamEvents){
    
  float eventID, trackID, pdgCode, parentID;
  float posX, posY, posZ, posT;
  float momX, momY, momZ;

  frontFaceTree->SetBranchAddress("EventID",&eventID);
  frontFaceTree->SetBranchAddress("TrackID",&trackID);
  frontFaceTree->SetBranchAddress("PDGid",&pdgCode);
  frontFaceTree->SetBranchAddress("ParentID",&parentID);
  frontFaceTree->SetBranchAddress("x",&posX);
  frontFaceTree->SetBranchAddress("y",&posY);
  frontFaceTree->SetBranchAddress("z",&posZ);
  frontFaceTree->SetBranchAddress("t",&posT);
  frontFaceTree->SetBranchAddress("Px",&momX);
  frontFaceTree->SetBranchAddress("Py",&momY);
  frontFaceTree->SetBranchAddress("Pz",&momZ);

  // Loop over all particles and group them by events
  for(unsigned int p = 0; p < frontFaceTree->GetEntries(); ++p){

    frontFaceTree->GetEntry(p);

    // Don't consider nuclei here
    const int intPdgCode = static_cast<int>(pdgCode);
    if(intPdgCode > 10000) continue;

    // If this particle is travelling backwards then it won't hit the detector
    if(momZ < 0) continue;

    // Convert to detector coordinate system
    ConvertCoordinates(posX,posY,posZ); //TODO

    // Keep only those particles that might reach the detector
    // Need for NP02?
    if(fReduceNP04frontArea){ //TODO
      if(posX < -500 || posX > 500) continue;
      if(posY < -150 || posY > 850) continue;
    }

    // Convert momentum
    ConvertMomentum(momX,momY,momZ); //TODO

    const int intEventID = static_cast<int>(eventID);
    const int intTrackID = static_cast<int>(trackID);
    const int intParentID= static_cast<int>(parentID);
    BeamParticle newParticle(intTrackID, intPdgCode, intParentID,
                             posX, posY, posZ, posT, momX, momY, momZ);

    std::map<int,BeamEvent>::iterator evIter = allBeamEvents.find(intEventID); //TODO
    if(evIter == allBeamEvents.end()){
      BeamEvent newBeamEvent(intEventID);
      newBeamEvent.AddParticle(newParticle);
      allBeamEvents.insert(std::make_pair(intEventID,newBeamEvent));
    }
    else{
      evIter->second.AddParticle(newParticle); 
    }

  }

  std::cout << "Found " << allBeamEvents.size() << " potential events" << std::endl;

  // Reset the branch addresses as the variables are going out of scope
  frontFaceTree->ResetBranchAddresses();
}

void evgen::ProtoDUNETriggeredBeamUtils::ConvertMomentum(float px, float py, float pz) {
  // Convert to GeV
  px = px / 1000.;
  py = py / 1000.;
  pz = pz / 1000.;
  
  TVector3 momVec(px,py,pz);

  if (!fIsNP02) {
    //NP04: If we want to rotate by changing theta and phi, do it here.
    momVec.SetTheta(momVec.Theta() + fBeamThetaShift);
    momVec.SetPhi(momVec.Phi() + fBeamPhiShift);

    px = momVec.X();
    py = momVec.Y();
    pz = momVec.Z();
  }
  else {
    //NP02: Beam ntuples are in the 'world' coordinate system
    //The beam is similar to the NP04 direction: ~ -8deg from Z
    //In our simulation, it needs to be at -135deg
    //Then we need to swap x and y
    momVec.RotateY(fNP02Rotation*TMath::Pi()/180.);
    px = (fNP02XDrift ? momVec.Y() : momVec.X());
    py = (fNP02XDrift ? -1.*momVec.X() : momVec.Y());
    pz = momVec.Z();
  }
}

void evgen::ProtoDUNETriggeredBeamUtils::FillInstrumentInformation(
    std::vector<int> &eventIDs, TTree *instrumentTree,
    std::map<int, BeamEvent> & allBeamEvents) {

  float eventID, trackID, pdgCode, parentID;
  float posX, posY, posZ, posT;
  float momX, momY, momZ;

  instrumentTree->SetBranchAddress("EventID",&eventID);
  instrumentTree->SetBranchAddress("TrackID",&trackID);
  instrumentTree->SetBranchAddress("PDGid",&pdgCode);
  instrumentTree->SetBranchAddress("ParentID",&parentID);
  instrumentTree->SetBranchAddress("x",&posX);
  instrumentTree->SetBranchAddress("y",&posY);
  instrumentTree->SetBranchAddress("z",&posZ);
  instrumentTree->SetBranchAddress("t",&posT);
  instrumentTree->SetBranchAddress("Px",&momX);
  instrumentTree->SetBranchAddress("Py",&momY);
  instrumentTree->SetBranchAddress("Pz",&momZ);

  // Buffer all of the tree entries for trigger events
  std::map<const int,std::vector<unsigned int>> triggerIndices; 
  std::map<const int,const bool> arePionDecays;
  std::map<const int,const int> trig1TrackIDs;
  std::map<const int,const int> trig2TrackIDs;
  std::map<const int,bool> foundTrackInEvent;
  for(const int &trigEventID : eventIDs){
    BeamEvent &event = allBeamEvents.at(trigEventID);
    const int trig1TrackID = event.fTriggeredParticleInfo.at("TRIG1").fTrackID;
    const int trig2TrackID = event.fTriggeredParticleInfo.at("TRIG2").fTrackID;
    const int trig1TrackPDG = event.fTriggeredParticleInfo.at("TRIG1").fPDG;
    const int trig2TrackPDG = event.fTriggeredParticleInfo.at("TRIG2").fPDG;
    const bool isPionDecay = ((trig1TrackPDG==211) && (trig2TrackPDG==-13)) || 
                             ((trig1TrackPDG==-211) && (trig2TrackPDG==13));

    arePionDecays.insert(std::make_pair(trigEventID,isPionDecay));
    trig1TrackIDs.insert(std::make_pair(trigEventID,trig1TrackID));
    trig2TrackIDs.insert(std::make_pair(trigEventID,trig2TrackID));
    foundTrackInEvent.insert(std::make_pair(trigEventID,false));

    triggerIndices.insert(std::make_pair(trigEventID,std::vector<unsigned int>()));
  }

  // Strip the first part of the tree name to get just the instrument name
  std::string treeName = instrumentTree->GetName();
  treeName = treeName.substr(treeName.find("/")+1);

  for(unsigned int p = 0; p < instrumentTree->GetEntries(); ++p){
    instrumentTree->GetEntry(p);
    // If this isn't a triggered event then move on
    const int thisEvent = static_cast<int>(eventID);
    if(std::find(eventIDs.begin(),eventIDs.end(),thisEvent)==eventIDs.end()) continue;

    triggerIndices.at(thisEvent).push_back(p);
    const int thisParticle = static_cast<int>(trackID);
    if(trig2TrackIDs.at(thisEvent) == thisParticle){
      BeamParticle particle(static_cast<int>(trackID), static_cast<int>(pdgCode), static_cast<int>(parentID),
                            posX, posY, posZ, posT, momX, momY, momZ);
      allBeamEvents.at(thisEvent).fTriggeredParticleInfo.insert(std::make_pair(treeName,particle));
      foundTrackInEvent.at(thisEvent) = true;
    }
  }


  for (auto it = eventIDs.begin(); it != eventIDs.end();) {
    const int ev = *it;
    if(foundTrackInEvent.at(ev)) {
      ++it;
      continue;
    }

    // If we didn't find TRIG2 particle, then look for the TRIG1 one
    for(const unsigned int index : triggerIndices.at(ev)){
      instrumentTree->GetEntry(index);
      const int thisEvent = static_cast<int>(eventID);
      const int thisParticle = static_cast<int>(trackID);
      if(trig1TrackIDs.at(thisEvent) == thisParticle){
        BeamParticle particle(static_cast<int>(trackID), static_cast<int>(pdgCode), static_cast<int>(parentID),
                            posX, posY, posZ, posT, momX, momY, momZ);
        allBeamEvents.at(thisEvent).fTriggeredParticleInfo.insert(std::make_pair(treeName,particle));
        foundTrackInEvent.at(thisEvent) = true;
        break;
      }
    }
  

    if(foundTrackInEvent.at(ev)) {
      ++it;
      continue;
    }
    // In the rare case that we still don't have the particle, try the TRIG1 parent
    const int parentTrack = allBeamEvents.at(ev).fTriggeredParticleInfo.at("TRIG1").fParentID;
    for(const unsigned int index : triggerIndices.at(ev)){
      instrumentTree->GetEntry(index);
      const int thisEvent = static_cast<int>(eventID);
      const int thisParticle = static_cast<int>(trackID);
      if(parentTrack == thisParticle){
        BeamParticle particle(static_cast<int>(trackID), static_cast<int>(pdgCode), static_cast<int>(parentID),
                            posX, posY, posZ, posT, momX, momY, momZ);
        allBeamEvents.at(thisEvent).fTriggeredParticleInfo.insert(std::make_pair(treeName,particle));
        foundTrackInEvent.at(thisEvent) = true;

        break;
      }

    }

    if(foundTrackInEvent.at(ev) == false){
      allBeamEvents.at(ev).fTriggerID = -999;
      // Remove this event from the input vector
      it = eventIDs.erase(it);
      std::cout << "Issue found with event " << ev << ". Removing it from the trigger list - " << eventIDs.size() << " remain" << std::endl;
      //std::cout << "We didn't find tracks " << trig2TrackIDs.at(ev) << " or " << trig1TrackIDs.at(ev) << " in " << treeName << std::endl;
      //std::cout << " - PDGs: 1 = " << allBeamEvents.at(ev).fTriggeredParticleInfo.at("TRIG1").fPDG
      //          << " and 2 = " << allBeamEvents.at(ev).fTriggeredParticleInfo.at("TRIG2").fPDG << std::endl;
      //for(const unsigned int index : triggerIndices.at(ev)){
      //  instrumentTree->GetEntry(index);
      //  std::cout << "- Choice = " << static_cast<int>(trackID) << " :: " << static_cast<int>(pdgCode) << std::endl;
      //}
    }
    else {
      ++it; 
    }
  }

  instrumentTree->ResetBranchAddresses();
}
