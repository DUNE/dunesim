#include <string>
#include <iostream>
#include <sstream>
#include <fstream>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"
#include "cetlib/search_path.h"

#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TF1.h"

#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

namespace evgen {
  class NuEScatterGen;
}

class evgen::NuEScatterGen : public art::EDProducer {
public:
  explicit NuEScatterGen(fhicl::ParameterSet const & p);
  virtual ~NuEScatterGen();

  void produce(art::Event & e)                    override;
  void beginJob()               		  override;
  void beginRun(art::Run & run) 		  override;
  void reconfigure(fhicl::ParameterSet const & p) ;

  std::vector<simb::MCParticle> GenerateEventKinematics(bool isNewNu);

private:

  // Energies are in GeV, throughout
  double fMinEnu;
  double fMaxEnu;

  //dla double fGF; // Fermi constant in GeV^-2
  double fWMA; // sin^2 (Weak mixing angle)

  // Detector coordinates, in cm
  double fMinX;
  double fMaxX;
  double fMinY;
  double fMaxY;
  double fMinZ;
  double fMaxZ;
  double fMinT;
  double fMaxT;

  std::string fEventRateFileName;

  bool fIsSupernova;
  int fNNu; // number of neutrinos to generate per supernova

  // Event rate distributions in energy for each flavor
  TF1 *fNueE;
  TF1 *fNumuE;
  TF1 *fNutauE;
  TF1 *fNuebarE;
  TF1 *fNumubarE;
  TF1 *fNutaubarE;

  // The neutrino direction, which may be the same for multiple neutrinos
  TVector3 fNuDir;

  // Handels the cross section calculation, as a function of outgoing e mom.
  // Depends on ga, gv, the electron mass, and neutrino energy
  TF1   *fdsigdT;

  TRandom3 *fRand;
};

//------------------------------------------------------------------------------
evgen::NuEScatterGen::NuEScatterGen(fhicl::ParameterSet const & p) : EDProducer{p}
{
  this->reconfigure(p);

  produces< std::vector<simb::MCTruth>   >();
  produces< sumdata::RunData, art::InRun >();

  // There's a prefactor, too, ignoring since I'll just be pulling from dist
  // [0] = gv   [1] = ga   [2] = Enu   [3] = eMass
  // ga and gv depend on neutrino flavor
  fdsigdT = new TF1("xsecform","TMath::Power([0] + [1],2) + TMath::Power([0] - [1],2) * TMath::Power(1-x/[2],2) - ([1]*[1] - [0]*[0] * TMath::Power([3]/[2],2))*x");
}

//------------------------------------------------------------------------------
evgen::NuEScatterGen::~NuEScatterGen()
{
}

//------------------------------------------------------------------------------
void evgen::NuEScatterGen::beginJob()
{
  TFile* eventrates = TFile::Open(fEventRateFileName.c_str());
  eventrates->cd();

  fNueE   = (TF1*)eventrates->Get("NueE");
  fNumuE  = (TF1*)eventrates->Get("NumuE");
  fNutauE = (TF1*)eventrates->Get("NutauE");
  fNuebarE   = (TF1*)eventrates->Get("NuebarE");
  fNumubarE  = (TF1*)eventrates->Get("NumubarE");
  fNutaubarE = (TF1*)eventrates->Get("NutaubarE");

  fRand = new TRandom3(0);
}

//------------------------------------------------------------------------------
void evgen::NuEScatterGen::beginRun(art::Run& run)
{
    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;
    std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));

    run.put(std::move(runcol));

    return;
  }

//------------------------------------------------------------------------------
void evgen::NuEScatterGen::produce(art::Event & e)
{
  // Set up a vector of ptrs to eventually put into the event
  std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);
  // And we'll fill one
  simb::MCTruth truth;

  bool isNewNu = (e.event()%fNNu == 0);
  std::vector<simb::MCParticle> parts = GenerateEventKinematics(isNewNu);

  assert(parts.size()==2);
  truth.Add(parts[0]);
  truth.Add(parts[1]);

  truthcol->push_back(truth);

  e.put(std::move(truthcol));

  return;
}

//------------------------------------------------------------------------------
void evgen::NuEScatterGen::reconfigure(fhicl::ParameterSet const & p)
{
  fMinEnu = p.get<double>("MinEnu");
  fMaxEnu = p.get<double>("MaxEnu");

  fWMA  = p.get<double>("WMA");

  fMinX = p.get<double>("MinX");
  fMaxX = p.get<double>("MaxX");
  fMinY = p.get<double>("MinY");
  fMaxY = p.get<double>("MaxY");
  fMinZ = p.get<double>("MinZ");
  fMaxZ = p.get<double>("MaxZ");
  fMinT = p.get<double>("MinT");
  fMaxT = p.get<double>("MaxT");

  fIsSupernova = p.get<bool>("IsSupernova");
  fNNu = p.get<int>("NNu");

  fEventRateFileName = p.get<std::string>("EventRateFileName");

  return;
}

std::vector<simb::MCParticle> evgen::NuEScatterGen::GenerateEventKinematics(bool isNewNu)
{
  // First, need to choose a flavor and energy for our interacting neutrino
  int flav = 0;
  double Enu = 0;

  double totNueE   = fNueE  ->Integral(fMinEnu,fMaxEnu,1e-2);
  double totNumuE  = fNumuE ->Integral(fMinEnu,fMaxEnu,1e-2);
  double totNutauE = fNutauE->Integral(fMinEnu,fMaxEnu,1e-2);
  double totNuebarE   = fNuebarE  ->Integral(fMinEnu,fMaxEnu,1e-2);
  double totNumubarE  = fNumubarE ->Integral(fMinEnu,fMaxEnu,1e-2);
  double totNutaubarE = fNutaubarE->Integral(fMinEnu,fMaxEnu,1e-2);
  double tot = totNueE + totNumuE + totNutauE +
               totNuebarE + totNumubarE + totNutaubarE;

  double randflav = fRand->Uniform(0,1);
  if (randflav < totNueE/tot){
    flav = 12;
    Enu = fNueE->GetRandom(fMinEnu,fMaxEnu);
  }
  else if (randflav < (totNueE+totNumuE)/tot){
    flav = 14;
    Enu = fNumuE->GetRandom(fMinEnu,fMaxEnu);
  }
  else if (randflav < (totNueE+totNumuE+totNutauE)/tot){
    flav = 16;
    Enu = fNutauE->GetRandom(fMinEnu,fMaxEnu);
  }
  else if (randflav < (totNueE+totNumuE+totNutauE+totNuebarE)/tot){
    flav = -12;
    Enu = fNuebarE->GetRandom(fMinEnu,fMaxEnu);
  }
  else if (randflav < (totNueE+totNumuE+totNutauE+totNuebarE+totNumubarE)/tot){
    flav = -14;
    Enu = fNumubarE->GetRandom(fMinEnu,fMaxEnu);
  }
  else{
    flav = -16;
    Enu = fNutaubarE->GetRandom(fMinEnu,fMaxEnu);
  }

  // Throw a random neutrino direction if we're not generating events for
  // a supernova.  Then, generate a neutrino direction only once / supernova
  if (fNuDir.Mag() == 0 || !fIsSupernova || isNewNu){
    double nuCos = fRand->Uniform(-1,1);
    double nuPhi = fRand->Uniform(0,2*TMath::Pi());
    fNuDir.SetX(sqrt(1-nuCos*nuCos)*sin(nuPhi));
    fNuDir.SetY(sqrt(1-nuCos*nuCos)*cos(nuPhi));
    fNuDir.SetZ(nuCos);
  }

  // Pull out the electron mass - needed for dif xsec formulat
  int ePDG = 11;
  static TDatabasePDG  pdg;
  TParticlePDG* pdgp = pdg.GetParticle(ePDG);
  double eMass = 0;
  if (pdgp) eMass = pdgp->Mass();

  // ga, gv depend on whether we're scattering nue-e or numu/nutau-e
  double ga = abs(flav)==12 ? 0.5 : -0.5;
  double gv = abs(flav)==12 ? 2*fWMA+0.5 : 2*fWMA-0.5;
  fdsigdT->SetParameter(0,gv);
  fdsigdT->SetParameter(1,ga);
  fdsigdT->SetParameter(2,Enu);
  fdsigdT->SetParameter(3,eMass);
  fdsigdT->SetMinimum(0);

  double tmax = 2*eMass / (TMath::Power(1+eMass/Enu,2) - 1);
  double eKE = fdsigdT->GetRandom(0,tmax);
  double eMom = sqrt(eKE*eKE+2*eMass*eKE);
  std::cout << "Throwing a neutrino of energy " << Enu << " with a maximum electron energy " << tmax << " from cross section " << fdsigdT->GetExpFormula() << ", and got " << eKE <<  std::endl;

  // Now, pull a direction for your electron
  TVector3 eDir;
  double ECos =
    TMath::Sqrt( (eKE*TMath::Power(1+(eMass/Enu),2) ) / ( 2*eMass + eKE));
  double EPhi = fRand->Uniform(0,2*TMath::Pi());
  eDir.SetX(sqrt(1-ECos*ECos)*sin(EPhi));
  eDir.SetY(sqrt(1-ECos*ECos)*cos(EPhi));
  eDir.SetZ(ECos);
  // But, that's relative to the z-axis move to relative the nu direction
  TVector3 zaxis(0,0,1);
  TVector3 rotationAxis = zaxis.Cross(fNuDir);
  double rotationAngle = zaxis.Angle(fNuDir);
  eDir.Rotate(rotationAngle,rotationAxis);

  TVector3 e = eMom * eDir;;
  TVector3 nu = Enu * fNuDir;;

  TLorentzVector e4d(e,eMass+eKE);
  TLorentzVector nu4d(nu,Enu);

  double vtxx = fRand->Uniform(fMinX,fMaxX);
  double vtxy = fRand->Uniform(fMinY,fMaxY);
  double vtxz = fRand->Uniform(fMinZ,fMaxZ);
  double vtxt = fRand->Uniform(fMinT,fMaxT);
  TLorentzVector vtx(vtxx,vtxy,vtxz,vtxt);

  // arguments are particle#, pdg code, process, mother, mass, status
  // status must be 1 to tel GEANT to propogate
  simb::MCParticle mcNu(0, flav, "primary", 0, 0, 1);
  mcNu.AddTrajectoryPoint(vtx, nu4d);
  simb::MCParticle mcE(1, 11, "primary", 0, eMass, 1);
  mcE.AddTrajectoryPoint(vtx, e4d);

  std::vector<simb::MCParticle> ret;
  ret.push_back(mcNu);
  ret.push_back(mcE);

  return ret;
}


DEFINE_ART_MODULE(evgen::NuEScatterGen)
