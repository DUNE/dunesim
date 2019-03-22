#include <string>
#include <iostream>
#include <fstream>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib_except/exception.h"

#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TDatabasePDG.h"
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

  std::vector<simb::MCParticle> GenerateEventKinematics();

private:

  // Energies are in GeV, throughout
  double fMinEnu;
  double fMaxEnu;

  double fGF; // Fermi constant in GeV^-2
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

  double fNueFluxFrac;

  // Handels the cross section calculation, as a function of outgoing e mom.
  // Depends on ga, gv, the electron mass, and neutrino energy
  TF1   *fdsigdT;

  TRandom3 *fRand;
};

//------------------------------------------------------------------------------
evgen::NuEScatterGen::NuEScatterGen(fhicl::ParameterSet const & p)
{
  this->reconfigure(p);

  produces< std::vector<simb::MCTruth>   >();
  produces< sumdata::RunData, art::InRun >();

  // There's a prefactor, too, ignoring since I'll just be pulling from dist
  // [0] = gv   [1] = ga   [2] = Enu   [3] = eMass
  // ga and gv depend on neutrino flavor
  fdsigdT = new TF1("xsecform","TMath::Power([0] + [1],2) + TMath::Power([0] - [1],2) * TMath::Power(1-x/[2],2) - ([1]*[1] - [0]*[0] * TMath::Power([3]/[2],2))*x");
  fRand = new TRandom3(0);
}

//------------------------------------------------------------------------------
evgen::NuEScatterGen::~NuEScatterGen()
{
}

//------------------------------------------------------------------------------
void evgen::NuEScatterGen::beginJob()
{
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

  std::vector<simb::MCParticle> parts = GenerateEventKinematics();

  assert(parts.size()==2);
  truth.Add(parts[0]);
  truth.Add(parts[1]);

  TVector3 vec0(parts[0].Px(),parts[0].Py(),parts[0].Pz());
  TVector3 vec1(parts[1].Px(),parts[1].Py(),parts[1].Pz());

  truthcol->push_back(truth);

  e.put(std::move(truthcol));

  return;
}

//------------------------------------------------------------------------------
void evgen::NuEScatterGen::reconfigure(fhicl::ParameterSet const & p)
{
  fMinEnu = p.get<double>("MinEnu");
  fMaxEnu = p.get<double>("MaxEnu");

  fMinX = p.get<double>("MinX");
  fMaxX = p.get<double>("MaxX");
  fMinY = p.get<double>("MinY");
  fMaxY = p.get<double>("MaxY");
  fMinZ = p.get<double>("MinZ");
  fMaxZ = p.get<double>("MaxZ");
  fMinT = p.get<double>("MinT");
  fMaxT = p.get<double>("MaxT");

  fNueFluxFrac = p.get<double>("NueFluxFrac");

  return;
}

std::vector<simb::MCParticle> evgen::NuEScatterGen::GenerateEventKinematics()
{
  double Enu = fRand->Uniform(fMinEnu,fMaxEnu);

  int ePDG = 11;
  static TDatabasePDG  pdg;
  TParticlePDG* pdgp = pdg.GetParticle(ePDG);
  double eMass = 0;
  if (pdgp) eMass = pdgp->Mass();

  // The total cross section for nue-e  is A(1-2*ssWMA+4/3*ssWMA*ssWMA)
  // The total cross section for nue-mu is A(0.25-ssWMA+4/3*ssWMA*ssWMA)
  // The average nue flux is ~1/3 that of the total flux
  // Wrap these into a random number pull that decides the type of neutrino
  assert(fNueFluxFrac >= 0 && fNueFluxFrac <=1);
  double fracE = fNueFluxFrac * (1-2*fWMA+4./3*fWMA*fWMA) / (fNueFluxFrac * (1-2*fWMA+4./3*fWMA*fWMA) + (1-fNueFluxFrac) * (0.5-fWMA+4./3*fWMA*fWMA));
  // xsec is same for nutau and numu, for convenience make all neutrinos numu
  int nuPDG = gRandom->Uniform(0,1) < fracE ? 12 : 14;

  // ga, gv depend on whether we're scattering nue-e or numu/nutau-e
  double ga = nuPDG==12 ? 0.5 : -0.5;
  double gv = nuPDG==12 ? 2*fWMA+0.5 : 2*fWMA-0.5;
  fdsigdT->SetParameter(0,gv);
  fdsigdT->SetParameter(1,ga);
  fdsigdT->SetParameter(2,Enu);
  fdsigdT->SetParameter(3,eMass);
  fdsigdT->SetMinimum(0);

  double tmax = 2*eMass / (TMath::Power(1+eMass/Enu,2) - 1);
  double eKE = fdsigdT->GetRandom(0,tmax);
  double eMom = sqrt(eKE*eKE+2*eMass*eKE);
  std::cout << "Throwing a neutrino of energy " << Enu << " with a maximum electron energy " << tmax << " from cross section " << fdsigdT->GetExpFormula() << ", and got " << eKE <<  std::endl;

  // Throw a random neutrino direction
  double nuCos = fRand->Uniform(-1,1);
  double nuPhi = gRandom->Uniform(0,2*TMath::Pi());
  TVector3 nu(sqrt(1-nuCos*nuCos)*sin(nuPhi),
              sqrt(1-nuCos*nuCos)*cos(nuPhi),
              nuCos);

  // It's a bit tricky translating between the electron momentum in the nu
  // and detector frame.  Do it by making a Gram-schmidt basis around the nu
  // direction, and applying projecting onto it
  TVector3 zprime = nu.Unit();
  TVector3 xprime(0,0,0);
  xprime.SetX(1-nu[0]*nu[0]);
  xprime.SetY( -nu[0]*nu[1]);
  xprime.SetZ( -nu[0]*nu[2]);
  xprime = xprime.Unit();
  TVector3 yprime(0,0,0);
  yprime.SetX( -nu[1]*nu[0]-xprime[1]*xprime[0]);
  yprime.SetY(1-nu[1]*nu[1]-xprime[1]*xprime[1]);
  yprime.SetZ( -nu[1]*nu[2]-xprime[1]*xprime[2]);
  yprime = yprime.Unit();

  double eCos =
    TMath::Sqrt( (eKE*TMath::Power(1+(eMass/Enu),2) ) / ( 2*eMass + eKE));
  double ePhi = fRand->Uniform(0,2*TMath::Pi());

  double eSin = sqrt(1-eCos*eCos);
  TVector3 e(0,0,0);
  e.SetX(eCos*zprime[0]+eSin*sin(ePhi)*xprime[0]+eSin*cos(ePhi)*yprime[0]);
  e.SetY(eCos*zprime[1]+eSin*sin(ePhi)*xprime[1]+eSin*cos(ePhi)*yprime[1]);
  e.SetZ(eCos*zprime[2]+eSin*sin(ePhi)*xprime[2]+eSin*cos(ePhi)*yprime[2]);

  e *= eMom;
  nu *= Enu;

  TLorentzVector e4d(e,eMass+eKE);
  TLorentzVector nu4d(nu,Enu);

  double vtxx = fRand->Uniform(fMinX,fMaxX);
  double vtxy = fRand->Uniform(fMinY,fMaxY);
  double vtxz = fRand->Uniform(fMinZ,fMaxZ);
  double vtxt = fRand->Uniform(fMinT,fMaxT);
  TLorentzVector vtx(vtxx,vtxy,vtxz,vtxt);

  // arguments are particle#, pdg code, process, mother, mass, status
  // status must be 1 to tel GEANT to propogate
  simb::MCParticle mcNu(0, nuPDG, "primary", 0, 0, 1);
  mcNu.AddTrajectoryPoint(vtx, e4d);
  simb::MCParticle mcE(1, 11, "primary", 0, eMass, 1);
  mcE.AddTrajectoryPoint(vtx, nu4d);

  std::vector<simb::MCParticle> ret;
  ret.push_back(mcNu);
  ret.push_back(mcE);

  return ret;
}


DEFINE_ART_MODULE(evgen::NuEScatterGen)
