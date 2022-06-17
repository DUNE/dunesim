// ShapedAndCoherentProtoDUNENoise_service.cc
//
// Pierre Lasorak, Babak Abi
// Dec 2020

#include "dunesim/DetSim/Service/ShapedCohProtoDUNENoiseService.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "nurandom/RandomUtils/NuRandomService.h"
#include "dunepdlegacy/Services/ChannelMap/PdspChannelMapService.h"

#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "nurandom/RandomUtils/NuRandomService.h"

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandGaussQ.h"

ShapedCohProtoDUNENoiseService::ShapedCohProtoDUNENoiseService(fhicl::ParameterSet const& pset) {
  
  art::ServiceHandle<rndm::NuRandomService> seedSvc;
  m_pran = new CLHEP::HepJamesRandom();
  seedSvc->registerEngine(rndm::NuRandomService::CLHEPengineSeeder(m_pran), "ShapedCohProtoDUNENoiseService_rand");
  
  std::cout << "Initialising ShapedCohProtoDUNENoiseService\n";
  double amplitude_col = pset.get<double>("amplitude_multiplicator_collection");
  double amplitude_ind = pset.get<double>("amplitude_multiplicator_induction");
  m_collection_plane_noise = pset.get<double>("collection_plane_noise") * amplitude_col;
  m_induction_plane_noise  = pset.get<double>("induction_plane_noise" ) * amplitude_ind;
  m_collection_plane_noise_rms = pset.get<double>("collection_plane_noise_rms") * amplitude_col;
  m_induction_plane_noise_rms  = pset.get<double>("induction_plane_noise_rms" ) * amplitude_ind;
  
  m_FEMBCo_Frq_nominal = pset.get<std::vector<std::vector<double>>>("FEMBCo_Frq");
  m_FEMBCo_Amp_nominal = pset.get<std::vector<std::vector<double>>>("FEMBCo_Amp");
  m_FEMBCo_Phs_nominal = pset.get<std::vector<std::vector<double>>>("FEMBCo_Phs");
  m_random_phase = pset.get<bool>("random_phase_noise");
  
  assert(m_FEMBCo_Frq_nominal.size() == m_n_femb);
  assert(m_FEMBCo_Amp_nominal.size() == m_n_femb);
  assert(m_FEMBCo_Phs_nominal.size() == m_n_femb);

  for (size_t femb=0; femb<m_n_femb; ++femb) {
    assert(m_FEMBCo_Frq_nominal[femb].size() <= m_n_max_coh_noise);
    assert(m_FEMBCo_Amp_nominal[femb].size() <= m_n_max_coh_noise);
    assert(m_FEMBCo_Phs_nominal[femb].size() <= m_n_max_coh_noise);
  }
  
  // HV noise
  m_HV1_Frq_nominal = pset.get<std::vector<double>>("HV1_Frq");
  m_HV1_Amp_nominal = pset.get<std::vector<double>>("HV1_Amp");
  m_HV1_Phs_nominal = pset.get<std::vector<double>>("HV1_Phs");
  assert(m_HV1_Frq_nominal.size() <= m_n_max_coh_noise);
  assert(m_HV1_Amp_nominal.size() <= m_n_max_coh_noise);
  assert(m_HV1_Phs_nominal.size() <= m_n_max_coh_noise);
  m_init = false;
}

void ShapedCohProtoDUNENoiseService::init() {
  art::ServiceHandle<dune::PdspChannelMapService> channel_map;
  
  std::map<size_t,std::vector<size_t>> femb_se;
  for (size_t chan=0; chan<m_n_wire_per_apa; ++chan) {
    // int n_wire = chan%m_n_wire_per_apa; // that is one little useless statement
    int femb = channel_map->FEMBFromOfflineChannel  (chan)-1; // so that it starts from 0
    int slot = channel_map->SlotIdFromOfflineChannel(chan); // if we ever need this...
    // int femb_asic = channel_map->FEMBChannelFromOfflineChannel(chan) / 32; or this
    femb_se[femb+slot*4].push_back(chan);
    m_channel_femb[chan] = femb+slot*4;
  }
  m_init = true;
  std::cout << "ShapedCohProtoDUNENoiseService initialised!\n";
}

ShapedCohProtoDUNENoiseService::ShapedCohProtoDUNENoiseService(fhicl::ParameterSet const& pset, art::ActivityRegistry&): ShapedCohProtoDUNENoiseService(pset) { }

ShapedCohProtoDUNENoiseService::~ShapedCohProtoDUNENoiseService() { }



int ShapedCohProtoDUNENoiseService::addNoise(detinfo::DetectorClocksData const& clock,
                                             detinfo::DetectorPropertiesData const&,
                                             Channel c, AdcSignalVector& adcs) const {

  assert(adcs.size() <= m_n_tick);
  int ret = addFEMBNoise(c, adcs, clock) + addHVNoise(c, adcs, clock) + addShapedNoise(c, adcs, clock);
  return ret;
}

int ShapedCohProtoDUNENoiseService::addShapedNoise(const Channel c, AdcSignalVector& adcs, detinfo::DetectorClocksData const& clock) const {
  
  art::ServiceHandle<geo::Geometry> geo;
  const geo::View_t view = geo->View(c);
  CLHEP::HepRandomEngine& engine = *m_pran;
  CLHEP::RandGauss g_ind(engine, m_induction_plane_noise , m_induction_plane_noise_rms );
  CLHEP::RandGauss g_col(engine, m_collection_plane_noise, m_collection_plane_noise_rms);
  
  double ampl=-1;
  while (ampl<0) {
    switch (view) {
    case geo::kU:
    case geo::kV:
      ampl = g_ind.fire();
      break;
    default:
      ampl = g_col.fire();
      break;
    }
  }
  CLHEP::RandGauss rgauss(engine, 0, ampl);
  
  std::vector<double> noise_vec;
  for (size_t i=0; i<m_n_tick; ++i)
    noise_vec.push_back(rgauss.fire());
  
  art::ServiceHandle<util::SignalShapingServiceDUNE> sss;
  sss->ConvoluteElectronicResponse(clock, c, noise_vec);
  
  std::transform(adcs.begin(), adcs.end(), noise_vec.begin(), adcs.begin(), std::plus<double>());

  return 0;
}

int ShapedCohProtoDUNENoiseService::addFEMBNoise(const Channel c, AdcSignalVector& adcs, detinfo::DetectorClocksData const& clock) const {
  const size_t femb = m_channel_femb.at(c%m_n_wire_per_apa);
  const size_t apa = c/m_n_wire_per_apa;

  for (size_t sample=0; sample<adcs.size(); ++sample)
    adcs[sample] += m_FEMBCo_Wfm_this_event_vec[apa][femb][sample];
  
  return 0;
}

int ShapedCohProtoDUNENoiseService::addHVNoise(const Channel c, AdcSignalVector& adcs, detinfo::DetectorClocksData const& clock) const {

  for (size_t sample=0; sample<adcs.size(); ++sample)
    adcs[sample] += m_HV1_Wfm_this_event_vec[sample];
  
  return 0;
}



void ShapedCohProtoDUNENoiseService::newEvent() {
  if (not m_init) init();

  CLHEP::RandGauss gaus(*m_pran);
  CLHEP::RandFlat flat(*m_pran);
  const double two_pi_tick_period = 6.2831853072/2000000;
  for (size_t apa=0; apa<m_n_apa; ++apa) {
    for (size_t femb=0; femb<m_FEMBCo_Frq_nominal.size(); ++femb) {
      for (size_t noise=0; noise<m_FEMBCo_Frq_nominal[femb].size(); ++noise) {
      
        double frq = -1;
        double amp = -1;
        double phs = -1;
        
        while (frq<0)
          frq = gaus.fire(m_FEMBCo_Frq_nominal[femb][noise], m_frequency_rms);
     
        while (amp<0)
          amp = gaus.fire(m_FEMBCo_Amp_nominal[femb][noise], m_amplitude_rms);

        if (m_random_phase) {
          while (phs<0)
            phs = flat.fire(0, 6.2831853072);
        } else {
          throw std::runtime_error("ShapedCohProtoDUNENoiseService::newEvent(): Non random phases are not implemented yet");
        }
        
        if (noise == 0) {
          for (size_t sample=0; sample<m_n_tick; ++sample){
            double two_pi_t_f = two_pi_tick_period*sample*frq;
            m_FEMBCo_Wfm_this_event_vec[apa][femb][sample] = amp * sin(phs+two_pi_t_f);
          }
        } else {
          for (size_t sample=0; sample<m_n_tick; ++sample){
            double two_pi_t_f = two_pi_tick_period*sample*frq;
            m_FEMBCo_Wfm_this_event_vec[apa][femb][sample] += amp * sin(phs+two_pi_t_f);
          }
        }
      }
    }
  }
  
  for (size_t noise=0; noise<m_HV1_Frq_nominal.size(); ++noise) {
    double frq = -1;
    double amp = -1;
    double phs = -1;
      
    while (frq<0)
      frq = gaus.fire(m_HV1_Frq_nominal[noise], m_frequency_rms);
     
    while (amp<0)
      amp = gaus.fire(m_HV1_Amp_nominal[noise], m_amplitude_rms);

    if (m_random_phase) {
      while (phs<0)
        phs = flat.fire(0, 6.2831853072);
    } else {
      throw std::runtime_error("ShaoedCohProtoDUNENoise::newEvent(): Non random phases are not implemented yet");
    }
        
    if (noise == 0) {
      for (size_t sample=0; sample<m_n_tick; ++sample){
        double two_pi_t_f = two_pi_tick_period*sample*frq;
        m_HV1_Wfm_this_event_vec[sample] = amp * sin(phs+two_pi_t_f);
      }
    } else {
      for (size_t sample=0; sample<m_n_tick; ++sample){
        double two_pi_t_f = two_pi_tick_period*sample*frq;
        m_HV1_Wfm_this_event_vec[sample] += amp * sin(phs+two_pi_t_f);
      }
    }
  }
}

std::ostream& ShapedCohProtoDUNENoiseService::print(std::ostream& o, std::string) const {
  return o;
}

DEFINE_ART_SERVICE_INTERFACE_IMPL(ShapedCohProtoDUNENoiseService, ChannelNoiseService)
