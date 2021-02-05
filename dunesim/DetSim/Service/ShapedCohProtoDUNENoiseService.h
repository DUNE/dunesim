// ShapedCohProtoDUNENoiseService
// 
// Pierre Lasorak, Babak Abi
// Dec 2020

#pragma once

#include "dune/DuneInterface/Service/ChannelNoiseService.h"

#include "CLHEP/Random/JamesRandom.h"
#include "dune/Utilities/SignalShapingServiceDUNE.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

class ShapedCohProtoDUNENoiseService : public ChannelNoiseService {
public:

  ShapedCohProtoDUNENoiseService(fhicl::ParameterSet const&);
  ShapedCohProtoDUNENoiseService(fhicl::ParameterSet const&, art::ActivityRegistry&);
  ~ShapedCohProtoDUNENoiseService();

  int addNoise(detinfo::DetectorClocksData     const&,
               detinfo::DetectorPropertiesData const&,
               Channel, AdcSignalVector&) const;
  
  void newEvent();
   
  std::ostream& print(std::ostream&, std::string) const;

protected:
  int addFEMBNoise  (const Channel, AdcSignalVector&, detinfo::DetectorClocksData const&) const;
  int addHVNoise    (const Channel, AdcSignalVector&, detinfo::DetectorClocksData const&) const;
  int addShapedNoise(const Channel, AdcSignalVector&, detinfo::DetectorClocksData const&) const;

private:
  void init();
  static const size_t m_n_apa           =   12;
  static const size_t m_n_tick          = 8192;
  static const size_t m_n_wire_per_apa  = 2560;
  static const size_t m_n_femb          =   20;
  static const size_t m_n_max_coh_noise =    9;
  
  CLHEP::HepRandomEngine* m_pran;
  art::ServiceHandle<geo::Geometry> m_geo;

  double m_collection_plane_noise;
  double m_induction_plane_noise;
  double m_collection_plane_noise_rms;
  double m_induction_plane_noise_rms;
  
  std::map<size_t, size_t> m_channel_femb;
  // Coherent noise
  bool   m_random_phase ;
  double m_frequency_rms;
  double m_amplitude_rms;

  std::vector<std::vector<double>> m_FEMBCo_Frq_nominal;
  std::vector<std::vector<double>> m_FEMBCo_Amp_nominal;
  std::vector<std::vector<double>> m_FEMBCo_Phs_nominal;

  // Fast access
  //dla double m_FEMBCo_Frq_this_event_vec[m_n_apa][m_n_femb][m_n_max_coh_noise];
  //dla double m_FEMBCo_Amp_this_event_vec[m_n_apa][m_n_femb][m_n_max_coh_noise];
  //dla double m_FEMBCo_Phs_this_event_vec[m_n_apa][m_n_femb][m_n_max_coh_noise];

  double m_FEMBCo_Wfm_this_event_vec[m_n_apa][m_n_femb][m_n_tick];
  
  // HV noise
  std::vector<double> m_HV1_Frq_nominal;
  std::vector<double> m_HV1_Amp_nominal;
  std::vector<double> m_HV1_Phs_nominal;

  // Fast access
  //dla double m_HV1_Frq_this_event_vec[m_n_max_coh_noise];
  //dla double m_HV1_Amp_this_event_vec[m_n_max_coh_noise];
  //dla double m_HV1_Phs_this_event_vec[m_n_max_coh_noise];
  double m_HV1_Wfm_this_event_vec[m_n_tick];

  //dla double previous_x = 1;
  //dla double previous_sine = 1;
  bool m_init = 0;
};


DECLARE_ART_SERVICE_INTERFACE_IMPL(ShapedCohProtoDUNENoiseService, ChannelNoiseService, LEGACY)
