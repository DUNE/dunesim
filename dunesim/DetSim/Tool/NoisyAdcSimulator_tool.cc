// NoisyAdcSimulator_tool.cc

#include "NoisyAdcSimulator.h"
using rndm::NuRandomService;
using CLHEP::HepJamesRandom;

//**********************************************************************

NoisyAdcSimulator::NoisyAdcSimulator(double vsen, unsigned int nbit, int noise)
  :m_noise(noise), m_vsen(vsen), m_vmax(0.0), m_adcmax(0), m_random_engine(nullptr) {

  art::ServiceHandle<NuRandomService> seedSvc;
  m_random_engine = new HepJamesRandom;
  seedSvc->registerEngine(NuRandomService::CLHEPengineSeeder(m_random_engine), "NoisyAdcSimulator");

  if ( nbit == 0 ) return;
  unsigned int nbitmax = 8*sizeof(Count);
  if ( nbit > nbitmax ) return;
  unsigned int bit = 1;
  for ( unsigned int ibit=0; ibit<nbit; ++ibit ) {
    m_adcmax |= bit;
    bit = bit << 1;
  }
  m_vmax = m_vsen*m_adcmax;
  
}

//**********************************************************************

NoisyAdcSimulator::NoisyAdcSimulator(const fhicl::ParameterSet& ps)
  : NoisyAdcSimulator(ps.get<double>("Vsen"), ps.get<unsigned int>("Nbit"), ps.get<int>("Noise")) {
  std::cout << "Digitisation noise : " << m_noise << "\n";

}

//**********************************************************************

AdcSimulator::Count
NoisyAdcSimulator::count(double vin, Channel, Tick) const {
  CLHEP::RandBinomial binom(*m_random_engine);

  //std::cout << "v_sen " << m_vsen <<  " m_vmax " << m_vmax << " vin " << vin << std::endl;
  double halfsen = 0.5*m_vsen;
  if ( m_vsen <= 0.0 ) return 0;
  if ( vin < halfsen ) return 0;
  if ( vin > m_vmax - halfsen ) return m_adcmax;
  Count count = vin/m_vsen + 0.5;
  if (m_noise>0) {
    int n = binom.fire(m_noise,0.5)- m_noise/2;
    count += n;
  }
  if ( count > m_adcmax ) return m_adcmax;
  return count;
}

//**********************************************************************
