// IdealAdcSimulator_tool.cc

#include "IdealAdcSimulator.h"

//**********************************************************************

IdealAdcSimulator::IdealAdcSimulator(double vsen, unsigned int nbit)
: m_vsen(vsen), m_vmax(0.0), m_adcmax(0) {
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

IdealAdcSimulator::IdealAdcSimulator(const fhicl::ParameterSet& ps)
: IdealAdcSimulator(ps.get<double>("Vsen"), ps.get<unsigned int>("Nbit")) { }

//**********************************************************************

AdcSimulator::Count
IdealAdcSimulator::count(double vin, Channel, Tick) const {
  //std::cout << "v_sen " << m_vsen <<  " m_vmax " << m_vmax << " vin " << vin << std::endl;
  double halfsen = 0.5*m_vsen;
  if ( m_vsen <= 0.0 ) return 0;
  if ( vin < halfsen ) return 0;
  if ( vin > m_vmax - halfsen ) return m_adcmax;
  Count count = vin/m_vsen + 0.5;
  assert( count <= m_adcmax );
  //std::cout << " count " << count << " m_adcmax " << m_adcmax << std::endl;
  if ( count > m_adcmax ) return m_adcmax;//
  return count;
}

//**********************************************************************
