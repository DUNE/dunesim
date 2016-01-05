// FixedZeroSuppressService.cxx

#include "dune/DetSim/Service/FixedZeroSuppressService.h"
#include <cmath>
#include "fhiclcpp/ParameterSet.h"
#include "dune/DetSim/AdcCodeHelper.h"

using std::string;
using std::ostream;
using std::endl;

//**********************************************************************

FixedZeroSuppressService::
FixedZeroSuppressService(const fhicl::ParameterSet&, art::ActivityRegistry&) {
}
  
//**********************************************************************

FixedZeroSuppressService::FixedZeroSuppressService() { }

//**********************************************************************

int FixedZeroSuppressService::
filter(const AdcCountVector& sigs, Channel, AdcPedestal&, AdcFilterVector& keep) const {
  keep.clear();
  keep.resize(sigs.size(), true);
  return 0;
}

//**********************************************************************

ostream& FixedZeroSuppressService::print(ostream& out, string prefix) const {
  out << prefix << "FixedZeroSuppressService" << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(FixedZeroSuppressService, AdcSuppressService)

//**********************************************************************
