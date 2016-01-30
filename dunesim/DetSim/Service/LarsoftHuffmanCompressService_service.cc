// LarsoftHuffmanCompressService.cxx

#include "dune/DetSim/Service/LarsoftHuffmanCompressService.h"
#include "fhiclcpp/ParameterSet.h"
#include "RawData/raw.h"

using std::string;
using std::ostream;
using std::endl;

//**********************************************************************

LarsoftHuffmanCompressService::LarsoftHuffmanCompressService() { }

//**********************************************************************

LarsoftHuffmanCompressService::
LarsoftHuffmanCompressService(const fhicl::ParameterSet&, art::ActivityRegistry&) { }
  
//**********************************************************************

int LarsoftHuffmanCompressService::
compress(AdcCountVector& sigs, const AdcFilterVector& keep, AdcCount offset,
         raw::Compress_t& comp) const {
  const string myname = "LarsoftHuffmanCompressService::compress: ";
  cout << myname << "Not yet implemented." << endl;
  abort();
  comp = raw::kHuffman;
  return 0;
}

//**********************************************************************

ostream& LarsoftHuffmanCompressService::print(ostream& out, string prefix) const {
  out << prefix << "LarsoftHuffmanCompressService:" << endl;
  return out;
}

//**********************************************************************

DEFINE_ART_SERVICE_INTERFACE_IMPL(LarsoftHuffmanCompressService, AdcCompressService)

//**********************************************************************
