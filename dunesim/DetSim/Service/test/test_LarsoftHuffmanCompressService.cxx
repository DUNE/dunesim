// test_LarsoftHuffmanCompressService.cxx

// David Adams
// February 2015
//
// Test LarsoftHuffmanCompressService.

#include "lardataobj/RawData/raw.h"
#include "../LarsoftHuffmanCompressService.h"
#include <string>
#include <iostream>
#include <sstream>

using std::string;
using std::cout;
using std::endl;
using std::istringstream;

typedef std::vector<short> AdcVector;

#undef NDEBUG
#include <cassert>

int test_LarsoftHuffmanCompressService(bool useBlock, bool useHuffman, int logLevel) {
  const string myname = "test_LarsoftHuffmanCompressService: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << line << endl;
  cout << myname << "Testing with UseBlock=" << useBlock
       << ", useHuffman=" << useHuffman
       << ", LogLevel=" << logLevel
       << endl;
  raw::Compress_t               ctypeexp = raw::kNone;
  if ( useBlock && useHuffman ) ctypeexp = raw::kZeroHuffman;
  else if ( useBlock )          ctypeexp = raw::kZeroSuppression;
  else if ( useHuffman )        ctypeexp = raw::kHuffman;
  raw::Compress_t ctypeout = raw::kNone;

  const unsigned int nadc = 3200;
  AdcVector adc0(nadc, 50);
  cout << myname << "     Initial ADC size: " << adc0.size() << endl;
  assert( adc0.size() == nadc );

  AdcFilterVector keep(nadc, true);
  AdcCount offset = 1800;

  // Compress.
  AdcVector adc1(adc0);
  LarsoftHuffmanCompressService cmp(useBlock, useHuffman, logLevel);
  cout << myname << "Compressing." << endl;
  int rstat = cmp.compress(adc1, keep, offset, ctypeout);
  cout << myname << "        Return status: " << rstat << endl;
  cout << myname << "  Compressed ADC size: " << adc1.size() << endl;
  cout << myname << "     Compression type: " << ctypeexp << endl;
  assert( rstat == 0 );
  assert( ctypeout == ctypeexp );
  assert( adc1.size() > 0 );
  //assert( adc1.size() <= adc0.size() );

  // Check result.
  AdcVector adc2(nadc, 0);
  cout << myname << "Uncompressing." << endl;
  int uoffset = offset;
  raw::Uncompress(adc1, adc2, uoffset, ctypeout);
  cout << myname << "Uncompressed ADC size: " << adc2.size() << endl;
  assert( adc2.size() ==  adc0.size() );
  int nbad = 0;
  for ( unsigned int iadc=0; iadc<adc0.size(); ++iadc ) {
    if ( adc0[iadc] != adc2[iadc] ) {
      cout << myname << "  Mismatch for entry " << iadc << ": " << adc0[iadc] << " != " << adc2[iadc] << endl;
      ++nbad;
      if ( nbad > 20 ) break;
    }
  }
  assert( nbad == 0 );
  assert( uoffset == offset );

  cout << myname << "Done." << endl;
  return 0;
}

int main(int argc, char* argv[]) {
  int logLevel = 1;
  if ( argc > 1 ) {
    istringstream ssarg(argv[1]);
    ssarg >> logLevel;
  }
  test_LarsoftHuffmanCompressService(false, false, logLevel);
  test_LarsoftHuffmanCompressService(false, true,  logLevel);
  test_LarsoftHuffmanCompressService(true,  false, logLevel);
  test_LarsoftHuffmanCompressService(true,  true,  logLevel);
  return 0;
}
