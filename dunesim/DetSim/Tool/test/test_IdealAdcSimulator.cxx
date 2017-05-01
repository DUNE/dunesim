// test_IdealAdcSimulator.cxx
//
// David Adams
// April 2017
//
// Test IdealAdcSimulator.

#include <string>
#include <iostream>
#include <fstream>
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "art/Utilities/make_tool.h"
#include "dune/DuneInterface/AdcSimulator.h"

#undef NDEBUG
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ofstream;
using fhicl::ParameterSet;
using ToolPtr = std::unique_ptr<AdcSimulator>;
using Count = AdcSimulator::Count;

//**********************************************************************

template<class T1, class T2>
void checkEqual(T1 x1, T2 x2, string msg ="Assert failed") {
  if ( x1 != x2 ) {
    cout << msg << ": " << x1 << " != " << x2 << endl;
    assert(false);
  }
}

void checkCount(const ToolPtr& ptool, double vin, Count countExpected) {
  Count countActual = ptool->count(vin);
  cout << "  ADC(" << vin << ") = " << countActual << endl;
  checkEqual(countActual, countExpected);
}

//**********************************************************************

int test_IdealAdcSimulator(bool useExistingFcl) {
  const string myname = "test_IdealAdcSimulator: ";
#ifdef NDEBUG
  cout << myname << "NDEBUG must be off." << endl;
  abort();
#endif
  string line = "-----------------------------";

  cout << myname << line << endl;
  string fclfile = "test_IdealAdcSimulator.fcl";
  if ( ! useExistingFcl ) {
    cout << myname << "Creating top-level FCL." << endl;
    ofstream fout(fclfile.c_str());
    fout << "mytool: {" << endl;
    fout << "  tool_type: IdealAdcSimulator" << endl;
    fout << "  Vsen: 2.0" << endl;
    fout << "  Nbit: 12" << endl;
    fout << "}" << endl;
    fout.close();
  } else {
    cout << myname << "Using existing top-level FCL." << endl;
  }

  // Create pset
  cout << myname << line << endl;
  cout << myname << "Creating parameter set." << endl;
  //putenv(const_cast<char*>("FHICL_FILE_PATH=./test:."));
  cet::filepath_lookup policy("FHICL_FILE_PATH");
  fhicl::intermediate_table tbl;
  fhicl::parse_document(fclfile, policy, tbl);
  ParameterSet psTop;
  fhicl::make_ParameterSet(tbl, psTop);

  cout << myname << line << endl;
  cout << myname << "Checking tool type." << endl;
  ParameterSet psTool = psTop.get<fhicl::ParameterSet>("mytool");
  cout << "Tool type: " << psTool.get<string>("tool_type") << endl;
  assert( psTool.get<string>("tool_type") == "IdealAdcSimulator" );

  cout << myname << line << endl;
  cout << myname << "Instantiate tool." << endl;
  std::unique_ptr<AdcSimulator> ptool = art::make_tool<AdcSimulator>(psTool);
  assert( ptool != nullptr );
  double vin = 1200;
  cout << "ADC count is " << ptool->count(vin) << endl;
  checkCount(ptool, -999.0, 0);
  checkCount(ptool, -999.0, 0 );
  checkCount(ptool, -0.01, 0 );
  checkCount(ptool, 0.0, 0 );
  checkCount(ptool, 0.2, 0 );
  checkCount(ptool, 1.0, 1);
  checkCount(ptool, 201.0, 101 );
  checkCount(ptool, 2001.0, 1001 );
  checkCount(ptool, 8187.0, 4094 );
  checkCount(ptool, 8189.0, 4095 );
  checkCount(ptool, 8190.0, 4095 );
  checkCount(ptool, 8191.0, 4095 );
  checkCount(ptool, 999999.0, 4095 );

  cout << myname << line << endl;
  cout << myname << "Done." << endl;
  return 0;
}

//**********************************************************************

int main(int argc, char* argv[]) {
  bool useExistingFcl = false;
  if ( argc > 1 ) {
    string sarg(argv[1]);
    if ( sarg == "-h" ) {
      cout << "Usage: " << argv[0] << " [ARG]" << endl;
      cout << "  If ARG = true, existing FCL file is used." << endl;
      return 0;
    }
    useExistingFcl = sarg == "true" || sarg == "1";
  }
  return test_IdealAdcSimulator(useExistingFcl);
}

//**********************************************************************
