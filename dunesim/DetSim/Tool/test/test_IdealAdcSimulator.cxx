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
  assert( ptool->count(vin) == 1234 );

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
