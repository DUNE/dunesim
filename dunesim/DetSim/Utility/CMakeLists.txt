# dunetpc/dune/DetSim/Utility/CMakeLists.txt

# David Adams
# January 2017

FILE( GLOB src_files *.cxx )

art_make_library(
  LIBRARY_NAME dunesim::DetSim
  SOURCE ${src_files}
  LIBRARIES
  dunecore::DuneInterface_Data
)

install_headers()
install_source()
