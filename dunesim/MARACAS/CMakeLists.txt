# MARACAS ProtoDUNE-VD electric-field and position-distortion providers,
# read from PyTables HDF5 maps, exposed as ElectricFieldService /
# PositionDistorterService implementations.

cet_make_library(LIBRARY_NAME MARACAS
  SOURCE
  VDMaracasGrid.cxx
  VDMaracasElectricField.cxx
  VDMaracasDistorter.cxx
  VDMaracasDistortion.cxx
  LIBRARIES
  PUBLIC
  lardataalg::DetectorInfo
  larcoreobj::SimpleTypesAndConstants
  ROOT::Physics
  PRIVATE
  cetlib::cetlib
  cetlib_except::cetlib_except
  fhiclcpp::fhiclcpp
  fhiclcpp::types
  HDF5::HDF5
)

cet_build_plugin(VDMaracasElectricFieldService art::service
  LIBRARIES PRIVATE
  MARACAS
  lardata::ElectricFieldService
  lardataalg::DetectorInfo
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)

# cet_build_plugin(VDMaracasPositionDistorterService art::service
#   LIBRARIES PRIVATE
#   MARACAS
#   lardata::PositionDistorterService
#   lardataalg::DetectorInfo
#   fhiclcpp::fhiclcpp
#   cetlib_except::cetlib_except
# )

# IDistortion chain link (art tool) for ChainedDistorter.
cet_build_plugin(VDMaracasDistortion art::tool
  LIBRARIES PRIVATE
  MARACAS
  lardataalg::DetectorInfo
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
)

install_headers()
install_fhicl()
install_source()
