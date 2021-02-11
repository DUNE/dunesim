# dunetpc/dune/DetSim/Tool/CMakeLists.txt

# David Adams
# April 2017

# Add plugin for each tool.

simple_plugin(IdealAdcSimulator   "tool"
  art_Utilities canvas
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  )

simple_plugin(NoisyAdcSimulator   "tool"
  nurandom_RandomUtils_NuRandomService_service
  ${CLHEP}
  art_Utilities canvas
  cetlib cetlib_except
  ${MF_MESSAGELOGGER}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${ROOT_BASIC_LIB_LIST}
  ${ART_FRAMEWORK_CORE}
  )

add_subdirectory(test)
