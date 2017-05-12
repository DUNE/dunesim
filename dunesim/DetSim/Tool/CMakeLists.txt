# dunetpc/dune/DetSim/Tool/CMakeLists.txt

# David Adams
# April 2017

# Add plugin for each tool.

simple_plugin(IdealAdcSimulator   "tool"
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )

add_subdirectory(test)
