# dunetpc/dune/larsim_DetSim/CMakeLists.txt

# David Adams
# December 2015

# We can drop the MODULE_LIBRARIES block after the old modules
# SimWireDUNE35t_module.cc and SimWireDUNE10kt_module.cc are
# removed.

art_make( BASENAME_ONLY
          LIBRARY_NAME larsim_DetSim
          LIB_LIBRARIES
          MODULE_LIBRARIES larcore_Geometry
	                dune_Geometry
			larcore_Geometry_Geometry_service
			larsim_Simulation
	                lardata_Utilities
			larevt_Filters
			lardata_RawData
			lardata_Utilities_TimeService_service
                        SignalShapingServiceDUNE10kt_service
			SignalShapingServiceDUNEDPhase_service
			SignalShapingServiceDUNE35t_service
                        ${SEEDSERVICE_SERVICE} # artextensions
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
			${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
        )

add_subdirectory(Service)
add_subdirectory(Module)

install_headers()
install_fhicl()
install_source()

