include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )


art_make_library( LIBRARY_NAME ProtoDUNETriggeredBeamUtils
		  SOURCE ProtoDUNETriggeredBeamUtils.cxx
		  LIBRARIES
                       fhiclcpp::fhiclcpp
		       ROOT::Core
		       ROOT::Physics
	        )

cet_make_exec( NAME runEventCounter
	       SOURCE runEventCounter.c++
	       LIBRARIES ProtoDUNETriggeredBeamUtils 
	    )
cet_make_exec( NAME reformat_g4bl
	       SOURCE reformat_g4bl.c++
	       LIBRARIES
            ${ROOT_LIBRARIES}
            dunecore_DuneObj
	    )
simple_plugin(H4BeamInput "source"
  art::Framework_Services_Registry
  dunesim::EventGenerator_Utils_H4BeamFileService_service
)

simple_plugin(H4BeamFileService "service"
              ${ART_PERSISTENCY_ROOTDB}
              art::Framework_Services_Registry
              ${persistency_lib}
              SQLite::SQLite3
              messagefacility::MF_MessageLogger
              ${CETLIB_LIBS}
	      ProtoDUNETriggeredBeamUtils
)

simple_plugin( H4BeamFileTester "module"
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities ${canvas}
                        messagefacility::MF_MessageLogger
                        ${hep_concurrency}
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib ${cetlib_except}
                        ${CLHEP}
                        ${ROOT_BASIC_LIB_LIST}
                        dunesim_EventGenerator_ProtoDUNEbeamDataProducts
                        ProtoDUNETriggeredBeamUtils
                        dunesim::EventGenerator_Utils_H4BeamFileService_service
)
install_scripts()
install_fhicl()
