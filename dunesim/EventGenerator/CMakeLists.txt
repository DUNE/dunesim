include_directories( $ENV{IFDHC_FQ_DIR}/inc ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible
 

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - avoid directory-scope functions: use target_link_libraries() with target semantics or target_include_directories() whenever possible


#art_make_library( LIBRARY_NAME ProtoDUNETriggeredBeamUtils
#		  SOURCE ProtoDUNETriggeredBeamUtils.cxx
#		  LIBRARIES
#	    )

art_make( BASENAME_ONLY
          MODULE_LIBRARIES 
                        larcoreobj::SummaryData
                        larcorealg::Geometry
                        larcore::Geometry_Geometry_service
                        nurandom::RandomUtils_NuRandomService_service
                        larsim::Simulation lardataobj::Simulation
                        lardataobj::RecoBase
                        lardata::Utilities
                        nusimdata::SimulationBase
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities canvas::canvas
                        messagefacility::MF_MessageLogger
                        hep_concurrency::hep_concurrency
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except::cetlib_except
			CLHEP::CLHEP
			SQLite::SQLite3
                        ROOT::EG
			ROOT::Tree
                        ifdh::ifdh
                        dunesim_EventGenerator_ProtoDUNEbeamDataProducts
                        dunecore::ArtSupport
                        ProtoDUNETriggeredBeamUtils
 ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types


install_headers()
install_fhicl()
install_source()
install_scripts()

add_subdirectory(CRY)
add_subdirectory(GENIE)
add_subdirectory(Radiological)
add_subdirectory(CORSIKAprotodunedp)
add_subdirectory(NEUT)
add_subdirectory(MUSUN)
add_subdirectory(NuE)
add_subdirectory(ProtoDUNEbeamDataProducts)
add_subdirectory(Utils)
