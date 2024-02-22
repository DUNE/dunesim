art_make(
          MODULE_LIBRARIES
                        larcoreobj::SummaryData
                        larcore::Geometry_Geometry_service
                        larcorealg::Geometry
                        larsim_Simulation lardataobj_Simulation
                        nusimdata::SimulationBase
                        ROOT::Core
                        art::Persistency_Common
                        art::Persistency_Provenance
                        art::Utilities
                        canvas::canvas
                        messagefacility::MF_MessageLogger
                        fhiclcpp::fhiclcpp
                        cetlib::cetlib cetlib_except
			CLHEP::CLHEP
                        ROOT::Geom
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

install_headers()
install_fhicl()
install_source()
