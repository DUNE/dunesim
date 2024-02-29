
art_make( BASENAME_ONLY
          MODULE_LIBRARIES larcorealg::Geometry
			larcore::Geometry_Geometry_service
			larsim::Simulation nug4::ParticleNavigation lardataobj::Simulation
			nusimdata::SimulationBase
			ROOT::Core
			ROOT::Tree
			art_root_io::tfile_support
			art::Persistency_Common
			art::Persistency_Provenance
			art::Utilities canvas::canvas
			messagefacility::MF_MessageLogger
	                fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except::cetlib_except
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

# install_headers()
install_fhicl()
install_source()


