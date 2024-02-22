
art_make( BASENAME_ONLY
          MODULE_LIBRARIES larcorealg_Geometry
			larcore_Geometry_Geometry_service
			larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			nusimdata::SimulationBase
			ROOT::Core
			ROOT::Tree
			art_root_io::tfile_support
			art::Persistency_Common canvas
			art::Persistency_Provenance canvas
			art::Utilities canvas
			messagefacility::MF_MessageLogger
	                fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
                        ROOT::Geom
                        ROOT::XMLIO
                        ROOT::Gdml
        ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.22.02) - deprecated: use art_make_library(), art_dictonary(), and cet_build_plugin() with explicit source lists and plugin base types

# install_headers()
install_fhicl()
install_source()


