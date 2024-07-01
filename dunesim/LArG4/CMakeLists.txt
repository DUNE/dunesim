
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
        )

# install_headers()
install_fhicl()
install_source()


