
art_make( BASENAME_ONLY
          MODULE_LIBRARIES larcorealg_Geometry
			larcore_Geometry_Geometry_service
			larsim_Simulation nug4::ParticleNavigation lardataobj_Simulation
			nusimdata::SimulationBase
			ART_FRAMEWORK_CORE
			ART_FRAMEWORK_PRINCIPAL
			ART_FRAMEWORK_SERVICES_REGISTRY
			ART_ROOT_IO_TFILE_SUPPORT
			ROOT::Core
			ART_ROOT_IO_TFILESERVICE_SERVICE
			art::Persistency_Common canvas
			art::Persistency_Provenance canvas
			art::Utilities canvas
			messagefacility::MF_MessageLogger
			
	                fhiclcpp::fhiclcpp
			cetlib::cetlib cetlib_except
                        ROOT_GEOM
                        ROOT_XMLIO
                        ROOT_GDML
			ROOT_BASIC_LIB_LIST
        )

# install_headers()
install_fhicl()
install_source()


