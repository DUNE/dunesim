cet_build_plugin(GenExtCounterFilter art::module
		      lardataobj_RawData
		      lardataobj_RecoBase
		      dune_Geometry
		        larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			nusimdata_SimulationBase
                        ART_FRAMEWORK_CORE
			ART_FRAMEWORK_PRINCIPAL
			ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
 			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			messagefacility::MF_MessageLogger
			
			cetlib cetlib_except
			ROOT_BASIC_LIB_LIST
              BASENAME_ONLY
)

cet_build_plugin(LArG4ExtCounterFilter art::module
		      lardataobj_RawData
		      lardataobj_RecoBase
		      dune_Geometry
		        larcorealg_Geometry
			larcore_Geometry_Geometry_service
			lardata_Utilities
			nusimdata_SimulationBase
                        ART_FRAMEWORK_CORE
			ART_FRAMEWORK_PRINCIPAL
			ART_FRAMEWORK_SERVICES_REGISTRY
                        ART_ROOT_IO_TFILE_SUPPORT
                        ROOT::Core
                        ART_ROOT_IO_TFILESERVICE_SERVICE
 			art_Persistency_Common canvas
			art_Persistency_Provenance canvas
			art_Utilities canvas
			messagefacility::MF_MessageLogger
			
			cetlib cetlib_except
			ROOT_BASIC_LIB_LIST
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
install_scripts()
