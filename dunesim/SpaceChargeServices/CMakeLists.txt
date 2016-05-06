simple_plugin(SpaceChargeServiceProtoDUNE "service"
              dune_SpaceCharge
	      ${ART_FRAMEWORK_CORE}
	      ${ART_FRAMEWORK_PRINCIPAL}
	      ${ART_FRAMEWORK_SERVICES_REGISTRY}
	      ${ART_PERSISTENCY_COMMON}
	      ${ART_PERSISTENCY_PROVENANCE}
	      ${ART_FRAMEWORK_SERVICES_BASIC}
	      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
	      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
	      ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
	      ${ART_UTILITIES}
	      ${MF_MESSAGELOGGER}
	      ${MF_UTILITIES}
  	      ${CETLIB}
   	      ${ROOT_BASIC_LIB_LIST}
              
)

install_headers()
install_fhicl()
install_source()
