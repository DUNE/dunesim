simple_plugin(RandomChannelStatusService "service"
                dune_DetSim
                larcorealg_Geometry
                larsim_Simulation nutools_ParticleNavigation lardataobj_Simulation
                lardata_Utilities
                SignalShapingServiceDUNE_service
                art_Utilities canvas
                cetlib cetlib_except
                ${ROOT_BASIC_LIB_LIST}
             )
