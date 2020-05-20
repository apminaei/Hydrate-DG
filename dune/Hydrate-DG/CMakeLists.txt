#install headers
install(FILES IncludesDUNE.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/Hydrate-DG)
add_subdirectory("operators")
add_subdirectory("akerbp2D_pockmark")
add_subdirectory("extras")
