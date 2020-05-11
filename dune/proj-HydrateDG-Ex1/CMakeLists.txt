#install headers
install(FILES IncludesDUNE.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/proj-HydrateDG-Ex1)
add_subdirectory("operators")
add_subdirectory("hydrate-DG-old")
