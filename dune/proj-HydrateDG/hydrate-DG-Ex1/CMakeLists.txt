#install headers
install(FILES characteristicValues.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/proj-HydrateDG/hydrate-DG-Ex1)
install(FILES driver.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/proj-HydrateDG/hydrate-DG-Ex1)
install(FILES IncludesProblem.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/proj-HydrateDG/hydrate-DG-Ex1)
install(FILES boundary_conditions.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/proj-HydrateDG/hydrate-DG-Ex1)
install(FILES initial_conditions.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/proj-HydrateDG/hydrate-DG-Ex1)
install(FILES parameters.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/proj-HydrateDG/hydrate-DG-Ex1)

add_subdirectory("properties")
add_subdirectory("grids")
add_subdirectory("inputs")
add_subdirectory("outputs")