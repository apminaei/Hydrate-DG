#install headers
install(FILES IncludesDUNE.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/Hydrate-DG)
add_subdirectory("operators")
add_subdirectory("BM_phasechange_0d")
#add_subdirectory("problem_phaseChange_p0")
add_subdirectory("extras")
