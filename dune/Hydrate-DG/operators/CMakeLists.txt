#install headers
install(FILES LocalOperator.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/Hydrate-DG/operators)
install(FILES TimeOperator.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/Hydrate-DG/operators)

# install(FILES LocalOperatorLinear.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/Hydrate-DG/operators)
# install(FILES TimeOperatorLinear.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/Hydrate-DG/operators)
install(FILES Initial.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/Hydrate-DG/operators)
install(FILES indices.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/Hydrate-DG/operators)

install(FILES postprocess.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/Hydrate-DG/operators)
# install(FILES TimeOperator_2comps.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/Hydrate-DG/operators)
# install(FILES LocalOperator_Sh1.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/Hydrate-DG/operators)
# install(FILES TimeOperator_Sh1.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/Hydrate-DG/operators)
# install(FILES LocalOperator_T1.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/Hydrate-DG/operators)
# install(FILES TimeOperator_T1.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/Hydrate-DG/operators)
install(FILES operator_l2projection.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/Hydrate-DG/operators)