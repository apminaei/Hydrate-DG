#ifndef BM_PHASECHANGE_0d_INCLUDE_PROBLEM_HH_
#define BM_PHASECHANGE_0d_INCLUDE_PROBLEM_HH_

// #define PARALLEL
#include "../extras/ParameterTraits.hh"
#include "../extras/Evaluation.hh"

#include "../operators/indices.hh"

/*******************************************/
// PROBLEM SPECIFICATION
#include"characteristic_values.hh"
#include"grids/grid.hh"
#define STATEINDEPENDENTPROPERTIES
#include"parameters.hh"
#include"properties/include_properties.hh"
#include"initial_conditions.hh"
#include"boundary_conditions.hh"
/*******************************************/
#ifdef PLOT_VELOCITIES
//#include "../operators/phase_velocity.hh"
//#include "../operators/operator_l2projection.hh"
#endif
#include "../operators/postprocess.hh"
#include "../operators/Initial.hh"
// #include "../operators/LocalOperator_Sh1.hh"
// #include "../operators/LocalOperator_T1.hh"
// #include "../operators/TimeOperator_T1.hh"
// #include "../operators/TimeOperator_Sh1.hh"
// #include "../operators/LocalOperator_2comps.hh"
// #include "../operators/TimeOperator_2comps.hh"
// #include "driver_Sh.hh"
#include "../operators/LocalOperator.hh"
#include "../operators/TimeOperator.hh"
#include "driver.hh"

#endif /* BM_PHASECHANGE_0d_INCLUDE_PROBLEM_HH_ */
