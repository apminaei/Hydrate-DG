#ifndef EX2_INCLUDE_PROBLEM_HH_
#define EX2_INCLUDE_PROBLEM_HH_

// #define PARALLEL
#include "../extras/ParameterTraits.hh"
#include "../extras/Evaluation.hh"
#include "../extras/limiter.hh"

#include "../operators/indices.hh"

#define PENG_ROBINSON_EOS
/*******************************************/
// PROBLEM SPECIFICATION
#include "characteristic_values.hh"
#include "grids/grid.hh"
#include "parameters.hh"
#include "properties/include_properties.hh"
#include "initial_conditions.hh"
#include "boundary_conditions.hh"
/*******************************************/
#include "../operators/Initial.hh"
#include "../operators/LocalOperator.hh"
#include "../operators/TimeOperator.hh"
#include "../operators/postprocess.hh"
#include "driver.hh"
#endif /* BM_PHASECHANGE_0d_INCLUDE_PROBLEM_HH_ */
