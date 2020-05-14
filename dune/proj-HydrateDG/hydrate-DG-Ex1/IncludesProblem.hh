/*
 * IncludesProblem.hh
 *
 *  Created on: Sept 22, 2016
 *      Author: shubhangi
 *
 *  All files specifically related to the problem are included here, i.e.,
 *  - Problem description and parameters
 *  - Local and Time operators
 *  - Main file
 *  - Report template
 */

#ifndef INCLUDESPROBLEM_HH_
#define INCLUDESPROBLEM_HH_


#include"indices.hh"
#include"./properties/characteristicValues.hh"
#include"./grids/grid.hh"

/*PARAMETERS AND PROPERTIES*/
#include"./properties/parameters.hh"
#include"./properties/methane.hh"
#include"./properties/water.hh"
#include"./properties/mixture.hh"
#include"./properties/eosPengRobinson.hh"
#include"./properties/hydraulicProperties.hh"
#include"./properties/hydrate.hh"
#include"./properties/soil.hh"
#include"./properties/reactionKinetics.hh"
#include"IncludeProperties.hh"

/*MAIN*/

#include"initial_conditions.hh"
#include"../operators/Initial.hh"
#include"boundary_conditions.hh"
#include"../operators/LocalOperator.hh"
#include"../operators/TimeOperator.hh"

#include"driver.hh"


#endif /* INCLUDESPROBLEM_HH_ */
