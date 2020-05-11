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


#include"./parameters/indices.hh"
#include"./parameters/characteristicValues.hh"
#include"./grids/grid.hh"
/*PROBLEM FILES*/
#include"./problemSpecs/problem_20160922_test1.hh"

#include"./problemSpecs/includeProblemSpecifications.hh"

/*PARAMETERS AND PROPERTIES*/
#include"./parameters/parameters.hh"
#include"./parameters/methane.hh"
#include"./parameters/water.hh"
#include"./parameters/mixture.hh"
#include"./parameters/eosPengRobinson.hh"
#include"./parameters/hydraulicProperties.hh"
#include"./parameters/hydrate.hh"
#include"./parameters/soil.hh"
#include"./parameters/reactionKinetics.hh"
#include"IncludeClasses.hh"

/*MAIN*/
#include"../operators/Initial.hh"
#include"../operators/LocalOperator.hh"
#include"../operators/TimeOperator.hh"
#include"driver.hh"


#endif /* INCLUDESPROBLEM_HH_ */
