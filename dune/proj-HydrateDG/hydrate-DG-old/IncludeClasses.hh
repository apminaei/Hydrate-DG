/*
 * IncludeClasses.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef INCLUDECLASSES_HH_
#define INCLUDECLASSES_HH_

class IncludeClasses{
public:
	  Indices index;
	  ModelSubSystems modelSubSystem;
	  IncludeModels includeModel;
	  CharacteristicValues characteristicValue;
	  TimeStepControl timeStepControl;
	  TimeSteppingStrategies timeSteppingStrategy;
	  GravityVector gravityVector;
	  /*************************************************************************/
	  /*				PROBLEM SPECS				                           */
	  typedef IncludeProblemSpecifications::ProblemSpecifications ProblemSpecs;
	  ProblemSpecs problemSpecs;
	  /*************************************************************************/
	  Soil soil;
	  Methane methane;
	  Water water;
	  Hydrate hydrate;
	  HydraulicProperties hydraulicProperty;
	  Mixture mixture;
	  ReactionKinetics reactionKinetics;
	  PengRobinson eos;

};

#endif /* INCLUDECLASSES_HH_ */
