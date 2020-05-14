/*
 * Properties.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef INCLUDECLASSES_HH_
#define INCLUDECLASSES_HH_
template<typename GV, typename PTree>
class Properties{

	private:
	  const GV& gv;
	  const PTree& ptree;

	public:
		Indices index;
		ModelSubSystems modelSubSystem;
		IncludeModels includeModel;
		CharacteristicValues characteristicValue;
		MeshParameters<PTree> mesh;
		Parameters<PTree> parameter;
		TimeStepControl timeStepControl;
		TimeSteppingStrategies timeSteppingStrategy;
		//GravityVector gravityVector;
		/*************************************************************************/
		/*				PROBLEM SPECS				                           */
		//typedef IncludeProblemSpecifications::ProblemSpecifications ProblemSpecs;
		//ProblemSpecs problemSpecs;
		/*************************************************************************/
		
		Soil<GV, PTree> soil;
		Methane methane;
		Water water;
		Hydrate hydrate;
		HydraulicProperties< GV, PTree> hydraulicProperty;
		Mixture mixture;
		ReactionKinetics< GV, PTree> reactionKinetics;
		PengRobinson eos;
	  
	  	//! construct from grid view
	  	Properties ( const GV& gv_, const PTree& ptree_ )
		: gv( gv_ ),
		  ptree(ptree_),
		  mesh(ptree_),
		  parameter(ptree_),
		  soil(gv_,ptree_),
		  reactionKinetics(gv_, ptree_),
		  hydraulicProperty(gv_, ptree_)
	  	{}

};

#endif /* INCLUDECLASSES_HH_ */
