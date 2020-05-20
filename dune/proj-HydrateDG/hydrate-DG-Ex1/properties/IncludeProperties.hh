/*
 * Properties.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#include"methane.hh"
#include"salt.hh"
#include"water.hh"
#include"mixture.hh"
#include"eosPengRobinson.hh"
#include"hydraulicProperties.hh"
#include"hydrate.hh"
#include"soil.hh"
#include"reactionKinetics.hh"


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
		TimeStepControl timeStepControl;
		TimeSteppingStrategies timeSteppingStrategy;
		
		MeshParameters<PTree> mesh;
		Parameters<PTree> parameter;
		Soil<GV, PTree> soil;
		Methane methane;
		Water water;
		Hydrate hydrate;
		HydraulicProperties< GV, PTree> hydraulicProperty;
		Mixture<PTree> mixture;
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
		  hydraulicProperty(gv_, ptree_),
		  mixture(ptree_)
	  	{}

};

//#endif /* INCLUDECLASSES_HH_ */
