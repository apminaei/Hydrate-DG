// CONSTITUTIVE AND MATERIAL PROPERTIES
#include"salt.hh"
#include"H2O.hh"
#include"CH4.hh"
#include"eosCH4.hh"
#include"hydrate.hh"
#include"mixture.hh"
#include"soil.hh"
#include"hydraulic_properties.hh"
#include"hydrate_phase_change.hh"

template<typename GV, typename PTree>
class Properties
{
private:
	  const GV& gv;
	  const PTree& ptree;
	  double *time;
	  double *dt;

	  const static int dim = GV::dimension;
	  constexpr static double eps = 1.e-6;
	  
public:

	//PARAMETERS AND PROPERTIES
  	Indices index;
  	CharacteristicValues characteristicValue;
  	MeshParameters<PTree> mesh;
  	Parameters<PTree> parameter;
  	Methane gas;
#ifdef FRAUENHOFER_MODEL
  	FrauenhoferFunction eos;
#elif defined(PENG_ROBINSON_EOS)
  	PengRobinson eos;
#else
  	BaseEoS eos;
#endif
  	Water water;
  	Salt salt;
  	Mixture mixture;
  	Hydrate hydrate;
  	Soil<GV,Parameters<PTree>> soil;
  	HydraulicProperties<GV,Parameters<PTree>> hydraulicProperty;
  	HydratePhaseChangeKinetics<GV,Parameters<PTree>> kinetics;
  	
  	//! construct from grid view
  	Properties ( const GV& gv_ , 
  				 const PTree& ptree_, 
				 double *time_,
				 double *dt_  )
	: gv( gv_ ),
	  ptree(ptree_),
	  time(time_),
	  dt(dt_),
	  mesh(ptree_),
	  parameter(ptree_,time_,dt_),
	  soil(gv_,parameter),
	  hydraulicProperty(gv_,parameter),
	  kinetics(gv_,parameter)
  	{}

	/******************************************************************************/

  	void ReportStatistics( std::string file_name,
  						   double time /*s*/,
						   double dt /*s*/,
						   int total_newton_iterations,
						   double clock_time_elapsed /*s*/ ) {

  		std::fstream result;

  		if(time == 0. ){
  			result.open(file_name, std::fstream::out | std::fstream::trunc);
  			result	<< "time [s]" << '\t'
  					<< "dt [s]"	<< '\t'
					<< "total no. of newton iterations" << '\t'
					<< "clock time [s]"
  					<< std::endl;
  			result.close();
  		}

  		result.open(file_name, std::fstream::app);
  		double t_new = time+dt;

		result	<< time	<< '\t'
				<< dt	<< '\t'
				<< total_newton_iterations << '\t'
				<< clock_time_elapsed
				<< std::endl;
		result.close();
  	}

	/******************************************************************************/



	/******************************************************************************/

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

};
