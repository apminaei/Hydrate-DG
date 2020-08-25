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
	//   double *time;
	//   double *dt;

	  const static int dim = GV::dimension;
	  constexpr static double eps = 1.e-6;
	  
public:

	//PARAMETERS AND PROPERTIES
  	Indices index;
  	CharacteristicValues characteristicValue;
  	MeshParameters<PTree> mesh;
  	Parameters<PTree> parameter;
  	Methane<PTree> gas;
#ifdef FRAUENHOFER_MODEL
  	FrauenhoferFunction eos;
#elif defined(PENG_ROBINSON_EOS)
  	PengRobinson eos;
#else
  	BaseEoS eos;
#endif
  	Water<PTree> water;
  	Salt salt;
  	Mixture<PTree> mixture;
  	Hydrate<PTree> hydrate;
  	Soil<GV,Parameters<PTree>> soil;
  	HydraulicProperties<GV,Parameters<PTree>> hydraulicProperty;
  	HydratePhaseChangeKinetics<GV,PTree> kinetics;
  	
  	//! construct from grid view
  	Properties ( const GV& gv_ , 
  				 const PTree& ptree_)
	: gv( gv_ ),
	  ptree(ptree_),
	  mesh(ptree_),
	  parameter(ptree_),
	  gas(ptree_),
	  water(ptree_),
	  mixture(ptree_),
	  hydrate(ptree_),
	  soil(gv_,parameter),
	  hydraulicProperty(gv_,parameter),
	  kinetics(gv_,ptree_)
  	{}
	
	double dt_initial = ptree.get("time.dt_initial",(double)1.);
	int time_interval = ptree.get("output.time_interval",(double)1);
	/******************************************************************************/

  	void ReportStatistics( std::string file_name,
  						   double time /*s*/,
						   double dt /*s*/,
						   int total_newton_iterations,
						   double clock_time_elapsed /*s*/) {

  		std::fstream result;

  		if(time == 0. ){
  			result.open(file_name, std::fstream::out | std::fstream::trunc);
  			result	<< "time [s]" << '\t'
  					<< "dt [s]"	<< '\t'
					<< "total no. of newton iterations" << '\t'
					<< "clock time [s]" << '\t'
  					<< mesh.X_cells << "*" << mesh.Z_cells << '\t'
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
  	void ReportParameters( std::string file_name,
  						   std::string method_g,
							std::string method_w,
							std::string method_T,
							std::string method_x,
							std::string method_y,
							double alpha_g , double alpha_w, double alpha_s, double alpha_T, double alpha_x, double alpha_y) {

  		std::fstream result;

  		
  			result.open(file_name, std::fstream::out | std::fstream::trunc);
  			result	<< "penalty coeff.  " << alpha_g << '\t' << alpha_w << '\t'<< alpha_s << '\t'<< alpha_T << '\t'<< alpha_x << '\t'<< alpha_y << '\t'
					<< " S=1, N=0, I=-1,  " << method_g << '\t'<< method_w << '\t'<< method_T << '\t'<< method_x << '\t'<< method_y 
					<< std::endl;
  			result.close();
  		

  	}


	/******************************************************************************/

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

};
