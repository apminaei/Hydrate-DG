/*
 * parameters.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef PARAMETERS_PARAMETERS_HH_
#define PARAMETERS_PARAMETERS_HH_

template<typename PTree>
class Parameters
{
private:
	const PTree& ptree;
	MeshParameters<PTree> mesh;
//	const double pi = boost::math::constants::pi<double>();
	const double pi = 3.14159265358979323846;
	const static int dim = MeshParameters<PTree>::dimension;

	double Sg_t0;
	double Pw_t0;
	double Pg_t0;
	double Sh_t0;
	double XCH4_t0;
	double YH2O_t0;
	double Sg_x0;
	double Pw_x0;
	double Sgin_x0;

	int numMaterials;
	int numProps;
	std::vector<std::vector<double> > prop;

	double ref_salinity;
	double ref_saltconcentration;
	double ref_temperature;

	bool gravity_flag;
	double g_magnitude;

public:

  //! constructor
  Parameters (const PTree& ptree_)
  :ptree(ptree_),
   mesh(ptree_)
  {
		Sg_t0 = ptree.get("initial.Sg",(double)0.001);
		Pw_t0 = ptree.get("initial.Pw",(double)1.e6);
		Pg_t0 = ptree.get("initial.Pg",(double)0.001);
		//Sw_t0 = ptree.get("initial.Sw",(double)1.e6);
		Sh_t0 = ptree.get("initial.Sh",(double)0.001);
		//Pc_t0 = ptree.get("initial.Pc",(double)1.e6);
		YH2O_t0 = ptree.get("initial.YH2O",(double)0.001);
		XCH4_t0 = ptree.get("initial.XCH4",(double)1.e6);

		Pw_x0 = ptree.get("boundary.Pw_at_left",(double)1.e6);
		Sgin_x0 = ptree.get("boundary.Sg_at_inlet",(double)0.001);
		
		numMaterials = ptree.get("sediment.number_of_materials",(int)1);

		numProps = 7;
		prop = std::vector<std::vector<double> > (numMaterials,std::vector<double>(numProps, 0.));
		for(int n_mat=0; n_mat<numMaterials; n_mat++ ){
			std::string name = "sediment.material"+std::to_string(n_mat);
			prop[n_mat][0] = ptree.get(name+".por",	(double)0.5);
			prop[n_mat][1] = ptree.get(name+".K",	(double)1.e-12);
			prop[n_mat][2] = ptree.get(name+".pentry",(double)1.e4);
			prop[n_mat][3] = ptree.get(name+".lambda",(double)1.2);
			prop[n_mat][4] = ptree.get(name+".swr",	(double)0.);
			prop[n_mat][5] = ptree.get(name+".sgr",	(double)0.);
			prop[n_mat][6] = ptree.get(name+".beta",(double)1.);
		}

		//reference state
		ref_salinity = ptree.get("reference_state.salinity",(double)0.);
		ref_saltconcentration = ref_salinity * (18.0/58.4); /*MolarMass_H2O/MolarMass_salt*/
		ref_temperature = 273.15+ptree.get("reference_state.temperature",(double)0.);

		//gravity
		gravity_flag = ptree.get("gravity.flag",(bool)true);
		g_magnitude = ptree.get("gravity.magnitude",(double)9.81);
  }

	/**********************************************************************
	 * INPUTS
	 **********
	 * z_domain : height of the computational domain [m]
	 * z_cells	: no. of cells along Z-axis
	 * 
	 *
	 *
	 *
	 **********************************************************************/

	//2. Initial Values

	double InitialSg(Dune::FieldVector< double,dim > xglobal) const {
		double Sg = Sg_t0;
		return Sg;
	}

	double InitialPw(Dune::FieldVector< double,dim > xglobal) const {
		double Pw = Pw_t0;
		return Pw;
	}

	double InitialPg(Dune::FieldVector< double,dim > xglobal) const {
		double Pg = Pg_t0;
		return Pg;
	}

	// double InitialSw(Dune::FieldVector< double,dim > xglobal) const {
	// 	double Sw = Sw_t0;
	// 	return Sw;
	// }
	double InitialSh(Dune::FieldVector< double,dim > xglobal) const {
		double Sh = Sh_t0;
		return Sh;
	}

	double InitialXCH4(Dune::FieldVector< double,dim > xglobal) const {
		double XCH4 = XCH4_t0;
		return XCH4;
	}
	double InitialYH2O(Dune::FieldVector< double,dim > xglobal) const {
		double YH2O = YH2O_t0;
		return YH2O;
	}

	//3. Boundary values

	double InletSg(Dune::FieldVector< double,dim > xglobal) const {
		double Sg = Sgin_x0;
		return Sg;
	}

	double LeftPw(Dune::FieldVector< double,dim > xglobal) const {
		double Pw = Pw_x0;
		return Pw;
	}

	//4. Material properties

	std::vector< std::vector<double> > layer_properties() const {
		return prop;
	}

	/**********************************************************************/
	/* REFERENCE STATE */
	double ReferenceSalinity() const {
		return ref_salinity; /*kg/kg*/
	}
	double ReferenceSaltConcentration() const {
		return ref_saltconcentration;
	}
	double ReferenceTemperature() const {
		return ref_temperature; /*K*/
	}

    /**********************************************************************/
	Dune::FieldVector<double,dim>
	SedimentVelocity ( double time, double dt ) const {

		Dune::FieldVector<double,dim> vs( 0. );
		vs[dim-1] = 0.;
		vs[0] = 0.;

		return vs; /*m/s*/
	}

	/* GRAVITY VECTOR */
	Dune::FieldVector<double,dim>
	g( ) const {
		Dune::FieldVector<double,dim> gravity( 0. );
		double g = 0.;
		if(gravity_flag) g = g_magnitude;
		gravity[dim-1] = g;
		gravity[0] = 0.;
		return gravity; /*N/kg*/
	}

	/**********************************************************************/

};


class TimeStepControl{
private:
public:

	constexpr static double maxdt 				= 1.	;
	constexpr static double mindt 				= 0.001		;
	const static bool adaptiveTimeStepControl	= true 	;

};

class ModelSubSystems{
public:
	const static bool FlowAndSoil 	= false	;
	const static bool OnlyFlow		= true	;
};

class TimeSteppingStrategies{
public:
	const static bool SequentialWithFPIterations			= false	;
	const static bool ImplicitMultiRate_interpolation		= false	;
};

class IncludeModels{
public:
	  const static bool realGasEoS					= false	;
	  const static bool gravity						= true	;
};


// class GravityVector{
// private:
// 	IncludeProblemSpecifications::ProblemSpecifications problemSpecs;
// 	constexpr static double pi = std::atan(1.)*4;
// public:
// 	  std::vector<double> g( ) const {
// 		  std::vector<double> gravity( problemSpecs.dimension, 0. );

// 		  if( IncludeModels::gravity ){
// 			  double g = 9.81;
// 			  gravity[ problemSpecs.dimension - 1 ] = -g;
// 			  gravity[0] = 0.;
// 		  }

// 		  return gravity;
// 	  }
// };

template<typename GV, typename RF>
struct ParameterTraits
{
  //! \brief the grid view
  typedef GV GridViewType;

  //! \brief Enum for domain dimension
  enum {
    //! \brief dimension of the domain
    dimDomain = GV::dimension
  };

  //! \brief Export type for domain field
  typedef typename GV::Grid::ctype DomainFieldType;

  //! \brief domain type
  typedef Dune::FieldVector<DomainFieldType,dimDomain> DomainType;

  //! \brief domain type
  typedef Dune::FieldVector<DomainFieldType,dimDomain-1> IntersectionDomainType;

  //! \brief Export type for range field
  typedef RF RangeFieldType;

  //! \brief range type
  typedef Dune::FieldVector<RF,GV::dimensionworld> RangeType;

  //! \brief permeability tensor type
  typedef Dune::FieldMatrix<RangeFieldType,dimDomain,dimDomain> PermTensorType;

  //! grid types
  typedef typename GV::Traits::template Codim<0>::Entity ElementType;
  typedef typename GV::Intersection IntersectionType;
};



#endif /* PARAMETERS_PARAMETERS_HH_ */
