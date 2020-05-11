/*
 * parameters.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef PARAMETERS_PARAMETERS_HH_
#define PARAMETERS_PARAMETERS_HH_

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


class GravityVector{
private:
	IncludeProblemSpecifications::ProblemSpecifications problemSpecs;
	constexpr static double pi = std::atan(1.)*4;
public:
	  std::vector<double> g( ) const {
		  std::vector<double> gravity( problemSpecs.dimension, 0. );

		  if( IncludeModels::gravity ){
			  double g = 9.81;
			  gravity[ problemSpecs.dimension - 1 ] = -g;
			  gravity[0] = 0.;
		  }

		  return gravity;
	  }
};

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
