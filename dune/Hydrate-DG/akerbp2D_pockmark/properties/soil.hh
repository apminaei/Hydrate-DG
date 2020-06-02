/* ALL PARAMETERS ARE NONDIMENSIONAL */
template<typename GV, typename Parameters>
class Soil
{
private:
	  const GV& gv;
	  
	  
	const Parameters& parameter;
	  CharacteristicValues characteristicValue;
	  const static int dim = GV::dimension;

public:
  //! construct from grid view
  Soil ( const GV& gv_ , const Parameters& parameter_ )
  : gv( gv_ ), parameter(parameter_)
  {}


	double SedimentPorosity
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		auto prop_L = parameter.layer_properties();

		double por = 0.;

		por = prop_L[0][0];

		// if( parameter.mesh.isLens(x) ){
		// 	por = prop_L[1][0];
		// }
		

		return por;
	}

	double SedimentPermeability
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);
		
		
		auto prop_L = parameter.layer_properties();

		double K = 0.; /*m^2*/

		K = prop_L[0][1];

		// if( parameter.mesh.isLens(x) ){
		// 	K = prop_L[1][1];
		// }
		

		return K/characteristicValue.permeability_c; /*ndim*/
	}

	// vector coefficient
	Dune::FieldMatrix<double,dim,dim>
	SedimentPermeabilityTensor
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal) const {

		double K_xx = SedimentPermeability(element,xlocal);
		double K_yy = K_xx;
		Dune::FieldMatrix<double,dim, dim> PermeabilityTensor;
		
		PermeabilityTensor[0][0] = K_xx ;
		PermeabilityTensor[0][1] = 0. ;
		PermeabilityTensor[1][0] = 0. ;
		PermeabilityTensor[1][1] = K_yy ;
		
		return PermeabilityTensor; /*ndim*/
	}

	double SoilGrainRadius
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal) const {

		auto x = element.geometry().global(xlocal);

		// Bear et at, 1972
		double por  = SedimentPorosity( element,xlocal );
		double perm = SedimentPermeability( element,xlocal )*characteristicValue.permeability_c;
		double rp = sqrt( 45.0 * perm * pow( 1- por , 2.0 )/pow( por,3.0) );
		return rp/characteristicValue.x_c; /*ndim*/
	}

	double Density() const {
		/* unit -> kg/m^3 */
		double rho = 2600.0;
		return rho/characteristicValue.density_c; /*ndim*/
	}

	double ThermalConductivity() const {
		/* unit -> W/mK */
		double kth = 3.0;
		return kth/characteristicValue.thermalconductivity_c; /*ndim*/
	}

	double Cp() const {
		/* unit -> J/kg.K */
		double Cp = 1000.0;
		return Cp/characteristicValue.specificheat_c; /*ndim*/
	}

	double Cv() const {
		/* unit -> W/kg.K */
		double Cv = Cp(); /*ndim*/
		return Cv;/*ndim*/
	}

	double Tortuosity( double porosity ) const {
		return porosity * porosity ;/*ndim*/
	}

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

};
