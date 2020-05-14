/*
 * soil.hh
 *
 *  Created on: Sep 22, 2016 
 *      Author: shubhangi
 */


#ifndef PARAMETERS_SOIL_HH_
#define PARAMETERS_SOIL_HH_



/* Porosity specified -> (Vw + Vg)/Vt => ((1 - Vh)/Vp )*(Vp/Vt) => ( 1 - Sh ) * basePorosity
					  -> basePorosity is used in all balance equations
*/


/* 
	IMPORTANT NOTICE	=====>	Local and global coordinates should be checked in all parameters
	
		ALL PARAMETERS ARE DEFINED IN globalPos
*/

template<typename GV, typename PTree>
class Soil
{
private:
	const GV& gv;
	const PTree& ptree;
	Parameters<PTree> parameter;
	MeshParameters<PTree> mesh;

	HydraulicProperties<GV, PTree> hydraulicProperties;
	CharacteristicValues characteristicValue;
	const static int dim = GV::dimension;

public:

  //! construct from grid view
  Soil (const GV& gv_, const PTree& ptree_)
  : gv( gv_ ), ptree(ptree_), parameter(ptree_), mesh(ptree_),
		  hydraulicProperties(gv_, ptree_)
  {}


	double SedimentPorosity
	( const Dune::FieldVector<double,dim>& globalPos) const {

		//Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		auto prop_L = parameter.layer_properties();

		double por = 0.;

		por = prop_L[0][0];

		if( mesh.isLens(globalPos) ){
			por = prop_L[1][0];
		}

		return por;
	}

	double SedimentPermeability
	( const Dune::FieldVector<double,dim>& globalPos) const {

		//Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		auto prop_L = parameter.layer_properties();

		double K = 0.; /*m^2*/

		K = prop_L[0][1];

		if( mesh.isLens(globalPos) ){
			K = prop_L[1][1];
		}

		return K/characteristicValue.permeability_c; /*ndim*/
	}

	// vector coefficient
	Dune::FieldMatrix<double,dim,dim>
	SedimentPermeabilityTensor
	( const Dune::FieldVector<double,dim>& globalPos) const {

		double K_xx = SedimentPermeability(globalPos);
		double K_yy = K_xx;
		Dune::FieldMatrix<double,dim, dim> PermeabilityTensor;
		
		PermeabilityTensor[0][0] = K_xx ;
		PermeabilityTensor[0][1] = 0. ;
		PermeabilityTensor[1][0] = 0. ;
		PermeabilityTensor[1][1] = K_yy ;

		return PermeabilityTensor; /*ndim*/
	}

	double SoilGrainRadius( const Dune::FieldVector< double, dim > globalPos ) const {
	  // Bear et at, 1972
	  double rp = sqrt( 45.0 * SedimentPermeability( globalPos ) * pow( 1- SedimentPorosity(globalPos ) , 2.0 )/pow( SedimentPorosity( globalPos ),3.0) );
// 	  std::cout<< "rp = " << rp << std::endl;
	      return rp;
	}


	double AreaFactor() const {
	   return 1.;
	}

	double ReactionRateParam_k0() const {

		double k0 = 360000.;		//defined in mol/m².Pa.s

		return k0;
	}

	double ReactionRateParam_deltaEabyR() const {
		return 9752.73;
	}

	double ReactionRateParam_kForm() const {
		//Englezos et al. (1987a)			reference article: Kinetic simulation of methane hydrate formation and dissociation in porous media
											//  ------------------   author: K Mohanty
		double kd_form ;		//defined in mol/m².Pa.s
		kd_form = 0.;//0.5875 * 1.e-11;

		return kd_form;
	}



	double basePorosity( Dune::FieldVector< double, dim > globalPos ) const
	{
		double porositySediment = SedimentPorosity(globalPos );
		double ShInit = parameter.InitialSh(globalPos)[Indices::PVId_Sh] ;

		double porosityBase = porositySediment/( 1.0 - ShInit ) ;

		return porosityBase ;
	}

	double baseIntrinsicPermeability( Dune::FieldVector< double, dim > globalPos ) const
	{
		double kSediment = SedimentPermeability( globalPos );
		double ShInit = parameter.InitialSh( globalPos )[Indices::PVId_Sh] ;

		double kBase ;
		kBase = kSediment
			/ (   hydraulicProperties.KSF1( ShInit )
			    * 1./*hydraulicProperties.permeabilitySF2( porosityInit, porosityInit )*/
			    );

		return kBase ;
	}


	double density() const {
		/* unit -> kg/m^3 */
		double rho = 2600.0;
		return rho/characteristicValue.density_c; /*ndim*/
	}


	double thermalConductivity( double T, double Ppore ) const
	{
		/* unit -> W/mK */
		return 3.;

	}

	double Cp( double T, double Ppore ) const
	{
		/* unit -> J/kg.K */
		return 1000.0;

	}

	double Cv( double T, double Ppore ) const
	{
		/* unit -> W/kg.K */
		return Cp( T, Ppore );

	}

	double tortuosity( double porosity ) const
	{
		return porosity * porosity ;
	}

	double grainRadius( Dune::FieldVector< double, dim > globalPos ) const
	{
		return SoilGrainRadius(globalPos ) ;
	}
	

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

};



// ==============================================================================================

#endif /* PARAMETERS_SOIL_HH_ */
