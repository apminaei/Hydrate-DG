/*
 * soil.hh
 *
 *  Created on: Sep 22, 2016 jkjk
 *      Author: shubhangi
 */


#ifndef PARAMETERS_SOIL_HH_
#define PARAMETERS_SOIL_HH_


class Soil{
private:
    Indices indices;
    IncludeProblemSpecifications::ProblemSpecifications problemSpecs;
    HydraulicProperties hydraulicProperties;

public:

/* Porosity specified -> (Vw + Vg)/Vt => ((1 - Vh)/Vp )*(Vp/Vt) => ( 1 - Sh ) * basePorosity
					  -> basePorosity is used in all balance equations
*/
	double basePorosity( Dune::FieldVector< double, IncludeProblemSpecifications::ProblemSpecifications::dimension > globalPos ) const
	{
		double porositySediment = problemSpecs.SedimentPorosity( globalPos );
		double ShInit = problemSpecs.ProblemICValues( globalPos)[Indices::PVId_Sh] ;

		double porosityBase = porositySediment/( 1.0 - ShInit ) ;

		return porosityBase ;
	}

	double baseIntrinsicPermeability( Dune::FieldVector< double, IncludeProblemSpecifications::ProblemSpecifications::dimension > globalPos ) const
	{
		double kSediment = problemSpecs.SedimentPermeability( globalPos );
		double ShInit = problemSpecs.ProblemICValues( globalPos )[Indices::PVId_Sh] ;

		double kBase ;
		kBase = kSediment
			/ (   hydraulicProperties.KSF1( ShInit )
			    * 1./*hydraulicProperties.permeabilitySF2( porosityInit, porosityInit )*/
			    );

		return kBase ;
	}

	constexpr static double grainDensity = 2100.0 ;

	double density( ) const {
		/* unit -> kg/m^3 */
		double rho_s = grainDensity ;

		return rho_s;

	}

	double thermalConductivity( double T, double Ppore ) const
	{
		/* unit -> W/mK */
		return 1.9;

	}

	double Cp( double T, double Ppore ) const
	{
		/* unit -> J/kg.K */
		return 800.0;

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

	double grainRadius( Dune::FieldVector< double, IncludeProblemSpecifications::ProblemSpecifications::dimension > globalPos ) const
	{
		return problemSpecs.SoilGrainRadius( globalPos ) ;
	}

};


#endif /* PARAMETERS_SOIL_HH_ */
