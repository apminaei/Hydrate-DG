/*
 * reactionKinetics.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef PARAMETERS_REACTIONKINETICS_HH_
#define PARAMETERS_REACTIONKINETICS_HH_
template<typename GV, typename PTree>
class ReactionKinetics
{
private:
	const GV& gv;
	const PTree& ptree;
	Hydrate hydrate;
	Methane methane;
	Water water;
	HydraulicProperties<GV, PTree> hydraulicProperties;
	Soil<GV, PTree> soil;
	
	Indices index;
	const static int dim = GV::dimension;
public:
	// EQULIBRIUM CONDITIONS FOR HYDRATE:
	//------------------------------------
	 //! construct from grid view
  	ReactionKinetics (const GV& gv_, const PTree& ptree_)
  	: gv( gv_ ), ptree(ptree_), soil(gv_, ptree_),
		  hydraulicProperties(gv_, ptree_)
  	{}


	//equilibrium pressure
	double eqbpressure(double T) const
	{
		// KAMATH and HOLDER correlation:
		double T_l = 273.15 + 0.0;
		double T_r = 273.15 + 0.05;
		double P_eq = 0.;

		int Case_r = 1;
		int Case_l = 1.;

		double C_r = 1.e3, A_r = 38.98, B_r = 8533.8;
		double C_l = 1.e3, A_l = 14.717, B_l = 1886.79;

		if (Case_r == 1)
		{
			A_r = 38.98;
			B_r = 8533.8;
			C_r = 1.e3; //1.28*1.e3; // used it for kiel triax example 1.28,1.381; 1.5635->4.5 MPa //
		}
		else if (Case_r == 2)
		{
			A_r = 49.3185;
			B_r = 9459.0;
			C_r = 1.21;
		}
		else if (Case_r == 3)
		{
			A_r = 49.3185;
			B_r = 9459.0;
			C_r = 1.81;
		}

		if (Case_l == 1)
		{
			A_l = 14.717;
			B_l = 1886.79;
			C_l = 1.e3;
		}

		if (T > (T_r))
		{
			P_eq = C_r * exp(A_r - B_r / (T)); // defined in Pascals
		}
		else if (T < (T_l))
		{
			P_eq = C_l * exp(A_l - B_l / (T)); // defined in Pascals
		}
		else
		{

			double Peq_l = C_l * exp(A_l - B_l / (T));
			double Peq_r = C_r * exp(A_r - B_r / (T));
			double slope = (Peq_r - Peq_l) / (T_r - T_l);
			P_eq = Peq_l + (T - T_l) * slope;
		}

		//		P_eq = C_l * exp( A_l - B_l/( T ) ); // defined in Pascals
		return P_eq;
	}

	// KIM BISHNOI KINETICS MODEL FOR HYDRATE DISSOCIATION:

	//rate constant for hydrate dissociation;
	double diss_rateconstant(double T) const
	{
		// Clarke and Bishnoi 2001 gave 'k0' and 'delta E_activation'
		// Kim et al 1987 gave the Arhenious type eqn for kinetic rate constant 'k'

		double k0 = soil.ReactionRateParam_k0();
		double deltaEbyR = soil.ReactionRateParam_deltaEabyR();
		double kd_diss = k0 * exp(-deltaEbyR / (T)); //defined in mol/m².Pa.s

		//		double kd_diss = soil.problemSpecs.ReactionRateParam_k0(); //used in kiel triax test

		return kd_diss;
	}

	//rate constant for hydrate reformation;
	double form_rateconstant(double T) const
	{
		//Englezos et al. (1987a)			reference article: Kinetic simulation of methane hydrate formation and dissociation in porous media
		//  ------------------   author: K Mohanty
		//		double kd_form = 0.5875 * 1.e-11;		//defined in mol/m².Pa.s

		double kd_form = soil.ReactionRateParam_kForm(); //defined in mol/m².Pa.s
		return kd_form;
	}

	double Tau_formFactor(double Sh, double Sw, double porosity) const
	{
		double Sg = 1.0 - Sw - Sh;

		// used for Kiel-triax
		//		double term = 0.;
		//		if( Sh <= 0.02){
		//			term = Sg * Sw ;
		//		}
		//		else if( Sh > 0.02 and Sh <= 0.1 ){
		//			double y_1 = (1.-0.37-0.02)*0.37;
		//			double y_2 = pow(0.1, double(1./3.)) * pow(0.32 , double(5./3.) );;//pow(0.1, double(1./4.)) * pow(0.32 , double(7./4.) );//pow( 0.32*0.1*(1.-0.32-0.1) , double(2./3.) );//
		//			double x_1 = 0.02;
		//			double x_2 = 0.1;
		//			term = (y_2 - y_1)/(x_2 - x_1) * (Sh - x_1) + y_1 ;
		//		}
		//		else{
		//			term = pow(Sh, double(1./3.)) * pow(Sw , double(5./3.) );//pow( Sw*Sh*Sg , double(2./3.) );//
		//		}

		//		double term = pow(Sh, double(1./2.)) * pow(Sw , double(3./2.) ) ;
		//		double term = pow(Sh, double(2./5.)) * pow(Sw*Sg , double(4./5.) ) ;
		double term = Sg * Sw;
		//		double term = pow( Sw*Sh*Sg , double(2./3.) );
		return term;

		//		if( Sg>0.05){
		//			double term1 = Sw*Sh;
		//			term1 = std::pow( term1, double(2./3.) );
		//			double term2 = Sg;
		//			term2 = std::pow( term2, double(2./3.) );
		//			term = term1*term2;
		//		}
		//		else{
		//			term = Sg * sqrt( Sw * Sh );
		//		}
		//		return term ;
	}

	double Tau_dissFactor(double Sh, double Sw, double porosity) const
	{
		double Sg = 1. - Sw - Sh;
		if (Sh < 1.e-6)
		{
			Sh = 0.;
		}
		if (Sw < 1.e-6)
		{
			Sw = 0.;
		}
		if (Sg < 1.e-6)
		{
			Sg = 0.;
		}

		//		double term = pow( Sw*Sh*Sg , 2./3. );
		//		double term = 0.35 * Sh; //used it for kiel triax example
		//		double term = Sh*Sh;
		double term = porosity * Sh;
		return term;
	}

	//specific reaction area of hydrate in the sediment:
	double sp_surfacearea(double Sh, double porosity, Dune::FieldVector<double, dim> globalPos) const
	{
		if (Sh < 1.e-6)
		{
			Sh = 0.;
		}
		double A_s;
		// Yousif et al 1991:

		double eff_porosity = porosity * (1. - Sh);
		double permeability = soil.SedimentPermeability(globalPos) * hydraulicProperties.KSF1(Sh) * hydraulicProperties.KSF2(porosity, soil.SedimentPorosity(globalPos));
		A_s = sqrt(pow(eff_porosity, 3.) / (2. * permeability)); //defined in (m²/m³)

		// 		double eff_porosity = porosity * (1. - Sh);
		// 		double pi = 22./7. ;
		// 		A_s = ( ( 1 - porosity ) / ( ( 4./3.) * pi * pow(soil.grainRadius( globalPos ),3.) ) ) * ( 4. * pi * pow(soil.grainRadius( globalPos ),2.)) * pow( Sh , (2./3.)) ;
		// 		A_s = ( (1 - porosity ) * 3.0 / soil.grainRadius( globalPos ) ) * pow( Sh , (2./3.));

		// 		double M_base = 2.0 * ( porosity )
		// 							 * ( hydraulicProperties.lambda() / ( hydraulicProperties.lambda() - 1. ) )
		// 							 * ( 1./soil.grainRadius( globalPos ) ) ;
		// 		double m = hydraulicProperties.m;
		// 		double M_SF = pow( (1.0 - Sh), ( m - 1 )/m );
		// 		A_s = M_base * M_SF ;

		//		A_s = ( 3.75 * 1.e5 ) * pow( (1.0 - Sh), ( m - 1 )/m );

		// 		std::cout<< "A_s = " << A_s << std::endl ;

		return A_s;
	}

	// rate of gas generation:
	double gasGenerationRate(double T, double Pg, double Sh, double Sw, double porosity,
							 Dune::FieldVector<double, dim> globalPos) const
	{
		double gas_gen = 0.0;

		//		// Use for Kiel triax
		//		double Salinity_0 = 35.;
		//		double Salinity = Salinity_0 / ( 1. + hydrate.hydrationNumber()
		//												* (water.molarMass()/hydrate.molarMass())
		//												* (hydrate.grainDensity/water.density(T,Pg))
		//												* ( (soil.problemSpecs.ProblemICValues(globalPos)[indices.PVId_Sh] - Sh)
		//														/soil.problemSpecs.ProblemICValues(globalPos)[indices.PVId_Sw]) );
		////		double Salinity = Salinity_0 * ( soil.problemSpecs.ProblemICValues(globalPos)[indices.PVId_Sw] / Sw ) ;
		//
		//		double Peq = (0.5/40.0) * Salinity + ( eqbpressure( T ) - 0.50 );

		double Peq = eqbpressure(T);

		double potential = Peq - Pg;

		/*		if( potential > 1.e-6 ){
			std::cout<<""<<std::endl;
			std::cout<< "gasPressure: " << Pg << std::endl ;
			std::cout<< "eqbpressure: " << eqbpressure( T ) << std::endl ;
			std::cout<< "potential: " << potential << std::endl ;
			std::cout<< "k0^d: " << diss_rateconstant( T ) << std::endl;
			std::cout<< "As: " << sp_surfacearea( Sh, porosity, globalPos ) << std::endl;
			std::cout<< "tau_diss: " << Tau_dissFactor( Sh, Sw ) <<std::endl;
		}*/

		double potential_A = 0. + 1.e-3;
		double potential_B = 0. - 1.e-3;

		if (potential > potential_A)
		{
			gas_gen = diss_rateconstant(T) * methane.molarMass() * sp_surfacearea(Sh, porosity, globalPos) * Tau_dissFactor(Sh, Sw, porosity) * soil.AreaFactor() * potential;
			/*			std::cout << "dissociation: " << 	diss_rateconstant( T )
										      * methane.molarMass()
										      * sp_surfacearea( Sh, porosity, globalPos )
										      * soil.problemSpecs.AreaFactor()
									      << std::endl;*/
		}
		else if (potential < potential_B)
		{
			gas_gen = form_rateconstant(T) * methane.molarMass() * sp_surfacearea(Sh, porosity, globalPos) * Tau_formFactor(Sh, Sw, porosity) * soil.AreaFactor() * potential;
			/*			std::cout << "formation: " << 	form_rateconstant( T )
										   * methane.molarMass()
										   * sp_surfacearea( Sh, porosity, globalPos )
										   * soil.problemSpecs.AreaFactor()
									   << std::endl;*/
		}
		else
		{
			double reactionRate_A = diss_rateconstant(T) * Tau_dissFactor(Sh, Sw, porosity);
			double reactionRate_B = form_rateconstant(T) * Tau_formFactor(Sh, Sw, porosity);
			double slope = (reactionRate_B - reactionRate_A) / (potential_B - potential_A);
			double reactionRate = reactionRate_A + slope * (potential - potential_A);
			gas_gen = reactionRate * methane.molarMass() * sp_surfacearea(Sh, porosity, globalPos) * soil.AreaFactor() * potential;
		}

		/* 		if (gas_gen < 0.) {
 			std::cout << "gas_gen: " << gas_gen << std::endl;
 			std::cout<<""<<std::endl;
 		}*/

		return gas_gen;
	}

	// rate of water generation:
	double waterGenerationRate(double T, double Pg, double Sh, double Sw, double porosity,
							   Dune::FieldVector<double, dim> globalPos) const
	{
		double water_gen = gasGenerationRate(T, Pg, Sh, Sw, porosity, globalPos) * hydrate.hydrationNumber() * (water.molarMass() / methane.molarMass());
		return water_gen; /*[kg/m³s]*/
	}

	// rate of hydrate dissociation:
	double hydrateDissociationRate(double T, double Pg, double Sh, double Sw, double porosity,
								   Dune::FieldVector<double, dim> globalPos) const
	{
		double hyd_decomp = -gasGenerationRate(T, Pg, Sh, Sw, porosity, globalPos) * (hydrate.molarMass() / methane.molarMass());
		return hyd_decomp; /*[kg/m³s]*/
	}

	// heat of hydrate dissociation:
	double heatOfDissociation(double T, double Pg, double Sh, double Sw, double porosity,
							  Dune::FieldVector<double, dim> globalPos) const
	{
		double Q_decomp /*[W/m³]*/ =
			-(gasGenerationRate(T, Pg, Sh, Sw, porosity, globalPos) / hydrate.molarMass()) * (56599.0 - 16.744 * (T));
		//        		0.;

		return Q_decomp;
	}
};

#endif /* PARAMETERS_REACTIONKINETICS_HH_ */
