/*
 * characteristicValues.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef PARAMETERS_CHARACTERISTICVALUES_HH_
#define PARAMETERS_CHARACTERISTICVALUES_HH_

class CharacteristicValues{

private:

public:

	constexpr static double x_c = 1.;
	constexpr static double t_c = 1.;
	constexpr static double P_c = 1.;//1.e6;
	constexpr static double T_c = 1.;//273.15;
	constexpr static double permeability_c = 1.;//1.e-12;
	constexpr static double density_c = 1.;//1000.;
	constexpr static double viscosity_c = 1.;//1.e-5;
	constexpr static double specificheat_c = 1.;//4200.;
	constexpr static double thermalconductivity_c = 1.;

	constexpr static double X_source_mass 		= t_c/density_c;
	constexpr static double X_convective_mass 	= ( permeability_c/(x_c*x_c) ) * ( P_c*t_c/viscosity_c );
	constexpr static double X_diffusive_mass 	= 1.;//
	constexpr static double X_gravity 			= density_c*x_c/P_c;
	constexpr static double X_source_heat 		= t_c/(density_c*specificheat_c*T_c);
	constexpr static double X_convective_heat 	= ( permeability_c/(x_c*x_c) ) * ( P_c*t_c/viscosity_c );
	constexpr static double X_diffusive_heat 	= ( thermalconductivity_c/(x_c*x_c) ) * ( t_c/(density_c*specificheat_c));

	

};



#endif /* PARAMETERS_CHARACTERISTICVALUES_HH_ */
