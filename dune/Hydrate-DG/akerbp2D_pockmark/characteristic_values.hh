class CharacteristicValues{

public:

	constexpr static double x_c = 1.; 
	constexpr static double t_c = 1.;//7.2e4; // t_end Characteristic value for time
	constexpr static double permeability_c = 1.e0;
	constexpr static double density_c = 1.;
	constexpr static double viscosity_c = 1.e0;
	constexpr static double specificheat_c = 1.;
	constexpr static double thermalconductivity_c = 1.;
	constexpr static double dispersivity_c = 1.e0;
	constexpr static double P_c = viscosity_c / (t_c * x_c * x_c * permeability_c);//5.e4; // Pentry Characteristic value for Pressure
	constexpr static double T_c = t_c / ( density_c * specificheat_c);//273.15e0; // T_ref Characteristic value for Temperature

	constexpr static double X_source_mass 		= t_c/density_c;// Correct
	constexpr static double X_convective_mass 	= ( permeability_c/(x_c*x_c) ) * ( P_c*t_c/viscosity_c );// Correct
	constexpr static double X_diffusive_mass 	= ( dispersivity_c * t_c ) / ( x_c*x_c ) ;// Correct
	constexpr static double X_gravity 			= density_c*x_c/P_c;
	constexpr static double X_solidvelocity 	= viscosity_c*x_c/(P_c*permeability_c);
	constexpr static double X_source_heat 		= t_c/(density_c*specificheat_c*T_c);// Correct
	constexpr static double X_convective_heat 	= ( permeability_c/(x_c*x_c) ) * ( P_c*t_c/viscosity_c );// Correct
	constexpr static double X_diffusive_heat 	= ( thermalconductivity_c/(x_c*x_c) ) * ( t_c/(density_c*specificheat_c));// Correct

};
