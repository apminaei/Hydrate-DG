class CharacteristicValues{

public:

	constexpr static double x_c = 1.; 
	constexpr static double t_c = 3.6e3;//7.2e4; // t_end Characteristic value for time
	constexpr static double permeability_c = 1.e-13;
	constexpr static double density_c = 1.e3;
	constexpr static double viscosity_c = 1.e-3;
	constexpr static double specificheat_c = 1.e3;
	constexpr static double thermalconductivity_c = 1.;
	constexpr static double dispersivity_c = 1.e-9;//x_c * x_c / t_c;
	constexpr static double P_c = 1.e6;//viscosity_c * x_c * x_c / (t_c  * permeability_c);//5.e4; // Pentry Characteristic value for Pressure
	constexpr static double T_c = 1.e2;//t_c/(density_c*specificheat_c);//1.e0; // T_ref Characteristic value for Temperature
	constexpr static double volumetricheat_c =  specificheat_c  ;

	constexpr static double X_source_mass 		= t_c/density_c;// Correct
	constexpr static double X_convective_mass 	= ( permeability_c/(x_c*x_c) ) * ( P_c*t_c/viscosity_c );// Correct
	constexpr static double X_diffusive_mass 	= ( dispersivity_c * t_c ) / ( x_c*x_c ) ;// Correct
	constexpr static double X_gravity 			= 1. * P_c/(density_c*x_c); // Correct
	constexpr static double X_solidvelocity 	= viscosity_c*x_c/(P_c*permeability_c);
	constexpr static double X_source_heat 		= t_c/(density_c*T_c*specificheat_c);// Correct
	constexpr static double X_convective_heat 	= ( permeability_c/(x_c*x_c) ) * ( P_c*t_c/viscosity_c ) ;// Correct
	constexpr static double X_diffusive_heat 	= ( thermalconductivity_c/(x_c*x_c) ) * ( t_c/(density_c * specificheat_c));// Correct

};

