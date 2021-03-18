class CharacteristicValues{

public:

	constexpr static double x_c = 1.e0; 
	constexpr static double t_c = 1.;// * 3.6 * 1.e3; // t_end Characteristic value for time
	constexpr static double density_c = 1.e0;
	constexpr static double viscosity_c = 1.e0;
	constexpr static double thermalconductivity_c = 1.e0;// 
	constexpr static double T_c = 1.e0;//t_c/(density_c*specificheat_c);//1.e0; // T_ref Characteristic value for Temperature
	constexpr static double P_c = 1.e0;//viscosity_c * x_c * x_c / (t_c  * permeability_c);//5.e4; // Pentry Characteristic value for Pressure
	constexpr static double permeability_c = viscosity_c * x_c * x_c / (t_c  * P_c);
	constexpr static double volumetricheat_c =  thermalconductivity_c * t_c/ density_c / x_c / x_c  ;;
	constexpr static double specificheat_c = volumetricheat_c ;//* x_c * x_c * viscosity_c /(permeability_c * P_c * t_c);
	constexpr static double dispersivity_c = x_c * x_c / t_c;

	constexpr static double X_source_mass 		= density_c / t_c;//t_c/density_c;// Correct
	constexpr static double X_convective_mass 	= 1.;//( permeability_c/(x_c*x_c) ) * ( P_c*t_c/viscosity_c );// Correct
	constexpr static double X_diffusive_mass 	= 1.;//( dispersivity_c * t_c ) / ( x_c*x_c ) ;// Correct
	constexpr static double X_gravity 			= 1. * P_c/ x_c/density_c; // Correct
	constexpr static double X_solidvelocity 	= 1.;//viscosity_c*x_c/(P_c*permeability_c);
	constexpr static double X_source_heat 		= density_c * T_c * volumetricheat_c / t_c;// Correct
	constexpr static double X_convective_heat 	= 1.;//( permeability_c/(x_c*x_c) ) * ( P_c*t_c/viscosity_c ) ;// Correct
	constexpr static double X_diffusive_heat 	= 1.;//( thermalconductivity_c/(x_c*x_c) ) * ( t_c/(density_c * specificheat_c));// Correct
	constexpr static double X_penalty_heat 		= 1.;//t_c/(density_c*specificheat_c); // Correct
	constexpr static double X_penalty_pressure 	= 1;//P_c * X_source_mass; // Correct
};

