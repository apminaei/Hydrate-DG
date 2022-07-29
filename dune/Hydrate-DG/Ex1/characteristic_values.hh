class CharacteristicValues
{

public:
	constexpr static double x_c = 1.e2;
	constexpr static double t_c = 1.e3 * 12. * 3. * 24. * 36.; // Characteristic value for time
	constexpr static double thermalconductivity_c = 1.e0;	   //
	constexpr static double T_c = 1.e2;						   //  Characteristic value for Temperature
	constexpr static double P_c = 1.e6;						   // Characteristic value for Pressure
	constexpr static double viscosity_c = 1.e-3;
	constexpr static double X_gravity = 1.; // Correct
	constexpr static double density_c = 1.e3;
	constexpr static double permeability_c = 1.e-13;
	constexpr static double volumetricheat_c = 1.e3;
	constexpr static double specificheat_c = volumetricheat_c;
	constexpr static double dispersivity_c = 1.e-9;

	constexpr static double X_source_mass = t_c / density_c;																			   // Correct
	constexpr static double X_convective_mass = 1. * (permeability_c / (x_c)) * (t_c / viscosity_c);									   // Correct
	constexpr static double X_diffusive_mass = 1. * (dispersivity_c * t_c) / (x_c * x_c);												   // Correct
	constexpr static double X_source_heat = t_c / (density_c * T_c * volumetricheat_c);													   // Correct
	constexpr static double X_convective_heat = 1. * (permeability_c / (x_c)) * (t_c / viscosity_c) * (specificheat_c / volumetricheat_c); // Correct
	constexpr static double X_diffusive_heat = 1. * (thermalconductivity_c / (x_c * x_c)) * (t_c / (density_c * volumetricheat_c));		   // Correct
};
