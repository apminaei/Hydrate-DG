class CharacteristicValues
{

public:
	constexpr static double x_c = 1.e0; // Characteristic value for x
	constexpr static double t_c = 1.e0; // Characteristic value for time
	constexpr static double density_c = 1.e0;
	constexpr static double viscosity_c = 1.e0;
	constexpr static double thermalconductivity_c = 1.e0; //
	constexpr static double T_c = 1.e2;					  // Characteristic value for Temperature
	constexpr static double P_c = 1.e6;					  // Characteristic value for Pressure
	constexpr static double permeability_c = 1.e0;
	constexpr static double volumetricheat_c = 1.e0;
	constexpr static double specificheat_c = volumetricheat_c;
	constexpr static double dispersivity_c = 1.e0;

	constexpr static double X_source_mass = t_c / density_c;									// Correct
	constexpr static double X_convective_mass = (permeability_c / (x_c)) * (t_c / viscosity_c); // Correct
	constexpr static double X_diffusive_mass = (dispersivity_c * t_c) / (x_c * x_c);			// Correct
	constexpr static double X_gravity = 1.;
	constexpr static double X_source_heat = t_c / (density_c * T_c * volumetricheat_c);												  // Correct
	constexpr static double X_convective_heat = (permeability_c / (x_c)) * (t_c / viscosity_c) * (specificheat_c / volumetricheat_c); // Correct
	constexpr static double X_diffusive_heat = (thermalconductivity_c / (x_c * x_c)) * (t_c / (density_c * volumetricheat_c));		  // Correct
};
