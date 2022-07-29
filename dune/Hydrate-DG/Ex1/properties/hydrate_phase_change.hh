template <typename GV, typename PTree>
class HydratePhaseChangeKinetics
{
private:
	const GV &gv;
	const PTree &ptree;

	const static int dim = GV::dimension;

	CharacteristicValues characteristicValue;

	Parameters<PTree> parameter;
	Methane<PTree> methane;
	Water<PTree> water;
	Hydrate<GV, PTree> hydrate;

public:
	//! construct from grid view
	HydratePhaseChangeKinetics(const GV &gv_,
							   const PTree &ptree_)
		: gv(gv_),
		  ptree(ptree_),
		  parameter(ptree_),
		  water(ptree_),
		  methane(ptree_),
		  hydrate(gv_, ptree_)
	{
	}

	// EQULIBRIUM CONDITIONS FOR HYDRATE:
	//------------------------------------
	// equilibrium pressure
	double EquilibriumPressure(double T /*K*/, double S) const
	{

		double A = 38.592, B = 8533.8, C = 16.32;
		double P_eq = 1.e3 * exp(A - B / (T) + C * S); // defined in Pascals
	}

	double DissociationRateConstant(double T) const
	{

		double kd_diss = parameter.HydrateDissociationRateConstant(); // defined in mol/m².Pa.s
		return kd_diss;
	}

	double FormationRateConstant_ingas(double T) const
	{

		double kd_form = parameter.HydrateFormationRateConstant(); // defined in mol/m².Pa.s
		return kd_form;
	}

	double FormationLimitFactor(double Sh, double Sw, double porosity) const
	{

		double term = Sw * (1. - Sw - Sh);
		return term;
	}

	double DissociationLimitFactor(double Sh, double Sw, double porosity) const
	{

		double term = Sh;
		return term;
	}

	// specific reaction area of hydrate in the sediment:
	double SpecificSurfaceArea(double Sh, double porosity, double permeability) const
	{

		double A_s;

		double M_base = 1.e5; // pow( 2.*permeability, -1/2. );
		double M_SF = pow(porosity * (1. - Sh), 3. / 2.);
		A_s = M_base * M_SF;

		return A_s;
	}

	// rate of gas generation:
	double GasGenerationRate(double T, double Pg, double Sh, double Sw, double XCH4,
							 double zCH4, double S, double porosity, double permeability) const
	{

		double gas_gen = 0.0;
		double A = 38.592, B = 8533.8, C = 16.32;
		double P_eq = 1.e3 * exp(A - B / (T) + C * S); // defined in Pascals
		double Peq = P_eq;

		double potential_P = (Peq / Pg) - 1.;

		if (potential_P > 0.)
		{
			gas_gen = DissociationRateConstant(T) * methane.MolarMass() * SpecificSurfaceArea(Sh, porosity, permeability) * DissociationLimitFactor(Sh, Sw, porosity) * Pg * potential_P;
		}
		else if (potential_P < 0.)
		{
			gas_gen = FormationRateConstant_ingas(T) * methane.MolarMass() * SpecificSurfaceArea(Sh, porosity, permeability) * FormationLimitFactor(Sh, Sw, porosity) * Pg * potential_P;
		}
		return gas_gen;
	}

	// rate of water generation:
	double WaterGenerationRate(double gasGenRate) const
	{
		double water_gen = gasGenRate * hydrate.HydrationNumber() * (water.MolarMass() / methane.MolarMass());
		return water_gen; /*[kg/m³s]*/
	}

	// rate of hydrate dissociation:
	double HydrateDissociationRate(double gasGenRate) const
	{
		double hyd_decomp = -gasGenRate * (hydrate.MolarMass() / methane.MolarMass());
		return hyd_decomp; /*[kg/m³s]*/
	}

	// heat of hydrate dissociation:
	double HeatOfDissociation(double gasGenRate, double T) const
	{
		double Q_decomp /*[W/m³]*/ = -(gasGenRate / methane.MolarMass()) * (56599.0 - 16.744 * (T));

		return Q_decomp;
	}
	//! get a reference to the grid view
	inline const GV &getGridView() { return gv; }
};
