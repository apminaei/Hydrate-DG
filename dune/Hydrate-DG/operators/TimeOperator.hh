using namespace Dune::PDELab;

template <class GV, typename Params>
class TimeOperator : public Dune::PDELab::NumericalJacobianApplyVolume<TimeOperator<GV, Params>>,
					 public Dune::PDELab::NumericalJacobianVolume<TimeOperator<GV, Params>>,
					 public Dune::PDELab::FullVolumePattern,
					 public Dune::PDELab::LocalOperatorDefaultFlags,
					 public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
	const GV &gv;
	const Params &property;

	constexpr static double eps = 1.0e-6;
	constexpr static double pi = atan(1.) * 4;
	unsigned int intorder;
	double Xc_K;
	double Xc_mu;
	double Xc_rho;
	double Xc_kth;
	double Xc_C;
	double Xc_P;
	double Xc_T;
	double Xc_X;
	double Xc_Y;
	double T_ref;

public:
	// pattern assembly flags
	enum
	{
		doPatternVolume = true
	};
	// residual assembly flags
	enum
	{
		doAlphaVolume = true
	};

	typedef typename GV::IndexSet IndexSet;

	// constructor remembers parameters
	TimeOperator(const GV &gv_, const Params &property_, unsigned int intorder_ = 4)
		: gv(gv_), property(property_), intorder(intorder_)
	{
		Xc_K = property.characteristicValue.permeability_c;
		Xc_mu = property.characteristicValue.viscosity_c;
		Xc_rho = property.characteristicValue.density_c;
		Xc_kth = property.characteristicValue.thermalconductivity_c;
		Xc_C = property.characteristicValue.specificheat_c;
		Xc_P = property.characteristicValue.P_c;
		Xc_T = property.characteristicValue.T_c;
		Xc_X = property.characteristicValue.x_c;
		Xc_Y = property.characteristicValue.x_c;
		T_ref = property.parameter.ReferenceTemperature() / Xc_T;
	}
	// volume integral depending on test and ansatz functions
	template <typename EG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume(const EG &eg, const LFSU &lfsu, const X &x, const LFSV &lfsv, R &r) const
	{

		// subspaces
		// Gas pressure
		const auto &lfsv_Pw = lfsv.template child<Indices::PVId_Pw>();
		const auto &lfsu_Pw = lfsu.template child<Indices::PVId_Pw>();

		// Water Saturation
		const auto &lfsv_Sg = lfsv.template child<Indices::PVId_Sg>();
		const auto &lfsu_Sg = lfsu.template child<Indices::PVId_Sg>();

		// Hydrate Saturation
		const auto &lfsv_Sh = lfsv.template child<Indices::PVId_Sh>();
		const auto &lfsu_Sh = lfsu.template child<Indices::PVId_Sh>();

		// Temperature
		const auto &lfsv_T = lfsv.template child<Indices::PVId_T>();
		const auto &lfsu_T = lfsu.template child<Indices::PVId_T>();

		// Methane mole fraction
		const auto &lfsv_XCH4 = lfsv.template child<Indices::PVId_XCH4>();
		const auto &lfsu_XCH4 = lfsu.template child<Indices::PVId_XCH4>();

		// Water mole fraction
		const auto &lfsv_YH2O = lfsv.template child<Indices::PVId_YH2O>();
		const auto &lfsu_YH2O = lfsu.template child<Indices::PVId_YH2O>();

		// Salt mole fraction
		const auto &lfsv_XC = lfsv.template child<Indices::PVId_C>();
		const auto &lfsu_XC = lfsu.template child<Indices::PVId_C>();

		// define types
		using RF = typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::FiniteElementType::
			Traits::LocalBasisType::Traits::RangeFieldType;
		using RangeType = typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::FiniteElementType::
			Traits::LocalBasisType::Traits::RangeType;
		using size_type = typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::SizeType;

		// Reference to cell
		const auto &cell = eg.entity();

		// Get geometry
		auto geo = eg.geometry();

		// Transformation matrix
		typename EG::Geometry::JacobianInverseTransposed jac;
		Dune::QuadratureType::Enum qt = Dune::QuadratureType::GaussLobatto;

		// loop over quadrature points
		for (const auto &ip : quadratureRule(geo, intorder, qt))
		{
			auto qp = ip.position();
			// evaluate test shape functions
			std::vector<RangeType> phi_Pw(lfsu_Pw.size());
			lfsu_Pw.finiteElement().localBasis().evaluateFunction(qp, phi_Pw);
			std::vector<RangeType> psi_Pw(lfsv_Pw.size());
			lfsv_Pw.finiteElement().localBasis().evaluateFunction(qp, psi_Pw);

			std::vector<RangeType> phi_Sg(lfsu_Sg.size());
			lfsu_Sg.finiteElement().localBasis().evaluateFunction(qp, phi_Sg);
			std::vector<RangeType> psi_Sg(lfsv_Sg.size());
			lfsv_Sg.finiteElement().localBasis().evaluateFunction(qp, psi_Sg);

			std::vector<RangeType> phi_Sh(lfsu_Sh.size());
			lfsu_Sh.finiteElement().localBasis().evaluateFunction(qp, phi_Sh);
			std::vector<RangeType> psi_Sh(lfsv_Sh.size());
			lfsv_Sh.finiteElement().localBasis().evaluateFunction(qp, psi_Sh);

			std::vector<RangeType> phi_T(lfsu_T.size());
			lfsu_T.finiteElement().localBasis().evaluateFunction(qp, phi_T);
			std::vector<RangeType> psi_T(lfsv_T.size());
			lfsv_T.finiteElement().localBasis().evaluateFunction(qp, psi_T);

			std::vector<RangeType> phi_XCH4(lfsu_XCH4.size());
			lfsu_XCH4.finiteElement().localBasis().evaluateFunction(qp, phi_XCH4);
			std::vector<RangeType> psi_XCH4(lfsv_XCH4.size());
			lfsv_XCH4.finiteElement().localBasis().evaluateFunction(qp, psi_XCH4);

			std::vector<RangeType> phi_YH2O(lfsu_YH2O.size());
			lfsu_YH2O.finiteElement().localBasis().evaluateFunction(qp, phi_YH2O);
			std::vector<RangeType> psi_YH2O(lfsv_YH2O.size());
			lfsv_YH2O.finiteElement().localBasis().evaluateFunction(qp, psi_YH2O);

			std::vector<RangeType> phi_XC(lfsu_XC.size());
			lfsu_XC.finiteElement().localBasis().evaluateFunction(qp, phi_XC);
			std::vector<RangeType> psi_XC(lfsv_XC.size());
			lfsv_XC.finiteElement().localBasis().evaluateFunction(qp, psi_XC);

			// evaluate Pw
			RF Pw = 0.0;
			for (size_type i = 0; i < lfsu_Pw.size(); i++)
				Pw += x(lfsu_Pw, i) * phi_Pw[i];

			// evaluate Sg
			RF Sg = 0.0;
			for (size_type i = 0; i < lfsu_Sg.size(); i++)
				Sg += x(lfsu_Sg, i) * phi_Sg[i];

			// evaluate Sh
			RF Sh = 0.0;
			for (size_type i = 0; i < lfsu_Sh.size(); i++)
				Sh += x(lfsu_Sh, i) * phi_Sh[i];

			// evaluate T
			RF T = 0.0;
			for (size_type i = 0; i < lfsu_T.size(); i++)
				T += x(lfsu_T, i) * phi_T[i];

			// evaluate XCH4
			RF XCH4 = 0.0;
			for (size_type i = 0; i < lfsu_XCH4.size(); i++)
				XCH4 += x(lfsu_XCH4, i) * phi_XCH4[i];

			// evaluate YH2O
			RF YH2O = 0.0;
			for (size_type i = 0; i < lfsu_YH2O.size(); i++)
				YH2O += x(lfsu_YH2O, i) * phi_YH2O[i];

			// evaluate XC
			RF XC = 0.0;
			for (size_type i = 0; i < lfsu_XC.size(); i++)
			{
				XC += x(lfsu_XC, i) * phi_XC[i];
			}

			RF Sw = (1. - Sg - Sh);

			// evaluate Pg
			auto por = property.soil.SedimentPorosity(cell, qp);
			auto Pc = property.hydraulicProperty.CapillaryPressure(cell, qp, Sw, Sh, por); /* ndim */

			RF Pg = Pw + Pc;
			RF Peff = (Pg * Sg + Pw * Sw) / (1. - Sh);

			auto Pw_dim = Pw * Xc_P;
			auto Pg_dim = Pg * Xc_P;
			auto T_dim = T * Xc_T;

			double S = XC * (property.salt.MolarMass() / property.gas.MolarMass());
			auto zCH4 = property.eos.EvaluateCompressibilityFactor(T_dim, Pg_dim);
			auto H_CH4_w = property.gas.SolubilityCoefficient(T_dim /*K*/, S);	   /*ndim */
			auto P_H_satu = property.water.SaturatedVaporPressure(T_dim /*K*/, S); /*ndim */
			auto rho_g = property.gas.Density(T_dim, Pg_dim, zCH4);
			auto rho_w = property.water.Density(T_dim, Pw_dim, S);
			auto rho_h = property.hydrate.Density();
			auto rho_s = property.soil.Density();
			auto Cv_g = property.gas.Cv(T_dim, Pg_dim, zCH4);
			auto Cv_w = property.water.Cv(T_dim, Pw_dim, S);
			auto Cv_h = property.hydrate.Cv(T_dim, Peff * Xc_P);
			auto Cv_s = property.soil.Cv();
			auto Cv_eff = (1. - por) * rho_s * Cv_s + por * (rho_g * Sg * Cv_g + rho_w * Sw * Cv_w + rho_h * Sh * Cv_h);

			auto H2O_g = rho_g * por * YH2O * Sg;
			auto CH4_w = rho_w * por * XCH4 * Sw;
			auto SALT_w = rho_w * por * XC * Sw;
			auto H2O_w = rho_w * por * (1. - XC - XCH4) * Sw;
			auto CH4_g = rho_g * por * (1. - YH2O) * Sg;

			// integrate (A grad u - bu)*grad phi_i + a*u*phi_i
			RF factor = ip.weight() * geo.integrationElement(qp);
			auto mass = por * (rho_g * Sg + rho_w * Sw);
			for (size_type i = 0; i < lfsv_Sg.size(); i++)
			{
				r.accumulate(lfsv_Sg, i, ((CH4_g + CH4_w) * psi_Sg[i]) * factor); //
			}
			for (size_type i = 0; i < lfsv_XC.size(); i++)
			{
				r.accumulate(lfsv_XC, i, (SALT_w * psi_XC[i]) * factor);
			}
			for (size_type i = 0; i < lfsv_Pw.size(); i++)
			{
				r.accumulate(lfsv_Pw, i, ((H2O_g + H2O_w) * psi_Pw[i]) * factor);
			}
			for (size_type i = 0; i < lfsv_Sh.size(); i++)
			{
				r.accumulate(lfsv_Sh, i, (rho_h * por * Sh * psi_Sh[i]) * factor);
			}
			for (size_type i = 0; i < lfsv_T.size(); i++)
			{
				r.accumulate(lfsv_T, i, Cv_eff * (T - T_ref) * psi_T[i] * factor);
			}

		} // End Quadrature Rule
	}	  // End of alpha volume
};