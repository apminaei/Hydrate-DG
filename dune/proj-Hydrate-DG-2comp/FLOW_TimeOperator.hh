/*
 * FLOW_TimeOperator.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef FLOW_TIMEOPERATOR_HH_
#define FLOW_TIMEOPERATOR_HH_

using namespace Dune::PDELab;

class FLOW_TimeOperator
	: public Dune::PDELab::NumericalJacobianApplyVolume<FLOW_TimeOperator>,
	  public Dune::PDELab::NumericalJacobianVolume<FLOW_TimeOperator>,
	  public Dune::PDELab::FullVolumePattern,
	  public Dune::PDELab::LocalOperatorDefaultFlags,
	  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
	IncludeClasses paramclass;
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

	// constructor remembers parameters
	FLOW_TimeOperator(unsigned int intorder_ = 4)
		: intorder(intorder_)
	{
		Xc_K = paramclass.characteristicValue.permeability_c;
		Xc_mu = paramclass.characteristicValue.viscosity_c;
		Xc_rho = paramclass.characteristicValue.density_c;
		Xc_kth = paramclass.characteristicValue.thermalconductivity_c;
		Xc_C = paramclass.characteristicValue.specificheat_c;
		Xc_P = paramclass.characteristicValue.P_c;
		Xc_T = paramclass.characteristicValue.T_c;
		Xc_X = paramclass.characteristicValue.x_c;
		Xc_Y = paramclass.characteristicValue.x_c;
	}

	// volume integral depending on test and ansatz functions
	template <typename EG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume(const EG &eg, const LFSU &lfsu, const X &x, const LFSV &lfsv, R &r) const
	{

		// subspaces
		//Gas pressure
		const auto &lfsv_Pg = lfsv.template child<Indices::PVId_Pg>();
		const auto &lfsu_Pg = lfsu.template child<Indices::PVId_Pg>();

		//Capillary Pressure
		const auto &lfsv_Pc = lfsv.template child<Indices::PVId_Pc>();
		const auto &lfsu_Pc = lfsu.template child<Indices::PVId_Pc>();

		//Water Saturation
		const auto &lfsv_Sw = lfsv.template child<Indices::PVId_Sw>();
		const auto &lfsu_Sw = lfsu.template child<Indices::PVId_Sw>();

		//Hydrate Saturation
		const auto &lfsv_Sh = lfsv.template child<Indices::PVId_Sh>();
		const auto &lfsu_Sh = lfsu.template child<Indices::PVId_Sh>();

		//Temperature
		const auto &lfsv_T = lfsv.template child<Indices::PVId_T>();
		const auto &lfsu_T = lfsu.template child<Indices::PVId_T>();

		//Hydrate mole fraction
    	const auto &lfsv_XCH4 = lfsv.template child<Indices::PVId_XCH4>();
    	const auto &lfsu_XCH4 = lfsu.template child<Indices::PVId_XCH4>();

    	//Water mole fraction
    	const auto &lfsv_YH2O = lfsv.template child<Indices::PVId_YH2O>();
    	const auto &lfsu_YH2O = lfsu.template child<Indices::PVId_YH2O>();

		// define types
		using RF = typename LFSU::template Child<Indices::PVId_Pg>::Type::Traits::FiniteElementType::
			Traits::LocalBasisType::Traits::RangeFieldType;
		using RangeType = typename LFSU::template Child<Indices::PVId_Pg>::Type::Traits::FiniteElementType::
			Traits::LocalBasisType::Traits::RangeType;
		using size_type = typename LFSU::template Child<Indices::PVId_Pg>::Type::Traits::SizeType;

		// Get geometry
		auto geo = eg.geometry();

		// Transformation matrix
		typename EG::Geometry::JacobianInverseTransposed jac;

		// loop over quadrature points
		for (const auto &ip : quadratureRule(geo, intorder))
		{
			// evaluate test shape functions
			std::vector<RangeType> phi_Pg(lfsu_Pg.size());
			lfsu_Pg.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Pg);
			std::vector<RangeType> psi_Pg(lfsv_Pg.size());
			lfsv_Pg.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Pg);

			std::vector<RangeType> phi_Pc(lfsu_Pc.size());
			lfsu_Pc.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Pc);
			std::vector<RangeType> psi_Pc(lfsv_Pc.size());
			lfsv_Pc.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Pc);

			std::vector<RangeType> phi_Sw(lfsu_Sw.size());
			lfsu_Sw.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Sw);
			std::vector<RangeType> psi_Sw(lfsv_Sw.size());
			lfsv_Sw.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Sw);

			std::vector<RangeType> phi_Sh(lfsu_Sh.size());
			lfsu_Sh.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Sh);
			std::vector<RangeType> psi_Sh(lfsv_Sh.size());
			lfsv_Sh.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Sh);

			std::vector<RangeType> phi_T(lfsu_T.size());
			lfsu_T.finiteElement().localBasis().evaluateFunction(ip.position(), phi_T);
			std::vector<RangeType> psi_T(lfsv_T.size());
			lfsv_T.finiteElement().localBasis().evaluateFunction(ip.position(), psi_T);

			std::vector<RangeType> phi_XCH4(lfsu_XCH4.size());
			lfsu_XCH4.finiteElement().localBasis().evaluateFunction(ip.position(), phi_XCH4);
			std::vector<RangeType> psi_XCH4(lfsv_XCH4.size());
			lfsv_XCH4.finiteElement().localBasis().evaluateFunction(ip.position(), psi_XCH4);

			std::vector<RangeType> phi_YH2O(lfsu_YH2O.size());
			lfsu_YH2O.finiteElement().localBasis().evaluateFunction(ip.position(), phi_YH2O);
			std::vector<RangeType> psi_YH2O(lfsv_YH2O.size());
			lfsv_YH2O.finiteElement().localBasis().evaluateFunction(ip.position(), psi_YH2O);

			auto ip_global = geo.global(ip.position());

			// evaluate Pg
			RF Pg = 0.0;
			for (size_type i = 0; i < lfsu_Pg.size(); i++)
				Pg += x(lfsu_Pg, i) * phi_Pg[i];

			// evaluate Sw
			RF Sw = 0.0;
			for (size_type i = 0; i < lfsu_Sw.size(); i++)
				Sw += x(lfsu_Sw, i) * phi_Sw[i];

			// evaluate Sh
			RF Sh = 0.0;
			for (size_type i = 0; i < lfsu_Sh.size(); i++)
				Sh += x(lfsu_Sh, i) * phi_Sh[i];

			// evaluate Pc
			RF Pc = 0.0;
			for (size_type i = 0; i < lfsu_Pc.size(); i++)
				Pc += x(lfsu_Pc, i) * phi_Pc[i];

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

			// evaluate Pw
			RF Pw = Pg - Pc;
			RF Peff = (Pg * (1. - Sw - Sh) + Pw * Sw) / (1. - Sh);

			auto por = paramclass.problemSpecs.SedimentPorosity(ip_global);
			auto rho_g = paramclass.methane.density(T * Xc_T, Pg * Xc_P, 1.) / Xc_rho;
			auto rho_w = paramclass.water.density(T * Xc_T, Pw * Xc_P) / Xc_rho;
			auto rho_h = paramclass.hydrate.density() / Xc_rho;
			auto rho_s = paramclass.soil.density() / Xc_rho;
			auto Cv_g = paramclass.methane.Cv(T * Xc_T, Pg * Xc_P, 1.) / Xc_C;
			auto Cv_w = paramclass.water.Cv(T * Xc_T, Pw * Xc_P) / Xc_C;
			auto Cv_h = paramclass.hydrate.Cv(T * Xc_T, Peff * Xc_P) / Xc_C;
			auto Cv_s = paramclass.soil.Cv(T * Xc_T, Peff * Xc_P) / Xc_C;
			auto Cv_eff = (1. - por) * rho_s * Cv_s + por * (rho_g * (1. - Sw - Sh) * Cv_g + rho_w * Sw * Cv_w + rho_h * Sh * Cv_h);

			//  adding terms regarding components
			auto Sg = 1. - Sw - Sh;
			auto YCH4 = paramclass.mixture.mole_y_CH4(T * Xc_T, Pg * Xc_P);
			auto XH2O = paramclass.mixture.mole_x_H2O(T * Xc_T, Pg * Xc_P);
			//  end of terms regarding components

			// integrate (A grad u - bu)*grad phi_i + a*u*phi_i
			RF factor = ip.weight() * geo.integrationElement(ip.position());
			for (size_type i = 0; i < lfsv_Pg.size(); i++)
			{
				r.accumulate(lfsv_Pg, i, ((rho_g * por * YCH4 * Sg + rho_w * por * XCH4 * Sw) * psi_Pg[i]) * factor);
			}
			for (size_type i = 0; i < lfsv_Pc.size(); i++)
			{
				r.accumulate(lfsv_Pc, i, ((rho_g * por * YH2O * Sg + rho_w * por * XH2O * Sw) * psi_Pc[i]) * factor);
			}
			for (size_type i = 0; i < lfsv_Sh.size(); i++)
			{
				r.accumulate(lfsv_Sh, i, (rho_h * por * Sh * psi_Sh[i]) * factor);
			}
			for (size_type i = 0; i < lfsv_T.size(); i++)
			{
				r.accumulate(lfsv_T, i, (Cv_eff * T * psi_T[i]) * factor);
			}

		} 	//End Quadrature Rule
	}	// End of alpha volume
};
#endif /* FLOW_TIMEOPERATOR_HH_ */
