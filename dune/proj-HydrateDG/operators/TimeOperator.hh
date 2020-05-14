/*
 * TimeOperator.hh
 *
 *  Created on: Sep 22, 2016
 *      Author: shubhangi
 */

#ifndef FLOW_TIMEOPERATOR_HH_
#define FLOW_TIMEOPERATOR_HH_
using namespace Dune::PDELab;

template <class GV, typename Params>
class TimeOperator
	: public Dune::PDELab::NumericalJacobianApplyVolume<TimeOperator<GV, Params>>,
	  public Dune::PDELab::NumericalJacobianVolume<TimeOperator<GV, Params>>,
	  public Dune::PDELab::FullVolumePattern,
	  public Dune::PDELab::LocalOperatorDefaultFlags,
	  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
	const GV &gv;
	const Params&	  property;
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
	TimeOperator(const GV &gv_, const Params&	 property_, unsigned int intorder_ = 4)
		:gv(gv_), property( property_ ),intorder(intorder_)
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
	}
	// volume integral depending on test and ansatz functions
	template <typename EG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume(const EG &eg, const LFSU &lfsu, const X &x, const LFSV &lfsv, R &r) const
	{

		// subspaces
		//Gas pressure
		const auto &lfsv_Pw = lfsv.template child<Indices::PVId_Pw>();
		const auto &lfsu_Pw = lfsu.template child<Indices::PVId_Pw>();

		//Capillary Pressure
		const auto &lfsv_Pc = lfsv.template child<Indices::PVId_Pc>();
		const auto &lfsu_Pc = lfsu.template child<Indices::PVId_Pc>();

		//Water Saturation
		const auto &lfsv_Sg = lfsv.template child<Indices::PVId_Sg>();
		const auto &lfsu_Sg = lfsu.template child<Indices::PVId_Sg>();

		//Hydrate Saturation
		const auto &lfsv_Sh = lfsv.template child<Indices::PVId_Sh>();
		const auto &lfsu_Sh = lfsu.template child<Indices::PVId_Sh>();

		//Temperature
		const auto &lfsv_T = lfsv.template child<Indices::PVId_T>();
		const auto &lfsu_T = lfsu.template child<Indices::PVId_T>();

		//Methane mole fraction
    	const auto &lfsv_XCH4 = lfsv.template child<Indices::PVId_XCH4>();
    	const auto &lfsu_XCH4 = lfsu.template child<Indices::PVId_XCH4>();

    	//Water mole fraction
    	const auto &lfsv_YH2O = lfsv.template child<Indices::PVId_YH2O>();
    	const auto &lfsu_YH2O = lfsu.template child<Indices::PVId_YH2O>();
		
		//Salt mole fraction
    	const auto &lfsv_XC = lfsv.template child<Indices::PVId_C>();
    	const auto &lfsu_XC = lfsu.template child<Indices::PVId_C>();

		// define types
		using RF = typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::FiniteElementType::
			Traits::LocalBasisType::Traits::RangeFieldType;
		using RangeType = typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::FiniteElementType::
			Traits::LocalBasisType::Traits::RangeType;
		using size_type = typename LFSU::template Child<Indices::PVId_Pw>::Type::Traits::SizeType;

		// Get geometry
		auto geo = eg.geometry();

		// Transformation matrix
		typename EG::Geometry::JacobianInverseTransposed jac;

		// loop over quadrature points
		for (const auto &ip : quadratureRule(geo, intorder))
		{
			// evaluate test shape functions
			std::vector<RangeType> phi_Pw(lfsu_Pw.size());
			lfsu_Pw.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Pw);
			std::vector<RangeType> psi_Pw(lfsv_Pw.size());
			lfsv_Pw.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Pw);

			std::vector<RangeType> phi_Pc(lfsu_Pc.size());
			lfsu_Pc.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Pc);
			std::vector<RangeType> psi_Pc(lfsv_Pc.size());
			lfsv_Pc.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Pc);

			std::vector<RangeType> phi_Sg(lfsu_Sg.size());
			lfsu_Sg.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Sg);
			std::vector<RangeType> psi_Sg(lfsv_Sg.size());
			lfsv_Sg.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Sg);

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

			std::vector<RangeType> phi_XC(lfsu_XC.size());
			lfsu_XC.finiteElement().localBasis().evaluateFunction(ip.position(), phi_XC);
			std::vector<RangeType> psi_XC(lfsv_XC.size());
			lfsv_XC.finiteElement().localBasis().evaluateFunction(ip.position(), psi_XC);

			auto ip_global = geo.global(ip.position());

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

			// evaluate XC
      		RF XC = 0.0;
      		for (size_type i = 0; i < lfsu_XC.size(); i++)
        		XC += x(lfsu_XC, i) * phi_XC[i];

			RF Sw = 1. - Sg - Sh;
			// evaluate Pw
			RF Pg = Pw + Pc;
			RF Peff = (Pg * Sg + Pw * Sw) / (1. - Sh);

			auto por = property.soil.SedimentPorosity(ip_global);
			auto rho_g = property.methane.density(T * Xc_T, Pg * Xc_P, 1.) / Xc_rho;
			auto rho_w = property.water.density(T * Xc_T, Pw * Xc_P) / Xc_rho;
			auto rho_h = property.hydrate.density() / Xc_rho;
			auto rho_s = property.soil.density() / Xc_rho;
			auto Cv_g = property.methane.Cv(T * Xc_T, Pg * Xc_P, 1.) / Xc_C;
			auto Cv_w = property.water.Cv(T * Xc_T, Pw * Xc_P) / Xc_C;
			auto Cv_h = property.hydrate.Cv(T * Xc_T, Peff * Xc_P) / Xc_C;
			auto Cv_s = property.soil.Cv(T * Xc_T, Peff * Xc_P) / Xc_C;
			auto Cv_eff = (1. - por) * rho_s * Cv_s + por * (rho_g * (1. - Sw - Sh) * Cv_g + rho_w * Sw * Cv_w + rho_h * Sh * Cv_h);

			//  adding terms regarding components
			auto YCH4 = property.mixture.mole_y_CH4(T * Xc_T, Pg * Xc_P);
			auto XH2O = property.mixture.mole_x_H2O(T * Xc_T, Pg * Xc_P);
			//  end of terms regarding components

			// integrate (A grad u - bu)*grad phi_i + a*u*phi_i
			RF factor = ip.weight() * geo.integrationElement(ip.position());
			for (size_type i = 0; i < lfsv_Pw.size(); i++)
			{
				r.accumulate(lfsv_Pw, i, ((rho_g * por * YCH4 * Sg + rho_w * por * XCH4 * Sw) * psi_Pw[i]) * factor);
			}
			for (size_type i = 0; i < lfsv_XC.size(); i++)
			{
				r.accumulate(lfsv_XC, i, (rho_w * por * XC * Sw * psi_XC[i]) * factor);
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
