/*
 * TimeOperator.hh
 *
 *  
 */

using namespace Dune::PDELab;

template <class GV, typename Params, typename U_Sh, class GFS_Sh,
          typename U, class GFS, typename U_T, class GFS_T, class FEM_S>
class TimeOperator_XC:
		public Dune::PDELab::NumericalJacobianApplyVolume<TimeOperator_XC<GV, Params, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, FEM_S>>,
		public Dune::PDELab::NumericalJacobianVolume<TimeOperator_XC<GV, Params, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, FEM_S>>,
		public Dune::PDELab::FullVolumePattern,
		public Dune::PDELab::LocalOperatorDefaultFlags,
		public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
	const GV &gv;
	const Params&	  property;
	U_Sh unew_Sh;
	GFS_Sh gfs_Sh;
	U unew;
	GFS gfs;
	U_T unew_T;
	GFS_T gfs_T;
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
	enum{doPatternVolume = true	};
	// residual assembly flags
	enum{doAlphaVolume = true	};

  	typedef typename GV::IndexSet IndexSet;

	using PathPw = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_Pw>>;
	using SUBGFS_Pw = Dune::PDELab::GridFunctionSubSpace<GFS,PathPw>;
	//SUBGFS_Pw    subgfs_Pw(gfs);

	using PathSg = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_Sg>>;
	using SUBGFS_Sg = Dune::PDELab::GridFunctionSubSpace<GFS,PathSg>;
	//SUBGFS_Sg    subgfs_Sg(gfs);
	
	using PathXCH4 = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_XCH4>>;
	using SUBGFS_XCH4 = Dune::PDELab::GridFunctionSubSpace<GFS,PathXCH4>;
	//SUBGFS_XCH4    subgfs_XCH4(gfs);

	using PathYH2O = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_YH2O>>;
  	using SUBGFS_YH2O = Dune::PDELab::GridFunctionSubSpace<GFS,PathYH2O>;

	using DGF_Sg = typename Dune::PDELab::DiscreteGridFunction<SUBGFS_Sg, U> ;
	using DGF_Sh = typename Dune::PDELab::DiscreteGridFunction<GFS_Sh, U_Sh> ;
	using DGF_T = typename Dune::PDELab::DiscreteGridFunction<GFS_T, U_T> ;
	using DGF_XCH4 = typename Dune::PDELab::DiscreteGridFunction<SUBGFS_XCH4, U> ;
	using DGF_YH2O = typename Dune::PDELab::DiscreteGridFunction<SUBGFS_YH2O, U> ;
	using DGF_Pw = typename Dune::PDELab::DiscreteGridFunction<SUBGFS_Pw, U> ;

	using LocalBasisType_XC = typename FEM_S::Traits::FiniteElementType::Traits::LocalBasisType;
  	using Cache_XC = typename Dune::PDELab::LocalBasisCache<LocalBasisType_XC>;
	using RF = typename LocalBasisType_XC::Traits::RangeFieldType;
  	using JacobianType = typename LocalBasisType_XC::Traits::JacobianType ;
	using RFT = typename Dune::FieldVector<double, 1>;
  	//using int = typename LocalBasisType_Sh::Traits::Type::Traits::SizeType;
  
  	// constructor remembers parameters
	TimeOperator_XC(const GV &gv_, const Params&	 property_, 
					const U_Sh &unew_Sh_, GFS_Sh gfs_Sh_, const U &unew_, GFS gfs_, 
                    const U_T &unew_T_, GFS_T gfs_T_,
                    unsigned int intorder_ = 4)
		:gv(gv_), property( property_ ),
		unew_Sh(unew_Sh_), gfs_Sh(gfs_Sh_), unew(unew_), gfs(gfs_),
        unew_T(unew_T_), gfs_T(gfs_T_), 
        intorder(intorder_)
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
		T_ref = property.parameter.ReferenceTemperature()/Xc_T;/* ndim*/
		#ifdef STATEINDEPENDENTPROPERTIES
		  		T_ref = property.parameter.RefT()/Xc_T;
		#endif
		
	}
	// volume integral depending on test and ansatz functions
	template <typename EG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume(const EG &eg, const LFSU &lfsu_XC, const X &x, const LFSV &lfsv_XC, R &r) const
	{

		
		SUBGFS_Pw gfs_Pw(gfs);
		SUBGFS_Sg gfs_Sg(gfs);
		SUBGFS_XCH4 gfs_XCH4(gfs);
		SUBGFS_YH2O gfs_YH2O(gfs);
		DGF_Sg dgf_Sg(gfs_Sg, unew);	
		DGF_Sh dgf_Sh(gfs_Sh, unew_Sh);
		DGF_T dgf_T(gfs_T, unew_T);
		DGF_XCH4 dgf_XCH4(gfs_XCH4, unew);
		DGF_YH2O dgf_YH2O(gfs_YH2O, unew);
		DGF_Pw dgf_Pw(gfs_Pw, unew);
		
		// Reference to cell
	  	const auto& cell = eg.entity();
		const IndexSet &indexSet = gv.indexSet();
		int cell_number = indexSet.index(cell);

		// Get geometry
		auto geo = eg.geometry();

		// Transformation matrix
		typename EG::Geometry::JacobianInverseTransposed jac;

		// loop over quadrature points
		for (const auto &ip : quadratureRule(geo, intorder))
		{
			// evaluate test shape functions
			// std::vector<RangeType> phi_Pw(lfsu_Pw.size());
			// lfsu_Pw.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Pw);
			// std::vector<RangeType> psi_Pw(lfsv_Pw.size());
			// lfsv_Pw.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Pw);

			std::vector<RFT> phi_XC(lfsu_XC.size());
			lfsu_XC.finiteElement().localBasis().evaluateFunction(ip.position(), phi_XC);
			std::vector<RFT> psi_XC(lfsv_XC.size());
			lfsv_XC.finiteElement().localBasis().evaluateFunction(ip.position(), psi_XC);

			// std::vector<RF> phi_Sh(lfsu_Sh.size());
			// lfsu_Sh.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Sh);
			// std::vector<RF> psi_Sh(lfsv_Sh.size());
			// lfsv_Sh.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Sh);

			// std::vector<RangeType> phi_T(lfsu_T.size());
			// lfsu_T.finiteElement().localBasis().evaluateFunction(ip.position(), phi_T);
			// std::vector<RangeType> psi_T(lfsv_T.size());
			// lfsv_T.finiteElement().localBasis().evaluateFunction(ip.position(), psi_T);

			// std::vector<RangeType> phi_XCH4(lfsu_XCH4.size());
			// lfsu_XCH4.finiteElement().localBasis().evaluateFunction(ip.position(), phi_XCH4);
			// std::vector<RangeType> psi_XCH4(lfsv_XCH4.size());
			// lfsv_XCH4.finiteElement().localBasis().evaluateFunction(ip.position(), psi_XCH4);

			// std::vector<RangeType> phi_YH2O(lfsu_YH2O.size());
			// lfsu_YH2O.finiteElement().localBasis().evaluateFunction(ip.position(), phi_YH2O);
			// std::vector<RangeType> psi_YH2O(lfsv_YH2O.size());
			// lfsv_YH2O.finiteElement().localBasis().evaluateFunction(ip.position(), psi_YH2O);

			// std::vector<RangeType> phi_XC(lfsu_XC.size());
			// lfsu_XC.finiteElement().localBasis().evaluateFunction(ip.position(), phi_XC);
			// std::vector<RangeType> psi_XC(lfsv_XC.size());
			// lfsv_XC.finiteElement().localBasis().evaluateFunction(ip.position(), psi_XC);

			auto ip_global = geo.global(ip.position());
			auto ip_local = geo.local(ip_global);
			
			// evaluate XC
			RF XC = 0.0;
			for (int i = 0; i < lfsu_XC.size(); i++){
			  XC += x(lfsu_XC, i) * phi_XC[i];
			}

			// evaluate Sg
			RFT Sg0 = 0.0;
			dgf_Sg.evaluate(cell, ip.position(), Sg0);
			RF Sg = Sg0[0];
						
			// evaluate Sh
			RFT Sh0 = 0.0;
			dgf_Sh.evaluate(cell, ip.position(), Sh0);
			RF Sh = Sh0[0];
			
			// evaluate T
			RFT T0 = 0.0;
			dgf_T.evaluate(cell, ip.position(), T0);
			RF T =T0[0];
			
			// evaluate XCH4
			RFT XCH40 = 0.0;
			dgf_XCH4.evaluate(cell, ip.position(), XCH40);
			RF XCH4 = XCH40[0];
			
			// evaluate YH2O
			RFT YH2O0 = 0.0;
			dgf_YH2O.evaluate(cell, ip.position(), YH2O0);
			RF YH2O = YH2O0[0] ;
			
			// evaluate Pw
			RFT Pw0 = 0.0;
			dgf_Pw.evaluate(cell, ip.position(), Pw0);
			RF Pw = Pw0[0];
			
			RF Sw = 1. - Sg - Sh;

			// evaluate Pg
      		auto por = property.soil.SedimentPorosity(cell, ip_local);
      		auto Pc = property.hydraulicProperty.CapillaryPressure(cell, ip_local, Sw, Sh, por) ; /* ndim */
      
			RF Pg = Pw + Pc ;
			RF Peff = (Pg * Sg + Pw * Sw) / (1. - Sh);

			double S = XC * (property.salt.MolarMass()/property.water.MolarMass());
      		auto zCH4 = property.eos.EvaluateCompressibilityFactor(T * Xc_T, Pg * Xc_P);
			
			auto rho_w = property.water.Density(T * Xc_T, Pw * Xc_P, S);
			
			// integrate (A grad u - bu)*grad phi_i + a*u*phi_i
			RF factor = ip.weight() * geo.integrationElement(ip.position());
			// for (int i = 0; i < lfsv_Sg.size(); i++)
			// {
			// 	r.accumulate(lfsv_Sg, i, ((rho_g * por * (1. - YH2O) * Sg + rho_w * por * XCH4 * Sw) * psi_Sg[i]) * factor);
			// }
			for (int i = 0; i < lfsv_XC.size(); i++)
			{
				r.accumulate(lfsv_XC, i, (rho_w * por * XC * Sw * psi_XC[i]) * factor);
			}
			// for (int i = 0; i < lfsv_Pw.size(); i++)
			// {
			// 	r.accumulate(lfsv_Pw, i, ((rho_g * por * YH2O * Sg + rho_w * por * (1. -XC - XCH4) * Sw) * psi_Pw[i]) * factor);
			// }
			// for (int i = 0; i < lfsv_Sh.size(); i++)
			// {
			// 	r.accumulate(lfsv_Sh, i, (rho_h * por * Sh * psi_Sh[i]) * factor);
			// }
			// for (int i = 0; i < lfsv_T.size(); i++)
			// {
			// 	r.accumulate(lfsv_T, i, Cv_eff * (T-T_ref) * psi_T[i] * factor);
			// }

		} 	//End Quadrature Rule
	}	// End of alpha volume
	
	
 	// jacobian contribution of volume term
	// template <typename EG, typename LFSU, typename X, typename LFSV, typename M>
	// void jacobian_volume(const EG &eg, const LFSU &lfsu_Sh, const X &x, const LFSV &lfsv_Sh, M& mat) const
	// {
	// 	DGF_Pw dgf_Pw(gfs_Pw, unew_Pw);	
	// 	DGF_Sg dgf_Sg(gfs_Sg, unew_Sg);
	// 	DGF_T dgf_T(gfs_T, unew_T);
	// 	DGF_XCH4 dgf_XCH4(gfs_XCH4, unew_XCH4);
	// 	DGF_YH2O dgf_YH2O(gfs_YH2O, unew_YH2O);
	// 	DGF_XC dgf_XC(gfs_XC, unew_XC);

	// 	// Reference to cell
	//   	const auto& cell = eg.entity();
	// 	const IndexSet &indexSet = gv.indexSet();
	// 	int cell_number = indexSet.index(cell);
		
	// 	// Get geometry
	// 	auto geo = eg.geometry();

	// // 	// Transformation matrix
	// // 	typename EG::Geometry::JacobianInverseTransposed jac;

	// 	// loop over quadrature points
	// 	for (const auto &ip : quadratureRule(geo, intorder))
	// 	{
	// // 		// evaluate test shape functions
	// // 		std::vector<RangeType> phi_Pw(lfsu_Pw.size());
	// // 		lfsu_Pw.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Pw);
	// // 		std::vector<RangeType> psi_Pw(lfsv_Pw.size());
	// // 		lfsv_Pw.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Pw);

	// // 		std::vector<RangeType> phi_Sg(lfsu_Sg.size());
	// // 		lfsu_Sg.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Sg);
	// // 		std::vector<RangeType> psi_Sg(lfsv_Sg.size());
	// // 		lfsv_Sg.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Sg);

	// 		std::vector<RF> phi_Sh(lfsu_Sh.size());
	// 		lfsu_Sh.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Sh);
	// 		std::vector<RF> psi_Sh(lfsv_Sh.size());
	// 		lfsv_Sh.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Sh);

	// // 		std::vector<RangeType> phi_T(lfsu_T.size());
	// // 		lfsu_T.finiteElement().localBasis().evaluateFunction(ip.position(), phi_T);
	// // 		std::vector<RangeType> psi_T(lfsv_T.size());
	// // 		lfsv_T.finiteElement().localBasis().evaluateFunction(ip.position(), psi_T);

	// // 		std::vector<RangeType> phi_XCH4(lfsu_XCH4.size());
	// // 		lfsu_XCH4.finiteElement().localBasis().evaluateFunction(ip.position(), phi_XCH4);
	// // 		std::vector<RangeType> psi_XCH4(lfsv_XCH4.size());
	// // 		lfsv_XCH4.finiteElement().localBasis().evaluateFunction(ip.position(), psi_XCH4);

	// // 		std::vector<RangeType> phi_YH2O(lfsu_YH2O.size());
	// // 		lfsu_YH2O.finiteElement().localBasis().evaluateFunction(ip.position(), phi_YH2O);
	// // 		std::vector<RangeType> psi_YH2O(lfsv_YH2O.size());
	// // 		lfsv_YH2O.finiteElement().localBasis().evaluateFunction(ip.position(), psi_YH2O);

	// // 		std::vector<RangeType> phi_XC(lfsu_XC.size());
	// // 		lfsu_XC.finiteElement().localBasis().evaluateFunction(ip.position(), phi_XC);
	// // 		std::vector<RangeType> psi_XC(lfsv_XC.size());
	// // 		lfsv_XC.finiteElement().localBasis().evaluateFunction(ip.position(), psi_XC);

	// 		auto ip_global = geo.global(ip.position());
	// 		auto ip_local = geo.local(ip_global);

	// 		// evaluate Pw
	// 		RF Pw = 0.0;
	// 		dgf_Pw.evaluate(cell, ip.position(), Pw);
	// 		// for (int i = 0; i < lfsu_Pw.size(); i++){
	// 		//   Pw += x(lfsu_Pw, i) * phi_Pw[i];
	// 		// }

	// 		// evaluate Sg
	// 		RF Sg = 0.1;
	// 		dgf_Sg.evaluate(cell, ip.position(), Sg);
	// 		// for (int i = 0; i < lfsu_Sg.size(); i++){
	// 		//   Sg += x(lfsu_Sg, i) * phi_Sg[i];
	// 		// }
			
	// 		// evaluate Sh
	// 		RF Sh = 0.0;
	// 		for (int i = 0; i < lfsu_Sh.size(); i++){
	// 			Sh += x(lfsu_Sh, i) * phi_Sh[i];
	// 		}

	// 		// evaluate T
	// 		RF T = 2.0;
	// 		dgf_T.evaluate(cell, ip.position(), T);
	// 		// for (int i = 0; i < lfsu_T.size(); i++)
	// 		//   T += x(lfsu_T, i) * phi_T[i];
			
	// 		// evaluate XCH4
	// 		RF XCH4 = 0.0;
	// 		dgf_XCH4.evaluate(cell, ip.position(), XCH4);
	// 		// for (int i = 0; i < lfsu_XCH4.size(); i++)
	// 		//   XCH4 += x(lfsu_XCH4, i) * phi_XCH4[i];

	// 		// evaluate YH2O
	// 		RF YH2O = 0.0;
	// 		dgf_YH2O.evaluate(cell, ip.position(), YH2O);
	// 		// for (int i = 0; i < lfsu_YH2O.size(); i++)
	// 		//   YH2O += x(lfsu_YH2O, i) * phi_YH2O[i];

	// 		// evaluate XC
	// 		RF XC = 0.0;
	// 		dgf_XC.evaluate(cell, ip.position(), XC);
	// 		// for (int i = 0; i < lfsu_XC.size(); i++)
	// 		//   XC += x(lfsu_XC, i) * phi_XC[i];

	// 		RF Sw = 1. - Sg - Sh;

	// 		// evaluate Pg
    //   		auto por = property.soil.SedimentPorosity(cell, ip_local);
    //   		auto Pc = property.hydraulicProperty.CapillaryPressure(cell, ip_local, Sw, Sh, por) ; /* ndim */
	// 		RF Pg = Pw + Pc;
	// 		RF Peff = (Pg * Sg + Pw * Sw) / (1. - Sh);
			
	// 		double S = XC * (property.salt.MolarMass()/property.gas.MolarMass());
    //   		auto zCH4 = property.eos.EvaluateCompressibilityFactor(T * Xc_T, Pg * Xc_P);
	// // 		auto rho_g = property.gas.Density(T * Xc_T, Pg * Xc_P, zCH4);
	// // 		auto rho_w = property.water.Density(T * Xc_T, Pw * Xc_P, S);
	// 		auto rho_h = property.hydrate.Density() ;
	// // 		auto rho_s = property.soil.Density() ;
	// // 		auto Cv_g = property.gas.Cv(T * Xc_T, Pg * Xc_P, zCH4) ;
	// // 		auto Cv_w = property.water.Cv(T * Xc_T, Pw * Xc_P, S) ;
	// // 		auto Cv_h = property.hydrate.Cv(T * Xc_T, Peff * Xc_P) ;
	// // 		auto Cv_s = property.soil.Cv();
	// // 		auto Cv_eff = (1. - por) * rho_s * Cv_s + por * (rho_g * (1. - Sw - Sh) * Cv_g + rho_w * Sw * Cv_w + rho_h * Sh * Cv_h);

	// // 		// integrate (A grad u - bu)*grad phi_i + a*u*phi_i
	// 		RF factor = ip.weight() * geo.integrationElement(ip.position());
	// // 		for (int i = 0; i < lfsv_Sg.size(); i++)
	// // 		{
	// // 			for (int j = 0; j < lfsu_XCH4.size(); j++)
	// // 			{
	// // 				mat.accumulate(lfsv_Sg, i, lfsu_XCH4, j , por * rho_w * Sw * phi_XCH4[j] * psi_Sg[i] * factor);
	// // 			}
	// // 			for (int j = 0; j < lfsu_Sg.size(); j++)
	// // 			{
	// // 				mat.accumulate(lfsv_Sg, i, lfsu_Sg, j , por * (-rho_w * XCH4 + rho_g * (1. - YH2O)) * phi_Sg[j] * psi_Sg[i] * factor);
	// // 			}
	// // 			for (int j = 0; j < lfsu_Sh.size(); j++)
	// // 			{
	// // 				mat.accumulate(lfsv_Sg, i, lfsu_Sh, j , por * rho_w * XCH4 * -phi_Sh[j] * psi_Sg[i] * factor);
	// // 			}
	// // 			for (int j = 0; j < lfsu_YH2O.size(); j++)
	// // 			{
	// // 				mat.accumulate(lfsv_Sg, i, lfsu_YH2O, j , por * rho_g * Sg * -phi_YH2O[j] * psi_Sg[i] * factor);
	// // 			}
	// // 		}

	// // 		for (int i = 0; i < lfsv_XC.size(); i++)
	// // 		{
	// // 			for (int j = 0; j < lfsu_XC.size(); j++)
	// // 			{
	// // 				mat.accumulate(lfsv_XC, i, lfsu_XC, j , por * rho_w * Sw * phi_XC[j] * psi_XC[i] * factor);
	// // 			}
	// // 			for (int j = 0; j < lfsu_Sg.size(); j++)
	// // 			{
	// // 				mat.accumulate(lfsv_XC, i, lfsu_Sg, j , por * rho_w * XC  * -phi_Sg[j] * psi_XC[i] * factor);
	// // 			}
	// // 			for (int j = 0; j < lfsu_Sh.size(); j++)
	// // 			{
	// // 				mat.accumulate(lfsv_XC, i, lfsu_Sh, j , por * rho_w * XC * -phi_Sh[j] * psi_XC[i] * factor);
	// // 			}
	// // 		}
	// // 		for (int i = 0; i < lfsv_Pw.size(); i++)
	// // 		{
	// // 			for (int j = 0; j < lfsu_XCH4.size(); j++)
	// // 			{
	// // 				mat.accumulate(lfsv_Pw, i, lfsu_XCH4, j , por * rho_w * Sw * -phi_XCH4[j] * psi_Pw[i] * factor);
	// // 			}
	// // 			for (int j = 0; j < lfsu_XC.size(); j++)
	// // 			{
	// // 				mat.accumulate(lfsv_Pw, i, lfsu_XC, j , por * rho_w * Sw * -phi_XC[j] * psi_Pw[i] * factor);
	// // 			}
	// // 			for (int j = 0; j < lfsu_Sg.size(); j++)
	// // 			{
	// // 				mat.accumulate(lfsv_Pw, i, lfsu_Sg, j , por * (-rho_w * (1. - XC - XCH4) + rho_g * YH2O) * phi_Sg[j] * psi_Pw[i] * factor);
	// // 			}
	// // 			for (int j = 0; j < lfsu_Sh.size(); j++)
	// // 			{
	// // 				mat.accumulate(lfsv_Pw, i, lfsu_Sh, j , por * rho_w * (1. - XC - XCH4) * -phi_Sh[j] * psi_Pw[i] * factor);
	// // 			}
	// // 			for (int j = 0; j < lfsu_YH2O.size(); j++)
	// // 			{
	// // 				mat.accumulate(lfsv_Pw, i, lfsu_YH2O, j , por * rho_g * Sg * phi_YH2O[j] * psi_Pw[i] * factor);
	// // 			}
	// // 		}
	// 		for (int i = 0; i < lfsv_Sh.size(); i++)
	// 		{
	// 			for (int j = 0; j < lfsu_Sh.size(); j++)
	// 			{
	// 				mat.accumulate(lfsv_Sh, i, lfsu_Sh, j , por * rho_h * phi_Sh[j] * psi_Sh[i] * factor);
	// 			}
	// 		}
	// // 		for (int i = 0; i < lfsv_T.size(); i++)
	// // 		{
	// // 			for (int j = 0; j < lfsu_T.size(); j++)
	// // 			{
	// // 				mat.accumulate(lfsv_T, i, lfsu_T, j , Cv_eff * phi_T[j] * psi_T[i] * factor);
	// // 			}
	// // 			for (int j = 0; j < lfsu_Sg.size(); j++)
	// // 			{
	// // 				mat.accumulate(lfsv_T, i, lfsu_Sg, j , por * (-rho_w * Cv_w + rho_g * Cv_g) * (T-T_ref) * phi_Sg[j] * psi_T[i] * factor);
	// // 			}
	// // 			for (int j = 0; j < lfsu_Sh.size(); j++)
	// // 			{
	// // 				mat.accumulate(lfsv_T, i, lfsu_Sh, j , por * (-rho_w * Cv_w + rho_h * Cv_h) * (T-T_ref) * phi_Sh[j] * psi_T[i] * factor);
	// // 			}
	// // 		}

	// 	} 	//End Quadrature Rule
	// }	// End of jacobian volume

};