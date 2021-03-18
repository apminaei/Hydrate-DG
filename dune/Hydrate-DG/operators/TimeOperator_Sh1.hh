/*
 * TimeOperator.hh
 *
 *  
 */

using namespace Dune::PDELab;

template <class GV, typename Params,
          typename U, class GFS, typename U_T, class GFS_T,
          class FEM_S>
class TimeOperator_Sh:
		// public Dune::PDELab::NumericalJacobianApplyVolume<TimeOperator_Sh<GV, Params, 
        //               U_Pw, GFS_Pw, U_Sg, GFS_Sg, U_T, GFS_T, U_XCH4, GFS_XCH4, U_YH2O, GFS_YH2O, U_XC, GFS_XC, FEM_S>>,
		// public Dune::PDELab::NumericalJacobianVolume<TimeOperator_Sh<GV, Params,
                    //   U_Pw, GFS_Pw, U_Sg, GFS_Sg, U_T, GFS_T, U_XCH4, GFS_XCH4, U_YH2O, GFS_YH2O, U_XC, GFS_XC, FEM_S>>,
		public Dune::PDELab::FullVolumePattern,
		public Dune::PDELab::LocalOperatorDefaultFlags,
		public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
	const GV &gv;
	const Params&	  property;
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

	  using PathXC = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_XC>>;
	using SUBGFS_XC = Dune::PDELab::GridFunctionSubSpace<GFS,PathXC>;
	//SUBGFS_XCH4    subgfs_XCH4(gfs);

	using DGF_Pw = typename Dune::PDELab::DiscreteGridFunction<SUBGFS_Pw, U> ;
  	using DGF_Sg = typename Dune::PDELab::DiscreteGridFunction<SUBGFS_Sg, U> ;
	using DGF_T = typename Dune::PDELab::DiscreteGridFunction<GFS_T, U_T> ;
	using DGF_XCH4 = typename Dune::PDELab::DiscreteGridFunction<SUBGFS_XCH4, U> ;
	using DGF_YH2O = typename Dune::PDELab::DiscreteGridFunction<SUBGFS_YH2O, U> ;
	using DGF_XC = typename Dune::PDELab::DiscreteGridFunction<SUBGFS_XC, U> ;

	using LocalBasisType_Sh = typename FEM_S::Traits::FiniteElementType::Traits::LocalBasisType;
  	using Cache_Sh = typename Dune::PDELab::LocalBasisCache<LocalBasisType_Sh>;
	using RF = typename LocalBasisType_Sh::Traits::RangeFieldType;
  	using JacobianType = typename LocalBasisType_Sh::Traits::JacobianType ;
	using RFT = typename Dune::FieldVector<double, 1>;
  	//using size_type = typename LocalBasisType_Sh::Traits::Type::Traits::SizeType;
  
  	// constructor remembers parameters
	TimeOperator_Sh(const GV &gv_, const Params&	 property_, 
					const U &unew_, GFS gfs_, 
                    const U_T &unew_T_, GFS_T gfs_T_, unsigned int intorder_ = 4)
		:gv(gv_), property( property_ ), unew(unew_), gfs(gfs_),
        unew_T(unew_T_), gfs_T(gfs_T_),  intorder(intorder_)
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
	void alpha_volume(const EG &eg, const LFSU &lfsu_Sh, const X &x, const LFSV &lfsv_Sh, R &r) const
	{
		
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

			std::vector<RFT> phi_Sh(lfsu_Sh.size());
			lfsu_Sh.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Sh);
			std::vector<RFT> psi_Sh(lfsv_Sh.size());
			lfsv_Sh.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Sh);

			auto ip_global = geo.global(ip.position());
			auto ip_local = geo.local(ip_global);
			
			
			// evaluate Sh
			RF Sh = 0.0;
			for (int i = 0; i < lfsu_Sh.size(); i++){
				Sh += x(lfsu_Sh, i) * phi_Sh[i];
			}

      		auto por = property.soil.SedimentPorosity(cell, ip_local);
			auto rho_h = property.hydrate.Density() ;

			// integrate (A grad u - bu)*grad phi_i + a*u*phi_i
			RF factor = ip.weight() * geo.integrationElement(ip.position());
			
			for (int i = 0; i < lfsv_Sh.size(); i++)
			{
				r.accumulate(lfsv_Sh, i, (rho_h * por * Sh * psi_Sh[i]) * factor);
			}
		} 	//End Quadrature Rule
	}	// End of alpha volume
	
	
 	// jacobian contribution of volume term
	template <typename EG, typename LFSU, typename X, typename LFSV, typename M>
	void jacobian_volume(const EG &eg, const LFSU &lfsu_Sh, const X &x, const LFSV &lfsv_Sh, M& mat) const
	{
		
		// Reference to cell
	  	const auto& cell = eg.entity();
		const IndexSet &indexSet = gv.indexSet();
		int cell_number = indexSet.index(cell);
		
		// Get geometry
		auto geo = eg.geometry();


		// loop over quadrature points
		for (const auto &ip : quadratureRule(geo, intorder))
		{

			std::vector<RFT> phi_Sh(lfsu_Sh.size());
			lfsu_Sh.finiteElement().localBasis().evaluateFunction(ip.position(), phi_Sh);
			std::vector<RFT> psi_Sh(lfsv_Sh.size());
			lfsv_Sh.finiteElement().localBasis().evaluateFunction(ip.position(), psi_Sh);

			auto ip_global = geo.global(ip.position());
			auto ip_local = geo.local(ip_global);

			
			// evaluate Sh
			RF Sh = 0.0;
			for (int i = 0; i < lfsu_Sh.size(); i++){
				Sh += x(lfsu_Sh, i) * phi_Sh[i];
			}

			
      		auto por = property.soil.SedimentPorosity(cell, ip_local);
			auto rho_h = property.hydrate.Density() ;

			// integrate (A grad u - bu)*grad phi_i + a*u*phi_i
			RF factor = ip.weight() * geo.integrationElement(ip.position());
	
			for (int i = 0; i < lfsv_Sh.size(); i++)
			{
				for (int j = 0; j < lfsu_Sh.size(); j++)
				{
					mat.accumulate(lfsv_Sh, i, lfsu_Sh, j , por * rho_h * phi_Sh[j] * psi_Sh[i] * factor);
				}
			}

		} 	//End Quadrature Rule
	}	// End of jacobian volume

};