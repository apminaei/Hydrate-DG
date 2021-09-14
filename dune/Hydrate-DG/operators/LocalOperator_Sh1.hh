/*
 * LocalOperator.hh
 *
 * 
 */



using namespace Dune::PDELab;


template <typename GV, typename Params, class BC, typename U_Sh, class GFS_Sh,
          typename U, class GFS, typename U_T, class GFS_T, class FEM_S>
class LocalOperator_Sh : 
                      public Dune::PDELab::NumericalJacobianApplyVolume<LocalOperator_Sh<GV, Params, BC, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, FEM_S>>,
                      public Dune::PDELab::NumericalJacobianVolume<LocalOperator_Sh<GV, Params, BC, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, FEM_S>>,
                      public Dune::PDELab::NumericalJacobianApplySkeleton<LocalOperator_Sh<GV, Params, BC, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, FEM_S>>,
                      public Dune::PDELab::NumericalJacobianSkeleton<LocalOperator_Sh<GV, Params, BC, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, FEM_S>>,  
                      public Dune::PDELab::NumericalJacobianApplyBoundary<LocalOperator_Sh<GV, Params, BC, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, FEM_S>>,
                      public Dune::PDELab::NumericalJacobianBoundary<LocalOperator_Sh<GV, Params, BC, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, FEM_S>>,
                      public Dune::PDELab::FullSkeletonPattern, // matrix entries skeleton
                      public Dune::PDELab::FullVolumePattern,
                      public Dune::PDELab::LocalOperatorDefaultFlags,
                      public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
  
	
  const GV &gv;
  const Params&	  property;
	const BC&	 	  bc;

  U_Sh *unew_Sh;
  GFS_Sh gfs_Sh;
  U unew;
  GFS gfs;
  U_T unew_T;
  GFS_T gfs_T;

  double *time;
  double *dt;
  double alpha_g;
  double alpha_w;
  double alpha_s;
  double alpha_T;
  double alpha_x;
  double alpha_y;

  double theta_g;
  double theta_w;
  double theta_T;
  double theta_x;
  double theta_y;
  constexpr static double eps = 1.0e-6;
  constexpr static double eps_ap	= 0.;
  constexpr static double pi = atan(1.) * 4;
  unsigned int intorder;
  double Xc_conv_m;
  double Xc_conv_h;
  double Xc_source_m;
  double Xc_source_h;
  double Xc_diff_m;
  double Xc_diff_h;
  double Xc_grav;
  double Xc_K;
  double Xc_mu;
  double Xc_rho;
  double Xc_kth;
  double Xc_C;
  double Xc_P;
  double Xc_T;
  double Xc_t;

public:
  // pattern assembly flags
	enum { doPatternVolume	= true };
	enum { doPatternSkeleton	= false };

	// residual assembly flags
	enum { doAlphaVolume  	= true };
	enum { doAlphaSkeleton	= false };
	enum { doAlphaBoundary	= false };
  
  typedef typename GV::IndexSet IndexSet;

  
  typedef Dune::PDELab::LocalFunctionSpace<GFS_Sh> LFS;
  typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
  typedef typename U_Sh::template LocalView<LFSCache> VectorView;
  
  using PathPw = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_Pw>>;
  using SUBGFS_Pw = Dune::PDELab::GridFunctionSubSpace<GFS,PathPw>;

	using PathSg = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_Sg>>;
  using SUBGFS_Sg = Dune::PDELab::GridFunctionSubSpace<GFS,PathSg>;
	
	using PathXCH4 = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_XCH4>>;
  using SUBGFS_XCH4 = Dune::PDELab::GridFunctionSubSpace<GFS,PathXCH4>;

	using PathYH2O = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_YH2O>>;
  using SUBGFS_YH2O = Dune::PDELab::GridFunctionSubSpace<GFS,PathYH2O>;

  using PathXC = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::VId_XC>>;
  using SUBGFS_XC = Dune::PDELab::GridFunctionSubSpace<GFS,PathXC>;
  
	
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

  // In theory it is possible that one and the same local operator is
  // called first with a finite element of one type and later with a
  // finite element of another type.  Since finite elements of different
  // type will usually produce different results for the same local
  // coordinate they cannot share a cache.  Here we use a vector of caches
  // to allow for different orders of the shape functions, which should be
  // enough to support p-adaptivity.  (Another likely candidate would be
  // differing geometry types, i.e. hybrid meshes.)

  std::vector<Cache_Sh> cache_Sh;

  // constructor stores parameters
  LocalOperator_Sh( const GV &gv_, const Params&	property_, const BC& bc_,
                    U_Sh *unew_Sh_, GFS_Sh gfs_Sh_, const U &unew_, GFS gfs_, 
                    const U_T &unew_T_, GFS_T gfs_T_,
                    double *time_,
                    double *dt_,
                    unsigned int intorder_ = 6)
      : gv(gv_), property( property_ ), bc( bc_ ),
        unew_Sh(unew_Sh_), gfs_Sh(gfs_Sh_), unew(unew_), gfs(gfs_),
        unew_T(unew_T_), gfs_T(gfs_T_), 
        time(time_),
        dt(dt_),
        intorder(intorder_), cache_Sh(20)
  {
    Xc_conv_m = property.characteristicValue.X_convective_mass;
    Xc_conv_h = property.characteristicValue.X_convective_heat;
    Xc_source_m = property.characteristicValue.X_source_mass;
    Xc_source_h = property.characteristicValue.X_source_heat;
    Xc_diff_m = property.characteristicValue.X_diffusive_mass;
    Xc_diff_h = property.characteristicValue.X_diffusive_heat;
    Xc_grav = property.characteristicValue.X_gravity;

    Xc_K = property.characteristicValue.permeability_c;
    Xc_mu = property.characteristicValue.viscosity_c;
    Xc_rho = property.characteristicValue.density_c;
    Xc_kth = property.characteristicValue.thermalconductivity_c;
    Xc_C = property.characteristicValue.specificheat_c;
    Xc_P = property.characteristicValue.P_c;
    Xc_T = property.characteristicValue.T_c;
    Xc_t = property.characteristicValue.t_c;
    
  }

  // volume integral depending on test and ansatz functions
  template <typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG &eg, const LFSU &lfsu_Sh, const X &x, const LFSV &lfsv_Sh, R &r) const
  {
    SUBGFS_Pw gfs_Pw(gfs);
    SUBGFS_Sg gfs_Sg(gfs);
    SUBGFS_XCH4 gfs_XCH4(gfs);
    SUBGFS_XC gfs_XC(gfs);
    SUBGFS_YH2O gfs_YH2O(gfs);
    DGF_Pw dgf_Pw(gfs_Pw, unew);	
    DGF_Sg dgf_Sg(gfs_Sg, unew);
    DGF_T dgf_T(gfs_T, unew_T);
    DGF_XCH4 dgf_XCH4(gfs_XCH4, unew);
    DGF_YH2O dgf_YH2O(gfs_YH2O, unew);
    DGF_XC dgf_XC(gfs_XC, unew);

    // dimensions
    const int dim = EG::Entity::dimension;
    const int order_s = std::max(lfsu_Sh.finiteElement().localBasis().order(),
                               lfsv_Sh.finiteElement().localBasis().order());

    // Reference to cell
	  const auto& cell = eg.entity();
		const IndexSet &indexSet = gv.indexSet();
		int cell_number = indexSet.index(cell);

    // Get geometry
    auto geo = eg.geometry();
    
    // loop over quadrature points
    //      auto intorder = intorderadd + quadrature_factor * order;
    for (const auto &ip : quadratureRule(geo, intorder))
    {
      auto &phi_Sh = cache_Sh[order_s].evaluateFunction(ip.position(), lfsu_Sh.finiteElement().localBasis());
      auto &psi_Sh = cache_Sh[order_s].evaluateFunction(ip.position(), lfsv_Sh.finiteElement().localBasis());

      
      auto ip_global = geo.global(ip.position());
      auto ip_local = geo.local(ip_global);

      // evaluate Sg
			RFT Sg0 = 0.0;
			dgf_Sg.evaluate(cell, ip.position(), Sg0);
			RF Sg = Sg0[0];
						
			// evaluate T
			RFT T0 = 0.0;
			dgf_T.evaluate(cell, ip.position(), T0);
			RF T = T0[0];
			
			// evaluate Pw
			RFT Pw0 = 0.0;
			dgf_Pw.evaluate(cell, ip.position(), Pw0);
			RF Pw =Pw0[0];
			
			// evaluate XCH4
			RFT XCH40 = 0.0;
			dgf_XCH4.evaluate(cell, ip.position(), XCH40);
			RF XCH4 = XCH40[0];
			
			// evaluate YH2O
			RFT YH2O0 = 0.0;
			dgf_YH2O.evaluate(cell, ip.position(), YH2O0);
			RF YH2O = YH2O0[0] ;
			
			// evaluate XC
			RFT XC0 = 0.0;
			dgf_XC.evaluate(cell, ip.position(), XC0);
			RF XC = XC0[0];
			
			
			// evaluate Sh
			RF Sh = 0.0;
			for (int i = 0; i < lfsu_Sh.size(); i++)
			  Sh += x(lfsu_Sh, i) * phi_Sh[i];
      
      // evaluate Sw
      RF Sw = 1. - Sg - Sh;

      // evaluate Pg
      auto por = property.soil.SedimentPorosity(cell, ip_local);
      auto Pc = property.hydraulicProperty.CapillaryPressure(cell, ip_local, Sw, Sh, por) ; /* ndim */
      RF Pg = Pw + Pc; /* ndim */
      
      
      auto permeability = property.soil.SedimentPermeability(cell, ip_local )/*ndim K from soil.hh*/
							  * property.hydraulicProperty.PermeabilityScalingFactor(cell, ip_local, Sh, por );

      double S = XC * (property.salt.MolarMass()/property.gas.MolarMass());
      auto zCH4 = property.eos.EvaluateCompressibilityFactor(T * Xc_T, Pg * Xc_P);
      
      // compute source terms
			auto q_g  = property.kinetics.GasGenerationRate( T*Xc_T,
														    Pg*Xc_P,
														    Sh,
														    Sw,
														    XCH4,
														    zCH4,
														    S,
														    por,
														    permeability*Xc_K); /*[kg/m³s]*/
			auto q_h  = property.kinetics.HydrateDissociationRate( q_g ); /*[kg/m³s]*/
			// std::cout << q_g <<"   "<< q_h << "   " << Xc_source_m<<std::endl;
      // exit(0);
      // integrate (A grad u - bu)*grad phi_i + a*u*phi_i
      RF factor = ip.weight() * geo.integrationElement(ip.position());
      
      for (int i = 0; i < lfsv_Sh.size(); i++)
      {
        r.accumulate(lfsv_Sh, i, (-Xc_source_m*q_h * psi_Sh[i]) * factor);
      }
    } //End Quadrature Rule
  }  // End of alpha_volume

};
