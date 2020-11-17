/*
 * LocalOperator.hh
 *
 * 
 */



using namespace Dune::PDELab;


template <typename GV, typename Params, class BC, typename U_Sh, class GFS_Sh,
          typename U, class GFS, typename U_T, class GFS_T,
          class FEM_S>
class LocalOperator_T : 
                      // public Dune::PDELab::NumericalJacobianApplyVolume<LocalOperator_T<GV, Params, BC, U_Sh, GFS_Sh,
                      // U, GFS, U_T, GFS_T, FEM_S>>,
                      // public Dune::PDELab::NumericalJacobianVolume<LocalOperator_T<GV, Params, BC, U_Sh, GFS_Sh,
                      // U, GFS, U_T, GFS_T, FEM_S>>,
                      public Dune::PDELab::NumericalJacobianApplySkeleton<LocalOperator_T<GV, Params, BC, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, FEM_S>>,
                      public Dune::PDELab::NumericalJacobianSkeleton<LocalOperator_T<GV, Params, BC, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, FEM_S>>,  
                      public Dune::PDELab::NumericalJacobianApplyBoundary<LocalOperator_T<GV, Params, BC, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, FEM_S>>,
                      public Dune::PDELab::NumericalJacobianBoundary<LocalOperator_T<GV, Params, BC, U_Sh, GFS_Sh,
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
  //typedef ProblemBoundaryConditions<GV,Params> BC ;

  U_Sh unew_Sh;
  GFS_Sh gfs_Sh;
  U unew;
  GFS gfs;
  U_T *unew_T;
  GFS_T gfs_T;

  double *time;
  double *dt;
  double alpha_t;
  double method_t;

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
  double T_ref;
  Dune::FieldVector<double,GV::dimension> gravity; /* ndim */

public:
  // pattern assembly flags
	  enum { doPatternVolume	= true };
	  enum { doPatternSkeleton	= true };

	  // residual assembly flags
	  enum { doAlphaVolume  	= true };
	  enum { doAlphaSkeleton	= true };
	  enum { doAlphaBoundary	= true };
  

  typedef typename GV::IndexSet IndexSet;

  
  typedef Dune::PDELab::LocalFunctionSpace<GFS_T> LFS;
  typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
  typedef typename U_T::template LocalView<LFSCache> VectorView;
  
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

  using DGF_Sg = typename Dune::PDELab::DiscreteGridFunction<SUBGFS_Sg, U> ;
  using DGF_Sh = typename Dune::PDELab::DiscreteGridFunction<GFS_Sh, U_Sh> ;
	using DGF_Pw = typename Dune::PDELab::DiscreteGridFunction<SUBGFS_Pw, U> ;
	using DGF_XCH4 = typename Dune::PDELab::DiscreteGridFunction<SUBGFS_XCH4, U> ;
	using DGF_YH2O = typename Dune::PDELab::DiscreteGridFunction<SUBGFS_YH2O, U> ;
	using DGF_XC = typename Dune::PDELab::DiscreteGridFunction<SUBGFS_XC, U> ;

 
  using LocalBasisType_T = typename FEM_S::Traits::FiniteElementType::Traits::LocalBasisType;
  using Cache_T = typename Dune::PDELab::LocalBasisCache<LocalBasisType_T>;
  
  using RF = typename LocalBasisType_T::Traits::RangeFieldType;
  using JacobianType = typename LocalBasisType_T::Traits::JacobianType ;
  using RFT = typename Dune::FieldVector<double, 1>;
  
  std::vector<Cache_T> cache_T;
  

  // constructor stores parameters
  LocalOperator_T( const GV &gv_, const Params&	property_, const BC& bc_,
                    const U_Sh &unew_Sh_, GFS_Sh gfs_Sh_, const U &unew_, GFS gfs_,
                    U_T *unew_T_, GFS_T gfs_T_, 
                    double *time_,
                    double *dt_,
                    unsigned int intorder_ = 6, double method_t_=0., double alpha_t_=1.)
      : gv(gv_), property( property_ ), bc( bc_ ),
        unew_Sh(unew_Sh_), gfs_Sh(gfs_Sh_), unew(unew_), gfs(gfs_),
        unew_T(unew_T_), gfs_T(gfs_T_), 
        time(time_),
        dt(dt_), method_t(method_t_),
        alpha_t(alpha_t_),
        intorder(intorder_), cache_T(20)
  {
    theta_T = method_t;
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
    T_ref = property.parameter.ReferenceTemperature()/Xc_T;/* ndim*/
    gravity = -property.parameter.g() / Xc_grav  ; /* ndim */
    #ifdef STATEINDEPENDENTPROPERTIES
      		T_ref = property.parameter.RefT()/Xc_T;
    #endif

    
  }

  // volume integral depending on test and ansatz functions
  template <typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG &eg, const LFSU &lfsu_T, const X &x, const LFSV &lfsv_T, R &r) const
  {
    SUBGFS_Pw gfs_Pw(gfs);
    SUBGFS_Sg gfs_Sg(gfs);
    SUBGFS_XCH4 gfs_XCH4(gfs);
    SUBGFS_YH2O gfs_YH2O(gfs);
    SUBGFS_XC gfs_XC(gfs);
    DGF_Sg dgf_Sg(gfs_Sg, unew);	
    DGF_Sh dgf_Sh(gfs_Sh, unew_Sh);
    DGF_Pw dgf_Pw(gfs_Pw, unew);
    DGF_XCH4 dgf_XCH4(gfs_XCH4, unew);
    DGF_YH2O dgf_YH2O(gfs_YH2O, unew);
    DGF_XC dgf_XC(gfs_XC, unew);
    
    // dimensions
    const int dim = EG::Entity::dimension;
    
    const int order_t = std::max(lfsu_T.finiteElement().localBasis().order(),
                               lfsv_T.finiteElement().localBasis().order());
   
    // Reference to cell
	  const auto& cell = eg.entity();
		const IndexSet &indexSet = gv.indexSet();
		int cell_number = indexSet.index(cell);

    // Get geometry
    auto geo = eg.geometry();

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T(lfsu_T.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T(lfsv_T.size());
    
    Dune::FieldVector<RF, dim> gradu_Pw(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pw(0.0);
    Dune::FieldVector<RF, dim> gradu_Sg(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sg(0.0);
    Dune::FieldVector<RF, dim> gradu_Sh(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sh(0.0);
    Dune::FieldVector<RF, dim> gradu_T(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O(0.0);
    Dune::FieldVector<RF, dim> gradu_XC(0.0);
    Dune::FieldVector<RF, dim> Kg(0.0);
    Dune::FieldVector<RF, dim> delta_x(0.0);
    Dune::FieldVector<RF, dim> delta_y(0.0);
    delta_x[0] = 1.e-3;
    delta_y[1] = 1.e-3; 
   
    // Transformation matrix
    typename EG::Geometry::JacobianInverseTransposed jac;

    // loop over quadrature points
    //      auto intorder = intorderadd + quadrature_factor * order;
    for (const auto &ip : quadratureRule(geo, intorder))
    {
      // evaluate basis functions
      auto &phi_T = cache_T[order_t].evaluateFunction(ip.position(), lfsu_T.finiteElement().localBasis());
      auto &psi_T = cache_T[order_t].evaluateFunction(ip.position(), lfsv_T.finiteElement().localBasis());

      auto qp_x = ip.position() + delta_x;
      auto qp_y = ip.position() + delta_y;
      auto ip_global = geo.global(ip.position());
      auto ip_local = geo.local(ip_global);
      
      
      RF T = 0.0;
      for (int i = 0; i < lfsu_T.size(); i++){
        T += x(lfsu_T, i) * phi_T[i];
      }
      
      // evaluate Sg
      RFT Sg0 = 0.0;
      dgf_Sg.evaluate(cell, ip.position(), Sg0);
      RF Sg = Sg0[0];
      RFT Sg_x0 = 0.0;
      dgf_Sg.evaluate(cell, qp_x, Sg_x0);
      RF Sg_x = Sg_x0[0];
      RFT Sg_y0 = 0.0;
      dgf_Sg.evaluate(cell, qp_y, Sg_y0);
      RF Sg_y = Sg_y0[0];
      
      
      // evaluate Sh
      RFT Sh0 = 0.0;
      dgf_Sh.evaluate(cell, ip.position(), Sh0);
      RF Sh = Sh0[0];
      RFT Sh_x0 = 0.0;
      dgf_Sh.evaluate(cell, qp_x, Sh_x0);
      RF Sh_x = Sh_x0[0];
      RFT Sh_y0 = 0.0;
      dgf_Sh.evaluate(cell, qp_y, Sh_y0);
      RF Sh_y = Sh_y0[0];

      // evaluate T
      RFT Pw0 = 0.0;
      dgf_Pw.evaluate(cell, ip.position(), Pw0);
      RF Pw =Pw0[0];
      RFT Pw_x0 = 0.0;
      dgf_Pw.evaluate(cell, qp_x, Pw_x0);
      RF Pw_x = Pw_x0[0];
      RFT Pw_y0 = 0.0;
      dgf_Pw.evaluate(cell, qp_y, Pw_y0);
      RF Pw_y = Pw_y0[0];
      
      // evaluate XCH4
      RFT XCH40 = 0.0;
      dgf_XCH4.evaluate(cell, ip.position(), XCH40);
      RF XCH4 = XCH40[0];
      RFT XCH4_x0 = 0.0;
      dgf_XCH4.evaluate(cell, qp_x, XCH4_x0);
      RF XCH4_x = XCH4_x0[0];
      RFT XCH4_y0 = 0.0;
      dgf_XCH4.evaluate(cell, qp_y, XCH4_y0);
      RF XCH4_y = XCH4_y0[0];

      // evaluate YH2O
      RFT YH2O0 = 0.0;
      dgf_YH2O.evaluate(cell, ip.position(), YH2O0);
      RF YH2O = YH2O0[0] ;
      RFT YH2O_x0 = 0.0;
      dgf_YH2O.evaluate(cell, qp_x, YH2O_x0);
      RF YH2O_x = YH2O_x0[0] ;
      RFT YH2O_y0 = 0.0;
      dgf_YH2O.evaluate(cell, qp_y, YH2O_y0);
      RF YH2O_y = YH2O_y0[0] ;

      // evaluate XC
      RFT XC0 = 0.0;
      dgf_XC.evaluate(cell, ip.position(), XC0);
      RF XC = XC0[0];
      RFT XC_x0 = 0.0;
      dgf_XC.evaluate(cell, qp_x, XC_x0);
      RF XC_x = XC_x0[0];
      RFT XC_y0 = 0.0;
      dgf_XC.evaluate(cell, qp_y, XC_y0);
      RF XC_y = XC_y0[0];

      // evaluate Sw
      RF Sw = 1. - Sg - Sh;

      // evaluate Pg
      auto BrooksCParams = property.hydraulicProperty.BrooksCoreyParameters(cell, ip_local);/*BrooksCParams[0] gives Pentry in Pa*/
      auto por = property.soil.SedimentPorosity(cell, ip_local);
      auto Pc = property.hydraulicProperty.CapillaryPressure(cell, ip_local, Sw, Sh, por) ; /* ndim */
      RF Pg = Pw + Pc; /* ndim */
      RF Peff = (Pg * Sg + Pw * Sw) / (1. - Sh); /* ndim */
      
      // evaluate gradient of basis functions
      auto &js_T = cache_T[order_t].evaluateJacobian(ip.position(), lfsu_T.finiteElement().localBasis());
      auto &js_v_T = cache_T[order_t].evaluateJacobian(ip.position(), lfsv_T.finiteElement().localBasis());

      // transform gradients of shape functions to real element
      jac = geo.jacobianInverseTransposed(ip.position());

      for (int i = 0; i < lfsu_T.size(); i++)
        jac.mv(js_T[i][0], gradphi_T[i]);
      for (int i = 0; i < lfsv_T.size(); i++)
        jac.mv(js_v_T[i][0], gradpsi_T[i]);

      // compute gradient of Sg
      gradu_Sg[0] = (Sg_x - Sg) / delta_x[0] ;
      gradu_Sg[1] = (Sg_y - Sg) / delta_y[1] ;


      // compute gradient of T
      gradu_T = 0.0;
      for (int i = 0; i < lfsu_T.size(); i++)
        gradu_T.axpy(x(lfsu_T, i), gradphi_T[i]);

      // exit(0);
      // compute gradient of Pw
      gradu_Pw[0] = (Pw_x - Pw) / delta_x[0] ;
      gradu_Pw[1] = (Pw_y - Pw) / delta_y[1] ;

      // compute gradient of Sh
      gradu_Sh[0] = (Sh_x - Sh) / delta_x[0] ;
      gradu_Sh[1] = (Sh_y - Sh) / delta_y[1] ;

      // compute gradient of XCH4
      // gradu_XCH4 = 0.0;
      gradu_XCH4[0] = (XCH4_x - XCH4) / delta_x[0] ;
      gradu_XCH4[1] = (XCH4_y - XCH4) / delta_y[1] ;
 
      // // compute gradient of YH2O
      // gradu_YH2O = 0.0;
      gradu_YH2O[0] = (YH2O_x - YH2O) / delta_x[0] ;
      gradu_YH2O[1] = (YH2O_y - YH2O) / delta_y[1] ;
      

      // // compute gradient of XCH4
      // gradu_XC = 0.0;
      gradu_XC[0] = (XC_x - XC) / delta_x[0] ;
      gradu_XC[1] = (XC_y - XC) / delta_y[1] ;


      auto K = property.soil.SedimentPermeabilityTensor(cell, ip_local)
                    * property.hydraulicProperty.PermeabilityScalingFactor(cell, ip_local, Sh, por ); /*ndim K from soil.hh*/
      K.mv(gravity, Kg);

      // compute K * gradient of Pw
      K.mv(gradu_Pw, Kgradu_Pw);

      // compute K * gradient of Sg
      K.mv(gradu_Sg, Kgradu_Sg);

      // compute K * gradient of Sh
      K.mv(gradu_Sh, Kgradu_Sh);
      
      auto permeability = property.soil.SedimentPermeability(cell, ip_local )/*ndim K from soil.hh*/
							  * property.hydraulicProperty.PermeabilityScalingFactor(cell, ip_local, Sh, por );

      auto tau = property.soil.Tortuosity(por);/*ndim tau from soil.hh*/
      auto DH2O_g = tau * por * property.mixture.DiffCoeffH2OInGas(T * Xc_T, Pg * Xc_P); /*ndim D from mixture.hh*/
      auto DCH4_w = tau * por * property.mixture.DiffCoeffCH4InLiquid(T * Xc_T, Pw * Xc_P); /*ndim D from mixture.hh*/
      auto DC_w = tau * por * property.salt.DiffCoeff(T * Xc_T, Pw * Xc_P); /*ndim D from salt.hh*/
      
      double S = XC * (property.salt.MolarMass()/property.water.MolarMass());
      auto zCH4 = property.eos.EvaluateCompressibilityFactor(T * Xc_T, Pg * Xc_P);
      auto YCH4 =  property.mixture.YCH4(XCH4, T * Xc_T, Pg * Xc_P, XC, zCH4);
      auto XH2O =  property.mixture.XH2O(YH2O, T * Xc_T, Pg * Xc_P, XC);
      
      auto Swe = property.hydraulicProperty.EffectiveSw(Sw,Sh,0.0,0.0);
      auto dPc_dSwe = property.hydraulicProperty.dPc_dSwe(Swe, BrooksCParams[0], BrooksCParams[1]); /*ndim */
      auto dSwe_dSw =  property.hydraulicProperty.dSwe_dSw(Sw,Sh,0.0,0.0);
      auto coeff_grad_Sw = dPc_dSwe * dSwe_dSw ;

      auto dPcSF1_dSh =  property.hydraulicProperty.dPcSF1_dSh( Sh, BrooksCParams[1], BrooksCParams[4]);
      auto dSwe_dSh = property.hydraulicProperty.dSwe_dSh(Sw,Sh,0.0,0.0);
      auto coeff_grad_Sh = dPcSF1_dSh + dPc_dSwe * dSwe_dSh ;

      auto rho_g = property.gas.Density(T * Xc_T, Pg * Xc_P, zCH4); /*ndim density from CH4.hh; the input arguments are dimensional   */
      auto rho_w = property.water.Density(T * Xc_T, Pw * Xc_P, S); /*ndim density from H2O.hh; the input arguments are dimensional*/
      
      auto krW = property.hydraulicProperty.krw(cell, ip_local, Sw, Sh) / (property.water.DynamicViscosity(T * Xc_T, Pw * Xc_P, S) );
      auto krN = property.hydraulicProperty.krg(cell, ip_local, Sw, Sh) / (property.gas.DynamicViscosity(T * Xc_T, Pg * Xc_P));
      
      
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
      auto Q = property.kinetics.HeatOfDissociation( q_g, T*Xc_T ); /*[W/m³]*/
      
      auto Cp_g = property.gas.Cp(T * Xc_T, Pg * Xc_P, zCH4); /* ndim */
      auto Cp_w = property.water.Cp(T * Xc_T, Pw * Xc_P, S); /* ndim */
      auto kth_g = property.gas.ThermalConductivity(T * Xc_T, Pg * Xc_P); /* ndim */
      auto kth_w = property.water.ThermalConductivity(T * Xc_T, Pw * Xc_P, S); /* ndim */
      auto kth_h = property.hydrate.ThermalConductivity(T * Xc_T, Peff * Xc_P); /* ndim */
      auto kth_s = property.soil.ThermalConductivity(); /* ndim */
      auto kth_eff = (1. - por) * kth_s + por * (Sg * kth_g + Sw * kth_w + Sh * kth_h); /* ndim */
      
      auto gradu_Pg = gradu_Pw  - coeff_grad_Sw * gradu_Sg + (coeff_grad_Sh - coeff_grad_Sw) * gradu_Sh;
      auto Kgradu_Pg = Kgradu_Pw - coeff_grad_Sw * Kgradu_Sg + (coeff_grad_Sh - coeff_grad_Sw) * Kgradu_Sh;

      auto convectiveflux_Heat_w = rho_w * Cp_w * (T - T_ref) * krW * (Kgradu_Pw - rho_w * Kg);
      auto convectiveflux_Heat_g = rho_g * Cp_g * (T - T_ref) * krN * (Kgradu_Pg - rho_g * Kg);


      auto convectiveflux_Heat = convectiveflux_Heat_g + convectiveflux_Heat_w;

      auto diffusiveflux_Heat = kth_eff * gradu_T;

      // integrate (A grad u - bu)*grad phi_i + a*u*phi_i
      RF factor = ip.weight() * geo.integrationElement(ip.position());
      for (int i = 0; i < lfsv_T.size(); i++)
      {
        r.accumulate(lfsv_T, i, ((Xc_conv_h * convectiveflux_Heat + Xc_diff_h * diffusiveflux_Heat ) * gradpsi_T[i] 
                                  - Xc_source_h * Q * psi_T[i]) * factor);
      }
    } //End Quadrature Rule
  }  // End of alpha_volume
  
  //volume integral depending on test and ansatz functions
  template <typename EG, typename LFSU, typename X, typename LFSV, typename M>
  void jacobian_volume(const EG &eg, const LFSU &lfsu_T, const X &x, const LFSV &lfsv_T, M &mat) const
  {

    SUBGFS_Pw gfs_Pw(gfs);
    SUBGFS_Sg gfs_Sg(gfs);
    SUBGFS_XCH4 gfs_XCH4(gfs);
    SUBGFS_YH2O gfs_YH2O(gfs);
    SUBGFS_XC gfs_XC(gfs);
    DGF_Sg dgf_Sg(gfs_Sg, unew);	
    DGF_Sh dgf_Sh(gfs_Sh, unew_Sh);
    DGF_Pw dgf_Pw(gfs_Pw, unew);
    DGF_XCH4 dgf_XCH4(gfs_XCH4, unew);
    DGF_YH2O dgf_YH2O(gfs_YH2O, unew);
    DGF_XC dgf_XC(gfs_XC, unew);

    // dimensions
    const int dim = EG::Entity::dimension;
    
    const int order_t = std::max(lfsu_T.finiteElement().localBasis().order(),
                               lfsv_T.finiteElement().localBasis().order());
   
    // Reference to cell
	  const auto& cell = eg.entity();
		const IndexSet &indexSet = gv.indexSet();
		int cell_number = indexSet.index(cell);

    // Get geometry
    auto geo = eg.geometry();

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T(lfsu_T.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T(lfsv_T.size());
    
    Dune::FieldVector<RF, dim> gradu_Pw(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pw(0.0);
    Dune::FieldVector<RF, dim> gradu_Sg(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sg(0.0);
    Dune::FieldVector<RF, dim> gradu_Sh(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sh(0.0);
    Dune::FieldVector<RF, dim> gradu_T(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O(0.0);
    Dune::FieldVector<RF, dim> gradu_XC(0.0);
    Dune::FieldVector<RF, dim> Kg(0.0);
    Dune::FieldVector<RF, dim> delta_x(0.0);
    Dune::FieldVector<RF, dim> delta_y(0.0);
    delta_x[0] = 1.e-3;
    delta_y[1] = 1.e-3; 
   
    // Transformation matrix
    typename EG::Geometry::JacobianInverseTransposed jac;

    // loop over quadrature points
    //      auto intorder = intorderadd + quadrature_factor * order;
    for (const auto &ip : quadratureRule(geo, intorder))
    {
      // evaluate basis functions
      auto &phi_T = cache_T[order_t].evaluateFunction(ip.position(), lfsu_T.finiteElement().localBasis());
      auto &psi_T = cache_T[order_t].evaluateFunction(ip.position(), lfsv_T.finiteElement().localBasis());

      auto qp_x = ip.position() + delta_x;
      auto qp_y = ip.position() + delta_y;
      auto ip_global = geo.global(ip.position());
      auto ip_local = geo.local(ip_global);
      
      
      RF T = 0.0;
      for (int i = 0; i < lfsu_T.size(); i++){
        T += x(lfsu_T, i) * phi_T[i];
      }
      RF T_x = 0.0;
      for (int i = 0; i < lfsu_T.size(); i++){
        T_x += (x(lfsu_T, i)+delta_x[0]) * phi_T[i];
      }
      
      // evaluate Sg
      RFT Sg0 = 0.0;
      dgf_Sg.evaluate(cell, ip.position(), Sg0);
      RF Sg = Sg0[0];
      RFT Sg_x0 = 0.0;
      dgf_Sg.evaluate(cell, qp_x, Sg_x0);
      RF Sg_x = Sg_x0[0];
      RFT Sg_y0 = 0.0;
      dgf_Sg.evaluate(cell, qp_y, Sg_y0);
      RF Sg_y = Sg_y0[0];
      
      
      // evaluate Sh
      RFT Sh0 = 0.0;
      dgf_Sh.evaluate(cell, ip.position(), Sh0);
      RF Sh = Sh0[0];
      RFT Sh_x0 = 0.0;
      dgf_Sh.evaluate(cell, qp_x, Sh_x0);
      RF Sh_x = Sh_x0[0];
      RFT Sh_y0 = 0.0;
      dgf_Sh.evaluate(cell, qp_y, Sh_y0);
      RF Sh_y = Sh_y0[0];

      // evaluate T
      RFT Pw0 = 0.0;
      dgf_Pw.evaluate(cell, ip.position(), Pw0);
      RF Pw =Pw0[0];
      RFT Pw_x0 = 0.0;
      dgf_Pw.evaluate(cell, qp_x, Pw_x0);
      RF Pw_x = Pw_x0[0];
      RFT Pw_y0 = 0.0;
      dgf_Pw.evaluate(cell, qp_y, Pw_y0);
      RF Pw_y = Pw_y0[0];
      
      // evaluate XCH4
      RFT XCH40 = 0.0;
      dgf_XCH4.evaluate(cell, ip.position(), XCH40);
      RF XCH4 = XCH40[0];
      RFT XCH4_x0 = 0.0;
      dgf_XCH4.evaluate(cell, qp_x, XCH4_x0);
      RF XCH4_x = XCH4_x0[0];
      RFT XCH4_y0 = 0.0;
      dgf_XCH4.evaluate(cell, qp_y, XCH4_y0);
      RF XCH4_y = XCH4_y0[0];

      // evaluate YH2O
      RFT YH2O0 = 0.0;
      dgf_YH2O.evaluate(cell, ip.position(), YH2O0);
      RF YH2O = YH2O0[0] ;
      RFT YH2O_x0 = 0.0;
      dgf_YH2O.evaluate(cell, qp_x, YH2O_x0);
      RF YH2O_x = YH2O_x0[0] ;
      RFT YH2O_y0 = 0.0;
      dgf_YH2O.evaluate(cell, qp_y, YH2O_y0);
      RF YH2O_y = YH2O_y0[0] ;

      // evaluate XC
      RFT XC0 = 0.0;
      dgf_XC.evaluate(cell, ip.position(), XC0);
      RF XC = XC0[0];
      RFT XC_x0 = 0.0;
      dgf_XC.evaluate(cell, qp_x, XC_x0);
      RF XC_x = XC_x0[0];
      RFT XC_y0 = 0.0;
      dgf_XC.evaluate(cell, qp_y, XC_y0);
      RF XC_y = XC_y0[0];

      RF Sw = 1. - Sg - Sh;
      // evaluate Pg
      auto BrooksCParams = property.hydraulicProperty.BrooksCoreyParameters(cell, ip_local);/*BrooksCParams[0] gives Pentry in Pa*/
      auto por = property.soil.SedimentPorosity(cell, ip_local);
      auto Pc = property.hydraulicProperty.CapillaryPressure(cell, ip_local, Sw, Sh, por) ; /* ndim */
      RF Pg = Pw + Pc; /* ndim */
      RF Peff = (Pg * Sg + Pw * Sw) / (1. - Sh); /* ndim */
      
      // evaluate gradient of basis functions
      auto &js_T = cache_T[order_t].evaluateJacobian(ip.position(), lfsu_T.finiteElement().localBasis());
      auto &js_v_T = cache_T[order_t].evaluateJacobian(ip.position(), lfsv_T.finiteElement().localBasis());

      // transform gradients of shape functions to real element
      jac = geo.jacobianInverseTransposed(ip.position());

      for (int i = 0; i < lfsu_T.size(); i++)
        jac.mv(js_T[i][0], gradphi_T[i]);
      for (int i = 0; i < lfsv_T.size(); i++)
        jac.mv(js_v_T[i][0], gradpsi_T[i]);

      // compute gradient of Sg
      gradu_Sg[0] = (Sg_x - Sg) / delta_x[0] ;
      gradu_Sg[1] = (Sg_y - Sg) / delta_y[1] ;


      // compute gradient of T
      gradu_T = 0.0;
      for (int i = 0; i < lfsu_T.size(); i++)
        gradu_T.axpy(x(lfsu_T, i), gradphi_T[i]);

      // std::cout << gradu_Pw << "   -- " << qp_x << "  " << qp_y <<  std::endl;
      // std::cout << "ip.position() = " << ip.position()  << "   Pw_y = " << Pw_y<< std::endl;

      // exit(0);
      // compute gradient of Pw
      gradu_Pw[0] = (Pw_x - Pw) / delta_x[0] ;
      gradu_Pw[1] = (Pw_y - Pw) / delta_y[1] ;

      // compute gradient of Sh
      gradu_Sh[0] = (Sh_x - Sh) / delta_x[0] ;
      gradu_Sh[1] = (Sh_y - Sh) / delta_y[1] ;

      // compute gradient of XCH4
      // gradu_XCH4 = 0.0;
      gradu_XCH4[0] = (XCH4_x - XCH4) / delta_x[0] ;
      gradu_XCH4[1] = (XCH4_y - XCH4) / delta_y[1] ;
 
      // // compute gradient of YH2O
      // gradu_YH2O = 0.0;
      gradu_YH2O[0] = (YH2O_x - YH2O) / delta_x[0] ;
      gradu_YH2O[1] = (YH2O_y - YH2O) / delta_y[1] ;
      

      // // compute gradient of XCH4
      // gradu_XC = 0.0;
      gradu_XC[0] = (XC_x - XC) / delta_x[0] ;
      gradu_XC[1] = (XC_y - XC) / delta_y[1] ;


      auto K = property.soil.SedimentPermeabilityTensor(cell, ip_local)
                    * property.hydraulicProperty.PermeabilityScalingFactor(cell, ip_local, Sh, por ); /*ndim K from soil.hh*/
      K.mv(gravity, Kg);

      // compute K * gradient of Pw
      K.mv(gradu_Pw, Kgradu_Pw);

      // compute K * gradient of Sg
      K.mv(gradu_Sg, Kgradu_Sg);

      // compute K * gradient of Sh
      K.mv(gradu_Sh, Kgradu_Sh);
      
      auto permeability = property.soil.SedimentPermeability(cell, ip_local )/*ndim K from soil.hh*/
							  * property.hydraulicProperty.PermeabilityScalingFactor(cell, ip_local, Sh, por );

      auto tau = property.soil.Tortuosity(por);/*ndim tau from soil.hh*/
      auto DH2O_g = tau * por * property.mixture.DiffCoeffH2OInGas(T * Xc_T, Pg * Xc_P); /*ndim D from mixture.hh*/
      auto DCH4_w = tau * por * property.mixture.DiffCoeffCH4InLiquid(T * Xc_T, Pw * Xc_P); /*ndim D from mixture.hh*/
      auto DC_w = tau * por * property.salt.DiffCoeff(T * Xc_T, Pw * Xc_P); /*ndim D from salt.hh*/
      
      double S = XC * (property.salt.MolarMass()/property.water.MolarMass());
      auto zCH4 = property.eos.EvaluateCompressibilityFactor(T * Xc_T, Pg * Xc_P);
      auto YCH4 =  property.mixture.YCH4(XCH4, T * Xc_T, Pg * Xc_P, XC, zCH4);
      auto XH2O =  property.mixture.XH2O(YH2O, T * Xc_T, Pg * Xc_P, XC);
      
      auto Swe = property.hydraulicProperty.EffectiveSw(Sw,Sh,0.0,0.0);
      auto dPc_dSwe = property.hydraulicProperty.dPc_dSwe(Swe, BrooksCParams[0], BrooksCParams[1]); /*ndim */
      auto dSwe_dSw =  property.hydraulicProperty.dSwe_dSw(Sw,Sh,0.0,0.0);
      auto coeff_grad_Sw = dPc_dSwe * dSwe_dSw ;

      auto dPcSF1_dSh =  property.hydraulicProperty.dPcSF1_dSh( Sh, BrooksCParams[1], BrooksCParams[4]);
      auto dSwe_dSh = property.hydraulicProperty.dSwe_dSh(Sw,Sh,0.0,0.0);
      auto coeff_grad_Sh = dPcSF1_dSh + dPc_dSwe * dSwe_dSh ;

      auto rho_g = property.gas.Density(T * Xc_T, Pg * Xc_P, zCH4); /*ndim density from CH4.hh; the input arguments are dimensional   */
      auto rho_w = property.water.Density(T * Xc_T, Pw * Xc_P, S); /*ndim density from H2O.hh; the input arguments are dimensional*/
      
      auto krW = property.hydraulicProperty.krw(cell, ip_local, Sw, Sh) / (property.water.DynamicViscosity(T * Xc_T, Pw * Xc_P, S) );
      auto krN = property.hydraulicProperty.krg(cell, ip_local, Sw, Sh) / (property.gas.DynamicViscosity(T * Xc_T, Pg * Xc_P));
      //auto P_eq = property.kinetics.EquilibriumPressure(T * Xc_T, S)/Xc_P;
      
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
      auto Q = property.kinetics.HeatOfDissociation( q_g, T*Xc_T ); /*[W/m³]*/
      // compute source terms
			auto q_g_x  = property.kinetics.GasGenerationRate( T_x*Xc_T,
														    Pg*Xc_P,
														    Sh,
														    Sw,
														    XCH4,
														    zCH4,
														    S,
														    por,
														    permeability*Xc_K); /*[kg/m³s]*/
      auto Q_x = property.kinetics.HeatOfDissociation( q_g_x, T_x*Xc_T ); /*[W/m³]*/
      
      auto Cp_g = property.gas.Cp(T * Xc_T, Pg * Xc_P, zCH4); /* ndim */
      auto Cp_w = property.water.Cp(T * Xc_T, Pw * Xc_P, S); /* ndim */
      auto kth_g = property.gas.ThermalConductivity(T * Xc_T, Pg * Xc_P); /* ndim */
      auto kth_w = property.water.ThermalConductivity(T * Xc_T, Pw * Xc_P, S); /* ndim */
      auto kth_h = property.hydrate.ThermalConductivity(T * Xc_T, Peff * Xc_P); /* ndim */
      auto kth_s = property.soil.ThermalConductivity(); /* ndim */
      auto kth_eff = (1. - por) * kth_s + por * (Sg * kth_g + Sw * kth_w + Sh * kth_h); /* ndim */
      
      auto gradu_Pg = gradu_Pw  - coeff_grad_Sw * gradu_Sg + (coeff_grad_Sh - coeff_grad_Sw) * gradu_Sh;
      auto Kgradu_Pg = Kgradu_Pw - coeff_grad_Sw * Kgradu_Sg + (coeff_grad_Sh - coeff_grad_Sw) * Kgradu_Sh;

      
      // integrate (A grad u - bu)*grad phi_i + a*u*phi_i
      RF factor = ip.weight() * geo.integrationElement(ip.position());
      
      
      for (int i = 0; i < lfsv_T.size(); i++)
      {
        for (int j = 0; j < lfsu_T.size(); j++)
        {
          mat.accumulate(lfsv_T, i, lfsu_T, j, Xc_conv_h * (rho_w * Cp_w * krW * (Kgradu_Pw - rho_w * Kg) + rho_g * Cp_g * krN * (Kgradu_Pg - rho_g * Kg)) * double(phi_T[j]) * gradpsi_T[i] * factor);
        }
        for (int j = 0; j < lfsu_T.size(); j++)
        {
          mat.accumulate(lfsv_T, i, lfsu_T, j, Xc_diff_h * kth_eff * gradphi_T[j] * gradpsi_T[i] * factor);
        }
      }
    } //End Quadrature Rule
  }  // End of jacobian_volume

  // //! apply local jacobian of the volume term -> nonlinear variant
  // template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  // void jacobian_apply_volume (const EG& eg, const LFSU& lfsu_T,
  //                             const X& x, const X& z, const LFSV& lfsv_T,
  //                             R& r) const
  // {

  //   alpha_volume(eg,lfsu_T,z,lfsv_T,r);
  // }

  // skeleton integral depending on test and ansatz functions
  // each face is only visited ONCE!
  template <typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_skeleton(const IG &ig,
                      const LFSU &lfsu_T_s, const X &x_s, const LFSV &lfsv_T_s,
                      const LFSU &lfsu_T_n, const X &x_n, const LFSV &lfsv_T_n,
                      R &r_s, R &r_n) const
  {
   
    SUBGFS_Pw gfs_Pw(gfs);
    SUBGFS_Sg gfs_Sg(gfs);
    SUBGFS_XCH4 gfs_XCH4(gfs);
    SUBGFS_YH2O gfs_YH2O(gfs);
    SUBGFS_XC gfs_XC(gfs);
    DGF_Sg dgf_Sg(gfs_Sg, unew);	
    DGF_Sh dgf_Sh(gfs_Sh, unew_Sh);
    DGF_Pw dgf_Pw(gfs_Pw, unew);
    DGF_XCH4 dgf_XCH4(gfs_XCH4, unew);
    DGF_YH2O dgf_YH2O(gfs_YH2O, unew);
    DGF_XC dgf_XC(gfs_XC, unew);

    // dimensions
    const int dim= IG::Entity::dimension;
    const int dimension = GV::dimension;

    const int order_t = std::max(lfsu_T_s.finiteElement().localBasis().order(),
                                lfsv_T_s.finiteElement().localBasis().order());


    // References to inside and outside cells
    const auto &cell_inside = ig.inside();
    const auto &cell_outside = ig.outside();
    
    // Get geometries
    auto geo = ig.geometry();
    //const auto dimension = geo.mydimension;
    auto geo_inside = cell_inside.geometry();
    auto geo_outside = cell_outside.geometry();

    // Get geometry of intersection in local coordinates of cell_inside and cell_outside
    auto geo_in_inside = ig.geometryInInside();
    auto geo_in_outside = ig.geometryInOutside();
    // cell geometries
    auto ref_el_inside 	= referenceElement(geo_inside);
	  auto ref_el_outside = referenceElement(geo_outside);
	  auto inside_cell_center_local 	= ref_el_inside.position(0,0);
	  auto outside_cell_center_local 	= ref_el_outside.position(0,0);
    // face diameter; this should be revised for anisotropic meshes?
    auto h_F = std::min(geo_inside.volume(), geo_outside.volume()) / geo.volume(); // Houston!

    // compute weights
    RF omega_s;
    RF omega_n;
    RF harmonic_average(0.0);
    omega_s = omega_n = 0.5;
    harmonic_average = 1.0;

    //      if (weights==ConvectionDiffusionDGWeights::weightsOn)
    //        {
    //          RF delta_s = (An_F_s*n_F);
    //          RF delta_n = (An_F_n*n_F);
    //          omega_s = delta_n/(delta_s+delta_n+1e-20);
    //          omega_n = delta_s/(delta_s+delta_n+1e-20);
    //          harmonic_average = 2.0*delta_s*delta_n/(delta_s+delta_n+1e-20);
    //        }
    //      else
    //        {
    //          omega_s = omega_n = 0.5;
    //          harmonic_average = 1.0;
    //        }

    // get polynomial degree
    auto order_i = lfsv_T_s.finiteElement().localBasis().order();
    auto order_o = lfsv_T_n.finiteElement().localBasis().order();
    auto degree = std::max(order_i, order_o);

    // penalty factor
    auto penalty_factor_t = (alpha_t / h_F) * harmonic_average * degree * (degree + dim - 1);

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T_s(lfsu_T_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T_s(lfsv_T_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T_n(lfsu_T_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T_n(lfsv_T_n.size());

    Dune::FieldVector<RF, dim> gradu_Pw_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pw_s(0.0);
    Dune::FieldVector<RF, dim> gradu_Sg_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sg_s(0.0);
    Dune::FieldVector<RF, dim> gradu_Sh_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sh_s(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4_s(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O_s(0.0);
    Dune::FieldVector<RF, dim> gradu_XC_s(0.0);
    Dune::FieldVector<RF, dim> gradu_T_s(0.0);

    Dune::FieldVector<RF, dim> gradu_Pw_n(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pw_n(0.0);
    Dune::FieldVector<RF, dim> gradu_Sg_n(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sg_n(0.0);
    Dune::FieldVector<RF, dim> gradu_Sh_n(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sh_n(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4_n(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O_n(0.0);
    Dune::FieldVector<RF, dim> gradu_XC_n(0.0);
    Dune::FieldVector<RF, dim> gradu_T_n(0.0);

    Dune::FieldVector<RF, dim> v_g(0.0);
    Dune::FieldVector<RF, dim> v_w(0.0);

    Dune::FieldVector<RF, dim> Kg_s(0.0);
    Dune::FieldVector<RF, dim> Kg_n(0.0);


    Dune::FieldVector<RF, dim> delta_x(0.0);
    Dune::FieldVector<RF, dim> delta_y(0.0);
    delta_x[0] = 1.e-3;
    delta_y[1] = 1.e-3; 

    // Transformation matrix
    typename IG::Entity::Geometry::JacobianInverseTransposed jac;

    // loop over quadrature points
    // auto intorder = intorderadd+quadrature_factor*order;
    for (const auto &ip : quadratureRule(geo, intorder))
    {
      // exact normal
      auto n_F_local = ig.unitOuterNormal(ip.position());

      // position of quadrature point in local coordinates of elements
      auto iplocal_s = geo_in_inside.global(ip.position());
      auto iplocal_n = geo_in_outside.global(ip.position());
      auto qp_x_s = iplocal_s + delta_x;
      auto qp_y_s = iplocal_s + delta_y;
      auto qp_x_n = iplocal_n + delta_x;
      auto qp_y_n = iplocal_n + delta_y;

      auto ip_global_s = geo_inside.global(iplocal_s);
      auto ip_global_n = geo_outside.global(iplocal_n);

      // evaluate basis functions
      auto &phi_T_s = cache_T[order_t].evaluateFunction(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &psi_T_s = cache_T[order_t].evaluateFunction(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &phi_T_n = cache_T[order_t].evaluateFunction(iplocal_n, lfsu_T_n.finiteElement().localBasis());
      auto &psi_T_n = cache_T[order_t].evaluateFunction(iplocal_n, lfsv_T_n.finiteElement().localBasis());
      
      
      // evaluate Sg
      RFT Sg_s0 = 0.0;
      dgf_Sg.evaluate(cell_inside, iplocal_s, Sg_s0);
      RF Sg_s = Sg_s0[0];
      RFT Sg_x_s0 = 0.0;
      dgf_Sg.evaluate(cell_inside, qp_x_s, Sg_x_s0);
      RF Sg_x_s = Sg_x_s0[0];
      RFT Sg_y_s0 = 0.0;
      dgf_Sg.evaluate(cell_inside, qp_y_s, Sg_y_s0);
      RF Sg_y_s = Sg_y_s0[0];
      
      RFT Sg_n0 = 0.0;
      dgf_Sg.evaluate(cell_outside, iplocal_n, Sg_n0);
      RF Sg_n = Sg_n0[0];
      RFT Sg_x_n0 = 0.0;
      dgf_Sg.evaluate(cell_outside, qp_x_n, Sg_x_n0);
      RF Sg_x_n = Sg_x_n0[0];
      RFT Sg_y_n0 = 0.0;
      dgf_Sg.evaluate(cell_outside, qp_y_n, Sg_y_n0);
      RF Sg_y_n = Sg_y_n0[0];

      // evaluate T
      RF T_s = 0.0;
      for (int i = 0; i < lfsu_T_s.size(); i++)
        T_s += x_s(lfsu_T_s, i) * phi_T_s[i];
      RF T_n = 0.0;
      for (int i = 0; i < lfsu_T_n.size(); i++)
        T_n += x_n(lfsu_T_n, i) * phi_T_n[i];

      // evaluate Sh
      RFT Sh_s0 = 0.0;
      dgf_Sh.evaluate(cell_inside, iplocal_s, Sh_s0);
      RF Sh_s = Sh_s0[0];
      RFT Sh_x_s0 = 0.0;
      dgf_Sh.evaluate(cell_inside, qp_x_s, Sh_x_s0);
      RF Sh_x_s = Sh_x_s0[0];
      RFT Sh_y_s0 = 0.0;
      dgf_Sh.evaluate(cell_inside, qp_y_s, Sh_y_s0);
      RF Sh_y_s = Sh_y_s0[0];
      
      RFT Sh_n0 = 0.0;
      dgf_Sh.evaluate(cell_outside, iplocal_n, Sh_n0);
      RF Sh_n = Sh_n0[0];
      RFT Sh_x_n0 = 0.0;
      dgf_Sh.evaluate(cell_outside, qp_x_n, Sh_x_n0);
      RF Sh_x_n = Sh_x_n0[0];
      RFT Sh_y_n0 = 0.0;
      dgf_Sh.evaluate(cell_outside, qp_y_n, Sh_y_n0);
      RF Sh_y_n = Sh_y_n0[0];

      // evaluate Pw
      RFT Pw_s0 = 0.0;
      dgf_Pw.evaluate(cell_inside, iplocal_s, Pw_s0);
      RF Pw_s = Pw_s0[0];
      RFT Pw_x_s0 = 0.0;
      dgf_Pw.evaluate(cell_inside, qp_x_s, Pw_x_s0);
      RF Pw_x_s = Pw_x_s0[0];
      RFT Pw_y_s0 = 0.0;
      dgf_Pw.evaluate(cell_inside, qp_y_s, Pw_y_s0);
      RF Pw_y_s = Pw_y_s0[0];
      
      RFT Pw_n0 = 0.0;
      dgf_Pw.evaluate(cell_outside, iplocal_n, Pw_n0);
      RF Pw_n = Pw_n0[0];
      RFT Pw_x_n0 = 0.0;
      dgf_Pw.evaluate(cell_outside, qp_x_n, Pw_x_n0);
      RF Pw_x_n = Pw_x_n0[0];
      RFT Pw_y_n0 = 0.0;
      dgf_Pw.evaluate(cell_outside, qp_y_n, Pw_y_n0);
      RF Pw_y_n = Pw_y_n0[0];
      

      // evaluate XCH4
      RFT XCH4_s0 = 0.0;
      dgf_XCH4.evaluate(cell_inside, iplocal_s, XCH4_s0);
      RF XCH4_s = XCH4_s0[0];
      RFT XCH4_x_s0 = 0.0;
      dgf_XCH4.evaluate(cell_inside, qp_x_s, XCH4_x_s0);
      RF XCH4_x_s = XCH4_x_s0[0];
      RFT XCH4_y_s0 = 0.0;
      dgf_XCH4.evaluate(cell_inside, qp_y_s, XCH4_y_s0);
      RF XCH4_y_s = XCH4_y_s0[0];
      
      RFT XCH4_n0 = 0.0;
      dgf_XCH4.evaluate(cell_outside, iplocal_n, XCH4_n0);
      RF XCH4_n = XCH4_n0[0];
      RFT XCH4_x_n0 = 0.0;
      dgf_XCH4.evaluate(cell_outside, qp_x_n, XCH4_x_n0);
      RF XCH4_x_n = XCH4_x_n0[0];
      RFT XCH4_y_n0 = 0.0;
      dgf_XCH4.evaluate(cell_outside, qp_y_n, XCH4_y_n0);
      RF XCH4_y_n = XCH4_y_n0[0];

      // evaluate YH2O
      RFT YH2O_s0 = 0.0;
      dgf_YH2O.evaluate(cell_inside, iplocal_s, YH2O_s0);
      RF YH2O_s = YH2O_s0[0];
      RFT YH2O_x_s0 = 0.0;
      dgf_YH2O.evaluate(cell_inside, qp_x_s, YH2O_x_s0);
      RF YH2O_x_s = YH2O_x_s0[0];
      RFT YH2O_y_s0 = 0.0;
      dgf_YH2O.evaluate(cell_inside, qp_y_s, YH2O_y_s0);
      RF YH2O_y_s = YH2O_y_s0[0];
      
      RFT YH2O_n0 = 0.0;
      dgf_YH2O.evaluate(cell_outside, iplocal_n, YH2O_n0);
      RF YH2O_n = YH2O_n0[0];
      RFT YH2O_x_n0 = 0.0;
      dgf_YH2O.evaluate(cell_outside, qp_x_n, YH2O_x_n0);
      RF YH2O_x_n = YH2O_x_n0[0];
      RFT YH2O_y_n0 = 0.0;
      dgf_YH2O.evaluate(cell_outside, qp_y_n, YH2O_y_n0);
      RF YH2O_y_n = YH2O_y_n0[0];

      // evaluate XC
      RFT XC_s0 = 0.0;
      dgf_XC.evaluate(cell_inside, iplocal_s, XC_s0);
      RF XC_s = XC_s0[0];
      RFT XC_x_s0 = 0.0;
      dgf_XC.evaluate(cell_inside, qp_x_s, XC_x_s0);
      RF XC_x_s = XC_x_s0[0];
      RFT XC_y_s0 = 0.0;
      dgf_XC.evaluate(cell_inside, qp_y_s, XC_y_s0);
      RF XC_y_s = XC_y_s0[0];
      
      RFT XC_n0 = 0.0;
      dgf_XC.evaluate(cell_outside, iplocal_n, XC_n0);
      RF XC_n = XC_n0[0];
      RFT XC_x_n0 = 0.0;
      dgf_XC.evaluate(cell_outside, qp_x_n, XC_x_n0);
      RF XC_x_n = XC_x_n0[0];
      RFT XC_y_n0 = 0.0;
      dgf_XC.evaluate(cell_outside, qp_y_n, XC_y_n0);
      RF XC_y_n = XC_y_n0[0];


      //
      RF Sw_s = 1. - Sg_s - Sh_s;
      RF Sw_n = 1. - Sg_n - Sh_n;

      // evaluate Pg
      auto BrooksCParams_s = property.hydraulicProperty.BrooksCoreyParameters(cell_inside, iplocal_s);/*BrooksCParams[0] gives Pentry in Pa*/
      auto BrooksCParams_n = property.hydraulicProperty.BrooksCoreyParameters(cell_outside, iplocal_n);/*BrooksCParams[0] gives Pentry in Pa*/
      auto por_s = property.soil.SedimentPorosity(cell_inside, iplocal_s);
      auto por_n = property.soil.SedimentPorosity(cell_outside, iplocal_n);
      auto Pc_s = property.hydraulicProperty.CapillaryPressure(cell_inside, iplocal_s, Sw_s, Sh_s, por_s) ; /* ndim */
      RF Pg_s = Pw_s + Pc_s;
      
      auto Pc_n = property.hydraulicProperty.CapillaryPressure(cell_outside, iplocal_n, Sw_n, Sh_n, por_n) ; /* ndim */
      RF Pg_n = Pw_n + Pc_n;/* ndim */

      RF Peff_s = (Pg_s * Sg_s + Pw_s * Sw_s) / (1. - Sh_s);
      RF Peff_n = (Pg_n * Sg_n + Pw_n * Sw_n) / (1. - Sh_n);

      // evaluate gradient of basis functions
      auto &js_T_s = cache_T[order_t].evaluateJacobian(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &js_v_T_s = cache_T[order_t].evaluateJacobian(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &js_T_n = cache_T[order_t].evaluateJacobian(iplocal_n, lfsu_T_n.finiteElement().localBasis());
      auto &js_v_T_n = cache_T[order_t].evaluateJacobian(iplocal_n, lfsv_T_n.finiteElement().localBasis());

      // transform gradients of shape functions to real element
      jac = geo_inside.jacobianInverseTransposed(iplocal_s);

      for (int i = 0; i < lfsu_T_s.size(); i++)
        jac.mv(js_T_s[i][0], gradphi_T_s[i]);
      for (int i = 0; i < lfsv_T_s.size(); i++)
        jac.mv(js_v_T_s[i][0], gradpsi_T_s[i]);

      jac = geo_outside.jacobianInverseTransposed(iplocal_n);

      for (int i = 0; i < lfsu_T_n.size(); i++)
        jac.mv(js_T_n[i][0], gradphi_T_n[i]);
      for (int i = 0; i < lfsv_T_n.size(); i++)
        jac.mv(js_v_T_n[i][0], gradpsi_T_n[i]);

      // compute gradient of Pw
      gradu_Sg_s = 0.0;
      gradu_Sg_s[0] = (Sg_x_s - Sg_s) / delta_x[0] ;
      gradu_Sg_s[1] = (Sg_y_s - Sg_s) / delta_y[1] ;
      gradu_Sg_n = 0.0;
      gradu_Sg_n[0] = (Sg_x_n - Sg_n) / delta_x[0] ;
      gradu_Sg_n[1] = (Sg_y_n - Sg_n) / delta_y[1] ;

      // compute gradient of Sg
      gradu_T_s = 0.0;
      for (int i = 0; i < lfsu_T_s.size(); i++)
        gradu_T_s.axpy(x_s(lfsu_T_s, i), gradphi_T_s[i]);
      gradu_T_n = 0.0;
      for (int i = 0; i < lfsu_T_n.size(); i++)
        gradu_T_n.axpy(x_n(lfsu_T_n, i), gradphi_T_n[i]);

      // compute gradient of Sh
      gradu_Sh_s = 0.0;
      gradu_Sh_s[0] = (Sh_x_s - Sh_s) / delta_x[0] ;
      gradu_Sh_s[1] = (Sh_y_s - Sh_s) / delta_y[1] ;
      gradu_Sh_n = 0.0;
      gradu_Sh_n[0] = (Sh_x_n - Sh_n) / delta_x[0] ;
      gradu_Sh_n[1] = (Sh_y_n - Sh_n) / delta_y[1] ;

      // compute gradient of XCH4
      gradu_XCH4_s = 0.0;
      gradu_XCH4_s[0] = (XCH4_x_s - XCH4_s) / delta_x[0] ;
      gradu_XCH4_s[1] = (XCH4_y_s - XCH4_s) / delta_y[1] ;
      gradu_XCH4_n = 0.0;
      gradu_XCH4_n[0] = (XCH4_x_n - XCH4_n) / delta_x[0] ;
      gradu_XCH4_n[1] = (XCH4_y_n - XCH4_n) / delta_y[1] ;

      // compute gradient of YH2O
      gradu_YH2O_s = 0.0;
      gradu_YH2O_s[0] = (YH2O_x_s - YH2O_s) / delta_x[0] ;
      gradu_YH2O_s[1] = (YH2O_y_s - YH2O_s) / delta_y[1] ;
      gradu_YH2O_n = 0.0;
      gradu_YH2O_n[0] = (YH2O_x_n - YH2O_n) / delta_x[0] ;
      gradu_YH2O_n[1] = (YH2O_y_n - YH2O_n) / delta_y[1] ;

      // compute gradient of Pw
      gradu_Pw_s = 0.0;
      gradu_Pw_s[0] = (Pw_x_s - Pw_s) / delta_x[0] ;
      gradu_Pw_s[1] = (Pw_y_s - Pw_s) / delta_y[1] ;
      gradu_Pw_n = 0.0;
      gradu_Pw_n[0] = (Pw_x_n - Pw_n) / delta_x[0] ;
      gradu_Pw_n[1] = (Pw_y_n - Pw_n) / delta_y[1] ;

      auto K_s = property.soil.SedimentPermeabilityTensor(cell_inside, iplocal_s)
        * property.hydraulicProperty.PermeabilityScalingFactor(cell_inside,iplocal_s, Sh_s, por_s); /* ndim */
      auto K_n = property.soil.SedimentPermeabilityTensor(cell_outside, iplocal_n)
        * property.hydraulicProperty.PermeabilityScalingFactor(cell_outside,iplocal_n, Sh_n, por_n); /* ndim */

      // compute K * gradient of Pw
      K_s.mv(gradu_Pw_s, Kgradu_Pw_s);
      K_n.mv(gradu_Pw_n, Kgradu_Pw_n);

      // compute K * gradient of Sg
      K_s.mv(gradu_Sg_s, Kgradu_Sg_s);
      K_n.mv(gradu_Sg_n, Kgradu_Sg_n);
      
      // compute K * gradient of Sh
      K_s.mv(gradu_Sh_s, Kgradu_Sh_s);
      K_n.mv(gradu_Sh_n, Kgradu_Sh_n);

      Dune::FieldVector<RF, dim> Kn_F_s;
      K_s.mv(n_F_local, Kn_F_s);
      Dune::FieldVector<RF, dim> Kn_F_n;
      K_n.mv(n_F_local, Kn_F_n);

      auto Swe_s = property.hydraulicProperty.EffectiveSw(Sw_s,Sh_s,0.0,0.0);
      auto dPc_dSwe_s =  property.hydraulicProperty.dPc_dSwe(Swe_s, BrooksCParams_s[0], BrooksCParams_s[1]);/* ndim */
      auto dSwe_dSw_s = property.hydraulicProperty.dSwe_dSw(Sw_s, Sh_s, 0.0, 0.0);
      auto coeff_grad_Sw_s = dPc_dSwe_s * dSwe_dSw_s ;

      auto dPcSF1_dSh_s =  property.hydraulicProperty.dPcSF1_dSh( Sh_s, BrooksCParams_s[1], BrooksCParams_s[4]);
      auto dSwe_dSh_s = property.hydraulicProperty.dSwe_dSh(Sw_s, Sh_s, 0.0, 0.0);
      auto coeff_grad_Sh_s = dPcSF1_dSh_s + dPc_dSwe_s * dSwe_dSh_s ;

      auto Swe_n = property.hydraulicProperty.EffectiveSw(Sw_n,Sh_n,0.0,0.0);
      auto dPc_dSwe_n =  property.hydraulicProperty.dPc_dSwe(Swe_n, BrooksCParams_n[0], BrooksCParams_n[1]);/* ndim */
      auto dSwe_dSw_n = property.hydraulicProperty.dSwe_dSw(Sw_n, Sh_n, 0.0, 0.0);
      auto coeff_grad_Sw_n = dPc_dSwe_n * dSwe_dSw_n ;

      auto dPcSF1_dSh_n =  property.hydraulicProperty.dPcSF1_dSh( Sh_n, BrooksCParams_n[1], BrooksCParams_n[4]);
      auto dSwe_dSh_n = property.hydraulicProperty.dSwe_dSh(Sw_n, Sh_n, 0.0, 0.0);
      auto coeff_grad_Sh_n = dPcSF1_dSh_n + dPc_dSwe_n * dSwe_dSh_n ;


      auto Kgradu_Pg_s = Kgradu_Pw_s - coeff_grad_Sw_s * Kgradu_Sg_s + (coeff_grad_Sh_s - coeff_grad_Sw_s) * Kgradu_Sh_s;
      auto gradu_Pg_s = gradu_Pw_s - coeff_grad_Sw_s * gradu_Sg_s + (coeff_grad_Sh_s - coeff_grad_Sw_s) * gradu_Sh_s;

      auto Kgradu_Pg_n = Kgradu_Pw_n - coeff_grad_Sw_n * Kgradu_Sg_n + (coeff_grad_Sh_n - coeff_grad_Sw_n) * Kgradu_Sh_n;
      auto gradu_Pg_n = gradu_Pw_n - coeff_grad_Sw_n * gradu_Sg_n + (coeff_grad_Sh_n - coeff_grad_Sw_n) * gradu_Sh_n;
      
      K_s.mv(gravity, Kg_s);
      K_n.mv(gravity, Kg_n);

      double S_s = XC_s * (property.salt.MolarMass()/property.water.MolarMass());
      auto krW_s = property.hydraulicProperty.krw(cell_inside, iplocal_s, Sw_s, Sh_s) / (property.water.DynamicViscosity(T_s * Xc_T, Pw_s * Xc_P, S_s) ); /* ndim */
      auto krN_s = property.hydraulicProperty.krg(cell_inside, iplocal_s, Sw_s, Sh_s) / (property.gas.DynamicViscosity(T_s * Xc_T, Pg_s * Xc_P) ); /* ndim */
      
      //  adding terms regarding components
      auto tau_s = property.soil.Tortuosity(por_s); 
      auto DH2O_g_s = tau_s * por_s * property.mixture.DiffCoeffH2OInGas(T_s * Xc_T, Pg_s * Xc_P); /* ndim */
      auto DCH4_w_s = tau_s * por_s * property.mixture.DiffCoeffCH4InLiquid(T_s * Xc_T, Pw_s * Xc_P); /* ndim */
      auto DC_w_s = tau_s * por_s * property.salt.DiffCoeff(T_s * Xc_T, Pw_s * Xc_P); /* ndim */
      auto zCH4_s = property.eos.EvaluateCompressibilityFactor(T_s * Xc_T, Pg_s * Xc_P); 

      auto rho_g_s = property.gas.Density(T_s * Xc_T, Pg_s * Xc_P, zCH4_s) ; /* ndim */
      auto rho_w_s = property.water.Density(T_s * Xc_T, Pw_s * Xc_P, S_s) ; /* ndim */

      auto Cp_g_s = property.gas.Cp(T_s * Xc_T, Pg_s * Xc_P, zCH4_s) ; /* ndim */
      auto Cp_w_s = property.water.Cp(T_s * Xc_T, Pw_s * Xc_P, S_s) ; /* ndim */
      auto kth_g_s = property.gas.ThermalConductivity(T_s * Xc_T, Pg_s * Xc_P) ; /* ndim */
      auto kth_w_s = property.water.ThermalConductivity(T_s * Xc_T, Pw_s * Xc_P, S_s) ; /* ndim */
      auto kth_h_s = property.hydrate.ThermalConductivity(T_s * Xc_T, Peff_s * Xc_P); /* ndim */
      auto kth_s_s = property.soil.ThermalConductivity() ; /* ndim */
      auto kth_eff_s = (1. - por_s) * kth_s_s + por_s * (Sg_s * kth_g_s + Sw_s * kth_w_s + Sh_s * kth_h_s); /* ndim */
     
      double S_n = XC_n * (property.salt.MolarMass()/property.water.MolarMass());
      auto krW_n = property.hydraulicProperty.krw(cell_outside, iplocal_n, Sw_n, Sh_n) / (property.water.DynamicViscosity(T_n * Xc_T, Pw_n * Xc_P, S_n) ); /* ndim */
      auto krN_n = property.hydraulicProperty.krg(cell_outside, iplocal_n, Sw_n, Sh_n) / (property.gas.DynamicViscosity(T_n * Xc_T, Pg_n * Xc_P) ); /* ndim */
      
      auto tau_n = property.soil.Tortuosity(por_n);
      auto DH2O_g_n = tau_n * por_n * property.mixture.DiffCoeffH2OInGas(T_n * Xc_T, Pg_n * Xc_P); /* ndim */
      auto DCH4_w_n = tau_n * por_n * property.mixture.DiffCoeffCH4InLiquid(T_n * Xc_T, Pw_n * Xc_P); /* ndim */
      auto DC_w_n = tau_n * por_n * property.salt.DiffCoeff(T_n * Xc_T, Pw_n * Xc_P); /* ndim */
      auto zCH4_n = property.eos.EvaluateCompressibilityFactor(T_n * Xc_T, Pg_n * Xc_P);
      
      auto rho_g_n = property.gas.Density(T_n * Xc_T, Pg_n * Xc_P, zCH4_n) ; /* ndim */
      auto rho_w_n = property.water.Density(T_n * Xc_T, Pw_n * Xc_P, S_n) ; /* ndim */
      
      auto Cp_g_n = property.gas.Cp(T_n * Xc_T, Pg_n * Xc_P, zCH4_n) ; /* ndim */
      auto Cp_w_n = property.water.Cp(T_n * Xc_T, Pw_n * Xc_P, S_n) ; /* ndim */
      auto kth_g_n = property.gas.ThermalConductivity(T_n * Xc_T, Pg_n * Xc_P) ; /* ndim */
      auto kth_w_n = property.water.ThermalConductivity(T_n * Xc_T, Pw_n * Xc_P, S_n) ; /* ndim */
      auto kth_h_n = property.hydrate.ThermalConductivity(T_n * Xc_T, Peff_n * Xc_P) ; /* ndim */
      auto kth_s_n = property.soil.ThermalConductivity() ; /* ndim */
      auto kth_eff_n = (1. - por_n) * kth_s_n + por_n * (Sg_n * kth_g_n + Sw_n * kth_w_n + Sh_n * kth_h_n);

			for(int i = 0;i<dim;i++){
        v_g[i] = ( omega_s * (Kgradu_Pg_s[i] - rho_g_s * Kg_s[i]) + omega_n * (Kgradu_Pg_n[i] - rho_g_n * Kg_n[i]));
        v_w[i] = ( omega_s * (Kgradu_Pw_s[i] - rho_w_s * Kg_s[i]) + omega_n * (Kgradu_Pw_n[i] - rho_w_n * Kg_n[i]));
      }
      double normalflux_g = -1.*(v_g*n_F_local);
      double normalflux_w = -1.*(v_w*n_F_local);
      double normalflux_T = (omega_s * gradu_T_s + omega_n * gradu_T_n) * n_F_local;

      // upwinding wrt gas-phase velocity
      RF omegaup_g_s, omegaup_g_n;
      if (normalflux_g>=0.0)
      {
        omegaup_g_s = 1.0;
        omegaup_g_n = 0.0;
      }
      else
      {
        omegaup_g_s = 0.0;
        omegaup_g_n = 1.0;
      }
      // upwinding wrt water-phase velocity
      RF omegaup_w_s, omegaup_w_n;
      if (normalflux_w>=0.0)
      {
        omegaup_w_s = 1.0;
        omegaup_w_n = 0.0;
      }
      else
      {
        omegaup_w_s = 0.0;
        omegaup_w_n = 1.0;
      }

      RF omegaup_T_s, omegaup_T_n;
      if (normalflux_T>0.0)
      {
        omegaup_T_s = 0.5;
        omegaup_T_n = 0.5;
      }
      else
      {
        omegaup_T_s = 0.5;
        omegaup_T_n = 0.5;
      }
      auto kth_eff = 2. * kth_eff_s * kth_eff_n / (kth_eff_s + kth_eff_n);
      // integration factor
      auto factor = ip.weight() * geo.integrationElement(ip.position());
      
      // fluxes and diff. flux
      auto convectiveflux_Heat_w_s = rho_w_s * Cp_w_s * (T_s - T_ref) * krW_s * (Kgradu_Pw_s - rho_w_s * Kg_s);
      auto convectiveflux_Heat_g_s = rho_g_s * Cp_g_s * (T_s - T_ref) * krN_s * (Kgradu_Pg_s - rho_g_s * Kg_s);

      auto convectiveflux_Heat_s = omegaup_g_s * convectiveflux_Heat_g_s + omegaup_w_s * convectiveflux_Heat_w_s;

      auto diffusiveflux_Heat_s = kth_eff_s * gradu_T_s;
      // *******************   //
      auto convectiveflux_Heat_w_n = rho_w_n * Cp_w_n * (T_n - T_ref) * krW_n * (Kgradu_Pw_n - rho_w_n * Kg_n);
      auto convectiveflux_Heat_g_n = rho_g_n * Cp_g_n * (T_n - T_ref) * krN_n * (Kgradu_Pg_n - rho_g_n * Kg_n);
 
      auto convectiveflux_Heat_n = omegaup_g_n * convectiveflux_Heat_g_n + omegaup_w_n * convectiveflux_Heat_w_n;
 
      auto diffusiveflux_Heat_n = kth_eff_n * gradu_T_n;

      auto convectiveflux_Heat = - ( convectiveflux_Heat_s + convectiveflux_Heat_n) * n_F_local;

      auto diffusiveflux_Heat = - (omegaup_T_s * diffusiveflux_Heat_s +  omegaup_T_n * diffusiveflux_Heat_n) * n_F_local;

      /*ACCCUMULATE RESIDUALS*/
			double tmp=0.;
      // ENERGY balance
      tmp = Xc_conv_h * convectiveflux_Heat + Xc_diff_h * diffusiveflux_Heat;
      double term_nipg_T = theta_T * (T_s - T_n);
      double term_penalty_T = penalty_factor_t * (T_s - T_n);
      // diffusion term
      for (int i = 0; i < lfsv_T_s.size(); i++)
      {
        r_s.accumulate(lfsv_T_s, i, tmp * psi_T_s[i] * factor);
      }
      for (int i = 0; i < lfsv_T_n.size(); i++)
      {
        r_n.accumulate(lfsv_T_n, i, tmp * -psi_T_n[i] * factor);
      }
      // (non-)symmetric IP term
      for (int i = 0; i < lfsv_T_s.size(); i++)
      {
        r_s.accumulate(lfsv_T_s, i, -omegaup_T_s * Xc_diff_h * term_nipg_T * kth_eff_s * n_F_local * gradpsi_T_s[i] * factor);
      }
      for (int i = 0; i < lfsv_T_n.size(); i++)
      {
        r_n.accumulate(lfsv_T_n, i, -omegaup_T_n * Xc_diff_h * term_nipg_T * kth_eff_n * n_F_local * gradpsi_T_n[i] * factor);
      }
      // standard IP term integral
      for (int i = 0; i < lfsv_T_s.size(); i++)
      {
        r_s.accumulate(lfsv_T_s, i, term_penalty_T * psi_T_s[i] * factor);
      }
      for (int i = 0; i < lfsv_T_n.size(); i++)
      {
        r_n.accumulate(lfsv_T_n, i, term_penalty_T * -psi_T_n[i] * factor);
      }
     
    } //End Quadrature Rule
  }//End of alpha_skeleton

  
  template <typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary(const IG &ig,
                      const LFSU &lfsu_T_s,
                      const X &x,
                      const LFSV &lfsv_T_s,
                      R &r) const
  {
    SUBGFS_Pw gfs_Pw(gfs);
    SUBGFS_Sg gfs_Sg(gfs);
    SUBGFS_XCH4 gfs_XCH4(gfs);
    SUBGFS_YH2O gfs_YH2O(gfs);
    SUBGFS_XC gfs_XC(gfs);
    DGF_Sg dgf_Sg(gfs_Sg, unew);	
    DGF_Sh dgf_Sh(gfs_Sh, unew_Sh);
    DGF_Pw dgf_Pw(gfs_Pw, unew);
    DGF_XCH4 dgf_XCH4(gfs_XCH4, unew);
    DGF_YH2O dgf_YH2O(gfs_YH2O, unew);
    DGF_XC dgf_XC(gfs_XC, unew);

    // dimensions
    const int dimension = GV::dimension;
    const int dim = IG::Entity::dimension;
    const int order_t = std::max(lfsu_T_s.finiteElement().localBasis().order(),
                               lfsv_T_s.finiteElement().localBasis().order());;

    // References to inside and outside cells
    const auto &cell_inside = ig.inside();

    // Get geometries
    auto geo = ig.geometry();
    //const auto dimension = geo.mydimension;
    auto geo_inside = cell_inside.geometry();
    
    // Get geometry of intersection in local coordinates of cell_inside and cell_outside
    auto geo_in_inside = ig.geometryInInside();

    // cell geometries
    auto ref_el_inside 	= referenceElement(geo_inside);
    auto inside_cell_center_local 	= ref_el_inside.position(0,0);
    auto inside_cell_center_global 	= geo_inside.center();

    // face geometry
    auto ref_el = referenceElement(geo);
    auto face_center_local = ref_el.position(0,0);
    auto face_center_global = geo.center();
    
    // face diameter; this should be revised for anisotropic meshes?
    auto h_F = geo_inside.volume() / geo.volume(); // Houston!

    // compute weights
    RF omega_s;
    RF omega_n;
    RF harmonic_average(0.0);
    harmonic_average = 1.0;

    // get polynomial degree
    auto degree = lfsv_T_s.finiteElement().localBasis().order();

    // penalty factor
    auto penalty_factor_t = (alpha_t / h_F) * harmonic_average * degree * (degree + dim - 1);

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_T_s(lfsu_T_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_T_s(lfsv_T_s.size());

    Dune::FieldVector<RF, dim> gradu_Pw_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pw_s(0.0);
    Dune::FieldVector<RF, dim> gradu_Sg_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sg_s(0.0);
    Dune::FieldVector<RF, dim> gradu_Sh_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Sh_s(0.0);
    Dune::FieldVector<RF, dim> gradu_XCH4_s(0.0);
    Dune::FieldVector<RF, dim> gradu_YH2O_s(0.0);
    Dune::FieldVector<RF, dim> gradu_XC_s(0.0);
    Dune::FieldVector<RF, dim> gradu_T_s(0.0);

    Dune::FieldVector<RF, dim> v_g(0.0);
    Dune::FieldVector<RF, dim> v_w(0.0);
    Dune::FieldVector<RF, dim> Kg(0.0);

    Dune::FieldVector<RF, dim> delta_x(0.0);
    Dune::FieldVector<RF, dim> delta_y(0.0);
    delta_x[0] = 1.e-3;
    delta_y[1] = 1.e-3; 

    // Transformation matrix
    typename IG::Entity::Geometry::JacobianInverseTransposed jac;

    // auto intorder = intorderadd+quadrature_factor*order;
    // loop over quadrature points
    for (const auto &ip : quadratureRule(geo, intorder))
    {
      // integration factor
      auto factor = ip.weight() * geo.integrationElement(ip.position());

      // exact normal
      auto n_F_local = ig.unitOuterNormal(ip.position());

      // position of quadrature point in local coordinates of elements
      auto iplocal_s = geo_in_inside.global(ip.position());
      auto ip_global_s = geo_inside.global(iplocal_s);
      auto qp_x_s = iplocal_s + delta_x;
      auto qp_y_s = iplocal_s + delta_y;
      
      // evaluate boundary condition types for {Pw,Sg} or {Fw,Fg} 
			auto bctype = bc.type(ig, ip.position(), (*time)*Xc_t, (*dt)*Xc_t) ;
      auto veltype = bc.velType(ig, ip.position(), (*time)*Xc_t, (*dt)*Xc_t) ;
			// evaluate boundary condition values for {Pw,Sg} or {Fw,Fg} 
			auto bcvalue = bc.value(ig, ip.position(), (*time)*Xc_t, (*dt)*Xc_t ) ;
      auto velvalue = bc.velValue(ig, ip.position(), (*time)*Xc_t, (*dt)*Xc_t) ;

      // evaluate basis functions at local quadrature points 
      auto &psi_T_s = cache_T[order_t].evaluateFunction(iplocal_s, lfsv_T_s.finiteElement().localBasis());
      auto &phi_T_s = cache_T[order_t].evaluateFunction(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      
      // evaluate Sg
      RFT Sg_s0 = 0.0;
      dgf_Sg.evaluate(cell_inside, iplocal_s, Sg_s0);
      RF Sg_s = Sg_s0[0];
      RFT Sg_x_s0 = 0.0;
      dgf_Sg.evaluate(cell_inside, qp_x_s, Sg_x_s0);
      RF Sg_x_s = Sg_x_s0[0];
      RFT Sg_y_s0 = 0.0;
      dgf_Sg.evaluate(cell_inside, qp_y_s, Sg_y_s0);
      RF Sg_y_s = Sg_y_s0[0];

      // evaluate T
      RF T_s = 0.0;
      for (int i = 0; i < lfsu_T_s.size(); i++)
        T_s += x(lfsu_T_s, i) * phi_T_s[i];

      // evaluate Sh
      RFT Sh_s0 = 0.0;
      dgf_Sh.evaluate(cell_inside, iplocal_s, Sh_s0);
      RF Sh_s = Sh_s0[0];
      RFT Sh_x_s0 = 0.0;
      dgf_Sh.evaluate(cell_inside, qp_x_s, Sh_x_s0);
      RF Sh_x_s = Sh_x_s0[0];
      RFT Sh_y_s0 = 0.0;
      dgf_Sh.evaluate(cell_inside, qp_y_s, Sh_y_s0);
      RF Sh_y_s = Sh_y_s0[0];

      // evaluate XCH4
      RFT XCH4_s0 = 0.0;
      dgf_XCH4.evaluate(cell_inside, iplocal_s, XCH4_s0);
      RF XCH4_s = XCH4_s0[0];
      RFT XCH4_x_s0 = 0.0;
      dgf_XCH4.evaluate(cell_inside, qp_x_s, XCH4_x_s0);
      RF XCH4_x_s = XCH4_x_s0[0];
      RFT XCH4_y_s0 = 0.0;
      dgf_XCH4.evaluate(cell_inside, qp_y_s, XCH4_y_s0);
      RF XCH4_y_s = XCH4_y_s0[0];

      // evaluate YH2O
      RFT YH2O_s0 = 0.0;
      dgf_YH2O.evaluate(cell_inside, iplocal_s, YH2O_s0);
      RF YH2O_s = YH2O_s0[0];
      RFT YH2O_x_s0 = 0.0;
      dgf_YH2O.evaluate(cell_inside, qp_x_s, YH2O_x_s0);
      RF YH2O_x_s = YH2O_x_s0[0];
      RFT YH2O_y_s0 = 0.0;
      dgf_YH2O.evaluate(cell_inside, qp_y_s, YH2O_y_s0);
      RF YH2O_y_s = YH2O_y_s0[0];
      

      // evaluate Pw
      RFT Pw_s0 = 0.0;
      dgf_Pw.evaluate(cell_inside, iplocal_s, Pw_s0);
      RF Pw_s = Pw_s0[0];
      RFT Pw_x_s0 = 0.0;
      dgf_Pw.evaluate(cell_inside, qp_x_s, Pw_x_s0);
      RF Pw_x_s = Pw_x_s0[0];
      RFT Pw_y_s0 = 0.0;
      dgf_Pw.evaluate(cell_inside, qp_y_s, Pw_y_s0);
      RF Pw_y_s = Pw_y_s0[0];
      
      // evaluate XC
      RFT XC_s0 = 0.0;
      dgf_XC.evaluate(cell_inside, iplocal_s, XC_s0);
      RF XC_s = XC_s0[0];
      RFT XC_x_s0 = 0.0;
      dgf_XC.evaluate(cell_inside, qp_x_s, XC_x_s0);
      RF XC_x_s = XC_x_s0[0];
      RFT XC_y_s0 = 0.0;
      dgf_XC.evaluate(cell_inside, qp_y_s, XC_y_s0);
      RF XC_y_s = XC_y_s0[0];
      
      RF Pw_n = Pw_s;
      // if (bctype[Indices::PVId_Pw] == Indices::BCId_dirichlet)
      // {
      //   Pw_n = bcvalue[Indices::PVId_Pw] ;
      // }
        
      RF T_n = T_s;
      if (bctype[Indices::PVId_T] == Indices::BCId_dirichlet)
      {
        T_n = bcvalue[Indices::PVId_T] ;
      }

      RF Sh_n = Sh_s ;

      RF Sg_n = Sg_s ;//* (1. - Sh_n);
      // if (bctype[Indices::PVId_Sg] == Indices::BCId_dirichlet)
      // {
      //   Sg_n = bcvalue[Indices::PVId_Sg] ;
      // }


      RF Sw_s = 1. - Sg_s - Sh_s;
      RF Sw_n = 1. - Sg_n - Sh_n;
      
      RF XC_n = XC_s ;
      // if (bctype[Indices::PVId_C] == Indices::BCId_dirichlet)
      // {
      //   XC_n = bcvalue[Indices::PVId_C] ;
      // }

      RF XCH4_n = XCH4_s;
      RF YH2O_n = YH2O_s;


      // evaluate gradient of basis functions
      auto &js_T_s = cache_T[order_t].evaluateJacobian(iplocal_s, lfsu_T_s.finiteElement().localBasis());
      auto &js_v_T_s = cache_T[order_t].evaluateJacobian(iplocal_s, lfsv_T_s.finiteElement().localBasis());

      // transform gradients of shape functions to real element
      jac = geo_inside.jacobianInverseTransposed(iplocal_s);

      for (int i = 0; i < lfsu_T_s.size(); i++)
        jac.mv(js_T_s[i][0], gradphi_T_s[i]);
      for (int i = 0; i < lfsv_T_s.size(); i++)
        jac.mv(js_v_T_s[i][0], gradpsi_T_s[i]);

      // compute gradient of Sg
      gradu_Sg_s = 0.0;
      gradu_Sg_s[0] = (Sg_x_s - Sg_s) / delta_x[0] ;
      gradu_Sg_s[1] = (Sg_y_s - Sg_s) / delta_y[1] ;
      

      // compute gradient of Pw
      gradu_T_s = 0.0;
      for (int i = 0; i < lfsu_T_s.size(); i++)
        gradu_T_s.axpy(x(lfsu_T_s, i), gradphi_T_s[i]);
      
      // compute gradient of Sh
      gradu_Sh_s = 0.0;
      gradu_Sh_s[0] = (Sh_x_s - Sh_s) / delta_x[0] ;
      gradu_Sh_s[1] = (Sh_y_s - Sh_s) / delta_y[1] ;
     
      // compute gradient of XCH4
      gradu_XCH4_s = 0.0;
      gradu_XCH4_s[0] = (XCH4_x_s - XCH4_s) / delta_x[0] ;
      gradu_XCH4_s[1] = (XCH4_y_s - XCH4_s) / delta_y[1] ;

      // compute gradient of YH2O
      gradu_YH2O_s = 0.0;
      gradu_YH2O_s[0] = (YH2O_x_s - YH2O_s) / delta_x[0] ;
      gradu_YH2O_s[1] = (YH2O_y_s - YH2O_s) / delta_y[1] ;

      // compute gradient of Pw
      gradu_Pw_s = 0.0;
      gradu_Pw_s[0] = (Pw_x_s - Pw_s) / delta_x[0] ;
      gradu_Pw_s[1] = (Pw_y_s - Pw_s) / delta_y[1] ;

      double S_s = XC_s * (property.salt.MolarMass()/property.water.MolarMass());
      auto normalgravity = gravity * n_F_local;
      auto BrooksCParams = property.hydraulicProperty.BrooksCoreyParameters(cell_inside, iplocal_s);/*BrooksCParams[0] gives Pentry in Pa*/
      auto por_s = property.soil.SedimentPorosity(cell_inside, iplocal_s);
      auto Pc_s = property.hydraulicProperty.CapillaryPressure(cell_inside, iplocal_s, Sw_s, Sh_s, por_s) ; /* ndim */
      
      RF Pg_s = Pw_s + Pc_s;
      auto por_n = property.soil.SedimentPorosity(cell_inside, iplocal_s);
      auto Pc_n = property.hydraulicProperty.CapillaryPressure(cell_inside, iplocal_s, Sw_n, Sh_n, por_n) ; /* ndim */
      
      RF Pg_n = Pw_n + Pc_n;
      RF Peff_s = (Pg_s * Sg_s + Pw_s * Sw_s) / (1. - Sh_s);
      RF Peff_n = (Pg_n * Sg_n + Pw_n * Sw_n) / (1. - Sh_n);

      
      
      auto K = property.soil.SedimentPermeability(cell_inside,  iplocal_s)
      * property.hydraulicProperty.PermeabilityScalingFactor(cell_inside,iplocal_s, Sh_s, por_s);
      
      auto Swe_s = property.hydraulicProperty.EffectiveSw(Sw_s,Sh_s,0.0,0.0);
      auto dPc_dSwe_s =  property.hydraulicProperty.dPc_dSwe(Swe_s, BrooksCParams[0], BrooksCParams[1]);/* ndim */
      auto dSwe_dSw_s = property.hydraulicProperty.dSwe_dSw(Sw_s, Sh_s, 0.0, 0.0);
      auto coeff_grad_Sw_s = dPc_dSwe_s * dSwe_dSw_s ;

      auto dPcSF1_dSh_s =  property.hydraulicProperty.dPcSF1_dSh( Sh_s, BrooksCParams[1], BrooksCParams[4]);
      auto dSwe_dSh_s = property.hydraulicProperty.dSwe_dSh(Sw_s, Sh_s, 0.0, 0.0);
      auto coeff_grad_Sh_s = dPcSF1_dSh_s + dPc_dSwe_s * dSwe_dSh_s ;

      auto Swe_n = property.hydraulicProperty.EffectiveSw(Sw_n,Sh_n,0.0,0.0);
      auto dPc_dSwe_n =  property.hydraulicProperty.dPc_dSwe(Swe_n, BrooksCParams[0], BrooksCParams[1]);/* ndim */
      auto dSwe_dSw_n = property.hydraulicProperty.dSwe_dSw(Sw_n, Sh_n, 0.0, 0.0);
      auto coeff_grad_Sw_n = dPc_dSwe_n * dSwe_dSw_s ;

      auto dPcSF1_dSh_n =  property.hydraulicProperty.dPcSF1_dSh( Sh_n, BrooksCParams[1], BrooksCParams[4]);
      auto dSwe_dSh_n = property.hydraulicProperty.dSwe_dSh(Sw_n, Sh_n, 0.0, 0.0);
      auto coeff_grad_Sh_n = dPcSF1_dSh_n + dPc_dSwe_n * dSwe_dSh_s ;

      auto krW_s = property.hydraulicProperty.krw(cell_inside, iplocal_s, Sw_s, Sh_s) / (property.water.DynamicViscosity(T_s * Xc_T, Pw_s * Xc_P, S_s));
      auto krN_s = property.hydraulicProperty.krg(cell_inside, iplocal_s, Sw_s, Sh_s) / (property.gas.DynamicViscosity(T_s * Xc_T, Pg_s * Xc_P) );
      
      //  adding terms regarding components
      auto tau_s = property.soil.Tortuosity(por_s);
      auto DH2O_g_s = tau_s * por_s * property.mixture.DiffCoeffH2OInGas(T_s * Xc_T, Pg_s * Xc_P);
      auto DCH4_w_s = tau_s * por_s * property.mixture.DiffCoeffCH4InLiquid(T_s * Xc_T, Pw_s * Xc_P);
      auto DC_w_s = tau_s * por_s * property.salt.DiffCoeff(T_s * Xc_T, Pw_s * Xc_P);
      auto zCH4_s = property.eos.EvaluateCompressibilityFactor(T_s * Xc_T, Pg_s * Xc_P);
      auto YCH4_s =  property.mixture.YCH4(XCH4_s, T_s * Xc_T, Pg_s * Xc_P, XC_s, zCH4_s);
      auto XH2O_s =  property.mixture.XH2O(YH2O_s, T_s * Xc_T, Pg_s * Xc_P, XC_s);
      
      auto rho_g_s = property.gas.Density(T_s * Xc_T, Pg_s * Xc_P, zCH4_s) ;
      auto rho_w_s = property.water.Density(T_s * Xc_T, Pw_s * Xc_P, S_s);

      auto Cp_g_s = property.gas.Cp(T_s * Xc_T, Pg_s * Xc_P, zCH4_s);
      auto Cp_w_s = property.water.Cp(T_s * Xc_T, Pw_s * Xc_P, S_s);
      auto kth_g_s = property.gas.ThermalConductivity(T_s * Xc_T, Pg_s * Xc_P) ;
      auto kth_w_s = property.water.ThermalConductivity(T_s * Xc_T, Pw_s * Xc_P, S_s);
      auto kth_h_s = property.hydrate.ThermalConductivity(T_s * Xc_T, Peff_s * Xc_P);
      auto kth_s_s = property.soil.ThermalConductivity() ;
      auto kth_eff_s = (1. - por_s) * kth_s_s + por_s * (Sg_s * kth_g_s + Sw_s * kth_w_s + Sh_s * kth_h_s);

      double S_n = XC_n * (property.salt.MolarMass()/property.water.MolarMass());
      auto krW_n = property.hydraulicProperty.krw(cell_inside, iplocal_s, Sw_n, Sh_n) / (property.water.DynamicViscosity(T_n * Xc_T, Pw_n * Xc_P, S_n));
      auto krN_n = property.hydraulicProperty.krg(cell_inside, iplocal_s, Sw_n, Sh_n) / (property.gas.DynamicViscosity(T_n * Xc_T, Pg_n * Xc_P) );
      
      auto tau_n = property.soil.Tortuosity(por_n);
      auto DH2O_g_n = tau_n * por_n * property.mixture.DiffCoeffH2OInGas(T_n * Xc_T, Pg_n * Xc_P);
      auto DCH4_w_n = tau_n * por_n * property.mixture.DiffCoeffCH4InLiquid(T_n * Xc_T, Pw_n * Xc_P);
      auto DC_w_n = tau_n * por_n * property.salt.DiffCoeff(T_n * Xc_T, Pw_n * Xc_P);
      auto zCH4_n = property.eos.EvaluateCompressibilityFactor(T_n * Xc_T, Pg_n * Xc_P);


      
      auto rho_g_n = property.gas.Density(T_n * Xc_T, Pg_n * Xc_P, zCH4_n) ;
      auto rho_w_n = property.water.Density(T_n * Xc_T, Pw_n * Xc_P, S_n);
      
      auto Cp_g_n = property.gas.Cp(T_n * Xc_T, Pg_n * Xc_P, zCH4_n);
      auto Cp_w_n = property.water.Cp(T_n * Xc_T, Pw_n * Xc_P, S_n);
      auto kth_g_n = property.gas.ThermalConductivity(T_n * Xc_T, Pg_n * Xc_P) ;
      auto kth_w_n = property.water.ThermalConductivity(T_n * Xc_T, Pw_n * Xc_P, S_n);
      auto kth_h_n = property.hydrate.ThermalConductivity(T_n * Xc_T, Peff_n * Xc_P);
      auto kth_s_n = property.soil.ThermalConductivity() ;
      auto kth_eff_n = (1. - por_n) * kth_s_n + por_n * (Sg_n * kth_g_n + Sw_n * kth_w_n + Sh_n * kth_h_n);
      auto kth_eff = 2. * kth_eff_s * kth_eff_n / (kth_eff_s + kth_eff_n);
      auto h_g_n =  Cp_g_n * (T_n-T_ref) ;
      auto h_w_n =  Cp_w_n * (T_n-T_ref) ;

      auto YCH4_n = property.mixture.YCH4(XCH4_n, T_n * Xc_T, Pg_n * Xc_P, XC_n, zCH4_n);
      auto XH2O_n = property.mixture.XH2O(YH2O_n, T_n * Xc_T, Pg_n * Xc_P, XC_n);

	    // evaluate normal flux of Pw i.e. grad_Pw.n
      RF grad_Pw_s = gradu_Pw_s * n_F_local;
      RF grad_Pw_n = grad_Pw_s;
      // if (bctype[Indices::PVId_Pw] == Indices::BCId_neumann)
      // {
      //   grad_Pw_n = (-1./(K*krW_n)) * velvalue[Indices::BCId_water] + rho_w_n * normalgravity;// NOTE: put the correct coefficients K krw and Muw instead of 1.
      // }
      
      // evaluate normal flux of Sh
      RF grad_Sh_s = gradu_Sh_s * n_F_local;
      RF grad_Sh_n = grad_Sh_s;
      
      // evaluate normal flux of Sg
      RF grad_Sg_s = gradu_Sg_s * n_F_local;
      RF grad_Sg_n = grad_Sg_s;
      
      // if (bctype[Indices::PVId_Sg] == Indices::BCId_neumann)
      // {
      //   //std::cout << coeff_grad_Sw_n << " " << dPc_dSwe_n << " " << dSwe_dSw_n << std::endl;
      //   grad_Sg_n = 0.0;
      //   if (krN_n > 0.){
      //     grad_Sg_n = ((1./(K*krN_n)) * velvalue[Indices::BCId_gas] + grad_Pw_n - rho_g_n * normalgravity 
      //     + (coeff_grad_Sh_n - coeff_grad_Sw_n) * grad_Sh_n) / coeff_grad_Sw_n;// NOTE: put the correct coefficients K krg and Mug instead of 1.
      //   }
      // }

      // evaluate normal flux of XCH4
      RF grad_XCH4_s = gradu_XCH4_s * n_F_local;
      RF grad_XCH4_n = grad_XCH4_s;

      // evaluate normal flux of YH2O
      RF grad_YH2O_s = gradu_YH2O_s * n_F_local;
      RF grad_YH2O_n = grad_YH2O_s;
     
      // evaluate normal flux of T
      RF grad_T_s = gradu_T_s * n_F_local;
      RF grad_T_n = grad_T_s;
      if (veltype[Indices::BCId_heat] == Indices::BCId_neumann)
      {
        grad_T_n =  velvalue[Indices::BCId_heat]; // NOTE: grad_T . n
      }
     
      // evaluate Pg
      auto grad_Pg_s = grad_Pw_s - coeff_grad_Sw_s * grad_Sg_s + (coeff_grad_Sh_s - coeff_grad_Sw_s) * grad_Sh_s;

      auto grad_Pg_n = grad_Pw_n - coeff_grad_Sw_n * grad_Sg_n + (coeff_grad_Sh_n - coeff_grad_Sw_n) * grad_Sh_n;
      
			double tmp = 0.;		


      omega_s = 0.5;
      omega_n = 0.5;

			
			
      auto normalvelocity_g_s = K * krN_s * (grad_Pg_s - rho_g_s * normalgravity);
      
      auto normalvelocity_w_s = K * krW_s * (grad_Pw_s - rho_w_s * normalgravity);
     
      auto normalvelocity_g_n = K * krN_n * (grad_Pg_n - rho_g_n * normalgravity);
      if (veltype[Indices::BCId_gas] = Indices::BCId_neumann ){
        normalvelocity_g_n = velvalue[Indices::BCId_gas];
      }
      
      auto normalvelocity_w_n = K * krW_n * (grad_Pw_n - rho_w_n * normalgravity);
      if (veltype[Indices::BCId_water] = Indices::BCId_neumann){
        normalvelocity_w_n = velvalue[Indices::BCId_water];
      }

      double normalflux_g = -1.*(omega_s * normalvelocity_g_s + omega_n * normalvelocity_g_n);
      double normalflux_w = -1.*(omega_s * normalvelocity_w_s + omega_n * normalvelocity_w_n);
      double normalflux_T = (omega_s * grad_T_s + omega_n * grad_T_n);

      // upwinding wrt gas-phase velocity
      RF omegaup_g_s, omegaup_g_n;
      if (normalflux_g>=0.0)
      {
        omegaup_g_s = 1.0;
        omegaup_g_n = 0.0;
      }
      else
      {
        omegaup_g_s = 0.0;
        omegaup_g_n = 1.0;
      }
      // upwinding wrt water-phase velocity
      RF omegaup_w_s, omegaup_w_n;
      if (normalflux_w>=0.0)
      {
        omegaup_w_s = 1.0;
        omegaup_w_n = 0.0;
      }
      else
      {
        omegaup_w_s = 0.0;
        omegaup_w_n = 1.0;
      }

      RF omegaup_T_s, omegaup_T_n;
      if (normalflux_T>0.0)
      {
        omegaup_T_s = 0.5;
        omegaup_T_n = 0.5;
      }
      else
      {
        omegaup_T_s = 0.5;
        omegaup_T_n = 0.5;
      }

      //   fluxes and diff. flux
      auto convectiveflux_Heat_w_s = rho_w_s * Cp_w_s * (T_s - T_ref) * normalvelocity_w_s;
      auto convectiveflux_Heat_g_s = rho_g_s * Cp_g_s * (T_s - T_ref) * normalvelocity_g_s;

      auto convectiveflux_Heat_s = omegaup_g_s * convectiveflux_Heat_g_s + omegaup_w_s * convectiveflux_Heat_w_s;

      auto diffusiveflux_Heat_s = kth_eff_s * grad_T_s;  

      //   *******************   //
      auto convectiveflux_Heat_w_n = rho_w_n * Cp_w_n * (T_n - T_ref) * normalvelocity_w_n;
      auto convectiveflux_Heat_g_n = rho_g_n * Cp_g_n * (T_n - T_ref) * normalvelocity_g_n;

      auto convectiveflux_Heat_n = omegaup_g_n * convectiveflux_Heat_g_n + omegaup_w_n * convectiveflux_Heat_w_n;

      auto diffusiveflux_Heat_n = kth_eff_n * grad_T_n; 

      
      auto convectiveflux_Heat_g = omegaup_g_s * convectiveflux_Heat_g_s + omegaup_g_n * convectiveflux_Heat_g_n;
      if (veltype[Indices::BCId_gas] == Indices::BCId_neumann){
        convectiveflux_Heat_g = 0.5 * ( rho_g_s * Cp_g_s * (T_s - T_ref) + rho_g_n * Cp_g_n * (T_n - T_ref)) * normalvelocity_g_n;
      }

      auto convectiveflux_Heat_w = omegaup_w_s * convectiveflux_Heat_w_s + omegaup_w_n * convectiveflux_Heat_w_n;
      if (veltype[Indices::BCId_water] == Indices::BCId_neumann){
        convectiveflux_Heat_w = 0.5 * ( rho_w_s * Cp_w_s * (T_s - T_ref) + rho_w_n * Cp_w_n * (T_n - T_ref)) * normalvelocity_w_n;
      }

      auto convectiveflux_Heat = - ( convectiveflux_Heat_g + convectiveflux_Heat_w);
      auto diffusiveflux_Heat = - (omegaup_T_s * diffusiveflux_Heat_s + omegaup_T_n * diffusiveflux_Heat_n);
      if (veltype[Indices::BCId_heat] == Indices::BCId_neumann){
        diffusiveflux_Heat = - (kth_eff_n * grad_T_n);
      }

      //  ACCCUMULATE RESIDUALS  //
			tmp=0.;
      
      // ENERGY balance
      tmp = Xc_conv_h * convectiveflux_Heat + Xc_diff_h * diffusiveflux_Heat;
      double term_nipg_T = theta_T * (T_s - T_n);
      double term_penalty_T = penalty_factor_t * (T_s - T_n);

      for (int i = 0; i < lfsv_T_s.size(); i++)
      {
        r.accumulate(lfsv_T_s, i, tmp * psi_T_s[i] * factor);
      }
      
      // (non-)symmetric IP term
      for (int i = 0; i < lfsv_T_s.size(); i++)
      {
        r.accumulate(lfsv_T_s, i, - Xc_diff_h * kth_eff_s * term_nipg_T * n_F_local * gradpsi_T_s[i] * factor); // in the run testAveragingXC-T there is no upwinding for sym terms
      }
      
      // standard IP term integral
      for (int i = 0; i < lfsv_T_s.size(); i++)
      {
        r.accumulate(lfsv_T_s, i, term_penalty_T * psi_T_s[i] * factor);
      }
    } // end of quadrature rule
  } // end of alpha_boundary
  
};
