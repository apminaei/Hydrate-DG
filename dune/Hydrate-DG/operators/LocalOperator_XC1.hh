/*
 * LocalOperator.hh
 *
 * 
 */



using namespace Dune::PDELab;


template <typename GV, typename Params, class BC, typename U_Sh, class GFS_Sh,
          typename U, class GFS, typename U_T, class GFS_T,
          typename U_XC, class GFS_XC, class FEM_S>
class LocalOperator_XC : 
                      public Dune::PDELab::NumericalJacobianApplyVolume<LocalOperator_XC<GV, Params, BC, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, U_XC, GFS_XC, FEM_S>>,
                      public Dune::PDELab::NumericalJacobianVolume<LocalOperator_XC<GV, Params, BC, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, U_XC, GFS_XC, FEM_S>>,
                      public Dune::PDELab::NumericalJacobianApplySkeleton<LocalOperator_XC<GV, Params, BC, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, U_XC, GFS_XC, FEM_S>>,
                      public Dune::PDELab::NumericalJacobianSkeleton<LocalOperator_XC<GV, Params, BC, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, U_XC, GFS_XC, FEM_S>>,  
                      public Dune::PDELab::NumericalJacobianApplyBoundary<LocalOperator_XC<GV, Params, BC, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, U_XC, GFS_XC, FEM_S>>,
                      public Dune::PDELab::NumericalJacobianBoundary<LocalOperator_XC<GV, Params, BC, U_Sh, GFS_Sh,
                      U, GFS, U_T, GFS_T, U_XC, GFS_XC, FEM_S>>,
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
  U_T unew_T;
  GFS_T gfs_T;
  U_XC *unew_XC;
  GFS_XC gfs_XC;

  double *time;
  double *dt;
  double alpha_x;
  double method_x;

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

  
  typedef Dune::PDELab::LocalFunctionSpace<GFS_XC> LFS;
  typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
  typedef typename U_XC::template LocalView<LFSCache> VectorView;
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
  
  std::vector<Cache_XC> cache_XC;
  

  // constructor stores parameters
  LocalOperator_XC( const GV &gv_, const Params&	property_, const BC& bc_,
                    const U_Sh &unew_Sh_, GFS_Sh gfs_Sh_, const U &unew_, GFS gfs_, 
                    const U_T &unew_T_, GFS_T gfs_T_, 
                    U_XC *unew_XC_, GFS_XC gfs_XC_,
                    double *time_,
                    double *dt_,
                    unsigned int intorder_ = 6, double method_x_=0., double alpha_x_=1.)
      : gv(gv_), property( property_ ), bc( bc_ ),
        unew_Sh(unew_Sh_), gfs_Sh(gfs_Sh_), unew(unew_), gfs(gfs_),
        unew_T(unew_T_), gfs_T(gfs_T_), 
        unew_XC(unew_XC_), gfs_XC(gfs_XC_), 
        time(time_),
        dt(dt_), method_x(method_x_),
        alpha_x(alpha_x_),
        intorder(intorder_), cache_XC(20)
  {
    theta_x = method_x;
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
    
    // dimensions
    const int dim = EG::Entity::dimension;
    
    const int order_x = std::max(lfsu_XC.finiteElement().localBasis().order(),
                               lfsv_XC.finiteElement().localBasis().order());
   
    // Reference to cell
	  const auto& cell = eg.entity();
		const IndexSet &indexSet = gv.indexSet();
		int cell_number = indexSet.index(cell);

    // Get geometry
    auto geo = eg.geometry();

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XC(lfsu_XC.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XC(lfsv_XC.size());
    
    Dune::FieldVector<RF, dim> gradu_Pw(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pw(0.0);
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
      auto &phi_XC = cache_XC[order_x].evaluateFunction(ip.position(), lfsu_XC.finiteElement().localBasis());
      auto &psi_XC = cache_XC[order_x].evaluateFunction(ip.position(), lfsv_XC.finiteElement().localBasis());

      auto qp_x = ip.position() + delta_x;
      auto qp_y = ip.position() + delta_y;
      auto ip_global = geo.global(ip.position());
      auto ip_local = geo.local(ip_global);
      
      
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
      RFT Pw_x0 = 0.0;
      dgf_Pw.evaluate(cell, qp_x, Pw_x0);
      RF Pw_x = Pw_x0[0];
      RFT Pw_y0 = 0.0;
      dgf_Pw.evaluate(cell, qp_y, Pw_y0);
      RF Pw_y = Pw_y0[0];

      // evaluate Sw
      RF Sw = 1. - Sg - Sh;

      // evaluate Pg
      auto BrooksCParams = property.hydraulicProperty.BrooksCoreyParameters(cell, ip_local);/*BrooksCParams[0] gives Pentry in Pa*/
      auto por = property.soil.SedimentPorosity(cell, ip_local);
      auto Pc = property.hydraulicProperty.CapillaryPressure(cell, ip_local, Sw, Sh, por) ; /* ndim */
      RF Pg = Pw + Pc; /* ndim */
      RF Peff = (Pg * Sg + Pw * Sw) / (1. - Sh); /* ndim */
      
      // evaluate gradient of basis functions
      auto &js_XC = cache_XC[order_x].evaluateJacobian(ip.position(), lfsu_XC.finiteElement().localBasis());
      auto &js_v_XC = cache_XC[order_x].evaluateJacobian(ip.position(), lfsv_XC.finiteElement().localBasis());

      // transform gradients of shape functions to real element
      jac = geo.jacobianInverseTransposed(ip.position());

      for (int i = 0; i < lfsu_XC.size(); i++)
        jac.mv(js_XC[i][0], gradphi_XC[i]);
      for (int i = 0; i < lfsv_XC.size(); i++)
        jac.mv(js_v_XC[i][0], gradpsi_XC[i]);


      // compute gradient of XC
      gradu_XC = 0.0;
      for (int i = 0; i < lfsu_XC.size(); i++)
        gradu_XC.axpy(x(lfsu_XC, i), gradphi_XC[i]);

      
      // compute gradient of Pw
      gradu_Pw[0] = (Pw_x - Pw) / delta_x[0] ;
      gradu_Pw[1] = (Pw_y - Pw) / delta_y[1] ;


      auto K = property.soil.SedimentPermeabilityTensor(cell, ip_local)
                    * property.hydraulicProperty.PermeabilityScalingFactor(cell, ip_local, Sh, por ); /*ndim K from soil.hh*/
      K.mv(gravity, Kg);

      // compute K * gradient of Pw
      K.mv(gradu_Pw, Kgradu_Pw);

      
      auto permeability = property.soil.SedimentPermeability(cell, ip_local )/*ndim K from soil.hh*/
							  * property.hydraulicProperty.PermeabilityScalingFactor(cell, ip_local, Sh, por );

      auto tau = property.soil.Tortuosity(por);/*ndim tau from soil.hh*/
      auto DC_w = tau * por * property.salt.DiffCoeff(T * Xc_T, Pw * Xc_P); /*ndim D from salt.hh*/
      
      double S = XC * (property.salt.MolarMass()/property.water.MolarMass());
      auto zCH4 = property.eos.EvaluateCompressibilityFactor(T * Xc_T, Pg * Xc_P);
      
      auto rho_w = property.water.Density(T * Xc_T, Pw * Xc_P, S); /*ndim density from H2O.hh; the input arguments are dimensional*/
      
      auto krW = property.hydraulicProperty.krw(cell, ip_local, Sw, Sh) / (property.water.DynamicViscosity(T * Xc_T, Pw * Xc_P, S) );
      		
			auto q_s = property.salt.Source(); /*kg/m³s*/

      auto convectiveflux_SALT_w = rho_w * (XC) * krW * (Kgradu_Pw - rho_w * Kg);
      
      auto j_SALT_w = rho_w * Sw * DC_w * gradu_XC;


      auto diffusiveflux_SALT = j_SALT_w;

      // integrate (A grad u - bu)*grad phi_i + a*u*phi_i
      RF factor = ip.weight() * geo.integrationElement(ip.position());
      for (int i = 0; i < lfsv_XC.size(); i++)
      {
        r.accumulate(lfsv_XC, i, ((Xc_conv_m * convectiveflux_SALT_w  
                                  - Xc_diff_m * diffusiveflux_SALT ) * gradpsi_XC[i]
                                  - Xc_source_m*q_s * psi_XC[i]) * factor);
      }

    } //End Quadrature Rule
  }  // End of alpha_volume

  // skeleton integral depending on test and ansatz functions
  // each face is only visited ONCE!
  template <typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_skeleton(const IG &ig,
                      const LFSU &lfsu_XC_s, const X &x_s, const LFSV &lfsv_XC_s,
                      const LFSU &lfsu_XC_n, const X &x_n, const LFSV &lfsv_XC_n,
                      R &r_s, R &r_n) const
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

    // dimensions
    const int dim= IG::Entity::dimension;
    const int dimension = GV::dimension;

    const int order_x = std::max(lfsu_XC_s.finiteElement().localBasis().order(),
                                lfsv_XC_s.finiteElement().localBasis().order());


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
    auto order_i = lfsv_XC_s.finiteElement().localBasis().order();
    auto order_o = lfsv_XC_n.finiteElement().localBasis().order();
    auto degree = std::max(order_i, order_o);

    // penalty factor
    auto penalty_factor_x = (alpha_x / h_F) * harmonic_average * degree * (degree + dim - 1);

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XC_s(lfsu_XC_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XC_s(lfsv_XC_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XC_n(lfsu_XC_n.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XC_n(lfsv_XC_n.size());

    Dune::FieldVector<RF, dim> gradu_Pw_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pw_s(0.0);
    Dune::FieldVector<RF, dim> gradu_XC_s(0.0);

    Dune::FieldVector<RF, dim> gradu_Pw_n(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pw_n(0.0);
    Dune::FieldVector<RF, dim> gradu_XC_n(0.0);

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
      auto &phi_XC_s = cache_XC[order_x].evaluateFunction(iplocal_s, lfsu_XC_s.finiteElement().localBasis());
      auto &psi_XC_s = cache_XC[order_x].evaluateFunction(iplocal_s, lfsv_XC_s.finiteElement().localBasis());
      auto &phi_XC_n = cache_XC[order_x].evaluateFunction(iplocal_n, lfsu_XC_n.finiteElement().localBasis());
      auto &psi_XC_n = cache_XC[order_x].evaluateFunction(iplocal_n, lfsv_XC_n.finiteElement().localBasis());
      
      
      // evaluate Sg
      RFT Sg_s0 = 0.0;
      dgf_Sg.evaluate(cell_inside, iplocal_s, Sg_s0);
      RF Sg_s = Sg_s0[0];
      
      RFT Sg_n0 = 0.0;
      dgf_Sg.evaluate(cell_outside, iplocal_n, Sg_n0);
      RF Sg_n = Sg_n0[0];


      // evaluate XC
      RF XC_s = 0.0;
      for (int i = 0; i < lfsu_XC_s.size(); i++)
        XC_s += x_s(lfsu_XC_s, i) * phi_XC_s[i];
      RF XC_n = 0.0;
      for (int i = 0; i < lfsu_XC_n.size(); i++)
        XC_n += x_n(lfsu_XC_n, i) * phi_XC_n[i];

      // evaluate Sh
      RFT Sh_s0 = 0.0;
      dgf_Sh.evaluate(cell_inside, iplocal_s, Sh_s0);
      RF Sh_s = Sh_s0[0];
      
      RFT Sh_n0 = 0.0;
      dgf_Sh.evaluate(cell_outside, iplocal_n, Sh_n0);
      RF Sh_n = Sh_n0[0];

      // evaluate T
      RFT T_s0 = 0.0;
      dgf_T.evaluate(cell_inside, iplocal_s, T_s0);
      RF T_s = T_s0[0];
      RFT T_n0 = 0.0;
      dgf_T.evaluate(cell_outside, iplocal_n, T_n0);
      RF T_n = T_n0[0];
      

      // evaluate XCH4
      RFT XCH4_s0 = 0.0;
      dgf_XCH4.evaluate(cell_inside, iplocal_s, XCH4_s0);
      RF XCH4_s = XCH4_s0[0];
      
      RFT XCH4_n0 = 0.0;
      dgf_XCH4.evaluate(cell_outside, iplocal_n, XCH4_n0);
      RF XCH4_n = XCH4_n0[0];

      // evaluate YH2O
      RFT YH2O_s0 = 0.0;
      dgf_YH2O.evaluate(cell_inside, iplocal_s, YH2O_s0);
      RF YH2O_s = YH2O_s0[0];
      
      RFT YH2O_n0 = 0.0;
      dgf_YH2O.evaluate(cell_outside, iplocal_n, YH2O_n0);
      RF YH2O_n = YH2O_n0[0];

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
      auto &js_XC_s = cache_XC[order_x].evaluateJacobian(iplocal_s, lfsu_XC_s.finiteElement().localBasis());
      auto &js_v_XC_s = cache_XC[order_x].evaluateJacobian(iplocal_s, lfsv_XC_s.finiteElement().localBasis());
      auto &js_XC_n = cache_XC[order_x].evaluateJacobian(iplocal_n, lfsu_XC_n.finiteElement().localBasis());
      auto &js_v_XC_n = cache_XC[order_x].evaluateJacobian(iplocal_n, lfsv_XC_n.finiteElement().localBasis());

      // transform gradients of shape functions to real element
      jac = geo_inside.jacobianInverseTransposed(iplocal_s);

      for (int i = 0; i < lfsu_XC_s.size(); i++)
        jac.mv(js_XC_s[i][0], gradphi_XC_s[i]);
      for (int i = 0; i < lfsv_XC_s.size(); i++)
        jac.mv(js_v_XC_s[i][0], gradpsi_XC_s[i]);

      jac = geo_outside.jacobianInverseTransposed(iplocal_n);

      for (int i = 0; i < lfsu_XC_n.size(); i++)
        jac.mv(js_XC_n[i][0], gradphi_XC_n[i]);
      for (int i = 0; i < lfsv_XC_n.size(); i++)
        jac.mv(js_v_XC_n[i][0], gradpsi_XC_n[i]);

      
      // compute gradient of XC
      gradu_XC_s = 0.0;
      for (int i = 0; i < lfsu_XC_s.size(); i++)
        gradu_XC_s.axpy(x_s(lfsu_XC_s, i), gradphi_XC_s[i]);
      gradu_XC_n = 0.0;
      for (int i = 0; i < lfsu_XC_n.size(); i++)
        gradu_XC_n.axpy(x_n(lfsu_XC_n, i), gradphi_XC_n[i]);

      
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

      Dune::FieldVector<RF, dim> Kn_F_s;
      K_s.mv(n_F_local, Kn_F_s);
      Dune::FieldVector<RF, dim> Kn_F_n;
      K_n.mv(n_F_local, Kn_F_n);

      
      K_s.mv(gravity, Kg_s);
      K_n.mv(gravity, Kg_n);

      double S_s = XC_s * (property.salt.MolarMass()/property.water.MolarMass());
      auto krW_s = property.hydraulicProperty.krw(cell_inside, iplocal_s, Sw_s, Sh_s) / (property.water.DynamicViscosity(T_s * Xc_T, Pw_s * Xc_P, S_s) ); /* ndim */
      
      //  adding terms regarding components
      auto tau_s = property.soil.Tortuosity(por_s); 
      auto DC_w_s = tau_s * por_s * property.salt.DiffCoeff(T_s * Xc_T, Pw_s * Xc_P); /* ndim */
      auto zCH4_s = property.eos.EvaluateCompressibilityFactor(T_s * Xc_T, Pg_s * Xc_P); 

      auto rho_w_s = property.water.Density(T_s * Xc_T, Pw_s * Xc_P, S_s) ; /* ndim */
     
      double S_n = XC_n * (property.salt.MolarMass()/property.water.MolarMass());
      auto krW_n = property.hydraulicProperty.krw(cell_outside, iplocal_n, Sw_n, Sh_n) / (property.water.DynamicViscosity(T_n * Xc_T, Pw_n * Xc_P, S_n) ); /* ndim */
      
      auto tau_n = property.soil.Tortuosity(por_n);
      auto DC_w_n = tau_n * por_n * property.salt.DiffCoeff(T_n * Xc_T, Pw_n * Xc_P); /* ndim */
      auto zCH4_n = property.eos.EvaluateCompressibilityFactor(T_n * Xc_T, Pg_n * Xc_P);
      
      auto rho_w_n = property.water.Density(T_n * Xc_T, Pw_n * Xc_P, S_n) ; /* ndim */
      
			for(int i = 0;i<dim;i++){
        v_w[i] = ( omega_s * (Kgradu_Pw_s[i] - rho_w_s * Kg_s[i]) + omega_n * (Kgradu_Pw_n[i] - rho_w_n * Kg_n[i]));
      }
      double normalflux_w = -1.*(v_w*n_F_local);
      double normalflux_x = (omega_s * gradu_XC_s + omega_n * gradu_XC_n) * n_F_local;

     
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

      RF omegaup_x_s, omegaup_x_n;
      if (normalflux_x>0.0)
      {
        omegaup_x_s = 0.5;
        omegaup_x_n = 0.5;
      }
      else
      {
        omegaup_x_s = 0.5;
        omegaup_x_n = 0.5;
      }

      // integration factor
      auto factor = ip.weight() * geo.integrationElement(ip.position());
      // fluxes and diff. flux
      auto convectiveflux_SALT_w_s = rho_w_s * (XC_s) * krW_s * (Kgradu_Pw_s - rho_w_s * Kg_s);

      auto j_SALT_w_s = rho_w_s * Sw_s * DC_w_s * gradu_XC_s;


      auto diffusiveflux_SALT_s = j_SALT_w_s;
      // *******************   //
      auto convectiveflux_SALT_w_n = rho_w_n * (XC_n) * krW_n * (Kgradu_Pw_n - rho_w_n * Kg_n);
      
      auto j_SALT_w_n = rho_w_n * Sw_n * DC_w_n * gradu_XC_n;
      
      
      auto diffusiveflux_SALT_n = j_SALT_w_n;
      
      auto convectiveflux_SALT = -(omegaup_w_s * convectiveflux_SALT_w_s + omegaup_w_n * convectiveflux_SALT_w_n) * n_F_local;
      
      auto diffusiveflux_SALT = (omegaup_x_s * diffusiveflux_SALT_s + omegaup_x_n * diffusiveflux_SALT_n) * n_F_local;
      
      /*ACCCUMULATE RESIDUALS*/
			double tmp=0.;
      // SALT-component-wise mass-balance
      tmp = Xc_conv_m * convectiveflux_SALT + Xc_diff_m * diffusiveflux_SALT ;
            
      double term_nipg_c_x = theta_x * (XC_s - XC_n);

      double term_penalty_c = penalty_factor_x * (XC_s - XC_n);
      // diffusion term
      for (int i = 0; i < lfsv_XC_s.size(); i++)
      {
        r_s.accumulate(lfsv_XC_s, i, tmp * psi_XC_s[i] * factor);
      }
      for (int i = 0; i < lfsv_XC_n.size(); i++)
      {
        r_n.accumulate(lfsv_XC_n, i, tmp * -psi_XC_n[i] * factor);
      }
      // (non-)symmetric IP term
      for (int i = 0; i < lfsv_XC_s.size(); i++)
      {
        r_s.accumulate(lfsv_XC_s, i,  omegaup_x_s * Xc_diff_m * term_nipg_c_x * rho_w_s 
                                      * Sw_s * DC_w_s * gradpsi_XC_s[i] * n_F_local * factor);
      }
      for (int i = 0; i < lfsv_XC_n.size(); i++)
      {
        r_n.accumulate(lfsv_XC_n, i,  omegaup_x_n * Xc_diff_m * term_nipg_c_x * rho_w_n 
                                      * Sw_n * DC_w_n * gradpsi_XC_n[i] * n_F_local * factor);
      }
      // standard IP term integral
      for (int i = 0; i < lfsv_XC_s.size(); i++)
      {
        r_s.accumulate(lfsv_XC_s, i, term_penalty_c * psi_XC_s[i] * factor);
      }
      for (int i = 0; i < lfsv_XC_n.size(); i++)
      {
        r_n.accumulate(lfsv_XC_n, i, term_penalty_c * -psi_XC_n[i] * factor);
      }
     
    } //End Quadrature Rule
  }//End of alpha_skeleton

  
  template <typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary(const IG &ig,
                      const LFSU &lfsu_XC_s,
                      const X &x,
                      const LFSV &lfsv_XC_s,
                      R &r) const
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

    // dimensions
    const int dimension = GV::dimension;
    const int dim = IG::Entity::dimension;
    const int order_x = std::max(lfsu_XC_s.finiteElement().localBasis().order(),
                               lfsv_XC_s.finiteElement().localBasis().order());;

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
    auto degree = lfsv_XC_s.finiteElement().localBasis().order();

    // penalty factor
    auto penalty_factor_x = (alpha_x / h_F) * harmonic_average * degree * (degree + dim - 1);

    // Initialize vectors outside for loop
    std::vector<Dune::FieldVector<RF, dim>> gradphi_XC_s(lfsu_XC_s.size());
    std::vector<Dune::FieldVector<RF, dim>> gradpsi_XC_s(lfsv_XC_s.size());

    Dune::FieldVector<RF, dim> gradu_Pw_s(0.0);
    Dune::FieldVector<RF, dim> Kgradu_Pw_s(0.0);
    Dune::FieldVector<RF, dim> gradu_XC_s(0.0);

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
      auto &psi_XC_s = cache_XC[order_x].evaluateFunction(iplocal_s, lfsv_XC_s.finiteElement().localBasis());
      auto &phi_XC_s = cache_XC[order_x].evaluateFunction(iplocal_s, lfsu_XC_s.finiteElement().localBasis());
      
      // evaluate Sg
      RFT Sg_s0 = 0.0;
      dgf_Sg.evaluate(cell_inside, iplocal_s, Sg_s0);
      RF Sg_s = Sg_s0[0];

      // evaluate XC
      RF XC_s = 0.0;
      for (int i = 0; i < lfsu_XC_s.size(); i++)
        XC_s += x(lfsu_XC_s, i) * phi_XC_s[i];

      // evaluate Sh
      RFT Sh_s0 = 0.0;
      dgf_Sh.evaluate(cell_inside, iplocal_s, Sh_s0);
      RF Sh_s = Sh_s0[0];

      // evaluate T
      RFT T_s0 = 0.0;
      dgf_T.evaluate(cell_inside, iplocal_s, T_s0);
      RF T_s = T_s0[0];
      

      // evaluate XCH4
      RFT XCH4_s0 = 0.0;
      dgf_XCH4.evaluate(cell_inside, iplocal_s, XCH4_s0);
      RF XCH4_s = XCH4_s0[0];

      // evaluate YH2O
      RFT YH2O_s0 = 0.0;
      dgf_YH2O.evaluate(cell_inside, iplocal_s, YH2O_s0);
      RF YH2O_s = YH2O_s0[0];
      

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
      
      
      RF Pw_n = Pw_s;
      // if (bctype[Indices::PVId_Pw] == Indices::BCId_dirichlet)
      // {
      //   Pw_n = bcvalue[Indices::PVId_Pw] ;
      // }
        
      RF T_n = T_s;
      // if (bctype[Indices::PVId_T] == Indices::BCId_dirichlet)
      // {
      //   T_n = bcvalue[Indices::PVId_T] ;
      // }

      RF Sh_n = Sh_s ;

      RF Sg_n = Sg_s ;//* (1. - Sh_n);
      // if (bctype[Indices::PVId_Sg] == Indices::BCId_dirichlet)
      // {
      //   Sg_n = bcvalue[Indices::PVId_Sg] ;
      // }


      RF Sw_s = 1. - Sg_s - Sh_s;
      RF Sw_n = 1. - Sg_n - Sh_n;
      
      RF XC_n = XC_s ;
      if (bctype[Indices::PVId_C] == Indices::BCId_dirichlet)
      {
        XC_n = bcvalue[Indices::PVId_C] ;
      }

      RF XCH4_n = XCH4_s;
      RF YH2O_n = YH2O_s;


      // evaluate gradient of basis functions
      auto &js_XC_s = cache_XC[order_x].evaluateJacobian(iplocal_s, lfsu_XC_s.finiteElement().localBasis());
      auto &js_v_XC_s = cache_XC[order_x].evaluateJacobian(iplocal_s, lfsv_XC_s.finiteElement().localBasis());

      // transform gradients of shape functions to real element
      jac = geo_inside.jacobianInverseTransposed(iplocal_s);

      for (int i = 0; i < lfsu_XC_s.size(); i++)
        jac.mv(js_XC_s[i][0], gradphi_XC_s[i]);
      for (int i = 0; i < lfsv_XC_s.size(); i++)
        jac.mv(js_v_XC_s[i][0], gradpsi_XC_s[i]);

      
      // compute gradient of XC
      gradu_XC_s = 0.0;
      for (int i = 0; i < lfsu_XC_s.size(); i++)
        gradu_XC_s.axpy(x(lfsu_XC_s, i), gradphi_XC_s[i]);
      
     
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
      
      
      auto krW_s = property.hydraulicProperty.krw(cell_inside, iplocal_s, Sw_s, Sh_s) / (property.water.DynamicViscosity(T_s * Xc_T, Pw_s * Xc_P, S_s));
      
      //  adding terms regarding components
      auto tau_s = property.soil.Tortuosity(por_s);
      auto DC_w_s = tau_s * por_s * property.salt.DiffCoeff(T_s * Xc_T, Pw_s * Xc_P);
      auto zCH4_s = property.eos.EvaluateCompressibilityFactor(T_s * Xc_T, Pg_s * Xc_P);
      
      auto rho_w_s = property.water.Density(T_s * Xc_T, Pw_s * Xc_P, S_s);
     

      double S_n = XC_n * (property.salt.MolarMass()/property.water.MolarMass());
      auto krW_n = property.hydraulicProperty.krw(cell_inside, iplocal_s, Sw_n, Sh_n) / (property.water.DynamicViscosity(T_n * Xc_T, Pw_n * Xc_P, S_n));
      
      auto tau_n = property.soil.Tortuosity(por_n);
      auto DC_w_n = tau_n * por_n * property.salt.DiffCoeff(T_n * Xc_T, Pw_n * Xc_P);
      auto zCH4_n = property.eos.EvaluateCompressibilityFactor(T_n * Xc_T, Pg_n * Xc_P);


      
      auto rho_w_n = property.water.Density(T_n * Xc_T, Pw_n * Xc_P, S_n);
      
      // evaluate normal flux of Pw i.e. grad_Pw.n
      RF grad_Pw_s = gradu_Pw_s * n_F_local;
      RF grad_Pw_n = grad_Pw_s;
      if (bctype[Indices::PVId_Pw] == Indices::BCId_neumann)
      {
        grad_Pw_n = (-1./(K*krW_n)) * velvalue[Indices::BCId_water] + rho_w_n * normalgravity;// NOTE: put the correct coefficients K krw and Muw instead of 1.
      }
      

      
      // evaluate normal flux of XC
      RF grad_XC_s = gradu_XC_s * n_F_local;
      RF grad_XC_n = grad_XC_s;
      if (veltype[Indices::BCId_salt] == Indices::BCId_neumann)
      {
        grad_XC_n = velvalue[Indices::BCId_salt];
      }
     
      
			double tmp = 0.;		


      omega_s = 0.5;
      omega_n = 0.5;

			
			
      auto normalvelocity_w_s = K * krW_s * (grad_Pw_s - rho_w_s * normalgravity);
     
      auto normalvelocity_w_n = K * krW_n * (grad_Pw_n - rho_w_n * normalgravity);
      if (veltype[Indices::BCId_water] = Indices::BCId_neumann){
        normalvelocity_w_n = velvalue[Indices::BCId_water];
      }

      double normalflux_w = -1.*(omega_s * normalvelocity_w_s + omega_n * normalvelocity_w_n);
      double normalflux_x = (omega_s * grad_XC_s + omega_n * grad_XC_n);

      
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

      RF omegaup_x_s, omegaup_x_n;
      if (normalflux_x>0.0)
      {
        omegaup_x_s = 0.5;
        omegaup_x_n = 0.5;
      }
      else
      {
        omegaup_x_s = 0.5;
        omegaup_x_n = 0.5;
      }

      //   fluxes and diff. flux

      auto convectiveflux_SALT_w_s = rho_w_s * (XC_s) * normalvelocity_w_s;
      auto j_SALT_w_s = rho_w_s * Sw_s * DC_w_s * grad_XC_s;
      
      
      auto diffusiveflux_SALT_s = j_SALT_w_s;
     
      //   *******************   //
      auto convectiveflux_SALT_w_n = rho_w_n * (XC_n) * normalvelocity_w_n;

      auto j_SALT_w_n = rho_w_n * Sw_n * DC_w_n * grad_XC_n;
      

      auto diffusiveflux_SALT_n = j_SALT_w_n;
      
      auto convectiveflux_SALT_w = omegaup_w_s * convectiveflux_SALT_w_s + omegaup_w_n * convectiveflux_SALT_w_n;
      if (veltype[Indices::BCId_water] == Indices::BCId_neumann){
        convectiveflux_SALT_w = 0.5 * ( rho_w_s * ( XC_s ) + rho_w_n * (XC_n)) * normalvelocity_w_n;
      }

      auto convectiveflux_SALT = -convectiveflux_SALT_w;//(omegaup_w_s * convectiveflux_SALT_w_s + omegaup_w_n * convectiveflux_SALT_w_n);
      auto diffusiveflux_SALT = omegaup_x_s * diffusiveflux_SALT_s + omegaup_x_n * diffusiveflux_SALT_n;
      if (veltype[Indices::BCId_salt] == Indices::BCId_neumann){
        diffusiveflux_SALT = 0.5 * (rho_w_s * Sw_s * DC_w_s + rho_w_n * Sw_n * DC_w_n) * grad_XC_n;
      }
     
      //  ACCCUMULATE RESIDUALS  //
			tmp=0.;
      
      // SALT-component-wise mass-balance
      tmp = Xc_conv_m * convectiveflux_SALT + Xc_diff_m * diffusiveflux_SALT ;
      double term_nipg_c_x = theta_x * (XC_s  - XC_n  );
      double term_penalty_c = penalty_factor_x * (XC_s  - XC_n);
      // diffusion term
      for (int i = 0; i < lfsv_XC_s.size(); i++)
      {
        r.accumulate(lfsv_XC_s, i, tmp * psi_XC_s[i] * factor);
      }
      
      // (non-)symmetric IP term
      for (int i = 0; i < lfsv_XC_s.size(); i++)
      {
        r.accumulate(lfsv_XC_s, i,  Xc_diff_m * term_nipg_c_x * 
            ( omegaup_x_s * rho_w_s * Sw_s * DC_w_s )* gradpsi_XC_s[i] * n_F_local * factor);//
      }
      // standard IP term integral
      for (int i = 0; i < lfsv_XC_s.size(); i++)
      {
        r.accumulate(lfsv_XC_s, i, term_penalty_c * psi_XC_s[i] * factor);
      }
    } // end of quadrature rule
  } // end of alpha_boundary
  
};
